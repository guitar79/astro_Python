# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 23:15:44 2018

@author: user
"""
#%%
from glob import glob
import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import pandas as pd 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.stats import sigma_clip, sigma_clipped_stats
from photutils import DAOStarFinder
from photutils import CircularAperture as CircAp
from photutils import CircularAnnulus as CircAn
from photutils import aperture_photometry as APPHOT
import Python_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))

base_dr = "../Post_processing/M35_Light_-_2018-10-31_-_TMB130ss_STF-8300M_-_1bin/"

### make all fits file list...
fullnames = Python_utilities.getFullnameListOfallFiles("{}/input".format(base_dr))
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))

c_method = 'median'
master_dr = "master_files/"
reduced_dr = "readuced_files/"
result_dr = "result_files/"

if not os.path.exists('{0}'.format("{}{}".format(base_dr, result_dr))):
    os.makedirs("{}{}".format(base_dr, result_dr))
    print("{}{}is created".format(base_dr, result_dr))

fullnames_light = [w for w in fullnames \
            if ("_bias_" not in w.lower()) \
                and ("_dark_" not in w.lower()) \
                    and ("_flat_" not in w.lower())]
print ("len(fullnames_light): {}".format(len(fullnames_light)))

#%%
#Returns magnitude from flux.
def mag_inst(flux, ferr):
    m_inst = -2.5 * np.log10(flux)
    merr   = 2.5/ np.log(10) * ferr / flux
    return m_inst, merr

def sky_fit(all_sky, method='mode', sky_nsigma=3, sky_iter=5, \
            mode_option='sex', med_factor=2.5, mean_factor=1.5):
    '''
    Estimate sky from given sky values.

    Parameters
    ----------
    all_sky : ~numpy.ndarray
        The sky values as numpy ndarray format. It MUST be 1-d for proper use.
    method : {"mean", "median", "mode"}, optional
        The method to estimate sky value. You can give options to "mode"
        case; see mode_option.
        "mode" is analogous to Mode Estimator Background of photutils.
    sky_nsigma : float, optinal
        The input parameter for sky sigma clipping.
    sky_iter : float, optinal
        The input parameter for sky sigma clipping.
    mode_option : {"sex", "IRAF", "MMM"}, optional.
        sex  == (med_factor, mean_factor) = (2.5, 1.5)
        IRAF == (med_factor, mean_factor) = (3, 2)
        MMM  == (med_factor, mean_factor) = (3, 2)

    Returns
    -------
    sky : float
        The estimated sky value within the all_sky data, after sigma clipping.
    std : float
        The sample standard deviation of sky value within the all_sky data,
        after sigma clipping.
    nsky : int
        The number of pixels which were used for sky estimation after the
        sigma clipping.
    nrej : int
        The number of pixels which are rejected after sigma clipping.
    -------

    '''
    sky = all_sky.copy()
    if method == 'mean':
        return np.mean(sky), np.std(sky, ddof=1)

    elif method == 'median':
        return np.median(sky), np.std(sky, ddof=1)

    elif method == 'mode':
        sky_clip   = sigma_clip(sky, sigma=sky_nsigma, iters=sky_iter)
        sky_clipped= sky[np.invert(sky_clip.mask)]
        nsky       = np.count_nonzero(sky_clipped)
        mean       = np.mean(sky_clipped)
        med        = np.median(sky_clipped)
        std        = np.std(sky_clipped, ddof=1)
        nrej       = len(all_sky) - len(sky_clipped)

        if nrej < 0:
            raise ValueError('nrej < 0: check the code')

        if nrej > nsky: # rejected > survived
            raise Warning('More than half of the pixels rejected.')

        if mode_option == 'IRAF':
            if (mean < med):
                sky = mean
            else:
                sky = 3 * med - 2 * mean

        elif mode_option == 'MMM':
            sky = 3 * med - 2 * mean

        elif mode_option == 'sex':
            if (mean - med) / std > 0.3:
                sky = med
            else:
                sky = (2.5 * med) - (1.5 * mean)
        else:
            raise ValueError('mode_option not understood')

        return sky, std, nsky, nrej

#%%
n = 0
for fullname in fullnames_light[10:11]:
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames_light), 
                                            (n/len(fullnames_light))*100, os.path.basename(__file__)))
    print ("Starting...\nfullname: {}".format(fullname))
    
    fullname_el = fullname.split("/")
    
    hdu = fits.open(fullname)
    
    #img = np.array(hdu[0].data/65536.0, dtype=np.float32)
    img = hdu[0].data
    print("img: {}".format(img))
    
    from photutils import detect_threshold

    #thresh = detect_threshold(data=img.data, snr=3)
    thresh = detect_threshold(data=img.data, nsigma=3)
    thresh = thresh[0][0]

    print('detect_threshold', thresh)
    
    from photutils import DAOStarFinder
    
    FWHM   = 2.5
    DAOfind = DAOStarFinder(threshold=thresh, 
                        fwhm=FWHM 
                        #sharplo=0.5, sharphi=2.0,  # default values: sharplo=0.2, sharphi=1.0,
                        #roundlo=0.0, roundhi=0.2,  # default values: roundlo=-1.0, roundhi=1.0,
                        #sigma_radius=1.5,          # default values: sigma_radius=1.5,
                        #ratio=0.5,                 # 1.0: circular gaussian:  ratio=1.0,
                        #sky=None, exclude_border=False)       # To exclude sources near edges : exclude_border=True
                        )
    # The DAOStarFinder object ("DAOfind") gets at least one input: the image.
    # Then it returns the astropy table which contains the aperture photometry results:
    DAOfound = DAOfind(img)
    print('{} star(s) founded by DAOStarFinder...'.format(len(DAOfound)))
    
    if len(DAOfound)==0 :
        print ('No star was founded by DAOStarFinder...\n'*3)
    else : 
        # Use the object "found" for aperture photometry:
        N_stars = len(DAOfound)
        print('{} star(s) founded by DAOStarFinder...'.format(N_stars))
        DAOfound.pprint(max_width=1800)
        
        # save XY coordinates:
        DAOfound.write("{}{}{}_DAOStarfinder.csv".\
                        format(base_dr, result_dr, fullname_el[-1][:-4]), 
                        overwrite = True,
                        format='ascii.fast_csv')
        
        print('type(DAOfound): {}'.format(type(DAOfound)))
        print('DAOfound: {}'.format(DAOfound))

        DAOcoord = (DAOfound['xcentroid'], DAOfound['ycentroid']) 
        print('type(DAOcoord): {}'.format(type(DAOcoord)))
        print('DAOcoord: {}'.format(DAOcoord))

        DAOcoord_zip = zip(DAOfound['xcentroid'], DAOfound['ycentroid'])  
        print('type(DAOcoord_zip): {}'.format(type(DAOcoord_zip)))
        print('DAOcoord_zip: {}'.format(DAOcoord_zip))
        print('DAOcoord_zip[0]: {}'.format(DAOcoord_zip[0]))
        print('DAOcoord_zip[1]: {}'.format(DAOcoord_zip[1]))

        # Save apertures as circular, 4 pixel radius, at each (X, Y)
        DAOapert = CircAp(DAOcoord_zip, r=4.)  
        print('type(DAOapert): {}'.format(type(DAOapert)))
        print('DAOapert: {}'.format(DAOapert))
        
        DAOimgXY = np.array(DAOcoord)
        print('type(DAOimgXY): {}'.format(type(DAOimgXY)))
        print('DAOimgXY: {}'.format(DAOimgXY))

        DAOannul = CircAn(DAOcoord_zip, r_in = 4*FWHM, r_out = 6*FWHM) 
        print('type(DAOannul): {}'.format(type(DAOannul)))
        print('DAOannul: {}'.format(DAOannul))
        
        plt.figure(figsize=(16,12))
        ax = plt.gca()
        im = plt.imshow(img, vmax=thresh*4, origin='lower')
        DAOannul.plot(color='red', lw=2., alpha=0.7)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        plt.colorbar(im, cax=cax)
        plt.savefig("{}{}{}_DAOstarfinder_Annulus_result_all_stars.png".\
                        format(base_dr, result_dr, fullname_el[-1][:-4]), 
                        overwrite=True)
        plt.show()