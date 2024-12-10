# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 23:15:44 2018

@author: user


"""
#%%
from cmath import log
from glob import glob
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from photutils import detect_threshold
from photutils import DAOStarFinder
from photutils import CircularAnnulus as CircAn
from photutils import CircularAperture as CircAp
from photutils import aperture_photometry as APPHOT
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from astropy.wcs import WCS
import _Python_utilities
import _astro_utilities
from astropy.stats import sigma_clip, sigma_clipped_stats

plt.rcParams.update({'figure.max_open_warning': 0})


#######################################################
# for log file

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))
#######################################################


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
#######################################################
# read all files in base directory for processing

BASEDIR = "../Post_processing/M35_Light_-_2018-10-31_-_TMB130ss_STF-8300M_-_1bin/"
BASEDIR = "../Rne_2022/AMPELLA_Light_-_2022-09-06_-_GSON300_STF-8300M_-_1bin/"
BASEDIR = "../Rne_2022/INTERAMNIA_Light_-_2022-09-21_-_GSON300_STF-8300M_-_1bin/"

fullnames = _Python_utilities.getFullnameListOfallFiles(BASEDIR)
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))
######################################################

Result_dir = "DAO_Annul_result/"

if not os.path.exists('{0}'.format("{}{}".format(BASEDIR, Result_dir))):
    os.makedirs("{}{}".format(BASEDIR, Result_dir))
    print("{}{}is created".format(BASEDIR, Result_dir))

#%%
fullnames_light = [w for w in fullnames \
            if ("_bias_" not in w.lower()) \
                and ("_dark_" not in w.lower()) \
                    and ("_flat_" not in w.lower())
                    and (w.endswith(".fit") or w.endswith(".fits"))]
print ("len(fullnames_light): {}".format(len(fullnames_light)))

#%%
n = 0
for fullname in fullnames_light[:]:
#fullname = fullnames_light[1]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames_light), 
                                            (n/len(fullnames_light))*100, os.path.basename(__file__)))
    print ("Starting...\nfullname: {}".format(fullname))
 
    fullname_el = fullname.split("/")
    hdul = fits.open(fullname)
    hdr = hdul[0].header
    data = hdul[0].data

    img = hdul[0].data
    print("img: {}".format(img))
    print("img.shape: {}".format(img.shape))

    # Set WCS and print for your information
    w = WCS(hdr)
    print("WCS: {}".format(w))

    thresh = detect_threshold(data=img, nsigma=3)
    thresh = thresh[0][0]
    print('detect_threshold', thresh)

    #%%
    from photutils import DAOStarFinder
    FWHM   = 4
    DAOfind = DAOStarFinder(threshold=thresh, fwhm=FWHM, 
                            sharplo=0.2, sharphi=1.0,  # default values: sharplo=0.2, sharphi=1.0,
                            roundlo=-1.0, roundhi=1.0,  # default values -1 and +1
                            sigma_radius=1.5,           # default values 1.5
                            ratio=1.0,                  # 1.0: circular gaussian
                            exclude_border=True         # To exclude sources near edges
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
        DAOfound.write("{}{}{}_DAOStarfinder_fwhm{}.csv".\
                        format(BASEDIR, Result_dir, fullname_el[-1][:-4], FWHM), 
                        overwrite = True,
                        format='ascii.fast_csv')

        print('type(DAOfound): {}'.format(type(DAOfound)))
        print('DAOfound: {}'.format(DAOfound))

        #DAOcoord = (DAOfound['xcentroid'], DAOfound['ycentroid']) 
        DAOcoord = (DAOfound['xcentroid'], DAOfound['ycentroid']) 
        print('type(DAOcoord): {}'.format(type(DAOcoord)))
        print('DAOcoord: {}'.format(DAOcoord))

        DAOcoord_zip = list(zip(DAOfound['xcentroid'], DAOfound['ycentroid']))
        print('type(DAOcoord_zip): {}'.format(type(DAOcoord_zip)))
        print('DAOcoord_zip: {}'.format(DAOcoord_zip))
        print('DAOcoord_zip[0]: {}'.format(DAOcoord_zip[0]))
        print('DAOcoord_zip[1]: {}'.format(DAOcoord_zip[1]))

        # Save apertures as circular, 4 pixel radius, at each (X, Y)
        
        DAOapert = CircAp((DAOcoord_zip), r=4.)  
        #DAOapert = CircAp(DAOcoord, r=4.)
        print('type(DAOapert): {}'.format(type(DAOapert)))
        print('DAOapert: {}'.format(DAOapert))
        print('dir(DAOapert): {}'.format(dir(DAOapert)))

        DAOimgXY = np.array(DAOcoord)
        print('type(DAOimgXY): {}'.format(type(DAOimgXY)))
        print('DAOimgXY: {}'.format(DAOimgXY))

        DAOannul = CircAn(positions = (DAOcoord_zip), r_in = 4*FWHM, r_out = 6*FWHM) 
        print('type(DAOannul): {}'.format(type(DAOannul)))
        print('DAOannul: {}'.format(DAOannul))


        plt.figure(figsize=(24,16))
        ax = plt.gca()
        im = plt.imshow(img, vmax=thresh*4, origin='lower')
        DAOannul.plot(color='red', lw=2., alpha=0.4)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        plt.colorbar(im, cax=cax)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)

        plt.colorbar(im, cax=cax)

        ###########################################################
        # input some text for explaination. 
        plt.title("Result of DAOStarfinder", fontsize = 28, 
            ha='center')
        plt.annotate('filename: {}'.format(fullname_el[-1]), fontsize=10,
            xy=(1, 0), xytext=(-500, -40), va='top', ha='left',
            xycoords='axes fraction', textcoords='offset points')
                    
        plt.annotate('FWHM: {}'.format(FWHM), fontsize=10,
            xy=(1, 0), xytext=(-1200, -30), va='top', ha='left',
            xycoords='axes fraction', textcoords='offset points')
            
        plt.annotate('Sky threshold: {:02f}'.format(thresh), fontsize=10,
            xy=(1, 0), xytext=(-1200, -40), va='top', ha='left',
            xycoords='axes fraction', textcoords='offset points')

        plt.annotate('Number of star(s): {}'.format(len(DAOfound)), fontsize=10,
            xy=(1, 0), xytext=(-1200, -50), va='top', ha='left',
            xycoords='axes fraction', textcoords='offset points')

        plt.savefig("{}{}{}_DAOStarfinder_fwhm{}.png".\
                        format(BASEDIR, Result_dir, fullname_el[-1][:-4], FWHM))
        #plt.show()
        plt.close()

        
        ####################################
        if not 'RDNOISE' in hdul[0].header : 
            ronoise = 10.0 #STF-8300M
        else :
            ronoise = hdul[0].header['RDNOISE']
        
        if 'EGAIN' in hdul[0].header : 
            gain = hdul[0].header['EGAIN']
        elif 'GAIN' in hdul[0].header : 
            gain = hdul[0].header['GAIN']
        else :
            gain = 0.5
        
        ####################################
        N_stars = len(DAOfound)
        print('Star ID    msky  sky_std  nsky nrej')

        for i in range(0, N_stars):
            mask_annul = (DAOannul.to_mask(method='center'))[i]
            sky_apply  = mask_annul.multiply(img)
            sky_non0   = np.nonzero(sky_apply)
            sky_pixel  = sky_apply[sky_non0]
            msky, sky_std, nsky, nrej = sky_fit(sky_pixel, method='mode', mode_option='sex')
            print('{0}: {1:.5f} {2:.5f} {3:4d} {4:3d}'.format(i, msky, sky_std, nsky, nrej))
            plt.errorbar([i], msky, yerr=sky_std, capsize=3, marker='o', color='b')

        plt.xlabel('Star ID')
        plt.ylabel('msky +- sky_std')
        plt.xticks(np.arange(0, N_stars+5, step=5))
        plt.grid(ls=':')
        #plt.show()

        #%%
        # change img value to 16 bit Integer in ADU
        img = np.array(img*65536.0, dtype=np.uint16)
        ronoise = 9.00  # electrons
        gain = 2.54     # e/ADU
        N_star = len(DAOfound)

        mag_ann  = np.zeros(N_star)
        merr_ann = np.zeros(N_star)

        # aperture sum
        apert_sum = APPHOT(img, DAOapert, method='exact')['aperture_sum']
        ap_area   = DAOapert.area()
        #print(apert_sum)

        apert_result = 'ID, Msky, sky_std, Sky count Pixel_N, Sky reject Pixel_N, mag_ann, merr_ann\n'
        for i in range(0, N_star):
            # sky estimation
            mask_annul = (DAOannul.to_mask(method='center'))[i]
            sky_apply  = mask_annul.multiply(img)
            sky_non0   = np.nonzero(sky_apply)
            sky_pixel  = sky_apply[sky_non0]
            msky, sky_std, nsky, nrej = sky_fit(sky_pixel, method='mode', mode_option='sex')
            
            flux_star = apert_sum[i] - msky * ap_area  # total - sky
            flux_err  = np.sqrt(apert_sum[i] * gain    # Poissonian (star + sky)
                                + ap_area * ronoise**2 # Gaussian
                                + (ap_area * (gain * sky_std))**2 / nsky ) 
            mag_ann[i], merr_ann[i] = mag_inst(flux_star, flux_err)
            #print('{0:7d}: {1:.5f} {2:.5f} {3:4d} {4:3d} {5:.3f} {6:.3f}'.format(i, msky, sky_std, nsky, nrej, mag_ann[i], merr_ann[i]))
            apert_result += '{0}, {1:.5f}, {2:.5f}, {3:4d}, {4:3d}, {5:.3f}, {6:.3f}\n'.format(i, msky, sky_std, nsky, nrej, mag_ann[i], merr_ann[i])
            
            plt.errorbar(i, mag_ann[i], yerr=merr_ann[i], marker='x', ms=10, capsize=3)

        plt.xlabel('Star ID')
        plt.ylabel('Instrumental mag')
        plt.xticks(np.arange(0, N_stars+5, step=5))
        plt.grid(ls=':')
        plt.show()

        #print and save result
        print(apert_result)
        with open(dir_name+f_name[:-5]+'_APP_Annumus_result.csv', 'w') as f:
            f.write(apert_result)