# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 23:15:44 2018

@author: user
"""
#%%
from glob import glob
import numpy as np
import os
import astropy.units as u
from ccdproc import CCDData, ccd_process
import _Python_utilities

from photutils import DAOStarFinder
from photutils import CircularAperture as CircAp
from photutils import CircularAnnulus as CircAn

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))

BASEDIR = "../Post_processing/M35_Light_-_2018-10-31_-_TMB130ss_STF-8300M_-_1bin/"

### make all fits file list...
fullnames = _Python_utilities.getFullnameListOfallFiles("{}/input".format(BASEDIR))
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))

c_method = 'median'
master_dir = "master_files/"
reduced_dir = "readuced_files/"
result_dir = "result_files/"

if not os.path.exists('{0}'.format("{}{}".format(BASEDIR, result_dir))):
    os.makedirs("{}{}".format(BASEDIR, result_dir))
    print("{}{}is created".format(BASEDIR, result_dir))

fullnames_light = [w for w in fullnames \
            if ("_bias_" not in w.lower()) \
                and ("_dark_" not in w.lower()) \
                    and ("_flat_" not in w.lower())]
print ("len(fullnames_light): {}".format(len(fullnames_light)))

#%%
###
fullname = fullnames_light[0]

fullname_el = fullname.split("/")

hdu = fits.open("{}".format(fullname))

img = hdu[0].data
print("img: {}".format(img))

#%%
from astropy.stats import mad_std, sigma_clipped_stats

sky_th = 5    # sky_th * sky_sigma will be used for detection lower limit : sky_th = 5
sky_s  = mad_std(img)
thresh = sky_s * sky_th
print('sky_s x sky_th = threshold')
print('{0:8.6f} x {1:4d}   =   {2:8.6f}\n'.format(sky_s, sky_th, thresh))

# What if we do "sigma-clip" than MAD?
sky_mean, sky_median, sky_std = sigma_clipped_stats(img,
                        sigma=3.0) # default is 3-sigma, 5 iters
print(sky_mean, sky_median, sky_std)

thresh_sc = sky_std * sky_th
print('3 sigma 5 iters clipped case:')
print('{0:8.6f} x {1:4d}   =   {2:8.6f}\n'.format(sky_std, sky_th, thresh_sc))

#%%
from photutils import detect_threshold

#thresh = detect_threshold(data=img.data, snr=3)
thresh = detect_threshold(data=img.data, nsigma=3)
thresh = thresh[0][0]

print('detect_threshold', thresh)

#%%
FWHM   = 2.5

DAOfind = DAOStarFinder(threshold=thresh, 
                        fwhm=FWHM 
                        #sharplo=0.5, sharphi=2.0,  # default values: sharplo=0.2, sharphi=1.0,
                        #roundlo=0.0, roundhi=0.2,  # default values: roundlo=-1.0, roundhi=1.0,
                        #sigma_radius=1.5,          # default values: sigma_radius=1.5,
                        #ratio=0.5,                 # 1.0: circular gaussian:  ratio=1.0,
                        #sky=None, exclude_border=False)       # To exclude sources near edges : exclude_border=True
                        )
# The DAOStarFinder object ("find") gets at least one input: the image.
# Then it returns the astropy table which contains the aperture photometry results:
DAOfound = DAOfind(img)
print('{} stars founded by DAOStarFinder...'.format(len(DAOfound)))

#%%
if len(DAOfound)==0 :
    print ('No star founded using DAOStarFinder\n'*3)
else : 
    # Use the object "found" for aperture photometry:
    # save XY coordinates:
    DAOfound.pprint(max_width=1800)
    
    DAOfound.write("{}{}{}_DAOStarfinder.csv".\
                        format(BASEDIR, result_dir, fullname_el[-1][:-4]), 
                        overwrite = True,
                        format='ascii.fast_csv')
    
    #DAOcoord = (DAOfound['xcentroid'], DAOfound['ycentroid']) 
    DAOcoord = zip(DAOfound['xcentroid'], DAOfound['ycentroid'])   

    # Save apertures as circular, 4 pixel radius, at each (X, Y)
    DAOapert = CircAp(DAOcoord, r=4.)  
    print('DAOapert\n ', DAOapert)
    
    DAOimgXY = np.array(DAOcoord)
    print('DAOimgXY \n', DAOimgXY)
    
    print('type(DAOimgXY): {}'.format(type(DAOimgXY)))
    print('DAOimgXY: {}'.format(DAOimgXY))

    DAOannul = CircAn(positions = (DAOcoord_zip), r_in = 4*FWHM, r_out = 6*FWHM) 
    print('type(DAOannul): {}'.format(type(DAOannul)))
    print('DAOannul: {}'.format(DAOannul))
    
    plt.figure(figsize=(24, 18))
    ax = plt.gca()
    im = plt.imshow(img, vmax=thresh*4, origin='lower')

    DAOapert.plot(color='red', lw=2., alpha=0.3)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    plt.colorbar(im, cax=cax)
    
    plt.show()
    
#%%
import numpy as np
import matplotlib.pyplot as plt
from photutils import IRAFStarFinder
from photutils import CircularAperture as CircAp

IRAFfind = IRAFStarFinder(threshold=thresh,
                            fwhm=FWHM
                            #sigma_radius=1.5, minsep_fwhm=1.0,  # default values: sigma_radius=1.5, minsep_fwhm=2.5,
                            #sharplo=0.5, sharphi=2.0,   # default values: sharplo=0.5, sharphi=2.0,
                            #roundlo=0.0, roundhi=0.2,   # default values: roundlo=0.0, roundhi=0.2,
                            #sky=None, exclude_border=False)  # default values: sky=None, exclude_border=False)
                            )

# The DAOStarFinder object ("find") gets at least one input: the image.
# Then it returns the astropy table which contains the aperture photometry results:
IRAFfound = IRAFfind(img)
print('{} stars founded by IRAFStarFinder...'.format(len(IRAFfound)))

#%%
if len(IRAFfound)==0 :
    print ('No star founded using IRAFStarFinder')
else : 
    IRAFfound.pprint(max_width=1800)
    IRAFfound.write("{}{}{}_IRAFStarfinder.csv".\
                        format(BASEDIR, result_dir, fullname_el[-1][:-4]), 
                        overwrite =True,
                        format='ascii.fast_csv')
    
    # Use the object "found" for aperture photometry:
    # save XY coordinates:
    print (len(IRAFfound), 'stars founded')
    
    #IRAFcoord = (IRAFfound['xcentroid'], IRAFfound['ycentroid']) 
    IRAFcoord = zip(IRAFfound['xcentroid'], IRAFfound['ycentroid'])   
    
    # Save apertures as circular, 4 pixel radius, at each (X, Y)
    IRAFapert = CircAp(IRAFcoord, r=4.)  
    #print('IRAFapert\n ', IRAFapert)
    
    IRAFimgXY = np.array(IRAFcoord)
    #print('IRAFimgXY \n', IRAFimgXY)
            
    plt.figure(figsize=(12,12))
    ax = plt.gca()
    im = plt.imshow(img, vmax=65536, origin='lower')
    IRAFapert.plot(color='red', lw=2., alpha=0.7)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    plt.colorbar(im, cax=cax)
    
    plt.show()