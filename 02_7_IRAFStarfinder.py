# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 23:15:44 2018

@author: user

이 파일은 IRAF starfinder로 별을 찾아 정리해 줍니다.
"""
#%%
from cmath import log
from glob import glob
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from photutils import detect_threshold
from photutils import IRAFStarFinder
from photutils import CircularAnnulus as CircAn
from photutils import CircularAperture as CircAp
from photutils import aperture_photometry as APPHOT
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from astropy.wcs import WCS
import Python_utilities
import astro_utilities

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
#######################################################
# read all files in base directory for processing

base_dir = "../Post_processing/M35_Light_-_2018-10-31_-_TMB130ss_STF-8300M_-_1bin/"
base_dir = "../Rne_2022/AMPELLA_Light_-_2022-09-06_-_GSON300_STF-8300M_-_1bin/"
base_dir = "../Rne_2022/INTERAMNIA_Light_-_2022-09-21_-_GSON300_STF-8300M_-_1bin/"
base_dir = "../Rne_2022/KLEOPATRA_Light_-_2022-10-11_-_GSON300_STF-8300M_-_1bin/result_PI/calibrated/Light_BIN-1_EXPOSURE-20.00s_FILTER-L_Mono_fit/"

fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))
######################################################

Result_dir = "IRAFfinder_result/"

if not os.path.exists('{0}'.format("{}{}".format(base_dir, Result_dir))):
    os.makedirs("{}{}".format(base_dir, Result_dir))
    print("{}{}is created".format(base_dir, Result_dir))

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
    from photutils import IRAFStarFinder
    FWHM   = 8
    IRAFfind = IRAFStarFinder(fwhm=FWHM, threshold=thresh,
                          sigma_radius=1.5, minsep_fwhm=10.5,  # default values: sigma_radius=1.5, minsep_fwhm=2.5,
                          sharplo=0.5, sharphi=2.0,   # default values: sharplo=0.5, sharphi=2.0,
                          roundlo=0.0, roundhi=0.2,   # default values: roundlo=0.0, roundhi=0.2,
                          sky=None, exclude_border=False)  # default values: sky=None, exclude_border=Fal02_7_IRAFStarfinder.pyse))
                        
    # The IRAFStarFinder object ("IRAFfind") gets at least one input: the image.
    # Then it returns the astropy table which contains the aperture photometry results:
    IRAFfound = IRAFfind(img)
    print('{} star(s) founded by IRAFStarFinder...'.format(len(IRAFfound)))

    if len(IRAFfound)==0 :
        print ('No star was founded by IRAFStarFinder...\n'*3)
    else : 

        # Use the object "found" for aperture photometry:
        N_stars = len(IRAFfound)
        print('{} star(s) founded by IRAFStarFinder...'.format(N_stars))
        IRAFfound.pprint(max_width=1800)

        # save XY coordinates:
        IRAFfound.write("{}{}{}_IRAFStarFinder_fwhm{}.csv".\
                        format(base_dir, Result_dir, fullname_el[-1][:-4], FWHM), 
                        overwrite = True,
                        format='ascii.fast_csv')

        print('type(IRAFfound): {}'.format(type(IRAFfound)))
        print('IRAFfound: {}'.format(IRAFfound))

        #IRAFcoord = (IRAFfound['xcentroid'], IRAFfound['ycentroid']) 
        IRAFcoord = (IRAFfound['xcentroid'], IRAFfound['ycentroid']) 
        print('type(IRAFcoord): {}'.format(type(IRAFcoord)))
        print('IRAFcoord: {}'.format(IRAFcoord))

        IRAFcoord_zip = list(zip(IRAFfound['xcentroid'], IRAFfound['ycentroid']))
        print('type(IRAFcoord_zip): {}'.format(type(IRAFcoord_zip)))
        print('IRAFcoord_zip: {}'.format(IRAFcoord_zip))
        print('IRAFcoord_zip[0]: {}'.format(IRAFcoord_zip[0]))
        print('IRAFcoord_zip[1]: {}'.format(IRAFcoord_zip[1]))

        # Save apertures as circular, 4 pixel radius, at each (X, Y)
        IRAFapert = CircAp((IRAFcoord_zip), r=4.)  
        print('type(IRAFapert): {}'.format(type(IRAFapert)))
        print('IRAFapert: {}'.format(IRAFapert))
        print('dir(IRAFapert): {}'.format(dir(IRAFapert)))

        IRAFimgXY = np.array(IRAFcoord)
        print('type(IRAFimgXY): {}'.format(type(IRAFimgXY)))
        print('IRAFimgXY: {}'.format(IRAFimgXY))

        IRAFannul = CircAn(positions = (IRAFcoord_zip), r_in = 4*FWHM, r_out = 6*FWHM) 
        print('type(IRAFannul): {}'.format(type(IRAFannul)))
        print('IRAFannul: {}'.format(IRAFannul))


        plt.figure(figsize=(24,16))
        ax = plt.gca()
        im = plt.imshow(img, vmax=thresh*4, origin='lower')
        IRAFannul.plot(color='red', lw=2., alpha=0.4)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        plt.colorbar(im, cax=cax)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)

        plt.colorbar(im, cax=cax)

        ###########################################################
        # input some text for explaination. 
        plt.title("Result of IRAFStarFinder", fontsize = 28, 
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

        plt.annotate('Number of star(s): {}'.format(len(IRAFfound)), fontsize=10,
            xy=(1, 0), xytext=(-1200, -50), va='top', ha='left',
            xycoords='axes fraction', textcoords='offset points')

        plt.savefig("{}{}{}_IRAFStarFinder_fwhm{}.png".\
                        format(base_dir, Result_dir, fullname_el[-1][:-4], FWHM))
        #plt.show()
        plt.close()
  