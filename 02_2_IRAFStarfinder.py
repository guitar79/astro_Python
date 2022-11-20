# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 23:15:44 2018

@author: user

이 파일은 IRAF starfinder로 별을 찾아 정리해 줍니다.
"""
#%%
import os
from glob import glob
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

import Python_utilities
import astro_utilities

from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy.wcs import WCS

from photutils import DAOStarFinder
from photutils import IRAFStarFinder
from photutils.aperture import CircularAperture as CAp
from photutils.aperture import CircularAnnulus as CAn
from photutils import detect_threshold
from photutils.centroids import centroid_com
#from photutils import aperture_photometry as apphot

import warnings

from ccdproc import CCDData, ccd_process

from astropy.time import Time
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord

plt.rcParams.update({'figure.max_open_warning': 0})

#%%
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
BASEDIR = "../RnE_2022/"
BASEDIR = "../RnE_2022/RiLA600_STX-16803_2bin/"

c_method = 'median'
master_dir = "master_files_ys"
reduced_dir = "reduced"
solved_dir = "solved"
STARfinder_result_dir = "IRAFfinder_result"

#%%
BASEDIRs = sorted(Python_utilities.getFullnameListOfsubDir(BASEDIR))
print ("BASEDIRs: {}".format(BASEDIRs))

for BASEDIR in BASEDIRs :
    print ("Starting...\n{}".format(BASEDIR))

    BASEDIR = Path(BASEDIR)
    RESULTDIR = BASEDIR / STARfinder_result_dir
    SOLVEDDIR = BASEDIR / solved_dir

    if not RESULTDIR.exists():
        os.makedirs("{}".format(str(RESULTDIR)))
        print("{} is created...".format(str(RESULTDIR)))

    summary = yfu.make_summary(SOLVEDDIR / "*.fits")
    #print(summary)
    #print("len(summary):", len(summary))
    #print(summary["file"][0])

    #%%
    if len(summary) != 0:
        n = 0
        for fname in summary["file"]:
            #fpath = summary["file"][1]
            n += 1
            print('#'*40,
                "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(summary["file"]), 
                                                    (n/len(summary["file"]))*100, os.path.basename(__file__)))
            print ("Starting...\nfpath: {}".format(fname))
            
            fpath = Path(fname)
            hdul = fits.open(fpath)
            hdr = hdul[0].header
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
            try:
                
                FWHM   = 6
                
                IRAFfind = IRAFStarFinder(
                                        fwhm = FWHM, 
                                        threshold = thresh,
                                        sigma_radius = 1.5, minsep_fwhm = 2.5,  # default values: sigma_radius=1.5, minsep_fwhm=2.5,
                                        sharplo = 0.5, sharphi = 2.0,   # default values: sharplo=0.5, sharphi=2.0,
                                        roundlo = 0.0, roundhi = 0.2,   # default values: roundlo=0.0, roundhi=0.2,
                                        sky = None,                     # default values: sky=None
                                        exclude_border = True  # default values: exclude_border=False
                                        )
                                    
                # The IRAFStarFinder object ("IRAFfind") gets at least one input: the image.
                # Then it returns the astropy table which contains the aperture photometry results:
                IRAFfound = IRAFfind(img)
                print('{} star(s) founded by IRAFStarFinder...'.format(len(IRAFfound)))

                #%%
                if len(IRAFfound)==0 :
                    print ('No star was founded by IRAFStarFinder...\n'*3)
                else : 

                    # Use the object "found" for aperture photometry:
                    N_stars = len(IRAFfound)
                    print('{} star(s) founded by IRAFStarFinder...'.format(N_stars))
                    IRAFfound.pprint(max_width=1800)

                    # save XY coordinates:
                    IRAFfound.write("{}/{}_IRAFStarfinder_fwhm{}.csv".\
                                    format(RESULTDIR, fpath.stem, FWHM), 
                                    overwrite = True,
                                    format='ascii.fast_csv')

                    print('type(IRAFfound): {}'.format(type(IRAFfound)))
                    print('IRAFfound: {}'.format(IRAFfound))

                    IRAFcoord = np.array([IRAFfound['xcentroid'], IRAFfound['ycentroid']]).T
                    print('type(IRAFcoord): {}'.format(type(IRAFcoord)))
                    print('IRAFcoord: {}'.format(IRAFcoord))

                    #%%
                    # Save apertures as circular, 4 pixel radius, at each (X, Y)
                    IRAFapert = CAp((IRAFcoord), r=4.)  
                    print('type(IRAFapert): {}'.format(type(IRAFapert)))
                    print('IRAFapert: {}'.format(IRAFapert))
                    print('dir(IRAFapert): {}'.format(dir(IRAFapert)))

                    IRAFannul = CAn(positions = (IRAFcoord), r_in = 4*FWHM, r_out = 6*FWHM) 
                    print('type(IRAFannul): {}'.format(type(IRAFannul)))
                    print('IRAFannul: {}'.format(IRAFannul))

                    #%%
                    plt.figure(figsize=(20,20))
                    ax = plt.gca()

                    ###########################################################
                    # input some text for explaination. 
                    plt.title("Result of IRAFStarfinder", fontsize = 28, 
                        ha='center')

                    plt.annotate('filename: {}'.format(fpath.stem), fontsize=10,
                        xy=(1, 0), xytext=(-500, -40), va='top', ha='left',
                        xycoords='axes fraction', textcoords='offset points')
                                
                    plt.annotate('FWHM: {}'.format(FWHM), fontsize=10,
                        xy=(1, 0), xytext=(-1100, -30), va='top', ha='left',
                        xycoords='axes fraction', textcoords='offset points')
                        
                    plt.annotate('Sky threshold: {:02f}'.format(thresh), fontsize=10,
                        xy=(1, 0), xytext=(-1100, -40), va='top', ha='left',
                        xycoords='axes fraction', textcoords='offset points')

                    plt.annotate('Number of star(s): {}'.format(len(IRAFfound)), fontsize=10,
                        xy=(1, 0), xytext=(-1100, -50), va='top', ha='left',
                        xycoords='axes fraction', textcoords='offset points')

                    im = plt.imshow(img, 
                                    vmin = thresh, 
                                    vmax = thresh * 3,
                                    #zscale=True,
                                    origin='lower'
                                    )

                    IRAFannul.plot(color='red', lw=2., alpha=0.4)
                    
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="3%", pad=0.05)
                    plt.colorbar(im, cax=cax)

                    plt.savefig(
                                "{}/{}_IRAFStarfinder_fwhm{}.png".\
                                    format(RESULTDIR, fpath.stem, FWHM)
                                )
                    print("{}/{}_IRAFStarfinder_fwhm{}.png is created...".\
                                    format(RESULTDIR, fpath.stem, FWHM))
                    #plt.show()
                    plt.close() 

            except Exception as err:
                print('{0} with {1} '.format(err, fpath.name))
    