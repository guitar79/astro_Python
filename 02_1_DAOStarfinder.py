# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 23:15:44 2018

@author: user

이 파일은 DAO starfinder로 별을 찾아 정리해 줍니다.
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
STARfinder_result_dir = "DAOfinder_result"

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

                DAOfind = DAOStarFinder(
                                        fwhm = FWHM, 
                                        threshold = thresh, 
                                        sharplo = 0.2, sharphi = 1.0,  # default values: sharplo=0.2, sharphi=1.0,
                                        roundlo = -1.0, roundhi = 1.0,  # default values -1 and +1
                                        sigma_radius = 1.5,           # default values 1.5
                                        ratio = 1.0,                  # 1.0: circular gaussian
                                        exclude_border = True         # To exclude sources near edges
                                        )
                # The DAOStarFinder object ("DAOfind") gets at least one input: the image.
                # Then it returns the astropy table which contains the aperture photometry results:
                DAOfound = DAOfind(img)
                print('{} star(s) founded by DAOStarFinder...'.format(len(DAOfound)))

                #%%
                if len(DAOfound)==0 :
                    print ('No star was founded by DAOStarFinder...\n'*3)
                else : 

                    # Use the object "found" for aperture photometry:
                    N_stars = len(DAOfound)
                    print('{} star(s) founded by DAOStarFinder...'.format(N_stars))
                    DAOfound.pprint(max_width=1800)

                    # save XY coordinates:
                    DAOfound.write("{}/{}_DAOStarfinder_fwhm{}.csv".\
                                    format(RESULTDIR, fpath.stem, FWHM), 
                                    overwrite = True,
                                    format='ascii.fast_csv')
                    #%%
                    print('type(DAOfound): {}'.format(type(DAOfound)))
                    print('DAOfound: {}'.format(DAOfound))

                    DAOcoord = np.array([DAOfound['xcentroid'], DAOfound['ycentroid']]).T
                    print('type(DAOcoord): {}'.format(type(DAOcoord)))
                    print('DAOcoord: {}'.format(DAOcoord))

                    #%%
                    # Save apertures as circular, 4 pixel radius, at each (X, Y)
                    DAOapert = CAp((DAOcoord), r=4.)  
                    print('type(DAOapert): {}'.format(type(DAOapert)))
                    print('DAOapert: {}'.format(DAOapert))
                    print('dir(DAOapert): {}'.format(dir(DAOapert)))

                    DAOannul = CAn(positions = (DAOcoord), r_in = 4*FWHM, r_out = 6*FWHM) 
                    print('type(DAOannul): {}'.format(type(DAOannul)))
                    print('DAOannul: {}'.format(DAOannul))
                    
                    #%%
                    plt.figure(figsize=(20,20))
                    ax = plt.gca()

                    ###########################################################
                    # input some text for explaination. 
                    plt.title("Result of DAOStarfinder", fontsize = 28, 
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

                    plt.annotate('Number of star(s): {}'.format(len(DAOfound)), fontsize=10,
                        xy=(1, 0), xytext=(-1100, -50), va='top', ha='left',
                        xycoords='axes fraction', textcoords='offset points')

                    im = plt.imshow(img, 
                                    vmin = thresh, 
                                    vmax = thresh * 3,
                                    #zscale=True,
                                    origin='lower'
                                    )

                    DAOannul.plot(color='red', lw=2., alpha=0.4)
                    
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="3%", pad=0.05)
                    plt.colorbar(im, cax=cax)

                    plt.savefig(
                                "{}/{}_DAOStarfinder_fwhm{}.png".\
                                    format(RESULTDIR, fpath.stem, FWHM)
                                )
                    print("{}/{}_DAOStarfinder_fwhm{}.png is created...".\
                                    format(RESULTDIR, fpath.stem, FWHM))
                    #plt.show()
                    plt.close() 

            except Exception as err:
                print('{0} with {1} '.format(err, fpath.name))