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

from astropy.nddata import CCDData
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip, sigma_clipped_stats

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

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
#import ysvisutilpy as yvu

import _Python_utilities
import _astro_utilities
import _tool_visualization as vis

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
BASEDIR = Path(r"r:\CCD_obs")
BASEDIR = Path("/mnt/Rdata/OBS_data") 
#BASEDIR = Path("/mnt/OBS_data") 
DOINGDIR = Path(BASEDIR/ "RnE_2022/GSON300_STF-8300M")
#DOINGDIR = Path(BASEDIR/ "A1_CCD_new_files1")

#DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(DOINGDIR))
DOINGDIRs = sorted([x for x in DOINGDIR.iterdir() if x.is_dir()])
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################

#%%
#####################################################################
# Observed location
LOCATION = dict(lon = 127.0, lat = 37.3, elevation = 130)

# It is used as a rough estimate, so no need to be accurate:
PIX2ARCSEC = 0.98 * u.arcsec

# Used for any `astropy.SkyCoord` object:
SKYC_KW = dict(unit = u.deg, frame = 'icrs')

# Initial guess of FWHM in pixel
FWHM_INIT = 6

# Photometry parameters
R_AP = 1.5 * FWHM_INIT # Aperture radius
R_IN = 4 * FWHM_INIT   # Inner radius of annulus
R_OUT = 6 * FWHM_INIT  # Outer radius of annulus
#######################################################
#%%
for DOINGDIR in DOINGDIRs[:] :
    #DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)
    fits_in_dir = sorted(list(DOINGDIR.glob('*.fit*')))
    #print("fits_in_dir", fits_in_dir)
    print("len(fits_in_dir)", len(fits_in_dir))

    if len(fits_in_dir) == 0 :
        print(f"There is no fits fils in {DOINGDIR}")
        pass
    else : 
        print(f"Starting: {str(DOINGDIR.parts[-1])}")

        MASTERDIR = DOINGDIR / _astro_utilities.master_dir
        REDUCEDDIR = DOINGDIR / _astro_utilities.reduced_dir2
        SOLVEDDIR = DOINGDIR / _astro_utilities.solved_dir2
        DAORESULTDIR = DOINGDIR / _astro_utilities.DAOfinder_result_dir
        
        if not DAORESULTDIR.exists():
            os.makedirs("{}".format(str(DAORESULTDIR)))
            print("{} is created...".format(str(DAORESULTDIR)))

        #summary = yfu.make_summary(BASEDIR/"*.fit*")
        summary = yfu.make_summary(SOLVEDDIR/"*.fit*")

        if summary.empty:
                print("The dataframe(summary) is empty")
                pass
        else:
            print("len(summary):", len(summary))
            print("summary:", summary)

            #%%
            df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
            df_light = df_light.reset_index(drop=True)
            print("df_light:\n{}".format(df_light))

            for _, row  in df_light.iterrows():
                fpath = Path(row["file"])
                print("type(fpath)", type(fpath))
                print("fpath", fpath)
                
                try:
                    # hdul = fits.open(fpath)
                    # hdr = hdul[0].header
                    # img = hdul[0].data
                    # print("img: {}".format(img))
                    # print("img.shape: {}".format(img.shape))
                    ccd = CCDData.read(fpath, unit="adu")
                    print("ccd.data: ", ccd.data)
                    print("ccd.data.shape: ", ccd.data.shape)

                    # Set WCS and print for your information
                    w = ccd.wcs
                    print("WCS: ", w)
                    FWHM = FWHM_INIT
                    avg, med, std = sigma_clipped_stats(ccd.data)  # by default, 3-sigma 5-iteration.
                    thresh = 5.*std
                    DAOfind = DAOStarFinder(
                                            fwhm = FWHM, 
                                            threshold = thresh, 
                                            sharplo = 0.2, sharphi = 1.0,  # default values: sharplo=0.2, sharphi=1.0,
                                            roundlo = -1.0, roundhi = 1.0,  # default values -1 and +1
                                            sigma_radius = 1.5,           # default values 1.5
                                            ratio = 1.0,                  # 1.0: circular gaussian
                                            exclude_border = True         # To exclude sources near edges
                                            )

                    DAOfound = DAOfind(ccd.data - med)
                    DAOfound['RADEC'] = w.pixel_to_world(DAOfound['xcentroid'], 
                                                        DAOfound['ycentroid'])  

                    print(f"Threshold = {DAOfind.threshold:.2f} counts")
                    print(f"Found {len(DAOfound)} sources")

                    # save XY coordinates:
                    DAOfound.write(f"{DAORESULTDIR}/{fpath.stem}_DAOStarfinder_fwhm_{FWHM}.csv",
                                    overwrite = True,
                                    #format='ascii.fast_csv'
                                    )
                    print(f"{DAORESULTDIR}/{fpath.stem}_DAOStarfinder_fwhm_{FWHM}.csv is created...")

                    pos = np.transpose((DAOfound['xcentroid'], DAOfound['ycentroid']))
                    aps1 = CAp(pos, r=R_IN)
                    aps2 = CAp(pos, r=R_OUT)
                    fig, axs = plt.subplots(1, 1, figsize=(10, 7), sharex=False, sharey=False, gridspec_kw=None)

                    im = vis.norm_imshow(axs, ccd.data, zscale=True)
                    aps1.plot(color='r', lw=1, alpha=0.5)
                    aps2.plot(color='cyan', lw=1, alpha=0.5)

                    ###########################################################
                    # input some text for explaination. 
                    plt.title("Result of DAOStarfinder", fontsize = 14, 
                        ha='center')

                    plt.annotate(f'filename: {fpath.stem}', fontsize=8,
                        xy=(1, 0), xytext=(-10, -30), va='top', ha='right',
                        xycoords='axes fraction', textcoords='offset points')
                                
                    plt.annotate(f'FWHM: {FWHM}', fontsize=8,
                        xy=(0, 0), xytext=(-10, -30), va='top', ha='left',
                        xycoords='axes fraction', textcoords='offset points')
                        
                    plt.annotate(f'Sky threshold: {thresh:.02f}', fontsize=8,
                        xy=(0, 0), xytext=(-10, -40), va='top', ha='left',
                        xycoords='axes fraction', textcoords='offset points')

                    plt.annotate(f'Number of star(s): {len(DAOfound)}', fontsize=8,
                        xy=(0, 0), xytext=(-10, -50), va='top', ha='left',
                        xycoords='axes fraction', textcoords='offset points')

                    divider = make_axes_locatable(axs)
                    cax = divider.append_axes("right", size="3%", pad=0.05)
                    plt.colorbar(im, cax=cax)
                    plt.tight_layout()
                    plt.savefig(f"{DAORESULTDIR}/{fpath.stem}_DAOStarfinder_fwhm_{FWHM}.png")
                    print(f"{DAORESULTDIR}/{fpath.stem}_DAOStarfinder_fwhm_{FWHM}.png is created...")

                    #plt.show()
                    plt.close()

                except Exception as err:
                    _Python_utilities.write_log(err_log_file, err)
                    #print('{0} with {1} '.format(err, fpath.name))