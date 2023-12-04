# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com
"""
#%%
from glob import glob
from pathlib import Path
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
#import ysvisutilpy as yvu

import _astro_utilities
import _Python_utilities

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
BASEDIR = Path("/mnt/Rdata/OBS_data") 

DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/STX-16803_1bin' )
#DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/STF-8300M_1bin/LIGHT_GSON300')
#DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/STL-11000M_1bin/' )
#DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/' )

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
remove = 'BIAS'
DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
remove = 'DARK'
DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
remove = 'FLAT'
DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################
#%%
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)
    
    REDUCEDDIR2 = DOINGDIR / _astro_utilities.reduced_dir2
    
    fits_in_dir = sorted(list(REDUCEDDIR2.glob('*.fit*')))
    #print("fits_in_dir", fits_in_dir)
    print("len(fits_in_dir)", len(fits_in_dir))

    if len(fits_in_dir) == 0 :
        print(f"There is no fits fils in {DOINGDIR}")
        pass
    else : 
        print(f"Starting: {str(DOINGDIR.parts[-1])}")

        summary = yfu.make_summary(DOINGDIR/"*.fit*")
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])
        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)
        print("df_light:\n{}".format(df_light))

        for _, row  in df_light.iterrows():

            fpath = Path(row["file"])
            hdul = fits.open(fpath)
            
            if 'PIXSCALE' in hdul[0].header:
                PIXc = hdul[0].header['PIXSCALE']
            else : 
                PIXc = _astro_utilities.calPixScale(hdul[0].header['FOCALLEN'], 
                                                    hdul[0].header['XPIXSZ'],
                                                    hdul[0].header['XBINNING'])
            print("PIXc : ", PIXc)
            hdul.close()

            SOLVE, ASTAP, LOCAL = _astro_utilities.checkPSolve(fpath)
            print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)

            if ASTAP :
                print(f"{fpath.name} is solved by ASTAP")
            else : 
                print(f"{fpath.name} is solving now by ASTAP")
                solved = _astro_utilities.ASTAPSolver(fpath, 
                                                        #str(SOLVEDDIR), 
                                                        downsample = 2,
                                                        pixscale = PIXc,
                                                                )

            SOLVE, ASTAP, LOCAL = _astro_utilities.checkPSolve(fpath)
            print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)

            if LOCAL :
                print(f"{fpath.name} is solved by LOCAL")
            else : 
                print(f"{fpath.name} is solving now by LOCAL")
                if 'PIXSCALE' in hdul[0].header:
                    PIXc = hdul[0].header['PIXSCALE']
                else : 
                    PIXc = _astro_utilities.calPixScale(hdul[0].header['FOCALLEN'], 
                                                    hdul[0].header['XPIXSZ'],
                                                    hdul[0].header['XBINNING'])
                print("PIXc : ", PIXc)

                solved = _astro_utilities.LOCALPSolver(fpath, 
                                                        #str(SOLVEDDIR), 
                                                        downsample = 2,
                                                        pixscale = PIXc,
>>>>>>> 83e7bdbe06de1d034dcbf1e69d4df772628a78a6
                                                                )