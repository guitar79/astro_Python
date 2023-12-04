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
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData

import ysfitsutilpy as yfu
import ysphotutilpy as ypu

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
DOINGDIR = Path(BASEDIR / "ccd_test_folder")
DOINGDIR = Path(BASEDIR/ "asteroid/RiLA600_STX-16803_-_1bin")

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(DOINGDIR))
DOINGDIRs = sorted([x for x in DOINGDIR.iterdir() if x.is_dir()])
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################

DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/')
BASEDIR = Path("/mnt/Rdata/OBS_data") 
DOINGDIR = Path(BASEDIR/ "asteroid/RiLA600_STX-16803_-_1bin")

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
    fits_in_dir = sorted(list(DOINGDIR.glob('*.fit')))
    #print("fits_in_dir", fits_in_dir)
    print("len(fits_in_dir)", len(fits_in_dir))

    if len(fits_in_dir) == 0 :
        print(f"There is no fits fils in {DOINGDIR}")
        pass
    else : 
        print(f"Starting: {str(DOINGDIR.parts[-1])}")
    
        summary = yfu.make_summary(DOINGDIR/"*.fit")
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])

        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)
        print("df_light:\n{}".format(df_light))

        if df_light.empty:
            print("The dataframe(df_light) is empty")
            pass
        else:
            df_light = df_light.reset_index(drop=True)
            print("df_light:\n{}".format(df_light))

            for _, row  in df_light.iterrows():
                fpath = Path(row["file"])
                print("fpath.name ;", fpath.name)
                ccd = yfu.load_ccd(str(fpath))
                
                if 'PIXSCALE' in ccd.header:
                    PIXc = ccd.header['PIXSCALE']
                else : 
                    PIXc = _astro_utilities.calPixScale(ccd.header['FOCALLEN'], 
                                            ccd.header['XPIXSZ'],
                                            ccd.header['XBINNING'])
                print("PIXc : ", PIXc)

                solved = _astro_utilities.KevinSolver((str(fpath)), 
                                                        #str(SOLVEDDIR), 
                                                        downsample = 2,
                                                        pixscale = PIXc,
                                                        )