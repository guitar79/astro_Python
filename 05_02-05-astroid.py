# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

이 파일은 관측 자료를 전처리 해준다.

"""
#%%
from glob import glob
from pathlib import Path
import os
import shutil
import numpy as np
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData

import ysfitsutilpy as yfu
import ysfitsutilpy as ypu

import _astro_utilities
import _asteroid_utilities
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
#%%
#######################################################
# read all files in base directory for processing
BASEDIR = Path("/mnt/Rdata/OBS_data") 
DOINGDIR = Path(BASEDIR/ "asteroid" / "RiLA600_STX-16803_-_1bin")
#DOINGDIR = Path(BASEDIR/ "asteroid" / "GSON300_STF-8300M_-_1bin")

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(DOINGDIR))

filter_str = '2023-11-2'
DOINGDIRs = [x for x in DOINGDIRs if filter_str in x]
# remove = 'BIAS'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'DARK'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'FLAT'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
print ("DOINGDIRs: ", DOINGDIRs)
print ("len(DOINGDIRs): ", 
       len(DOINGDIRs))
#######################################################

#%%
#mas1 = Path(DOINGDIRs[0])
#mas1 = mas1 /_astro_utilities.master_dir
# if str(DOINGDIR.parts[-2]) == "RiLA600_STX-16803_-_1bin" :
#     DOINGDIR = DOINGDIR / _astro_utilities.REDUC_nightsky_dir
# if str(DOINGDIR.parts[-2]) == "GSON300_STF-8300M_-_1bin" :
#     DOINGDIR = DOINGDIR / _astro_utilities.reduced_dir
mas1 = BASEDIR / "asteroid" / "GSON300_STF-8300M_-_1bin" / "master_files_ys"
mas1 = DOINGDIR / "master_files_ys"
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)
    fits_in_dir = sorted(list(DOINGDIR.glob('*.fit*')))
    #print("fits_in_dir", fits_in_dir)
    print("len(fits_in_dir)", len(fits_in_dir))

    if len(fits_in_dir) == 0 :
        print(f"There is no fits fils in {DOINGDIR}")
        pass
    else : 
        print(f"Starting: {str(DOINGDIR.parts[-1])}")

        summary = yfu.make_summary(DOINGDIR/"*.fit*")
        #print(summary)
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])

        MASTERDIR = DOINGDIR / _astro_utilities.master_dir
        REDUCEDDIR = DOINGDIR / _astro_utilities.reduced_dir
        REDUC_nightsky = DOINGDIR / _astro_utilities.REDUC_nightsky_dir
        ASTRESULTDIR = DOINGDIR / _astro_utilities.Asteroid_result_dir

        if not MASTERDIR.exists():
            shutil.copytree(mas1, MASTERDIR)

        if not REDUCEDDIR.exists():
            os.makedirs(str(REDUCEDDIR))
            print("{} is created...".format(str(REDUCEDDIR)))

        if not REDUC_nightsky.exists():
            os.makedirs("{}".format(str(REDUC_nightsky)))
            print("{} is created...".format(str(REDUC_nightsky)))
    
        if not ASTRESULTDIR.exists():
            os.makedirs("{}".format(str(ASTRESULTDIR)))
            print("{} is created...".format(str(ASTRESULTDIR)))

        if str(DOINGDIR.parts[-2]) == "RiLA600_STX-16803_-_1bin" :
            _asteroid_utilities.makeNightskyflatReduceLightFrame(DOINGDIR, summary)
        if str(DOINGDIR.parts[-2]) == "GSON300_STF-8300M_-_1bin" :
            _asteroid_utilities.reduceLightFrame(DOINGDIR, summary)
        #_asteroid_utilities.solvingLightFrame(DOINGDIR, summary)
        #_asteroid_utilities.checkAsteroids(DOINGDIR, summary)

        