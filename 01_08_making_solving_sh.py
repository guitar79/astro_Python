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
import ysvisutilpy as yvu

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
BASEDIR = Path("/mnt/Rdata/CCD_obs") 
#BASEDIR = Path("/mnt/OBS_data") 
DOINGDIR = Path(BASEDIR/ "RnE_2022/GSON300_STF-8300M")
#DOINGDIR = Path(BASEDIR/ "CCD_new_files1")

#DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(DOINGDIR))
DOINGDIRs = sorted([x for x in DOINGDIR.iterdir() if x.is_dir()])
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#%%
for DOINGDIR in DOINGDIRs[:] :
    result = ''
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
        REDUCEDDIR = DOINGDIR / _astro_utilities.reduced_dir
        SOLVEDDIR = DOINGDIR / _astro_utilities.solved_dir2

        if not SOLVEDDIR.exists():
            os.makedirs("{}".format(str(SOLVEDDIR)))
            print("{} is created...".format(str(SOLVEDDIR))) 

        summary = yfu.make_summary(DOINGDIR/"*Light*.fit*")
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

            for _, row  in df_light.iterrows():

                result0 = _astro_utilities.makingAstrometrySH(row["file"], 
                                                            0.93,)
                result += result0

    with open(f'{DOINGDIR.stem}_solve.sh', 'w') as f:
        f.write(result)

