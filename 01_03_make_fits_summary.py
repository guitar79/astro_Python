# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

"""
#%%
from glob import glob
from pathlib import Path
import shutil
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData

import ysfitsutilpy as yfu
#import ysphotutilpy as ypu
import ysvisutilpy as yvu

import astro_utilities
import Python_utilities

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
BASEDIR = Path("/mnt/Rdata/CCD_obs") 
#BASEDIR = Path("/mnt/OBS_data") 
DOINGDIR = BASEDIR/ astro_utilities.CCD_obs_raw_dir

DOINGDIRs = sorted(Python_utilities.getFullnameListOfsubDirs(DOINGDIR))
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    DOINGSUBDIRs = sorted(Python_utilities.getFullnameListOfallsubDirs(str(DOINGDIR)))
    ccd_fpath = Path(f"{DOINGDIR/DOINGDIR.parts[-1]}.csv")
    print("ccd_fpath", ccd_fpath)

    if ccd_fpath.exists() or False:
        print(f"{str(ccd_fpath)} is already exist...")

    else : 
        summary_all = pd.DataFrame()

        for DOINGSUBDIR in DOINGSUBDIRs[:] :
            print (f"Starting...\n{DOINGSUBDIR}")
            DOINGSUBDIR = Path(DOINGSUBDIR)
            #print("DOINGSUBDIR", DOINGSUBDIR)
            fits_in_dir = sorted(list(DOINGSUBDIR.glob('*.fit*')))
            #print("fits_in_dir", fits_in_dir)
            print("len(fits_in_dir)", len(fits_in_dir))
            save_fpath2 = DOINGSUBDIR / f"{DOINGSUBDIR.parts[-1]}.csv"
            if save_fpath2.exists():
                os.remove(str(save_fpath2))
                print (f"{str(save_fpath2)} is deleted...")
        
            if len(fits_in_dir) == 0 :
                print(f"There is no fits fils in {DOINGSUBDIR}")
                pass
            else : 
                save_fpath = DOINGSUBDIR / f"summary_{DOINGSUBDIR.parts[-1]}.csv"
                if save_fpath.exists():
                    print(f"{str(save_fpath)} is already exist...")
                    summary = pd.read_csv(str(save_fpath))
                
                else : 
                    summary = yfu.make_summary(DOINGSUBDIR/"*.fit*",
                                output = save_fpath,
                                verbose = False
                                )
                    print(f"{save_fpath} is created...")
                summary_all = pd.concat([summary_all, summary], axis = 0)
        
        summary_all.reset_index(inplace=True)
        summary_all.to_csv(str(ccd_fpath))
        print(f"{ccd_fpath} is created...")