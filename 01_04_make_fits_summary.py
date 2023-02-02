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
import ysphotutilpy as ypu
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
BASEDIR = astro_utilities.base_dir
BASEDIR = Path(astro_utilities.CCD_obs_raw_dir) 

ccd_dirs = ["QSI683ws_1bin", "QSI683ws_2bin", "QSI683ws_3bin", 
            "ST-8300M_1bin", "ST-8300M_2bin", 
            "STF-8300M_1bin", "STF-8300M_2bin", 
            "STL-11000M_1bin", "STL-11000M_2bin",
            "STX-16803_1bin", "STX-16803_2bin" ]

for ccd_dir in ccd_dirs :
    #BASEDIR = BASEDIR / ccd_dirs[0]

    #BASEDIRs = sorted(Python_utilities.getFullnameListOfsubDir(BASEDIR))
    BASEDIRs = sorted(Python_utilities.getFullnameListOfallsubDirs(BASEDIR/ccd_dir))
    #print ("BASEDIRs: {}".format(BASEDIRs))
    print ("len(BASEDIRs): {}".format(len(BASEDIRs)))

    #%%
    ccd_fpath = Path(f"{BASEDIR/ccd_dir}.csv")
    if ccd_fpath.exists():
        print(f"{BASEDIR/ccd_dir}.csv is already exist...")
    else : 
        summary_all = pd.DataFrame()
        for fpath in BASEDIRs[:] :
            #BASEDIR = Path(BASEDIRs[0])
            fpath = Path(fpath)
            save_fpath2 = fpath/f"{fpath.parts[-1]}.csv"
            save_fpath = fpath/f"summary_{fpath.parts[-1]}.csv"
            print (f"Starting...\n{fpath.name}")
            if save_fpath2.exists():
                os.remove(str(save_fpath2))
                print (f"{str(save_fpath2)} is deleted...")
            
            if save_fpath.exists():
                print(f"{str(save_fpath)} is already exist...")
                summary = pd.read_csv(str(save_fpath))
            
            else : 
                summary = yfu.make_summary(fpath/"*.fit*",
                            output = save_fpath,
                            verbose = False
                            )
                print(f"{save_fpath} is created...")
            summary_all = pd.concat([summary_all, summary], axis = 0)
        
        summary_all.reset_index(inplace=True)
        summary_all.to_csv(f"{BASEDIR/ccd_dir}.csv")
        print(f"{BASEDIR/ccd_dir}.csv is created...")
