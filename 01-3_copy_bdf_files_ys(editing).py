# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

#first time
cd ~/Downloads/ && git clone https://github.com/ysBach/ysvisutilpy && cd ysvisutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/ysfitsutilpy && cd ysfitsutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/ysphotutilpy && cd ysphotutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/SNUO1Mpy && cd SNUO1Mpy && git pull && pip install -e . && cd ..

# second time...
cd ~/Downloads/ysvisutilpy && git pull && pip install -e . 
cd ~/Downloads/ysfitsutilpy && git pull && pip install -e . 
cd ~/Downloads/ysphouutilpy && git pull && pip install -e . 
cd ~/Downloads/SNUO1Mpy && git pull && pip install -e . 

"""
#%%
from glob import glob
from pathlib import Path
from datetime import datetime, timedelta, date
import pandas as pd
import os
import numpy as np
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

import Python_utilities
import astro_utilities

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
BASEDIR = "../RnE_2022/GSON300_STF-8300M/"

#%%
BASEDIRs = sorted(Python_utilities.getFullnameListOfsubDir(BASEDIR))
print ("BASEDIRs: {}".format(BASEDIRs))
print ("len(BASEDIRs): {}".format(len(BASEDIRs)))


for BASEDIR in BASEDIRs[:2] :
    print ("Starting...\n{}".format(BASEDIR))

    BASEDIR = Path(BASEDIR)
    
    MASTERDIR = BASEDIR / astro_utilities.master_dir
    REDUCEDDIR = BASEDIR / astro_utilities.reduced_dir
    SOLVEDDIR = BASEDIR / astro_utilities.solved_dir
    RESULTDIR = BASEDIR / astro_utilities.DAOfinder_result_dir
    OBSRAWDIR = BASEDIR / astro_utilities.CCD_obs_dir

    #%%
    summary = yfu.make_summary(BASEDIR/"*.fit*")
    if summary.empty:
        pass
    else:
        print(summary)
        print("len(summary):", len(summary))
        print(summary["file"][0])
        
        #%%
        light_fits = summary[summary["IMAGETYP"] == "LIGHT"]
        if light_fits.empty :
            pass
        else:
            print("len(light_fits):", len(light_fits))
            print("light_fits", light_fits)
            #%%
            light_fits["DATE-LOC-DT"] = pd.to_datetime(light_fits["DATE-LOC"])
            print("light_fits", light_fits)
            obs_date = light_fits["DATE-LOC-DT"].median().to_pydatetime().date()
            obs_date.strftime("%Y-%m-%d")
            print("BASEDIR.parts:", BASEDIR.parts)
            "{}_{}".format(BASEDIR.parts[-1].split("_")[5], BASEDIR.parts[-1].split("_")[6])
            
            calBDDIR = Path("../CCD_obs_raw/") /  \
                                "{}_{}".format(BASEDIR.parts[-1].split("_")[6], 
                                            BASEDIR.parts[-1].split("_")[8]) / "cal" 
            print("calBDDIR: ", calBDDIR)
            calBDDIRs = sorted(Python_utilities.getFullnameListOfallsubDirs(calBDDIR))
            print("calBDDIRs: ", calBDDIRs)

            #CCDOBSDIRs = sorted(Python_utilities.getFullnameListOfallsubDirs(CCD_obs_dir))


       
