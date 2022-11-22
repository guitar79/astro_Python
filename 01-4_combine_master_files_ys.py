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
cd ~/Downloads/ysphotutilpy && git pull && pip install -e . 
cd ~/Downloads/SNUO1Mpy && git pull && pip install -e . 

"""
#%%
from glob import glob
from pathlib import Path
import os
import numpy as np
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

import astro_utilities
import Python_utilities

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
BASEDIR = astro_utilities.base_dir

#%%
BASEDIRs = sorted(Python_utilities.getFullnameListOfsubDir(BASEDIR))
print ("BASEDIRs: {}".format(BASEDIRs))
print ("len(BASEDIRs): {}".format(len(BASEDIRs)))
#%%
for BASEDIR in BASEDIRs[5:] :
    print ("Starting...\n{}".format(BASEDIR))

    BASEDIR = Path(BASEDIR)
    
    MASTERDIR = BASEDIR / astro_utilities.master_dir

    if not MASTERDIR.exists():
        os.makedirs("{}".format(str(MASTERDIR)))
        print("{} is created...".format(str(MASTERDIR)))

    #%%
    summary = yfu.make_summary(BASEDIR/"*.fit*")
    #print(summary)
    print("len(summary):", len(summary))
    print("summary:", summary)
    #print(summary["file"][0])

    #%%
    try: 
        bias_fits = summary[summary["IMAGETYP"] == "BIAS"]["file"]
        bias_comb = yfu.group_combine(
                        bias_fits.tolist(),
                        type_key = ["IMAGETYP"],
                        type_val = ["BIAS"],
                        group_key = ["EXPTIME"],
                        fmt = "master_bias.fits",  # output file name format
                        outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                        verbose = True
                    )
    except Exception as err :
        print("X"*60)
        print('{0}'.format(err))

    try: 
        dark_fits = summary[summary["IMAGETYP"] == "DARK"]["file"]
        # Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
        dark_comb = yfu.group_combine(
                        dark_fits.tolist(),
                        type_key = ["IMAGETYP"],
                        type_val = ["DARK"],
                        group_key = ["EXPTIME"],
                        fmt = "master_dark_{:.0f}sec.fits",  # output file name format
                        outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                        verbose=True
                    )
    except Exception as err :
        print("X"*60)
        print('{0}'.format(err))


    try:
        flat_fits = summary[summary["IMAGETYP"] == "FLAT"]["file"] 
        # Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
        flat_comb_norm = yfu.group_combine(
                        flat_fits.tolist(),
                        type_key = ["IMAGETYP"],
                        type_val = ["FLAT"],
                        group_key = ["FILTER"],
                        fmt = "master_flat_{:s}_norm.fits",  # output file name format
                        scale="med_sc", #norm
                        scale_to_0th=False, #norm
                        outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                        verbose=True
                    )

        # Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
        flat_comb_norm = yfu.group_combine(
                        flat_fits.tolist(),
                        type_key = ["IMAGETYP"],
                        type_val = ["FLAT"],
                        group_key = ["FILTER"],
                        fmt = "master_flat_{:s}.fits",  # output file name format
                        #scale="med_sc", #norm
                        #scale_to_0th=False, #norm
                        outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                        verbose=True
                    )
    except Exception as err :
        print("X"*60)
        print('{0}'.format(err))
    
# %%
