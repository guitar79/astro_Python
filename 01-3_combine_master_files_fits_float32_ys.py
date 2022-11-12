# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user


#first time
cd ~/Downloads/ && git clone https://github.com/ysBach/ysvisutilpy && cd ysvisutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/ysfitsutilpy && cd ysfitsutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/ysphotutilpy && cd ysphotutilpy && git pull && pip install -e . && cd ..

# second time...
cd ~/Downloads/ysvisutilpy && git pull && pip install -e . 
cd ~/Downloads/ysfitsutilpy && git pull && pip install -e . 
cd ~/Downloads/ysphouutilpy && git pull && pip install -e . 

"""
#%%
from glob import glob
import numpy as np
import os
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData
import Python_utilities
import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

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
c_method = 'median'
master_dir = "master_files/"

base_dir = "../Rne_2022/KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin/"
base_dir = "../Rne_2022/"
base_dir = "../Post_process/M13_Light_-_2021-04_-_TEC140x75_STL-11000M_-_1bin/"

#base_dir = "../CCD_obs_raw/"

base_dirs = Python_utilities.getFullnameListOfsubDir(base_dir)
print ("base_dirs1: {}".format(base_dirs))
base_dirs = [w for w in base_dirs if not (w.endswith(master_dir) \
        or w.endswith("fits"))]
print ("base_dirs2: {}".format(base_dirs))

#%%
for base_dir in base_dirs :
    print ("Starting...\n{}".format(base_dir))
    ######################################################
    #%%
    try : 
        summary = yfu.make_summary(
                    "{}/*.fit".format(base_dir),
                    keywords = ["DATE-OBS", "FILTER", "OBJECT"],  # header keywords; actually it is case-insensitive
                    #fname_option = 'name',  # 'file' column will contain only the name of the file (not full path)
                    sort_by = "DATE-OBS",  # 'file' column will be sorted based on "DATE-OBS" value in the header
                    output = "{}summary.csv".format(base_dir)
                )
        print("summary:\n {}".format(summary))
    
    except Exception as err :
        print("X"*60)
        print('{0}'.format(err))


    #%%
    try: 
        bias_comb = yfu.group_combine(
                        "{}/*Bias*.fit".format(base_dir),
                        type_key = ["IMAGETYP"],
                        type_val = ["BIAS"],
                        group_key = ["EXPTIME"],
                        fmt = "master_bias.fits",  # output file name format
                        outdir = "{}{}".format(base_dir, master_dir)  # output directory (will automatically be made if not exist)
                    )
    except Exception as err :
        print("X"*60)
        print('{0}'.format(err))

    #%%
    try:
        bias_comb = yfu.group_combine(
                        "{}/*Bias*.fit".format(base_dir),
                        type_key = ["IMAGETYP"],
                        type_val = ["bias"],
                        group_key = ["EXPTIME"],
                        fmt = "master_bias.fits",  # output file name format
                        outdir = "{}{}".format(base_dir, master_dir)  # output directory (will automatically be made if not exist)
                    )
    except Exception as err :
        print("X"*60)
        print('{0}'.format(err))


    #%%
    try: 
        # Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
        dark_comb = yfu.group_combine(
                        "{}/*Dark*.fit".format(base_dir),
                        type_key = ["IMAGETYP"],
                        type_val = ["DARK"],
                        group_key = ["EXPTIME"],
                        fmt = "master_dark_{:.1f}sec.fits",  # output file name format
                        outdir = "{}{}".format(base_dir, master_dir)  # output directory (will automatically be made if not exist)
                    )
    except Exception as err :
        print("X"*60)
        print('{0}'.format(err))


    #%%
    try: 
        # Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
        dark_comb = yfu.group_combine(
                        "{}/*Dark*.fit".format(base_dir),
                        type_key = ["IMAGETYP"],
                        type_val = ["dark"],
                        group_key = ["EXPTIME"],
                        fmt = "master_dark_{:.1f}sec.fits",  # output file name format
                        outdir = "{}{}".format(base_dir, master_dir)  # output directory (will automatically be made if not exist)
                    )
    except Exception as err :
        print("X"*60)
        print('{0}'.format(err))


    #%%
    try: 
        # Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
        flat_comb = yfu.group_combine(
                        "{}/*Flat*.fit".format(base_dir),
                        type_key = ["IMAGETYP"],
                        type_val = ["FLAT"],
                        group_key = ["FILTER"],
                        fmt = "master_flat_{:s}.fits",  # output file name format
                        outdir = "{}{}".format(base_dir, master_dir)  # output directory (will automatically be made if not exist)
                    )
    except Exception as err :
        print("X"*60)
        print('{0}'.format(err))

