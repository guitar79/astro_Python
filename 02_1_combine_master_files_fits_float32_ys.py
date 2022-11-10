# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

git clone https://github.com/ysBach/ysvisutilpy && cd ysvisutilpy && git pull && pip install -e . && cd ..

git clone https://github.com/ysBach/ysfitsutilpy && cd ysfitsutilpy && git pull && pip install -e . && cd ..

git clone https://github.com/ysBach/ysphotutilpy && cd ysphotutilpy && git pull && pip install -e . && cd ..

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

#base_dir = "../Post_processing/M35_Light_-_2018-10-31_-_TMB130ss_STF-8300M_-_1bin/"
base_dir = "../RnE_2022/KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin/"

#base_dir = "../RnE_2022/"

#base_dirs = Python_utilities.getFullnameListOfsubDir(base_dir)
#print ("base_dirs: {}".format(base_dirs))


#%%
c_method = 'median'
master_dir = "master_files/"

summary = yfu.make_summary(
            "{}/*.fit".format(base_dir),
            keywords = ["DATE-OBS", "FILTER", "OBJECT"],  # header keywords; actually it is case-insensitive
            #fname_option = 'name',  # 'file' column will contain only the name of the file (not full path)
            sort_by = "DATE-OBS",  # 'file' column will be sorted based on "DATE-OBS" value in the header
            output = "{}summary.csv".format(base_dir)
        )

summary


#%%
bias_comb = yfu.group_combine(
                "{}/*Bias*.fit".format(base_dir),
                type_key = ["IMAGETYP"],
                type_val = ["BIAS"],
                group_key = ["EXPTIME"],
                fmt = "master_bias.fits",  # output file name format
                outdir = "{}{}".format(base_dir, master_dir)  # output directory (will automatically be made if not exist)
            )

#%%
bias_comb = yfu.group_combine(
                "{}/*Bias*.fit".format(base_dir),
                type_key = ["IMAGETYP"],
                type_val = ["bias"],
                group_key = ["EXPTIME"],
                fmt = "master_bias.fits",  # output file name format
                outdir = "{}{}".format(base_dir, master_dir)  # output directory (will automatically be made if not exist)
            )


#%%
# Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
dark_comb = yfu.group_combine(
                "{}/*Dark*.fit".format(base_dir),
                type_key = ["IMAGETYP"],
                type_val = ["DARK"],
                group_key = ["EXPTIME"],
                fmt = "master_dark_{:.1f}sec.fits",  # output file name format
                outdir = "{}{}".format(base_dir, master_dir)  # output directory (will automatically be made if not exist)
            )



#%%
# Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
dark_comb = yfu.group_combine(
                "{}/*Dark*.fit".format(base_dir),
                type_key = ["IMAGETYP"],
                type_val = ["dark"],
                group_key = ["EXPTIME"],
                fmt = "master_dark_{:.1f}sec.fits",  # output file name format
                outdir = "{}{}".format(base_dir, master_dir)  # output directory (will automatically be made if not exist)
            )


#%%
# Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
flat_comb = yfu.group_combine(
                "{}/*Flat*.fit".format(base_dir),
                type_key = ["IMAGETYP"],
                type_val = ["FLAT"],
                group_key = ["FILTER"],
                fmt = "master_flat_{:s}.fits",  # output file name format
                outdir = "{}{}".format(base_dir, master_dir)  # output directory (will automatically be made if not exist)
            )

