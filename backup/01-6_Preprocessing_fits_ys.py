# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user


#first time
cd ~/Downloads/ && git clone https://github.com/ysBach/ysvisutilpy && cd ysvisutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/ysfitsutilpy && cd ysfitsutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/ysphotutilpy && cd ysphotutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/SNUO1Mpy && cdSNUO1Mpy && git pull && pip install -e . && cd ..

# second time...
cd ~/Downloads/ysvisutilpy && git pull && pip install -e . 
cd ~/Downloads/ysfitsutilpy && git pull && pip install -e . 
cd ~/Downloads/ysphouutilpy && git pull && pip install -e . 
cd ~/Downloads/SNUO1Mpy && git pull && pip install -e . 

"""
#%%
from glob import glob
import numpy as np
import os
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData
import _Python_utilities
import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu


from pathlib import Path
from snuo1mpy import Preprocessor
import numpy as np

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
BASEDIR = "../Post_process/M13_Light_-_2021-04_-_TEC140x75_STL-11000M_-_1bin/"
BASEDIR = "../RnE_2022/KLEOPATRA_Light_-_2022-11-08_-_RiLA600_STX-16803_-_2bin_work/"
#BASEDIR = "../RnE_2022/"
#BASEDIR = "../CCD_obs_raw/"

master_dir = "master_files/"

#%%
#######################################################
# At the current version of SNUO1Mpy, the following bias_kw and dark_kw are 
# identical to the defaults. I just explicitly wrote them for clarity.
bias_kw = dict(bias_type_key=["IMAGETYP"], bias_type_val=["BIAS"], bias_group_key=[])
dark_kw = dict(dark_type_key=["IMAGETYP"], dark_type_val=["DARK"], dark_group_key=["EXPTIME"])
flat_kw = dict(flat_type_key=["IMAGETYP"], flat_type_val=["FLAT"], flat_group_key=["FILTER"])

#%%
TOPDIR = Path("{}".format(BASEDIR))
RAWDIR = TOPDIR
ARCHIVE = TOPDIR/"archive"
CALIBDIR = TOPDIR/"calib"

p = Preprocessor(topdir=TOPDIR, rawdir=RAWDIR, 
                instrument = "STX-16803",
                #overwrite=True,
                **bias_kw, **dark_kw, **flat_kw) 

#%%
p.organize_raw(archive_dir=ARCHIVE, 
                verbose=True)

# %%
