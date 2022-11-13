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
import numpy as np
import os
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData
from astropy.io import fits
import Python_utilities
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
base_dir = "../RnE_2022/KLEOPATRA_Light_-_2022-11-08_-_RiLA600_STX-16803_-_2bin/"
base_dir = "../RnE_2022/"

master_dir = "master_files/"

#%%
base_dirs = sorted(Python_utilities.getFullnameListOfsubDir(base_dir))
base_dirs = [w for w in base_dirs if not (w.endswith(master_dir) \
                or w.endswith(".fits"))]
print ("base_dirs: {}".format(base_dirs))


#%%
for base_dir in base_dirs :
    print ("Starting...\n{}".format(base_dir))

    try : 
        summary = yfu.make_summary("{}/*.fit".format(base_dir))
                            #keywords = ["DATE-OBS", "FILTER", "OBJECT", "IMAGETYP"],  # header keywords; actually it is case-insensitive
                            #fname_option = 'name',  # 'file' column will contain only the name of the file (not full path)
                            #output = "{}summary.csv".format(base_dir),
                            #sort_by = "DATE-OBS"  # 'file' column will be sorted based on "DATE-OBS" value in the header
                            #)
        print("summary:\n {}".format(summary))

    except Exception as err :
        print("X"*60)
        print('{0}'.format(err))

    #%%
    try:
        df_light = summary[summary["IMAGETYP"] == "LIGHT"]
        print ("df_light: {}".format(df_light))
        print ("len(df_light): {}".format(len(df_light)))

        #%%
        for _, row in df_light.iterrows():
            print (row["file"])
            yfu.wcsremove(row["file"], 
                        output=row["file"], 
                        overwrite=True)
    except Exception as err :
        print("X"*60)
        print('{0}'.format(err))
# %%
