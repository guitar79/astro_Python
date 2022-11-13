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
import numpy as np
import os
import astropy.units as u
from ccdproc import CCDData, ccd_process
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
base_dir = "../RnE_2022/"

c_method = 'median'
master_dir = "master_files_ys"
reduced_dir = "reduced"
solved_dir = "solved"

#%%
base_dirs = sorted(Python_utilities.getFullnameListOfsubDir(base_dir))
base_dirs = [w for w in base_dirs \
        if not (w.endswith("{}/".format(master_dir)) \
            or w.endswith("{}/".format(reduced_dir)) \
            or w.endswith("{}/".format(solved_dir)) \
                or w.endswith(".fits"))]

print ("base_dirs: {}".format(base_dirs))

#%%
for base_dir in base_dirs :
    print ("Starting...\n{}".format(base_dir))

    base_dir = Path(base_dir)

    if not (base_dir/reduced_dir).exists():
        os.makedirs(str((base_dir/reduced_dir)))

    #%%
    summary = yfu.make_summary(base_dir/"*.fit")

    df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
    df_light = df_light.reset_index(drop=True)

    # %%
    for _, row in df_light.iterrows():
        try:
            fpath = Path(row["file"])
            ccd = yfu.load_ccd(fpath)
            filt = ccd.header["FILTER"]
            expt = ccd.header["EXPTIME"]
            red = yfu.ccdred(
                ccd,
                output=Path(f"{base_dir}/reduced/{fpath.stem}.fits"),
                mdarkpath=str(base_dir/"master_files_ys/master_dark_{:.0f}sec.fits".format(expt)),
                mflatpath=str(base_dir/"master_files_ys/master_flat_{}_norm.fits".format(filt.upper())),
                # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                overwrite=True
            )
        except Exception as err: 
            print ('Error messgae .......')
            print (err)
   
# %%
