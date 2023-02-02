# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

이 파일은 fits file의 header에 있는 plate solving 관련 keyword를 삭제해 준다.
필요할때만 사용하면 된다.
"""
#%%
import os, shutil
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.stats import sigma_clip
from astropy.io import fits
from ccdproc import combine, ccd_process, CCDData

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

from snuo1mpy import Preprocessor

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
#%%
BASEDIR = Path(astro_utilities.CCD_obs_raw_dir)

fpaths = sorted(list(BASEDIR.glob("*.csv")))
print(fpaths)

for summary_fpath in fpaths[:-1] :
    summary_all = pd.read_csv(str(summary_fpath))
    print(summary_all)

    for _, row in summary_all.iterrows():
        # 파일명 출력
        print (row["file"])
        fpath = Path(row["file"])
        new_fpath = Path(f"{fpath.parents[0]}/{fpath.stem}_clean.fit")
        # fits hedaer 에 있는 wcs 정보를 지운다
        yfu.wcsremove(fpath, 
                        additional_keys=["COMMENT"],
                        verbose=True,
                        output=new_fpath,
                        ccddata=True,
                        overwrite=True)
        if new_fpath.exists() \
                and fpath.exists():
                print("rename", f"{str(new_fpath)}", f"{str(fpath)}")
                #os.rename(f"{str(new_fpath)}", f"{str(fpath)}")
                shutil.move(f"{str(new_fpath)}", f"{str(fpath)}")
                
        # except Exception as err :
        #     print("X"*60)
        #     print('{0}'.format(err))
