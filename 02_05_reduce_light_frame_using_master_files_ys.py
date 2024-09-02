# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

이 파일은 관측 자료를 전처리 해준다.

"""
#%%
from glob import glob
from pathlib import Path
import os
import shutil
import numpy as np
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData

import ysfitsutilpy as yfu
import ysfitsutilpy as yfu
import ysphotutilpy as ypu
#import ysvisutilpy as yvu

import _astro_utilities
import _Python_utilities

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
BASEDIR = Path("/mnt/Rdata/OBS_data") 
PROJECDIR = Path("/mnt/Rdata/OBS_data/2024-EXO")
TODODIR = PROJECDIR / "_-_-_2024-05_-_GSON300_STF-8300M_-_1bin"
TODODIR = PROJECDIR / "_-_-_2024-06_-_GSON300_STF-8300M_-_1bin"

# PROJECDIR = Path("/mnt/Rdata/OBS_data/2022-Asteroid")
# TODODIR = PROJECDIR / "GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "RiLA600_STX-16803_-_2bin"

# PROJECDIR = Path("/mnt/Rdata/OBS_data/2023-Asteroid")
# TODODIR = PROJECDIR / "GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "RiLA600_STX-16803_-_2bin"

# PROJECDIR = Path("/mnt/Rdata/OBS_data/2016-Variable")
# TODODIR = PROJECDIR / "-_-_-_2016-_-_RiLA600_STX-16803_-_2bin"

# PROJECDIR = Path("/mnt/Rdata/OBS_data/2017-Variable")
# TODODIR = PROJECDIR / "-_-_-_2017-_-_RiLA600_STX-16803_-_2bin"

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(TODODIR))
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

MASTERDIR = [x for x in DOINGDIRs if "CAL-BDF" in str(x)]
MASTERDIR = Path(MASTERDIR[0]) / _astro_utilities.master_dir
print ("MASTERDIR: ", format(MASTERDIR))

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_" in str(x)])
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

filter_str = '2024-08-30'
DOINGDIRs = [x for x in DOINGDIRs if filter_str in x]
# remove = 'BIAS'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'DARK'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'FLAT'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
print ("DOINGDIRs: ", DOINGDIRs)
print ("len(DOINGDIRs): ", len(DOINGDIRs))
#######################################################
#%%
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    print(f"Starting: {str(DOINGDIR.parts[-1])}")
    
    sMASTERDIR = DOINGDIR / _astro_utilities.master_dir
    REDUCEDDIR = DOINGDIR / _astro_utilities.reduced_dir

    if not sMASTERDIR.exists():
        os.makedirs(str(sMASTERDIR))
        print("{} is created...".format(str(sMASTERDIR)))

    if not REDUCEDDIR.exists():
        os.makedirs(str(REDUCEDDIR))
        print("{} is created...".format(str(REDUCEDDIR)))

    summary = yfu.make_summary(DOINGDIR/"*.fit*")
    if summary is not None :
        #print(summary)
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])

        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)

        for _, row in df_light.iterrows():

            fpath = Path(row["file"])
            ccd = yfu.load_ccd(fpath)
            filt = ccd.header["FILTER"]
            expt = ccd.header["EXPTIME"]
            try : 
                red = yfu.ccdred(
                    ccd,
                    output=Path(f"{REDUCEDDIR/ fpath.name}"),
                    mdarkpath=str(MASTERDIR / "master_dark_{:.0f}sec.fits".format(expt)),
                    mflatpath=str(MASTERDIR / "master_flat_{}_norm.fits".format(filt.upper())),
                    # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                    overwrite=True
                    )
            except FileNotFoundError: 
                _Python_utilities.write_log(err_log_file, "FileNotFoundError")