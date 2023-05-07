# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

이 파일은 관측 자료를 전처리 해준다.

"""
#%%
from glob import glob
from pathlib import Path
import numpy as np
import os
import astropy.units as u
from ccdproc import CCDData, ccd_process

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

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
# read all files in base directory for processing
BASEDIR = Path(r"r:\CCD_obs")
BASEDIR = Path("/mnt/Rdata/CCD_obs") 
#BASEDIR = Path("/mnt/OBS_data") 
DOINGDIR = Path(BASEDIR/ "RnE_2022/GSON300_STF-8300M")
DOINGDIR = Path(BASEDIR/ "RnE_2022/RiLA600_STX-16803_1bin")
#DOINGDIR = Path(BASEDIR/ "CCD_new_files1")

#DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(DOINGDIR))
DOINGDIRs = sorted([x for x in DOINGDIR.iterdir() if x.is_dir()])
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################

#%%
for DOINGDIR in DOINGDIRs[:] :
    #DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)
    fits_in_dir = sorted(list(DOINGDIR.glob('*.fit*')))
    #print("fits_in_dir", fits_in_dir)
    print("len(fits_in_dir)", len(fits_in_dir))

    if len(fits_in_dir) == 0 :
        print(f"There is no fits fils in {DOINGDIR}")
        pass
    else : 
        print(f"Starting: {str(DOINGDIR.parts[-1])}")

        MASTERDIR = DOINGDIR / _astro_utilities.master_dir
        REDUCEDDIR = DOINGDIR / _astro_utilities.reduced_dir

        if not REDUCEDDIR.exists():
            os.makedirs(str(REDUCEDDIR))
            print("{} is created...".format(str(REDUCEDDIR)))

        summary = yfu.make_summary(DOINGDIR/"*.fit*")
        #print(summary)
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])

        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)

        for _, row in df_light.iterrows():
            try:
                fpath = Path(row["file"])
                ccd = yfu.load_ccd(fpath)
                filt = ccd.header["FILTER"]
                expt = ccd.header["EXPTIME"]
                red = yfu.ccdred(
                    ccd,
                    output=Path(f"{REDUCEDDIR}/{fpath.stem}.fit"),
                    mdarkpath=str(MASTERDIR / "master_dark_{:.0f}sec.fits".format(expt)),
                    mflatpath=str(MASTERDIR / "master_flat_{}_norm.fits".format(filt.upper())),
                    # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                    overwrite=True
                )
            except FileNotFoundError: 
                _Python_utilities.write_log(err_log_file, "FileNotFoundError")
 