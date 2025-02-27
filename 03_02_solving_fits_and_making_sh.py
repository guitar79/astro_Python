# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from glob import glob
from pathlib import Path
from datetime import datetime, timedelta
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData

import ysfitsutilpy as yfu

import _astro_utilities
import _Python_utilities

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
BASEDIR = Path("/mnt/Rdata/ASTRO_data") 
PROJECDIR = Path("/mnt/Rdata/ASTRO_data/2024-EXO")
TODODIR = PROJECDIR / "_-_-_2024-05_-_GSON300_STF-8300M_-_1bin"
TODODIR = PROJECDIR / "_-_-_2024-06_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "RiLA600_STX-16803_-_2bin"

PROJECDIR = Path("/mnt/Rdata/ASTRO_data/2022-Asteroid")
TODODIR = PROJECDIR / "GSON300_STF-8300M_-_1bin"
TODODIR = PROJECDIR / "RiLA600_STX-16803_-_1bin"
TODODIR = PROJECDIR / "RiLA600_STX-16803_-_2bin"

PROJECDIR = Path("/mnt/Rdata/ASTRO_data/2023-Asteroid")
TODODIR = PROJECDIR / "GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "RiLA600_STX-16803_-_2bin"

# PROJECDIR = Path("/mnt/Rdata/ASTRO_data/2016-Variable")
# TODODIR = PROJECDIR / "-_-_-_2016-_-_RiLA600_STX-16803_-_2bin"

# PROJECDIR = Path("/mnt/Rdata/ASTRO_data/2017-Variable")
# TODODIR = PROJECDIR / "-_-_-_2017-_-_RiLA600_STX-16803_-_2bin"

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(TODODIR))
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

BDFDIR = [x for x in DOINGDIRs if "CAL-BDF" in str(x)]
MASTERDIR = Path(BDFDIR[0]) / _astro_utilities.master_dir
if not MASTERDIR.exists():
    os.makedirs("{}".format(str(MASTERDIR)))
    print("{} is created...".format(str(MASTERDIR)))

print ("MASTERDIR: ", format(MASTERDIR))

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_" in str(x)])
# print ("DOINGDIRs: ", format(DOINGDIRs))
# print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

# filter_str = '2024-08-29'
# DOINGDIRs = [x for x in DOINGDIRs if filter_str in x]
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
for DOINGDIR in DOINGDIRs :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)

    
    SOLVINGDIR = DOINGDIR / _astro_utilities.reduced_dir
    # SOLVINGDIR = DOINGDIR / _astro_utilities.reduced_nightsky_dir
    # SOLVINGDIR = DOINGDIR

    summary = yfu.make_summary(SOLVINGDIR/"*.fit*",
                                    verify_fix=True,
                                    ignore_missing_simple=True,
                                    )
    if summary is not None :
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])  
        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)
        print("df_light:\n{}".format(df_light))
    df_light

    for _, row  in df_light.iterrows():

        fpath = Path(row["file"])
        # fpath = Path(df_light["file"][1])
        print("fpath :" ,fpath)
        hdul = fits.open(fpath)

        if 'PIXSCALE' in hdul[0].header:
            PIXc = hdul[0].header['PIXSCALE']
        else : 
            PIXc = _astro_utilities.calPixScale(hdul[0].header['FOCALLEN'], 
                                                hdul[0].header['XPIXSZ'],
                                                hdul[0].header['XBINNING'])
        print("PIXc : ", PIXc)
        hdul.close()

        try : 
            solved = _astro_utilities.KevinSolver(fpath, 
                                                #str(SOLVEDDIR), 
                                                # downsample = 2,
                                                # pixscale = PIXc,
                                                tryASTAP = True, 
                                                tryLOCAL = False,
                                                makeLOCALsh = True,
                                                tryASTROMETRYNET = False, 
                                                )
        except Exception as err :
            print("X"*60)
            _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
        
# os.popen(f"sh __{datetime.now().strftime('%Y%m%d')}_todo_astrometry_solve.sh")
