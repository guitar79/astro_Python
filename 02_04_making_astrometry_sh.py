# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from glob import glob
from pathlib import Path
from datetime import datetime
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
BASEDIR = Path("/mnt/Rdata/OBS_data") 
PROJECDIR = Path("/mnt/Rdata/OBS_data/2024-EXO")
TODODIR = PROJECDIR / "_-_-_2024-05_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "_-_-_2024-06_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "RiLA600_STX-16803_-_2bin"

PROJECDIR = Path("/mnt/Rdata/OBS_data/2022-Asteroid")
TODODIR = PROJECDIR / "GSON300_STF-8300M_-_1bin"
TODODIR = PROJECDIR / "RiLA600_STX-16803_-_1bin"
TODODIR = PROJECDIR / "RiLA600_STX-16803_-_2bin"

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

CALDIR = [x for x in DOINGDIRs if "CAL-BDF" in str(x)]
MASTERDIR = Path(CALDIR[0]) / _astro_utilities.master_dir
if not MASTERDIR.exists():
    os.makedirs("{}".format(str(MASTERDIR)))
    print("{} is created...".format(str(MASTERDIR)))

print ("MASTERDIR: ", format(MASTERDIR))

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_" in str(x)])
# print ("DOINGDIRs: ", format(DOINGDIRs))
# print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

# filter_str = '2023-12-'
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
solve_sh = ""
# rename_sh = ""
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    SOLVINGDIR = DOINGDIR / _astro_utilities.reduced_dir
    SOLVINGDIR = DOINGDIR / _astro_utilities.reduced_nightsky_dir
    # SOLVINGDIR = DOINGDI
    # R
    
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

        # solve_sh = ""
        # rename_sh = ""
        for _, row  in df_light.iterrows():

            fpath = Path(row["file"])
            hdul = fits.open(fpath)
            
            if 'PIXSCALE' in hdul[0].header:
                pixscale = hdul[0].header['PIXSCALE']
            else : 
                pixscale = _astro_utilities.calPixScale(hdul[0].header['FOCALLEN'], 
                                            hdul[0].header['XPIXSZ'],
                                            hdul[0].header['XBINNING'])
            hdul.close()
            print(f"pixscale: {pixscale:.03f}, L: {pixscale*0.97:.03f}, U: {pixscale*1.03:.03f}")
            # fpath = Path(df_light["file"][1])
            print("fpath :" ,fpath)
            SOLVE, ASTAP, LOCAL = _astro_utilities.checkPSolve(fpath)
            print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)
            if not SOLVE : 
                #solve-field -O -g --cpulimit 15 --nsigma 15 --downsample 4 -u app -L 0.6 -U 0.63 --no-plots
                solve_sh += f"solve-field -O -g --cpulimit 15 --nsigma 15 --downsample 4 -u app  -L 0.6 --no-plots {str(fpath)}\n"
                solve_sh += f"mv {fpath.parent/fpath.stem}.new {str(fpath)}\n"
                solve_sh += f"rm {fpath.parent/fpath.stem}-indx.xyls\n"
                solve_sh += f"rm {fpath.parent/fpath.stem}.rdls\n"
                solve_sh += f"rm {fpath.parent/fpath.stem}.corr\n"
                solve_sh += f"rm {fpath.parent/fpath.stem}.solved\n"
                solve_sh += f"rm {fpath.parent/fpath.stem}.match\n"
                solve_sh += f"rm {fpath.parent/fpath.stem}.axy\n"
               
print("solve_sh:", solve_sh)
    
with open(f"_{datetime.now().strftime('%Y%m%d-%H%m%S')}_{TODODIR.stem}_astrometry_solve.sh", 'w') as f:
    f.write(solve_sh)        
                
        # print("solve_sh:", solve_sh)
            
        # with open(f"_{datetime.now().strftime('%Y%m%d-%H%m%S')}_{DOINGDIR.stem}_astrometry_solve.sh", 'w') as f:
        #     f.write(solve_sh)
        # with open(f"{DOINGDIR.stem}_astrometry_rename.sh", 'w') as f:
        #     f.write(rename_sh)
