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
import shutil
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

import warnings
warnings.filterwarnings('ignore')

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

PROJECDIR = BASEDIR / "C1-Variable"
TODODIR = PROJECDIR / "-_-_-_2016-_-_RiLA600_STX-16803_-_2bin"
TODODIR = PROJECDIR / "-_-_-_2017-01_-_RiLA600_STX-16803_-_2bin"
TODODIR = PROJECDIR / "-_-_-_2017-03_-_RiLA600_STX-16803_-_2bin"
TODODIR = PROJECDIR / "-_-_-_2017-05_-_RiLA600_STX-16803_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2017-06_-_RiLA600_STX-16803_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2021-10_-_RiLA600_STX-16803_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2022-01_-_RiLA600_STX-16803_-_2bin"

PROJECDIR = BASEDIR / "C2-Asteroid"
TODODIR = PROJECDIR / "-_-_-_2022-_-_GSON300_STF-8300M_-_1bin"
TODODIR = PROJECDIR / "-_-_-_2022-_-_RiLA600_STX-16803_-_1bin"
TODODIR = PROJECDIR / "-_-_-_2022-_-_RiLA600_STX-16803_-_2bin"
TODODIR = PROJECDIR / "-_-_-_2023-_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2023-_-_RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2023-_-_RiLA600_STX-16803_-_2bin"

# PROJECDIR = BASEDIR / "C3-EXO"
# TODODIR = PROJECDIR / "-_-_-_2024-05_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2024-05_-_RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2024-06_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2024-06_-_RiLA600_STX-16803_-_2bin"

# PROJECDIR = BASEDIR / "C4-Spectra"
# TODODIR = PROJECDIR / "-_-_-_2024-05_TEC140_ASI183MMPro_-_1bin"

PROJECDIR = BASEDIR / "C5-Test"
TODODIR = PROJECDIR / "-_-_-_-_GSON300_STF-8300M_-_1bin"

PROJECDIR = BASEDIR / "A3_CCD_obs_raw/STF-8300M_1bin/"
TODODIR = PROJECDIR / "LIGHT_FS60CB/"
TODODIR = PROJECDIR / "LIGHT_FSQ106ED-x73"
TODODIR = PROJECDIR / "LIGHT_GSON300"
TODODIR = PROJECDIR / "LIGHT_OON300"
TODODIR = PROJECDIR / "LIGHT_RILA600"
TODODIR = PROJECDIR / "LIGHT_SVX80T-x80"
TODODIR = PROJECDIR / "LIGHT_TEC140"
TODODIR = PROJECDIR / "LIGHT_TMB130ss"
TODODIR = PROJECDIR / "LIGHT_TMB130ss-x75"

PROJECDIR = BASEDIR / "A3_CCD_obs_raw/QSI683ws_1bin/"
TODODIR = PROJECDIR / "LIGHT_FSQ106ED-x73"
TODODIR = PROJECDIR / "LIGHT_GSON300"
TODODIR = PROJECDIR / "LIGHT_RILA600"
TODODIR = PROJECDIR / "LIGHT_SVX80T"
TODODIR = PROJECDIR / "LIGHT_SVX80T-x80"
TODODIR = PROJECDIR / "LIGHT_TMB130ss-x75"

PROJECDIR = BASEDIR / "A3_CCD_obs_raw/STL-11000M_1bin/"
TODODIR = PROJECDIR / "LIGHT_FSQ106ED"
TODODIR = PROJECDIR / "LIGHT_FSQ106ED-x72"
TODODIR = PROJECDIR / "LIGHT_TEC140-x75"
TODODIR = PROJECDIR / "LIGHT_TMB130ss"
TODODIR = PROJECDIR / "LIGHT_TMB130ss-x75"

PROJECDIR = BASEDIR / "A3_CCD_obs_raw/STX-16803_1bin/"
TODODIR = PROJECDIR / "LIGHT_RILA600"

PROJECDIR = BASEDIR / "A3_CCD_obs_raw/STX-16803_2bin/"
TODODIR = PROJECDIR / "LIGHT_RILA600"
               
DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(TODODIR))
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

try : 
    BDFDIR = [x for x in DOINGDIRs if "CAL-BDF" in str(x)]
    print ("BDFDIR: ", format(BDFDIR))
    MASTERDIR = Path(BDFDIR[0]) / _astro_utilities.master_dir
    if not MASTERDIR.exists():
        os.makedirs("{}".format(str(MASTERDIR)))
        print("{} is created...".format(str(MASTERDIR)))
    print ("MASTERDIR: ", format(MASTERDIR))
except : 
    pass

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_" in str(x)])
# print ("DOINGDIRs: ", format(DOINGDIRs))
# print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

# filter_str = '127JOHANNA_LIGHT_-_2023-11-17_-_GSON300_STF-8300M_-_1bin'
# DOINGDIRs = [x for x in DOINGDIRs if filter_str in str(x)]
# remove = 'BIAS'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'DARK'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'FLAT'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
print ("DOINGDIRs: ", DOINGDIRs)
print ("len(DOINGDIRs): ", len(DOINGDIRs))

#%%
tryagain = False
trynightsky = True

file_age = 80
downsample = 4

#%%
for DOINGDIR in DOINGDIRs[:] :
    _astro_utilities.solving_fits_file(DOINGDIR,
                downsample = downsample,
                count_stars= False,
                tryagain = tryagain,
                file_age = 80,
                tryASTROMETRYNET = False,
                )    

# for DOINGDIR in DOINGDIRs[:] :
#     DOINGDIR = Path(DOINGDIR)
#     print("DOINGDIR", DOINGDIR)
#     SOLVINGDIR = DOINGDIR / _astro_utilities.reduced_dir
#     SOLVINGDIR = DOINGDIR / _astro_utilities.reduced_nightsky_dir
#     SOLVINGDIR = DOINGDIR
#     BADFITSDIR = DOINGDIR / _astro_utilities.bad_fits_dir

#     if not BADFITSDIR.exists() :
#         os.mkdir(str(BADFITSDIR))
#         print(f"{str(BADFITSDIR)} is created...")

#     summary = yfu.make_summary(SOLVINGDIR/"*.fit*",
#                                     verify_fix=True,
#                                     ignore_missing_simple=True,
#                                     )
#     if summary is not None :
#         print("len(summary):", len(summary))
#         print("summary:", summary)
#         #print(summary["file"][0])  
#         df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
#         df_light = df_light.reset_index(drop=True)
#         print("df_light:\n{}".format(df_light))
#     df_light

#     for _, row  in df_light.iterrows():

#         fpath = Path(row["file"])
#         # fpath = Path(df_light["file"][1])
#         print("fpath :" ,fpath)
#         # hdul = fits.open(fpath)
#         Num_stars = _astro_utilities.count_Num_stars(fpath, 
#                     FWHM = 6,
#                     )
#         print("Num_stars :", Num_stars)
#         if Num_stars < 1 :
#             shutil.move(str(fpath), str(BADFITSDIR / fpath.name))
#             print(f"{str(fpath.name)} is moved to Bad_fits...")
#         else : 

#             try : 
#                 _astro_utilities.KevinSolver(fpath, 
#                                                     # solved_dir = None,
#                                                     # downsample = 4,
#                                                     # pixscale = None ,
#                                                     # cpulimit = 15,
#                                                     # tryASTAP = True, 
#                                                     # tryLOCAL = True,
#                                                     tryASTROMETRYNET = False, 
#                                                     makeLOCALsh = True,
#                                                     )
#             except Exception as err :
#                 print("X"*60)

#             try :
#                 fpath = DOINGDIR / _astro_utilities.reduced_dir / fpath.name
#                 print(f"Starting {fpath}") 
#                 _astro_utilities.KevinSolver(fpath, 
#                                                     # solved_dir = None,
#                                                     # downsample = 4,
#                                                     # pixscale = None ,
#                                                     # cpulimit = 15,
#                                                     # tryASTAP = True, 
#                                                     # tryLOCAL = True,
#                                                     tryASTROMETRYNET = True, 
#                                                     # makeLOCALsh = True,
#                                                     )
#             except Exception as err :
#                 print("X"*60)

#             try :
#                 fpath = DOINGDIR / _astro_utilities.reduced_nightsky_dir / fpath.name
#                 print(f"Starting {fpath}") 
#                 _astro_utilities.KevinSolver(fpath, 
#                                                     # solved_dir = None,
#                                                     # downsample = 4,
#                                                     # pixscale = None ,
#                                                     # cpulimit = 15,
#                                                     # tryASTAP = True, 
#                                                     # tryLOCAL = True,
#                                                     tryASTROMETRYNET = True, 
#                                                     # makeLOCALsh = True,
#                                                     )
                
#             except Exception as err :
#                 print("X"*60)
#                 # _Python_utilities.write_log(err_log_file, err)
        
# # os.popen(f"sh __{datetime.now().strftime('%Y%m%d')}_todo_astrometry_solve.sh")