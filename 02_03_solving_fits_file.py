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
import matplotlib.pyplot as plt
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
verbose = True # False
tryagain = False
trynightsky = False
file_age = 365
downsample = 4

#################################################
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
# TODODIR = PROJECDIR / "LIGHT_OON300"
# TODODIR = PROJECDIR / "LIGHT_RILA600"
# TODODIR = PROJECDIR / "LIGHT_SVX80T-x80"
# TODODIR = PROJECDIR / "LIGHT_TEC140"
# TODODIR = PROJECDIR / "LIGHT_TMB130ss"
# TODODIR = PROJECDIR / "LIGHT_TMB130ss-x75"

# PROJECDIR = BASEDIR / "A3_CCD_obs_raw/QSI683ws_1bin/"
# TODODIR = PROJECDIR / "LIGHT_FSQ106ED-x73"
# TODODIR = PROJECDIR / "LIGHT_GSON300"
# TODODIR = PROJECDIR / "LIGHT_RILA600"
# TODODIR = PROJECDIR / "LIGHT_SVX80T"
# TODODIR = PROJECDIR / "LIGHT_SVX80T-x80"
# TODODIR = PROJECDIR / "LIGHT_TMB130ss-x75"

# PROJECDIR = BASEDIR / "A3_CCD_obs_raw/STL-11000M_1bin/"
# TODODIR = PROJECDIR / "LIGHT_FSQ106ED"
# TODODIR = PROJECDIR / "LIGHT_FSQ106ED-x72"
# TODODIR = PROJECDIR / "LIGHT_TEC140-x75"
# TODODIR = PROJECDIR / "LIGHT_TMB130ss"
# TODODIR = PROJECDIR / "LIGHT_TMB130ss-x75"

# PROJECDIR = BASEDIR / "A3_CCD_obs_raw/STX-16803_1bin/"
# TODODIR = PROJECDIR / "LIGHT_RILA600"

PROJECDIR = BASEDIR / "A3_CCD_obs_raw/STX-16803_2bin/"
PROJECDIR = BASEDIR / "A3_CCD_obs_raw/ASI6200MMPro_3bin/"
TODODIR = PROJECDIR / "LIGHT_RILA600"
               
DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(TODODIR))
if verbose == True :
    print ("DOINGDIRs: ", format(DOINGDIRs))
    print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

try : 
    BDFDIR = [x for x in DOINGDIRs if "CAL-BDF" in str(x)]
    if verbose == True :
        print ("BDFDIR: ", format(BDFDIR))
    MASTERDIR = Path(BDFDIR[0]) / _astro_utilities.master_dir
    if not MASTERDIR.exists():
        os.makedirs("{}".format(str(MASTERDIR)))
        if verbose == True :
            print("{} is created...".format(str(MASTERDIR)))
    if verbose == True :
        print ("MASTERDIR: ", format(MASTERDIR))
except Exception as err :
    print("X"*60)
    _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
    pass

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_" in str(x)])
# print ("DOINGDIRs: ", format(DOINGDIRs))
# print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

filter_str = '2025-01-1'
DOINGDIRs = [x for x in DOINGDIRs if filter_str in str(x)]
# remove = 'BIAS'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'DARK'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'FLAT'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
if verbose == True :
    print ("DOINGDIRs: ", DOINGDIRs)
    print ("len(DOINGDIRs): ", len(DOINGDIRs))

#%%
for DOINGDIR in DOINGDIRs[:] :
    try : 
        DOINGDIR = Path(DOINGDIR)
        if verbose == True :
            print(f"Starting: {str(DOINGDIR.parts[-1])}")
        _astro_utilities.solving_fits_file(DOINGDIR,
                # downsample = downsample,
                count_stars= False,
                tryagain = tryagain,
                file_age = 80,
                # tryASTAP = False, # default True
                # tryLOCAL = False,  # default True
                # tryASTROMETRYNET = True,   # default False
                verbose = verbose
                )    
    except Exception as err :
        print("X"*60)
        _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
        pass
