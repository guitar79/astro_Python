# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from glob import glob
from pathlib import Path
import os
import matplotlib.pyplot as plt
from astropy.io import fits
import ysfitsutilpy as yfu

import _astro_utilities
import _Python_utilities

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
verbose = False     
Owrite = False  
tryagain = False 
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

# PROJECDIR = BASEDIR / "C2-Asteroid"
# TODODIR = PROJECDIR / "-_-_-_2022-_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2022-_-_RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2022-_-_RiLA600_STX-16803_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2023-_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2023-_-_RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2023-_-_RiLA600_STX-16803_-_2bin"

# PROJECDIR = BASEDIR / "C3-EXO"
# TODODIR = PROJECDIR / "-_-_-_2024-05_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2024-05_-_RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2024-06_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2024-06_-_RiLA600_STX-16803_-_2bin"

# PROJECDIR = BASEDIR / "C4-Spectra"
# TODODIR = PROJECDIR / "-_-_-_2024-05_TEC140_ASI183MMPro_-_1bin"

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
except : 
    pass

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_" in str(x)])
# print ("DOINGDIRs: ", format(DOINGDIRs))
# print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

# filter_str = 'BL'
# DOINGDIRs = [x for x in DOINGDIRs if filter_str in str(x)]
# remove = 'BIAS'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'DARK'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'FLAT'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
if verbose == True :
    print ("DOINGDIRs: ", DOINGDIRs)
    print ("len(DOINGDIRs): ", len(DOINGDIRs))
#######################################################
#%%
DOINGDIR = Path(BDFDIR[0])
if verbose == True :
    print(f"Starting: {str(DOINGDIR.parts[-1])}")

MASTERDIR = DOINGDIR / _astro_utilities.master_dir

summary = yfu.make_summary(DOINGDIR/"*.fit*", 
                                verify_fix=True,
                                ignore_missing_simple=True,
                            )
if summary is None :
    if verbose == True :
        print(f"summary is None...")
else :
    if verbose == True :
        #print(summary)
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])

# check_filters = summary['FILTER'].drop_duplicates()
# check_filters = check_filters.reset_index(drop=True)
# print("check_filters", check_filters)

# EXPkeys = ['EXPSORE', 'EXPTIME']
# for EXPkey in EXPkeys :
#     if EXPkey in summary :
#         check_exptimes = summary[EXPkey].drop_duplicates()
#         check_exptimes = check_exptimes.reset_index(drop=True)
# print("check_exptimes", check_exptimes)

if (MASTERDIR / "master_bias.fits").exists() and not tryagain:
    print("bias file is already exist....")
else :
    #bias_fits = summary[summary["IMAGETYP"] == "BIAS"]["file"]
    summary_bias = summary.loc[summary["IMAGETYP"] == "BIAS"].copy()
    summary_bias.reset_index(inplace=True)
    # print("summary_bias", summary_bias)

    bias_fits = summary_bias["file"]
    # print("type(bias_fits)", type(bias_fits))
    print("len(bias_fits)", len(bias_fits))
    # print("bias_fits", bias_fits)

    bias_comb = yfu.group_combine(
                    bias_fits.tolist(),
                    type_key = ["IMAGETYP"],
                    type_val = ["BIAS"],
                    group_key = ["EXPTIME"],
                    fmt = "master_bias.fits",  # output file name format
                    outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                    combine = "med",
                    memlimit = 2.e+10,
                    verbose = True,
                )

#dark_fits = summary[summary["IMAGETYP"] == "DARK"]["file"]
summary_dark = summary.loc[summary["IMAGETYP"] == "DARK"].copy()
summary_dark.reset_index(inplace=True)
# print("summary_dark", summary_dark)

# EXPkeys = ['EXPSORE', 'EXPTIME']
# for EXPkey in EXPkeys :
#     if EXPkey in summary_dark :
#         check_exptimes = summary_dark[EXPkey].drop_duplicates()
#         check_exptimes = check_exptimes.reset_index(drop=True)

if 'EXPTIME' in summary_dark :
    check_exptimes = summary_dark['EXPTIME'].drop_duplicates()
    check_exptimes = check_exptimes.reset_index(drop=True)
    print("check_exptimes", check_exptimes)

    for exptime in check_exptimes :
        if (MASTERDIR / f"master_dark_{exptime:.0f}sec.fits" ).exists() and not tryagain :
            print(f"master_dark_{exptime:.0f}sec.fits already exist....")
        else :
            summary_dark_each = summary_dark.loc[summary_dark['EXPTIME'] == exptime]
            dark_fits = summary_dark_each['file']
            # print("type(dark_fits)", type(dark_fits))
            print("len(dark_fits)", len(dark_fits))
            # print("dark_fits", dark_fits)

            # Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
            dark_comb = yfu.group_combine(
                        dark_fits.tolist(),
                        type_key = ["IMAGETYP"],
                        type_val = ["DARK"],
                        group_key = ["EXPTIME"],
                        fmt = "master_dark_{:.0f}sec.fits",  # output file name format
                        outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                        combine = "med",
                        memlimit = 2.e+10,
                        verbose = True,
                    )
            

summary_flat = summary.loc[summary["IMAGETYP"] == "FLAT"].copy()
summary_flat.reset_index(inplace=True)

if 'FILTER' in summary_flat :
    check_filters = summary_flat['FILTER'].drop_duplicates()
    check_filters = check_filters.reset_index(drop=True)
    print("check_filters", check_filters)

    for filter in check_filters :
        if (MASTERDIR / f"master_flat_{filter:s}_norm.fits" ).exists() and not tryagain :
            print(f"master_flat_{filter:s}_norm.fits already exist....")
        else : 
            summary_flat_each = summary_flat.loc[summary_flat['FILTER'] == filter]
            flat_fits = summary_flat_each['file']
            # print("type(flat_fits)", type(flat_fits))
            print("len(flat_fits)", len(flat_fits))
            # print("flat_fits", flat_fits)

            # Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
            flat_comb_norm = yfu.group_combine(
                            flat_fits.tolist(),
                            type_key = ["IMAGETYP"],
                            type_val = ["FLAT"],
                            group_key = ["FILTER"],
                            fmt = "master_flat_{:s}_norm.fits",  # output file name format
                            scale="med_sc", #norm
                            scale_to_0th=False, #norm
                            outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                            combine = "med",
                            memlimit = 2.e+10,
                            verbose=True,
                        )

            # Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
            flat_comb = yfu.group_combine(
                            flat_fits.tolist(),
                            type_key = ["IMAGETYP"],
                            type_val = ["FLAT"],
                            group_key = ["FILTER"],
                            fmt = "master_flat_{:s}.fits",  # output file name format
                            #scale="med_sc", #norm
                            #scale_to_0th=False, #norm
                            outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                            combine = "med",
                            memlimit = 2.e+10,
                            verbose=True,
                        )
