# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

"""
#%%
from glob import glob
from pathlib import Path
import os
import numpy as np
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData
from astropy.io import fits
import matplotlib.pyplot as plt

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
# read all files in base directory for processing
BASEDIR = Path("/mnt/Rdata/OBS_data") 
DOINGDIR = Path(BASEDIR/ "asteroid" / "RiLA600_STX-16803_-_1bin")
DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(DOINGDIR))
DOINGDIRs = sorted([x for x in DOINGDIR.iterdir() if x.is_dir()])
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################

#%%
for DOINGDIR in DOINGDIRs[:1] :
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

        if not MASTERDIR.exists():
            os.makedirs("{}".format(str(MASTERDIR)))
            print("{} is created...".format(str(MASTERDIR)))

        
        summary = yfu.make_summary(DOINGDIR/"*.fit*")
        #print(summary)
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])

        #%%
        if (MASTERDIR / "master_bias.fits").exists():
            print("is exist")
        else :
            try: 
                #bias_fits = summary[summary["IMAGETYP"] == "BIAS"]["file"]
                bias_fits = summary.loc[summary["IMAGETYP"] == "BIAS"].copy()
                bias_fits.reset_index(inplace=True)
                bias_fits = bias_fits["file"]
                print(type(bias_fits))
                print(len(bias_fits))
                print(bias_fits)
                bias_comb = yfu.group_combine(
                                bias_fits.tolist(),
                                type_key = ["IMAGETYP"],
                                type_val = ["BIAS"],
                                group_key = ["EXPTIME"],
                                fmt = "master_bias.fits",  # output file name format
                                outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                                combine = "med",
                                memlimit = 2.e+10,
                                verbose = True
                            )
            except Exception as err :
                print("X"*60)
                _Python_utilities.write_log(err_log_file, err)

        try: 
            #dark_fits = summary[summary["IMAGETYP"] == "DARK"]["file"]
            dark_fits = summary.loc[summary["IMAGETYP"] == "DARK"].copy()
            dark_fits.reset_index(inplace=True)
            dark_fits = dark_fits["file"]
            print(type(dark_fits))
            print(len(dark_fits))
            print(dark_fits)
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
                            verbose = True
                        )
        except Exception as err :
            print("X"*60)
            _Python_utilities.write_log(err_log_file, err)

        try:
            flat_fits = summary[summary["IMAGETYP"] == "FLAT"]["file"] 
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
                            verbose=True
                        )

            # Say dark frames have header OBJECT = "calib" && "IMAGE-TYP" = "DARK"
            flat_comb_norm = yfu.group_combine(
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
                            verbose=True
                        )
        except Exception as err :
            print("X"*60)
            _Python_utilities.write_log(err_log_file, err)