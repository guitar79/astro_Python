# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from glob import glob
from pathlib import Path
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData
from astroquery.astrometry_net import AstrometryNet

import ysfitsutilpy as yfu

import _astro_utilities
import _Python_utilities
import _tool_visualization

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
TODODIR = PROJECDIR / "_-_-_2024-06_-_GSON300_STF-8300M_-_1bin"
TODODIR = PROJECDIR / "RiLA600_STX-16803_-_1bin"
TODODIR = PROJECDIR / "RiLA600_STX-16803_-_2bin"

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

CALDIR = [x for x in DOINGDIRs if "CAL-BDF" in str(x)]
MASTERDIR = Path(CALDIR[0]) / _astro_utilities.master_dir
if not MASTERDIR.exists():
    os.makedirs("{}".format(str(MASTERDIR)))
    print("{} is created...".format(str(MASTERDIR)))

print ("MASTERDIR: ", format(MASTERDIR))

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_" in str(x)])
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

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
ast = AstrometryNet()

# ger from nova.astrometry.net
ast.api_key = 'bldvwzzuvktnwfph' #must changed...

#%%
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)

    # if "RiLA600_STX-16803" in str(DOINGDIR.parts[-2]) :
    READINGDIR = DOINGDIR / _astro_utilities.reduced_nightsky_dir
    # if "GSON300_STF-8300M_-_1bin" in str(DOINGDIR.parts[-2]) :
    # DOINGDIR = DOINGDIR / _astro_utilities.reduced_dir
    
    summary = yfu.make_summary(READINGDIR/"*.fit*")
    if summary is not None :
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])  
        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)
        print("df_light:\n{}".format(df_light))

        for _, row  in df_light.iterrows():
            fpath = Path(row["file"])
            print(fpath)
            hdul = fits.open(fpath)

            submission_id = None
            solve_timeout = 600

            if 'PIXSCALE' in hdul[0].header:
                PIXc = hdul[0].header['PIXSCALE']
            else : 
                PIXc = _astro_utilities.calPixScale(hdul[0].header['FOCALLEN'], 
                                                    hdul[0].header['XPIXSZ'],
                                                    hdul[0].header['XBINNING'])
            print("PIXc : ", PIXc)
            hdul.close()

            SOLVE, ASTAP, LOCAL = _astro_utilities.checkPSolve(fpath)
            print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)

            if SOLVE :
                print(f"{fpath.name} is already solved...")
            else :             
                try_again = True                

                try : 
                    
                    while try_again:
                        try:
                            if not submission_id:
                                wcs_header = ast.solve_from_image(str(fpath),
                                                    force_image_upload=True,
                                                    solve_timeout = solve_timeout,
                                                    submission_id=submission_id)
                            else:
                                wcs_header = ast.monitor_submission(submission_id,
                                                                    solve_timeout = solve_timeout)
                        except TimeoutError as e:
                            submission_id = e.args[1]
                        else:
                            # got a result, so terminate
                            try_again = False

                    if not wcs_header:
                        # Code to execute when solve fails
                        print("fits file solving failure...")

                    else:
                        # Code to execute when solve succeeds
                        print("fits file solved successfully...")

                        with fits.open(str(fpath), mode='update') as hdul:
                            for card in wcs_header :
                                try: 
                                    print(card, wcs_header[card], wcs_header.comments[card])
                                    hdul[0].header.set(card, wcs_header[card], wcs_header.comments[card])
                                except : 
                                    print(card)
                            hdul.flush

                        print(str(fpath)+" is created...")
                
                except Exception as err: 
                    print("Err :", err)
                    continue