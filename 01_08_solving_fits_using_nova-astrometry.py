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
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData

import ysfitsutilpy as yfu
import ysphotutilpy as ypu

import _astro_utilities
import _Python_utilities

import shutil
from astroquery.astrometry_net import AstrometryNet

from astropy.io import fits

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
DOINGDIR = Path(BASEDIR / "ccd_test_folder")

#DOINGDIRs = sorted(_Python_utilities.getFullnameListOfallfilessubDirs(DOINGDIR))
#print ("len(DOINGDIRs): ", len(DOINGDIRs))
#######################################################
fpaths = sorted(list(DOINGDIR.glob("*.fit")))
print("fpaths;", fpaths)

ast = AstrometryNet()
ast.api_key = 'bldvwzzuvktnwfph'
#%%
for fpath in fpaths[:] :
    #fpath = Path(fpaths[0])
    #fpath = fpaths[5]
    print(fpath)

    submission_id = None
    solve_timeout = 600

    try:

        hdul = fits.open(str(fpath))
        if not 'B_1_1' in hdul[0].header :
            print("it's not solved file...")
            try_again = True    

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

        if wcs_header:
            # Code to execute when solve succeeds
            print("fits file solved successfully...")
            shutil.copy(str(fpath), str(fpath.parents[0] / f"{fpath.stem}.tmp"))

            with fits.open(str(fpath.parents[0] / f"{fpath.stem}.tmp"), mode='update') as filehandle:
                print("filehandle[0].header :", filehandle[0].header)
                for card in wcs_header :
                    try: 
                        print(card, wcs_header[card], wcs_header.comments[card])
                        filehandle[0].header.set(card, wcs_header[card], wcs_header.comments[card])
                    except : 
                        print(card)
                filehandle.flush

            shutil.move(str(fpath.parents[0] / f"{fpath.stem}.tmp"), str(fpath.parents[0] / f"{fpath.stem}.fits"))
            print(str(fpath.parents[0] / f'{fpath.stem}.fits')+" is created...")
        else:
            # Code to execute when solve fails
            print("fits file solving failure...")
 
    except Exception as err:
        print("X"*30, f'\n{err}')
        # _Python_utilities.write_log(err_log_file, err)
        pass