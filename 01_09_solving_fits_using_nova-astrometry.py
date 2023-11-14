<<<<<<< HEAD:01_09_solving_fits_using_nova-astrometry.py
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from glob import glob
from pathlib import Path
import os

import _astro_utilities
import _Python_utilities
import ysfitsutilpy as yfu

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
DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/')
DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/STX-16803_1bin' )
#DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/STX-16803_2bin' )
#DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/QSI683ws_1bin/LIGHT_RiLA600')
DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/STF-8300M_2bin')
DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/STF-8300M_1bin')
DOINGDIR = Path('/mnt/Rdata/OBS_data/ccd_test_folder')


DOINGDIRs = sorted(_Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
print ("len(DOINGDIRs): ", len(DOINGDIRs))
remove = '/Cal'
DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
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
#######################################################
ast = AstrometryNet()
ast.api_key = 'bldvwzzuvktnwfph'


#%%
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)
    fits_in_dir = sorted(list(DOINGDIR.glob('*.fit')))
    #print("fits_in_dir", fits_in_dir)
    print("len(fits_in_dir)", len(fits_in_dir))

    if len(fits_in_dir) == 0 :
        print(f"There is no fits fils in {DOINGDIR}")
        pass
    else : 
        print(f"Starting: {str(DOINGDIR.parts[-1])}")
        summary = yfu.make_summary(DOINGDIR/"*.fit")
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])

        for _, row in summary.iterrows():
            # 파일명 출력
            fpath = Path(row["file"])
            print(fpath)
            SOLVE, ASTAP, LOCAL = _astro_utilities.checkPSolve(fpath)
            print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)
            if SOLVE : 
                print(f"{fpath.name} is already solved...")
            else : 
                try: 
                    submission_id = None
                    solve_timeout = 600

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
                except Exception as err :
                    print(err)
                    continue
=======
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from glob import glob
from pathlib import Path
import os

import _astro_utilities
import _Python_utilities
import ysfitsutilpy as yfu

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
DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/')
DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/STX-16803_1bin' )
#DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/STX-16803_2bin' )
#DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/QSI683ws_1bin/LIGHT_RiLA600')
DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/STF-8300M_2bin')
DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw/STF-8300M_1bin')
DOINGDIR = Path('/mnt/Rdata/OBS_data/ccd_test_folder')


DOINGDIRs = sorted(_Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
print ("len(DOINGDIRs): ", len(DOINGDIRs))
remove = '/Cal'
DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
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
#######################################################
ast = AstrometryNet()
ast.api_key = 'bldvwzzuvktnwfph'


#%%
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)
    fits_in_dir = sorted(list(DOINGDIR.glob('*.fit')))
    #print("fits_in_dir", fits_in_dir)
    print("len(fits_in_dir)", len(fits_in_dir))

    if len(fits_in_dir) == 0 :
        print(f"There is no fits fils in {DOINGDIR}")
        pass
    else : 
        print(f"Starting: {str(DOINGDIR.parts[-1])}")
        summary = yfu.make_summary(DOINGDIR/"*.fit")
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])

        for _, row in summary.iterrows():
            # 파일명 출력
            fpath = Path(row["file"])
            print(fpath)
            SOLVE, ASTAP, LOCAL = _astro_utilities.checkPSolve(fpath)
            print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)
            if SOLVE : 
                print(f"{fpath.name} is already solved...")
            else : 
                try: 
                    submission_id = None
                    solve_timeout = 600

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
                except Exception as err :
                    print(err)
                    continue
>>>>>>> 83e7bdbe06de1d034dcbf1e69d4df772628a78a6:01_08_solving_fits_using_nova-astrometry.py
