# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

이 파일은 fits file의 header에 있는 plate solving 관련 keyword를 삭제해 준다.
필요할때만 사용하면 된다.
"""
#%%
import os, shutil
from glob import glob
from pathlib import Path
from datetime import datetime
import numpy as np

from astropy.io import fits

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

import astro_utilities
import Python_utilities

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
BASEDIR = Path("/mnt/OBS_data") 
BASEDIR = Path("/mnt/Rdata/CCD_obs") 
#BASEDIR = Path("/mnt/OBS_data") 
DOINGDIR = BASEDIR/ astro_utilities.CCD_obs_raw_dir
DOINGDIR = BASEDIR/ astro_utilities.CCD_NEW_dir 
#DOINGDIR = Path(BASEDIR/ "CCD_new_files1")

DOINGDIRs = sorted(Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

#%%
for DOINGDIR in DOINGDIRs :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)
    fits_in_dir = sorted(list(DOINGDIR.glob('*.fit*')))
    #print("fits_in_dir", fits_in_dir)
    print("len(fits_in_dir)", len(fits_in_dir))

    if len(fits_in_dir) == 0 :
        print(f"There is no fits fils in {DOINGDIR}")
        pass
    else : 
        print(f"Starting: {str(DOINGDIR.parts[-1])}")
        summary = None 
        summary = yfu.make_summary(DOINGDIR/"*.fit*",
                    #output = save_fpath,
                    verbose = True
                    )
        print("summary: ", summary)
        print("type(summary): ", type(summary))

        for _, row in summary.iterrows():
            # 파일명 출력
            print (row["file"])
            fpath = Path(row["file"])
            new_fpath = Path(f"{fpath.parents[0]}/{fpath.stem}_clean.fit")
            # fits hedaer 에 있는 wcs 정보를 지운다
            try:
                yfu.wcsremove(fpath, 
                            additional_keys=["COMMENT"],
                            verbose=True,
                            output=new_fpath,
                            ccddata=False,
                            overwrite=True)
                if new_fpath.exists() \
                    and fpath.exists():
                    print("rename", f"{str(new_fpath)}", f"{str(fpath)}")
                    #os.rename(f"{str(new_fpath)}", f"{str(fpath)}")
                    shutil.move(f"{str(new_fpath)}", f"{str(fpath)}")
                    
            except Exception as err :
                print("X"*60)
                with open(err_log_file, 'a') as f:
                    f.write(f'{datetime.now()} ::: {str(fpath)}, {err}\n')