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

import numpy as np

import astropy.units as u
from astropy.stats import sigma_clip
from astropy.io import fits
from ccdproc import combine, ccd_process, CCDData

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

from snuo1mpy import Preprocessor

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
DOINGDIR = Path( BASEDIR/ astro_utilities.CCD_obs_raw_dir)

DOINGDIRs = sorted(Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

#%%
for DOINGDIR in DOINGDIRs :
    # 디렉토리 하나씩 loop...
    print ("Starting...\n{}".format(DOINGDIR))

    DOINGDIR = Path(DOINGDIR)
    
    RESULTDIR = DOINGDIR / astro_utilities.DAOfinder_result_dir
    SOLVEDDIR = DOINGDIR / astro_utilities.solved_dir
    MASTERDIR = DOINGDIR / astro_utilities.master_dir
    REDUCEDDIR = DOINGDIR / astro_utilities.reduced_dir
    MASTERDIR = DOINGDIR / astro_utilities.master_dir

    #try : 
    #파일 목록을 DataFrame으로 만들어서 이용하자
    #summary = yfu.make_summary(REDUCEDDIR / "*.fit*")
    summary = yfu.make_summary(DOINGDIR / "*.fit*")
                        #keywords = ["DATE-OBS", "FILTER", "OBJECT", "IMAGETYP"],  # header keywords; actually it is case-insensitive
                        #fname_option = 'name',  # 'file' column will contain only the name of the file (not full path)
                        #output = "{}summary.csv".format(DOINGDIR),
                        #sort_by = "DATE-OBS"  # 'file' column will be sorted based on "DATE-OBS" value in the header
                        #)
    print("summary:\n {}".format(summary))

    try:
        # light frame  만 선택
        df_light = summary[summary["IMAGETYP"] == "LIGHT"]
        #print ("df_light: {}".format(df_light))
        #print ("len(df_light): {}".format(len(df_light)))

        for _, row in df_light.iterrows():
            # 파일명 출력
            print (row["file"])
            fpath = Path(row["file"])
            new_fpath = Path(f"{fpath.parents[0]}/{fpath.stem}_clean.fit")
            # fits hedaer 에 있는 wcs 정보를 지운다
            yfu.wcsremove(fpath, 
                        additional_keys=["COMMENT"],
                        verbose=True,
                        output=new_fpath,
                        ccddata=True,
                        overwrite=True)
            if new_fpath.exists() \
                and fpath.exists():
                print("rename", f"{str(new_fpath)}", f"{str(fpath)}")
                #os.rename(f"{str(new_fpath)}", f"{str(fpath)}")
                shutil.move(f"{str(new_fpath)}", f"{str(fpath)}")
                
    except Exception as err :
        print("X"*60)
        print('{0}'.format(err))

