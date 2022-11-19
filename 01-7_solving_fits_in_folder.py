# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

#first time
cd ~/Downloads/ && git clone https://github.com/ysBach/ysvisutilpy && cd ysvisutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/ysfitsutilpy && cd ysfitsutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/ysphotutilpy && cd ysphotutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/SNUO1Mpy && cd SNUO1Mpy && git pull && pip install -e . && cd ..

# second time...
cd ~/Downloads/ysvisutilpy && git pull && pip install -e . 
cd ~/Downloads/ysfitsutilpy && git pull && pip install -e . 
cd ~/Downloads/ysphouutilpy && git pull && pip install -e . 
cd ~/Downloads/SNUO1Mpy && git pull && pip install -e . 

이 파일은 BASEDIR 폴더 안에 있는 모든 fit 파일에 대해서 
plste solving을 수행합니다.
이미 solving이 완료된 파일은 건너뛰고, 
먼저 ASTAP로 시도하고, 실패할 경우 Astrometry로 시도합니다. 

"""
#%%
import os
import subprocess
from datetime import datetime
from astropy.io import fits
from pathlib import Path
import shutil 
import Python_utilities
import astro_utilities

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

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
BASEDIR = "../RnE_2022/"
BASEDIR = "../RnE_2022/RiLA600_STX-16803_2bin/"

c_method = "median"
master_dir = "master_files_ys"
reduced_dir = "reduced"
solved_dir = "solved"
DAOfinder_result = "DAOfinder_result"

#%%
BASEDIRs = sorted(Python_utilities.getFullnameListOfsubDir(BASEDIR))
print ("BASEDIRs: {}".format(BASEDIRs))

#%%
for BASEDIR in BASEDIRs :
    print ("Starting...\n{}".format(BASEDIR))

    BASEDIR = Path(BASEDIR)
    
    RESULTDIR = BASEDIR / DAOfinder_result
    SOLVEDDIR = BASEDIR / solved_dir
    MASTERDIR = BASEDIR / master_dir
    REDUCEDDIR = BASEDIR / reduced_dir
    MASTERDIR = BASEDIR / master_dir

    if not SOLVEDDIR.exists():
        os.makedirs("{}".format(str(SOLVEDDIR)))
        print("{} is created...".format(str(SOLVEDDIR)))

    summary = yfu.make_summary(REDUCEDDIR / "*.fits")

    df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
    df_light = df_light.reset_index(drop=True)
    print("df_light:\n{}".format(df_light))

    #%%
    n = 0
    for _, row  in df_light.iterrows():

        #fullname = fullnames[5]
        n += 1
        print('#'*40,
            "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(df_light), (n/len(df_light))*100, os.path.basename(__file__)))
        print ("Starting...\nfullname: {}".format(row["file"]))

        #astro_utilities.KevinSolver(row["file"], solved_dir)
        #astro_utilities.AstrometrySolver(df_light["file"][0], BASEDIR/solved_dir)
        solved = astro_utilities.AstrometrySolver(row["file"], str(SOLVEDDIR))
        #astro_utilities.AstrometrySolver(row["file"], "../{}".format(solved_dir))
                
    # #############################################################################
    # #Check existence tmp file and rename ...
    # #############################################################################
    # summary_new = yfu.make_summary(BASEDIR/solved_dir/"*.new")
    # print ("summary_new: {}".format(summary_new))

    # #%%
    # n = 0
    # try:
    #     for _, row in summary_new.iterrows():
    #         n += 1
    #         print('#'*40,
    #             "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(summary_new), (n/len(summary_new))*100, os.path.basename(__file__)))
    #         print ("Starting...\nfullname: {}".format(row["file"]))
    #         shutil.move(r"{}".format(row["file"]), \
    #                     r"{}.fits".format(row["file"][:-4]))

    # except Exception as err:
    #     Python_utilities.write_log(err_log_file,
    #             '{2} ::: {0} There is no {1} '.format(err, row["file"], datetime.now())) 
