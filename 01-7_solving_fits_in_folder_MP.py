# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

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
plate solving을 수행합니다.
이미 solving이 완료된 파일은 건너뛰고, 
먼조 ASTAP로 시도하고, 실패할 경우 Astrometry로 시도합니다.

다만 Multiprocessing을 적용하여 여러개의 core를 사용합니다.
(터짐 주의)
"""
#%%
from glob import glob
from pathlib import Path
import os
import numpy as np
import shutil
from datetime import datetime 
from astropy.io import fits

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

import Python_utilities
import astro_utilities


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
from multiprocessing import Process, Queue
class Multiprocessor():
    def __init__(self):
        self.processes = []
        self.queue = Queue()

    @staticmethod
    def _wrapper(func, queue, args, kwargs):
        ret = func(*args, **kwargs)
        queue.put(ret)

    def restart(self):
        self.processes = []
        self.queue = Queue()

    def run(self, func, *args, **kwargs):
        args2 = [func, self.queue, args, kwargs]
        p = Process(target=self._wrapper, args=args2)
        self.processes.append(p)
        p.start()

    def wait(self):
        rets = []
        for p in self.processes:
            ret = self.queue.get()
            rets.append(ret)
        for p in self.processes:
            p.join()
        return rets
#######################################################

#%%
#######################################################
# read all files in base directory for processing
BASEDIR = "../RnE_2022/"
BASEDIR = "../RnE_2022/RiLA600_STX-16803_2bin/"

#%%
BASEDIRs = sorted(Python_utilities.getFullnameListOfsubDir(BASEDIR))
print ("BASEDIRs: {}".format(BASEDIRs))
print ("len(BASEDIRs): {}".format(len(BASEDIRs)))
#BASEDIRs = ["{}{}".format(w, reduced_dir) for w in BASEDIRs]
#print ("BASEDIRs2: {}".format(BASEDIRs))

#%%
#BASEDIR = Path("../RnE_2022/KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin/")

#%%
for BASEDIR in BASEDIRs [:]:
    print ("Starting...\n{}".format(BASEDIR))

    BASEDIR = Path(BASEDIR)

    MASTERDIR = BASEDIR / astro_utilities.master_dir
    REDUCEDDIR = BASEDIR / astro_utilities.reduced_dir
    SOLVEDDIR = BASEDIR / astro_utilities.solved_dir    
    RESULTDIR = BASEDIR / astro_utilities.DAOfinder_result_dir
    REDUCEDDIR = BASEDIR / astro_utilities.reduced_dir2
    SOLVEDDIR = BASEDIR / astro_utilities.solved_dir2    

    if not (SOLVEDDIR).exists():
        os.makedirs(str(SOLVEDDIR))

    summary = yfu.make_summary(REDUCEDDIR/"*.fits")

    df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
    df_light = df_light.reset_index(drop=True)
    print("df_light:\n{}".format(df_light))

    myMP = Multiprocessor()
    num_cpu = 6
    values = []
    fullnames = df_light["file"].tolist()
    num_batches = len(fullnames) // num_cpu + 1

    for batch in range(num_batches):
        myMP.restart()
        for fullname in fullnames[batch*num_batches:(batch+1)*num_batches]:
            #myMP.run(astro_utilities.KevinSolver, fullname, solved_dir)
            myMP.run(astro_utilities.AstrometrySolver, fullname, str(SOLVEDDIR))

        print("Batch " + str(batch))
        #myMP.wait()

        values.append(myMP.wait())
        print("OK batch" + str(batch))  
