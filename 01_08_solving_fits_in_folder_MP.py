# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

이 파일은 BASEDIR 폴더 안에 있는 모든 fit 파일에 대해서 
Astrometry plate solving을 수행합니다.
이미 solving이 완료된 파일은 건너뛰고, 

다만 Multiprocessing을 적용하여 여러개의 core를 사용합니다.
(터짐 주의)
"""
#%%
import os
import subprocess
from datetime import datetime
from astropy.io import fits
from pathlib import Path
import shutil 
import _Python_utilities
import _astro_utilities

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

#######################################################
# read all files in base directory for processing
BASEDIR = Path(r"r:\CCD_obs")
BASEDIR = Path("/mnt/Rdata/OBS_data") 
#BASEDIR = Path("/mnt/OBS_data") 
DOINGDIR = Path(BASEDIR/ "RnE_2022/GSON300_STF-8300M")
#DOINGDIR = Path(BASEDIR/ "CCD_new_files1")

#DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(DOINGDIR))
DOINGDIRs = sorted([x for x in DOINGDIR.iterdir() if x.is_dir()])
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################
#%%
for DOINGDIR in DOINGDIRs[:] :
    result = ''
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
        REDUCEDDIR = DOINGDIR / _astro_utilities.reduced_dir2
        SOLVEDDIR = DOINGDIR / _astro_utilities.solved_dir2

        if not SOLVEDDIR.exists():
            os.makedirs("{}".format(str(SOLVEDDIR)))
            print("{} is created...".format(str(SOLVEDDIR))) 

        summary = yfu.make_summary(REDUCEDDIR/"*Light*.fit*")
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])

        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)
        print("df_light:\n{}".format(df_light))

        if df_light.empty:
            print("The dataframe(df_light) is empty")
            pass
        else:
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
                    #myMP.run(_astro_utilities.KevinSolver, fullname, solved_dir)
                    myMP.run(_astro_utilities.AstrometrySolver, fullname, 
                                                                str(SOLVEDDIR), 
                                                                4,
                                                                0.93,)

                print("Batch " + str(batch))
                #myMP.wait()

                values.append(myMP.wait())
                print("OK batch" + str(batch))