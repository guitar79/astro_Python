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

이 파일은 base_dir 폴더 안에 있는 모든 fit 파일에 대해서 
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
base_dir = "../RnE_2022/"

c_method = 'median'
master_dir = "master_files_ys"
reduced_dir = "reduced"
solved_dir = "solved"

#%%
base_dirs = sorted(Python_utilities.getFullnameListOfallsubDirs(base_dir))
print ("base_dirs1: {}".format(base_dirs))
base_dirs = [w for w in base_dirs \
        if not (w.endswith("{}/".format(master_dir)) \
            or w.endswith("{}/".format(reduced_dir)) \
            or w.endswith("{}/".format(solved_dir)) \
                or w.endswith(".fits"))]
print ("base_dirs2: {}".format(base_dirs))

#%%
#base_dir = Path("../RnE_2022/KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin/")


#%%
for base_dir in base_dirs :
    print ("Starting...\n{}".format(base_dir))

    base_dir = Path(base_dir)

    if not (base_dir/reduced_dir/solved_dir).exists():
        os.makedirs(str((base_dir/reduced_dir/solved_dir)))

    summary = yfu.make_summary(base_dir/reduced_dir/"*.fits")

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
        for fullname  in fullnames[batch*num_batches:(batch+1)*num_batches]:
            myMP.run(astro_utilities.KevinSolver, fullname, solved_dir)
        print("Batch " + str(batch))
        #myMP.wait()

        values.append(myMP.wait())
        print("OK batch" + str(batch))  

    #############################################################################
    #Check existence tmp file and rename ...
    #############################################################################
    summary_new = yfu.make_summary(base_dir/reduced_dir/solved_dir/"*.new")
    print ("summary_new: {}".format(summary_new))

    #%%
    n = 0
    try:
        for _, row in summary_new.iterrows():
            n += 1
            print('#'*40,
                "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(summary_new), (n/len(summary_new))*100, os.path.basename(__file__)))
            print ("Starting...\nfullname: {}".format(row["file"]))

        
            shutil.move(r"{}".format(row["file"]), \
                            r"{}.fits".format(row["file"][:-4]))

    except Exception as err:
        Python_utilities.write_log(err_log_file,
                    '{2} ::: {0} There is no {1} '.format(err, row["file"], datetime.now())) 