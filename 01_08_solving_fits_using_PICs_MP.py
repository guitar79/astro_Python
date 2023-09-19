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
from datetime import datetime
from astropy.io import fits
from pathlib import Path
import _Python_utilities
import _astro_utilities

import ysfitsutilpy as yfu

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

DOINGDIR = Path('/mnt/Rdata/OBS_data/CCD_obs_raw')

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
remove = 'BIAS'
DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
remove = 'DARK'
DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
remove = 'FLAT'
DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################
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
                    fpath = Path(fullname)
                    print("fpath.name ;", fpath.name)
                    ccd = yfu.load_ccd(str(fpath))
                    
                    if 'PIXSCALE' in ccd.header:
                        PIXc = ccd.header['PIXSCALE']
                    else : 
                        PIXc = _astro_utilities.calPixScale(ccd.header['FOCALLEN'], ccd.header['XPIXSZ'])
                    PIXc = _astro_utilities.calPixScale(ccd.header['FOCALLEN'], ccd.header['XPIXSZ'])
                    print("PIXc : ", PIXc)

                    myMP.run(_astro_utilities.KevinSolver, (str(fpath)), 
                                                    #str(SOLVEDDIR), 
                                                    downsample = 4,
                                                    pixscale = PIXc,)

                print("Batch " + str(batch))
                #myMP.wait()

                values.append(myMP.wait())
                print("OK batch" + str(batch))