# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

이 파일은 관측 자료를 전처리 해준다.

"""
#%%
from glob import glob
from pathlib import Path
import numpy as np
import os
import astropy.units as u
from ccdproc import CCDData, ccd_process

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
OBSRAWDIR = astro_utilities.CCD_obs_dir
BASEDIR = astro_utilities.base_dir
BASEDIRs = sorted(Python_utilities.getFullnameListOfsubDir(BASEDIR))
print ("BASEDIRs: {}".format(BASEDIRs))
print ("len(BASEDIRs): {}".format(len(BASEDIRs)))

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
def f(fpath):
    fpath = Path(row["file"])
    ccd = yfu.load_ccd(fpath)
    filt = ccd.header["FILTER"]
    expt = ccd.header["EXPTIME"]
    red = yfu.ccdred(
        ccd,
        output=Path(f"{REDUCEDDIR}/{fpath.stem}.fits"),
        mdarkpath=str(MASTERDIR / "master_dark_{:.0f}sec.fits".format(expt)),
        mflatpath=str(MASTERDIR / "master_flat_{}_norm.fits".format(filt.upper())),
        # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
        overwrite=True
    )
    return 0
    
#%%
for BASEDIR in BASEDIRs[4:] :
    print ("Starting...\n{}".format(BASEDIR))

    BASEDIR = Path(BASEDIR)
    
    MASTERDIR = BASEDIR / astro_utilities.master_dir
    REDUCEDDIR = BASEDIR / astro_utilities.reduced_dir

    if not REDUCEDDIR.exists():
        os.makedirs(str(REDUCEDDIR))
        print("{} is created...".format(str(REDUCEDDIR)))

    #%%
    summary = yfu.make_summary(BASEDIR/"*.fit")
    #print(summary)
    print("len(summary):", len(summary))

    df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()

    if df_light.empty:
            print("The dataframe(df_light) is empty")
            pass
    else:
        df_light = df_light.reset_index(drop=True)
        print("df_light:\n{}".format(df_light))

    # %%
    for _, row in df_light.iterrows():
        try:
            fpath = Path(row["file"])
            ccd = yfu.load_ccd(fpath)
            filt = ccd.header["FILTER"]
            expt = ccd.header["EXPTIME"]
            red = yfu.ccdred(
                ccd,
                output=Path(f"{REDUCEDDIR}/{fpath.stem}.fits"),
                mdarkpath=str(MASTERDIR / "master_dark_{:.0f}sec.fits".format(expt)),
                mflatpath=str(MASTERDIR / "master_flat_{}_norm.fits".format(filt.upper())),
                # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                overwrite=True
            )
        except Exception as err: 
            print ('Error messgae .......')
            print (err)
 