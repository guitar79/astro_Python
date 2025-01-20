# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
import os
from pathlib import Path
from astropy.io import fits

import ysfitsutilpy as yfu

import _astro_utilities
import _Python_utilities

import warnings
warnings.filterwarnings('ignore')

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
########################################################

#%%
verbose = True
tryagain = False
trynightsky = True
#######################################################
BASEDIR = Path("/mnt/Rdata/ASTRO_data")  

PROJECDIR = BASEDIR / "A3_CCD_obs_raw"
TODODIR = PROJECDIR / "ASI2600MC_1bin" / "LIGHT_OON300"
TODODIR = PROJECDIR / "STF-8300M_1bin" / "LIGHT_GSON300"
TODODIR = PROJECDIR / "ASI6200MMPro_3bin" / "LIGHT_RiLA600"
TODODIR = PROJECDIR / "QSI683ws_1bin" / "LIGHT_GSON300"
# TODODIR = PROJECDIR / "QSI683ws_1bin" / "LIGHT_FSQ106ED-x73"

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(TODODIR))
if verbose == True :
    print ("DOINGDIRs: ", format(DOINGDIRs))
    print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

try : 
    BDFDIR = [x for x in DOINGDIRs if "CAL-BDF" in str(x)]
    if verbose == True :
        print ("BDFDIR: ", format(BDFDIR))
    BDFDIR = Path(BDFDIR[0])    
except : 
    BDFDIR = TODODIR
    pass

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_" in str(x)])

# filter_str = 'IC443'
# DOINGDIRs = [x for x in DOINGDIRs if filter_str in str(x)]

if verbose == True :
    print ("DOINGDIRs: ", DOINGDIRs)
    print ("len(DOINGDIRs): ", len(DOINGDIRs))

#######################################################

for DOINGDIR in DOINGDIRs[157:] :
    DOINGDIR = Path(DOINGDIR)
    if verbose == True :
        print("DOINGDIR", DOINGDIR)
    SOLVINGDIR = DOINGDIR

    summary = yfu.make_summary(SOLVINGDIR/"*.fit*")
    if summary is not None :
        if verbose == True :
            print("len(summary):", len(summary))
            print("summary:", summary)
            #print(summary["file"][0])  
        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)
        if verbose == True :
            print("df_light:\n{}".format(df_light))
    df_light

    for _, row  in df_light.iterrows():

        myMP = Multiprocessor()
        num_cpu = 14
        values = []
        fullnames = df_light["file"].tolist()
        num_batches = len(fullnames) // num_cpu + 1

        for batch in range(num_batches):
            myMP.restart()
            for fullname in fullnames[batch*num_batches:(batch+1)*num_batches]:
                fpath = Path(fullname)
                if verbose == True :
                    print("fpath :" ,fpath)
                try :
                    myMP.run(_astro_utilities.KevinSolver, fpath, 
                                                #str(SOLVEDDIR), 
                                                # downsample = 2,
                                                # pixscale = PIXc,
                                                tryASTAP = True,  
                                                tryLOCAL = True,
                                                # makeLOCALsh = True,
                                                # tryASTROMETRYNET = False, 
                                                verbose = verbose,
                                                )
                except Exception as err :
                    print("X"*60)
                    _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
                    pass
            print("+"*10+"\nBatch " + str(batch) + str(DOINGDIR.stem))
            # myMP.wait()

            values.append(myMP.wait())
            print("+"*20+"\nOK batch" + str(batch) + str(DOINGDIR.stem))

    _Python_utilities.write_log(err_log_file, f"{str(DOINGDIR)} is finighed..", verbose=verbose)