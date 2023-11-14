<<<<<<< HEAD
"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from glob import glob
from pathlib import Path, PosixPath, WindowsPath
import os
import shutil

import ysfitsutilpy as yfu

import _Python_utilities
import _astro_utilities
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
BASEDIR = Path(r"r:\CCD_obs")
BASEDIR = Path("/mnt/Rdata/OBS_data/CCD_keyward_uncompleted_files") 
SAVEDIR = Path("/mnt/Rdata/OBS_data/CCD_new_files")
DOINGDIRs = sorted(_Python_utilities.getFullnameListOfallsubDirs(BASEDIR))
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################
dirnameKEYs = ["OBJECT", "IMAGETYP", "-", "-", 
            "-", "OPTIC", "INSTRUME", "-", "XBINNING"]
fnameKEYs = ["OBJECT", "IMAGETYP", "FILTER", "DATE-OBS", "TIME-OBS",
            "EXPOSURE", "OPTIC", "INSTRUME", "CCD-TEMP", "XBINNING"]
#%%
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)
    fits_in_dir = sorted(list(DOINGDIR.glob('*.fit*')))
    #print("fits_in_dir", fits_in_dir)
    print("len(fits_in_dir)", len(fits_in_dir))

    if len(fits_in_dir) == 0 :
        print(f"There is no fits fils in {DOINGDIR}")
        pass
    else : 
        try:
            print(f"Starting: {str(DOINGDIR.parts[-1])}")
            summary = None 
            summary = yfu.make_summary(DOINGDIR/"*.fit*",
                        #output = save_fpath,
                        verbose = False
                        )
            #print("summary: ", summary)
            print("len(summary)", len(summary))

            for _, row in summary.iterrows():
                # 파일명 출력
                fpath = Path(row["file"])
                print(fpath)
                dirname = ""
                for dirnameKEY in dirnameKEYs :
                    #print(dirnameKEY)
                    try:
                        dirname += row[dirnameKEY].replace("_","-")+"_"
                    except :
                        dirname += "-_"
                print("dirname :", dirname)
                fname = ""
                for fnameKEY in fnameKEYs :
                    #print(dirnameKEY)
                    try:
                        fname += row[fnameKEY].replace("_","-")+"_"
                    except :
                        fname += "-_"
                fname += ".fit"
                print("fname :", fname)
                new_fpath = SAVEDIR / dirname / fname
                if not new_fpath.parents[0].exists() :
                    os.mkdir(str(new_fpath.parents[0]))
                    print(str(new_fpath.parents[0]), "is created...")
                print(str(fpath), str(new_fpath))
                shutil.move(str(fpath), str(new_fpath))

        except Exception as err:
            print("X"*60)
            print(err)
            pass

#%%   
#############################################################################
#Check and delete empty folder....
#############################################################################
for i in range(4):
    BASEDIR = Path("/mnt/Rdata/OBS_data/CCD_keyward_uncompleted_files") 
    DOINGDIR = BASEDIR
    DOINGDIRs = sorted(_Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
    #print ("DOINGDIRs: ", format(DOINGDIRs))
    print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

    for DOINGDIR in DOINGDIRs :    
        if len(os.listdir(str(DOINGDIR))) == 0 :
            shutil.rmtree(f"{str(DOINGDIR)}") # Delete..
            print (f"rmtree {str(DOINGDIR)}")
        else : 
            fpaths = _Python_utilities.getFullnameListOfallFiles(str(DOINGDIR))
            print("fpaths", fpaths)

            for fpath in fpaths[:]:
                print("fpath", fpath)

                if fpath[-4:].lower() in [".txt", "xisf", ".zip", ".png", ".log",
                                            "seal", "tiff", ".axy", "atch", "lved",
                                            "rdls", "xyls", "corr", "xosm", ".ini",
                                            ".wcs", ".csv"] \
                                        and os.path.isfile(fpath):
                    os.remove("{}".format(fpath))
=======
"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from glob import glob
from pathlib import Path, PosixPath, WindowsPath
import os
import shutil

import ysfitsutilpy as yfu

import _Python_utilities
import _astro_utilities
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
BASEDIR = Path(r"r:\CCD_obs")
BASEDIR = Path("/mnt/Rdata/OBS_data/CCD_keyward_uncompleted_files") 
SAVEDIR = Path("/mnt/Rdata/OBS_data/CCD_new_files")
DOINGDIRs = sorted(_Python_utilities.getFullnameListOfallsubDirs(BASEDIR))
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################
dirnameKEYs = ["OBJECT", "IMAGETYP", "-", "-", 
            "-", "OPTIC", "INSTRUME", "-", "XBINNING"]
fnameKEYs = ["OBJECT", "IMAGETYP", "FILTER", "DATE-OBS", "TIME-OBS",
            "EXPOSURE", "OPTIC", "INSTRUME", "CCD-TEMP", "XBINNING"]
#%%
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)
    fits_in_dir = sorted(list(DOINGDIR.glob('*.fit*')))
    #print("fits_in_dir", fits_in_dir)
    print("len(fits_in_dir)", len(fits_in_dir))

    if len(fits_in_dir) == 0 :
        print(f"There is no fits fils in {DOINGDIR}")
        pass
    else : 
        try:
            print(f"Starting: {str(DOINGDIR.parts[-1])}")
            summary = None 
            summary = yfu.make_summary(DOINGDIR/"*.fit*",
                        #output = save_fpath,
                        verbose = False
                        )
            #print("summary: ", summary)
            print("len(summary)", len(summary))

            for _, row in summary.iterrows():
                # 파일명 출력
                fpath = Path(row["file"])
                print(fpath)
                dirname = ""
                for dirnameKEY in dirnameKEYs :
                    #print(dirnameKEY)
                    try:
                        dirname += row[dirnameKEY].replace("_","-")+"_"
                    except :
                        dirname += "-_"
                print("dirname :", dirname)
                fname = ""
                for fnameKEY in fnameKEYs :
                    #print(dirnameKEY)
                    try:
                        fname += row[fnameKEY].replace("_","-")+"_"
                    except :
                        fname += "-_"
                fname += ".fit"
                print("fname :", fname)
                new_fpath = SAVEDIR / dirname / fname
                if not new_fpath.parents[0].exists() :
                    os.mkdir(str(new_fpath.parents[0]))
                    print(str(new_fpath.parents[0]), "is created...")
                print(str(fpath), str(new_fpath))
                shutil.move(str(fpath), str(new_fpath))

        except Exception as err:
            print("X"*60)
            print(err)
            pass

#%%   
#############################################################################
#Check and delete empty folder....
#############################################################################
for i in range(4):
    BASEDIR = Path("/mnt/Rdata/OBS_data/CCD_keyward_uncompleted_files") 
    DOINGDIR = BASEDIR
    DOINGDIRs = sorted(_Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
    #print ("DOINGDIRs: ", format(DOINGDIRs))
    print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

    for DOINGDIR in DOINGDIRs :    
        if len(os.listdir(str(DOINGDIR))) == 0 :
            shutil.rmtree(f"{str(DOINGDIR)}") # Delete..
            print (f"rmtree {str(DOINGDIR)}")
        else : 
            fpaths = _Python_utilities.getFullnameListOfallFiles(str(DOINGDIR))
            print("fpaths", fpaths)

            for fpath in fpaths[:]:
                print("fpath", fpath)

                if fpath[-4:].lower() in [".txt", "xisf", ".zip", ".png", ".log",
                                            "seal", "tiff", ".axy", "atch", "lved",
                                            "rdls", "xyls", "corr", "xosm", ".ini",
                                            ".wcs", ".csv"] \
                                        and os.path.isfile(fpath):
                    os.remove("{}".format(fpath))
>>>>>>> 83e7bdbe06de1d034dcbf1e69d4df772628a78a6
                    print("{} is removed...".format(fpath)) 