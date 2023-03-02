"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

이 파일은 BASEDIR 폴더 안에 있는 모든 fit 파일에 대해서 
fits header에 있는 정보를 바탕으로  
destination_BASEDIR_name 안에 규칙적으로 폴더를 만들어서 저장합니다. 

"""
#%%
from glob import glob
from pathlib import Path, PosixPath, WindowsPath
import os
from astropy.time import Time
from datetime import datetime, timedelta
from astropy.io import fits
import shutil 

import ysfitsutilpy as yfu
#import ysphotutilpy as ypu
import ysvisutilpy as yvu

import astro_utilities
import Python_utilities

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

#######################################################
# read all files in base directory for processing
BASEDIR = Path(r"r:\CCD_obs") 
BASEDIR = Path("/mnt/Rdata/CCD_obs") 
#BASEDIR = Path("/mnt/OBS_data") 
DOINGDIR = BASEDIR / astro_utilities.CCD_NEW_dir
#DOINGDIR = BASEDIR / "CCD_new_files1"
               
DOINGDIRs = sorted(Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

fnameKEYs = ["OBJECT", "IMAGETYP", "FILTER", "DATE-OBS", 
            "EXPOSURE", "OPTIC", "CCDNAME", "CCD-TEMP", "XBINNING"]

#%%
for DOINGDIR in DOINGDIRs[:] : 
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR: ", DOINGDIR)
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
            #print (row["file"])
            fpath = Path(row["file"])            
            hdul = fits.open(str(fpath))
            print("fpath: ", fpath)
            try:
                
                for fnameKEY in fnameKEYs: 
                    print(f"{fnameKEY}: ", hdul[0].header[fnameKEY])

                try :
                    ccdtemp = str(int(hdul[0].header["CCD-TEMP"]))
                except : 
                    ccdtemp = "N"
                print("ccdtemp: ", ccdtemp)

                new_fname = hdul[0].header["OBJECT"]+"_"+hdul[0].header["IMAGETYP"]+"_"+hdul[0].header["FILTER"]+"_"
                new_fname += hdul[0].header["DATE-OBS"][:19].replace("T","-").replace(":","-")+"_"
                new_fname += str(int(hdul[0].header["EXPOSURE"]))+"sec_"
                new_fname += hdul[0].header["OPTIC"]+"_"+hdul[0].header["CCDNAME"]+"_"       
                new_fname += ccdtemp+"c_"+str(int(hdul[0].header["XBINNING"]))+"bin.fit"
                print("new_fname: ", new_fname)
                hdul.close()
                new_folder = astro_utilities.get_new_foldername_from_filename(new_fname)
                #print("new_folder: ", new_folder)

                new_fpath =  BASEDIR /astro_utilities.CCD_obs_raw_dir / new_folder / new_fname
                print("new_fpath: ", new_fpath)

                if not new_fpath.parent.exists():
                    os.makedirs(str(new_fpath.parent))
                    print(f'{str(new_fpath.parts[-2])} is created')  
                #os.makedirs(str(new_fpath.parent), exist_ok=True)
                if new_fpath.exists():
                    duplicate_fpath = BASEDIR / astro_utilities.CCD_duplicate_dir / new_fpath.name
                    #os.rename(str(new_fpath), str(duplicate_fpath))
                    shutil.move(str(new_fpath), str(duplicate_fpath))
                    print (f"move duplicate file to {str(duplicate_fpath)}")

                #os.rename(str(fpath), str(new_fpath))
                shutil.move(str(fpath), str(new_fpath))
                print(f"move {str(fpath.name)} to {str(new_fpath)}")
            except Exception as err :
                print("X"*60)
                with open(err_log_file, 'a') as f:
                    f.write(f'{datetime.now()} ::: {str(fpath)}, {err}\n')
                pass
#%%   
#############################################################################
#Check and delete empty folder....
#############################################################################
for i in range(4):
    DOINGDIR = Path( BASEDIR/ astro_utilities.CCD_NEW_dir)           
    DOINGDIRs = sorted(Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
    #print ("DOINGDIRs: ", format(DOINGDIRs))
    print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

    for DOINGDIR in DOINGDIRs :    
        if len(os.listdir(str(DOINGDIR))) == 0 :
            shutil.rmtree(f"{str(DOINGDIR)}") # Delete..
            print (f"rmtree {str(DOINGDIR)}")
        else : 
            fpaths = Python_utilities.getFullnameListOfallFiles(str(DOINGDIR))
            print("fpaths", fpaths)

            for fpath in fpaths[:]:
                print("fpath", fpath)

                if fpath[-4:].lower() in [".txt", "xisf", ".zip", ".png", ".log",
                                            "seal", "tiff", ".axy", "atch", "lved",
                                            "rdls", "xyls", "corr", "xosm", ".ini",
                                            ".wcs", ".csv"] \
                                        and os.path.isfile(fpath):
                    os.remove("{}".format(fpath))
                    print("{} is removed...".format(fpath)) 
