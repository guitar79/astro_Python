"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from glob import glob
from pathlib import Path
import os
import shutil 
import ysfitsutilpy as yfu

import _astro_utilities
import _Python_utilities
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
# Set directory variables.
#######################################################

# BASEDIR = Path(r"r:\OBS_data")   # for windows
BASEDIR = Path("/mnt/Rdata/ASTRO_data")  # for ubuntu
# BASEDIR = Path("/Volumes/OBS_data")  # for mac OS
 
DOINGDIR = BASEDIR/ _astro_utilities.CCD_NEWUP_dir

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################   

#%%
fnameKEYs = ["OBJECT", "IMAGETYP", "FILTER", "DATE-OBS", 
            "EXPOSURE", "OPTIC", "CCDNAME", "CCD-TEMP", "XBINNING"]
Owrite = False

#%%
for DOINGDIR in DOINGDIRs[:] : 
    DOINGDIR = Path(DOINGDIR)
    print(f"Starting: {str(DOINGDIR.parts[-1])}")
    try: 
        summary = yfu.make_summary(DOINGDIR/"*.fit*",)
        if summary is not None : 
            print("summary: ", summary)
            print("len(summary)", len(summary))
        
            for _, row in summary.iterrows():
                fpath = Path(row["file"])
                print (f"starting {fpath.name}...")
                new_fname = ""
                suffix = ".fit"
                try:
                    for KEY in fnameKEYs :
                        if KEY in ["OBJECT", "IMAGETYP", "FILTER", 
                            "OPTIC", "CCDNAME"] :
                            new_fname += str(row[KEY])+"_"
                        
                        if KEY == "DATE-OBS" : 
                            new_fname += row[KEY][:19].replace("T","-").replace(":","-")+"_"

                        if KEY == "EXPOSURE" : 
                            new_fname += str(int(row[KEY]))+"sec_"

                        if KEY == "CCD-TEMP" : 
                            try:
                                new_fname += str(int(row[KEY]))+"c_"
                            except:
                                new_fname += (row[KEY])+"c_"
                        if KEY == "XBINNING" : 
                            new_fname += str(row[KEY])+"bin"+suffix
                
                    print(new_fname)                      
                    new_folder = _astro_utilities.get_new_foldername_from_filename(new_fname)
                    #print("new_folder: ", new_folder)
                    new_fpath =  BASEDIR /_astro_utilities.A3_CCD_obs_raw_dir / new_folder / new_fname
                    print("new_fpath: ", new_fpath)

                    if not new_fpath.parents[0].exists():
                        os.makedirs(f'{new_fpath.parents[0]}')
                        print(f'{new_fpath.parts[-2]} is created')  
                
                    if new_fpath.exists() :
                        print(f'{new_fpath} is already exist')
                        duplicate_fpath = BASEDIR / _astro_utilities.CCD_duplicate_dir / new_fpath.name
                        if Owrite == False:
                            shutil.move(fpath, duplicate_fpath)
                            print(f'{fpath.parts[-1]} is move to duplicate folder...')
                        else :
                            shutil.move(str(fpath), str(new_fpath))
                            print(f"move {str(fpath.name)} to {str(new_fpath)}")
                    else : 
                        shutil.move(str(fpath), str(new_fpath))
                        print(f"move {str(fpath.name)} to {str(new_fpath)}")
                    
                except Exception as err:
                    print("X"*30, f'\n{err}')
                    #_Python_utilities.write_log(err_log_file, err)
                    pass
    except Exception as err:
        print("X"*30, f'\n{err}')
        #_Python_utilities.write_log(err_log_file, err)
        pass
#%%   
#############################################################################
#Check and delete empty folder....
#############################################################################
for i in range(4):
    DOINGDIR = ( BASEDIR/ _astro_utilities.CCD_NEWUP_dir)         
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
                    print("{} is removed...".format(fpath)) 