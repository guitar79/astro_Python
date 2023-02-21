"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

이 파일은 fits 파일의 헤더에 누락된 정보를 넣어주는 것입니다.
BASEDIR 폴더 안에 있는 모든 fit 파일에 대해서 바로 상위 디렉토리 명을 참조하여 
OPTIC, CCDNAME 등의 정보를 줍니다.
"""
#%%
from glob import glob
from pathlib import Path, PosixPath, WindowsPath
import os
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
# read all files in base directory for processing
BASEDIR = Path(r"r:\CCD_obs") 
BASEDIR = Path( BASEDIR/ astro_utilities.CCD_NEW_dir)

DOINGDIRs = sorted(Python_utilities.getFullnameListOfallsubDirs(BASEDIR))
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################

#%%
#######################################################
checkKEYs = ["OBJECT", "TELESCOP", "OPTIC", "CCDNAME", 
            "GAIN", "EGAIN", "RDNOISE", "FOCALLEN", "PIXSCALE",
            "XBINNING", "YBINNING", "FLIPSTAT"]
#%%
for DOINGDIR in DOINGDIRs[:] :
    #fpath = Path(DOINGDIRs[0])
    DOINGDIR = Path(DOINGDIR)
    #print(f"Starting: {str(fpath.parts[-1])}")
    #save_fpath = fpath/f"summary_{fpath.parts[-1]}.csv"
    try: 
        summary = yfu.make_summary(DOINGDIR/"*.fit*",
                    #output = save_fpath,
                    verbose = False
                    )
        print("summary", summary)
        print("len(summary)", len(summary))

        for _, row in summary.iterrows():
            # 파일명 출력
            print (row["file"])
            fpath = Path(row["file"])

            foldername_el = fpath.parts[-2].split('_')
            print("foldername_el", foldername_el)
            object_name = foldername_el[0]
            optic_name = foldername_el[5]
            ccd_name = foldername_el[6]
            print("object_name", object_name)
            print("ccd_name", ccd_name)

            new_fpath = Path(f"{fpath.parents[0]}/{fpath.stem}_new.fit")

            #hdul = fits.open(str(fpath))

            try: 
                with fits.open(str(fpath), mode="append") as hdul :
                    for checkKEY in checkKEYs: 
                        if not checkKEY in hdul[0].header :
                            hdul[0].header.append(checkKEY, 
                                            '', 
                                            f"The keyword '{checkKEY}' is added") 
                            print(f"The keyword '{checkKEY}' is added...")
                    hdul.flush()  # changes are written back to original.fits
                
                for checkKEY in checkKEYs: 
                    print(f"{checkKEY}: ", hdul[0].header[checkKEY])

                # Change something in hdul.
                with fits.open(str(fpath), mode="update") as hdul :
                    try : 
                        if 'qsi' in hdul[0].header['INSTRUME'].lower() :     
                            CCDNAME = 'QSI683ws'
                        elif  'st-8300' in hdul[0].header['INSTRUME'].lower() : 
                            CCDNAME = 'ST-8300M'
                        elif  'stf-8300' in hdul[0].header['INSTRUME'].lower() : 
                            CCDNAME = 'STF-8300M'
                        elif  '11000' in hdul[0].header['INSTRUME'].lower() : 
                            CCDNAME = 'STL-11000M'
                        elif  '16803' in hdul[0].header['INSTRUME'].lower() : 
                            CCDNAME = 'STX-16803'
                    except :
                            CCDNAME = ccd_name
                    print("CCDNAME", CCDNAME)

                    if object_name != hdul[0].header["OBJECT"] : 
                        hdul[0].header["OBJECT"] = object_name
                        print(f"The 'OBJECT' is set {object_name}")

                    if not "CCDNAME" in hdul[0].header\
                        or hdul[0].header["CCDNAME"] != CCDNAME : 
                        hdul[0].header["CCDNAME"] = CCDNAME
                        print(f"The 'CCDNAME' is set {CCDNAME}")

                    if optic_name != hdul[0].header["TELESCOP"] :
                        hdul[0].header["OPTIC"] = optic_name
                        print(f"The 'OPTIC' is set {optic_name}")
                    else:
                        hdul[0].header["OPTIC"] = hdul[0].header["TELESCOP"] 
                        print(f"The 'OPTIC' is set {hdul[0].header['TELESCOP']}")
                        
                    hdul[0].header['GAIN'] = astro_utilities.GAINDIC[CCDNAME]
                    hdul[0].header['EGAIN'] = astro_utilities.GAINDIC[CCDNAME]
                    hdul[0].header['RDNOISE'] = astro_utilities.RDNOISEDIC[CCDNAME]
                    print(f"The 'GAIN' is set {astro_utilities.GAINDIC[CCDNAME]}...")
                    print(f"The 'RDNOISE' is set {astro_utilities.RDNOISEDIC[CCDNAME]}...")
                    hdul.flush()  # changes are written back to original.fits
                    print('*'*30)
                    print(f"The header of {fpath.name} is updated..")

            except Exception as err :
                print("X"*60)
                print(err)

    except:
        pass

        



    