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
BASEDIR = Path("/mnt/Rdata/CCD_obs") 
#BASEDIR = Path("/mnt/OBS_data") 
DOINGDIR = BASEDIR / _astro_utilities.CCD_NEW_dir
DOINGDIR = BASEDIR / _astro_utilities.CCD_obs_raw_dir

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
#print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################

#%%
#######################################################
checkKEYs = ["OBJECT", "TELESCOP", "OPTIC", "CCDNAME", 'FILTER',
            "GAIN", "EGAIN", "RDNOISE", "FOCALLEN", "FOCRATIO", "PIXSCALE", "CCD-TEMP",
            "XBINNING", "YBINNING", "FLIPSTAT", "EXPTIME", "EXPOSURE"]
#%%
for DOINGDIR in DOINGDIRs[:] :
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
            
            fpath = Path(fpath)

            foldername_el = fpath.parts[-2].split('_')
            fname_el = fpath.parts[-1].split('_')
            print("foldername_el", foldername_el)
            print("fname_el", fname_el)
            object_name = foldername_el[0]
            image_type = foldername_el[1]
            filter_name = fname_el[2]
            optic_name = foldername_el[5]
            ccd_name = foldername_el[6]
            print("object_name", object_name)
            print("filter_name", filter_name)
            print("optic_name", optic_name)
            print("ccd_name", ccd_name)

            with fits.open(str(fpath), mode="append") as hdul :
                for checkKEY in checkKEYs: 
                    if not checkKEY in hdul[0].header :
                        hdul[0].header.append(checkKEY, 
                                        '', 
                                        f"The keyword '{checkKEY}' is added.") 
                    print(f"{checkKEY}: ", hdul[0].header[checkKEY])

                hdul.flush()  # changes are written back to original.fits

            # Change something in hdul.
            with fits.open(str(fpath), mode="update") as hdul :
                
                #if object_name != hdul[0].header["OBJECT"] : 
                hdul[0].header["OBJECT"] = object_name.upper()
                print(f"The 'OBJECT' is set {object_name.upper()}")

                if 'INSTRUME' in hdul[0].header :
                    if 'qsi' in hdul[0].header['INSTRUME'].lower() :     
                        CCDNAME = 'QSI683ws'
                    elif 'st-8300' in hdul[0].header['INSTRUME'].lower() : 
                        CCDNAME = 'ST-8300M'
                    elif 'qhy8' in hdul[0].header['INSTRUME'].lower() : 
                        CCDNAME = 'QHY8'
                    elif 'stf-8300' in hdul[0].header['INSTRUME'].lower() : 
                        CCDNAME = 'STF-8300M'
                    elif '11000' in hdul[0].header['INSTRUME'] : 
                        CCDNAME = 'STL-11000M'
                    elif '16803' in hdul[0].header['INSTRUME'] : 
                        CCDNAME = 'STX-16803'
                    elif "SBIG" in hdul[0].header['INSTRUME'] :
                        if hdul[0].header['XPIXSZ'] == 5.4 \
                                or hdul[0].header['XPIXSZ'] == 10.8 :
                            CCDNAME = 'STF-8300M'
                        elif hdul[0].header['XPIXSZ'] == 9.0 \
                                or hdul[0].header['XPIXSZ'] == 18.0 :
                            if hdul[0].header['NAXIS1'] == 2048 \
                                    or  hdul[0].header['NAXIS1'] == 4096 :
                                CCDNAME = 'STX-16803'
                            elif hdul[0].header['NAXIS1'] == 4008 \
                                    or  hdul[0].header['NAXIS1'] == 2672 \
                                    or  hdul[0].header['NAXIS1'] == 2004 \
                                    or  hdul[0].header['NAXIS1'] == 1336 :
                                CCDNAME = 'STL-11000M'
                            else:
                                CCDNAME = ccd_name
                        else :
                            CCDNAME = ccd_name
                    else :
                        CCDNAME = ccd_name
                else :
                    CCDNAME = ccd_name    
                print("CCDNAME", CCDNAME)

                hdul[0].header["CCDNAME"] = CCDNAME
                print(f"The 'CCDNAME' is set {hdul[0].header['CCDNAME']}...")

                if len(hdul[0].header['DATE-OBS']) == 10 \
                    and 'TIME-OBS' in hdul[0].header : 
                    hdul[0].header['DATE-OBS'] += 'T' + hdul[0].header['TIME-OBS']
                    print(f"The 'DATE-OBS' is set {hdul[0].header['DATE-OBS']}")
                
                if not "IMAGETYP" in hdul[0].header :
                    hdul[0].header["IMAGETYP"] = image_type
                    
                if "ze" in hdul[0].header["IMAGETYP"].lower() \
                        or "bi" in hdul[0].header["IMAGETYP"].lower() :
                    hdul[0].header["IMAGETYP"] = "BIAS"
                    print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")
                elif "da" in hdul[0].header["IMAGETYP"].lower() :
                    hdul[0].header["IMAGETYP"] = "DARK"
                    print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")
                elif "fl" in hdul[0].header["IMAGETYP"].lower() :
                    hdul[0].header["IMAGETYP"] = "FLAT"
                    print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")
                elif "da" in hdul[0].header["IMAGETYP"].lower() \
                        or "lig" in hdul[0].header["IMAGETYP"].lower() :
                    hdul[0].header["IMAGETYP"] = "LIGHT"
                    print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")
                
                if "BIAS" in hdul[0].header["IMAGETYP"] \
                    or "DARK" in hdul[0].header["IMAGETYP"] :
                    for _KEY in ['FILTER', 'OPTIC', 'FOCALLEN', 'FOCRATIO', 'PIXSCALE'] :
                        hdul[0].header[_KEY] = "-"
                        print(f"The '{_KEY}' is set {hdul[0].header[_KEY]}")

                if "FLAT" in hdul[0].header["IMAGETYP"] \
                    or "LIGHT" in hdul[0].header["IMAGETYP"] :
                    if not "FILTER" in hdul[0].header :
                        if  hdul[0].header["FILTER"] != filter_name.upper() :
                            hdul[0].header["FILTER"] = filter_name.upper()
                            print(f"The timestamp = datetime.new().strftime(r'%Y-%m-%d %H:%M:%S')'FILTER' is set {hdul[0].header['FILTER']}")
                    if not "OPTIC" in hdul[0].header :
                        hdul[0].header["OPTIC"] = optic_name
                        print(f"The 'OPTIC' is set {hdul[0].header['OPTIC']}")
                    elif  hdul[0].header["OPTIC"] != optic_name :
                        hdul[0].header["OPTIC"] = optic_name
                        print(f"The 'OPTIC' is set {hdul[0].header['OPTIC']}")

                    #hdul[0].header['FOCALLEN'] = FOCALLENDIC[hdul[0].header['OPTIC']]
                    hdul[0].header['FOCALLEN'] = OPTICDIC[hdul[0].header['OPTIC']]["FOCALLEN"]
                    print(f"The 'FOCALLEN' is set {hdul[0].header['FOCALLEN']}...")
                    hdul[0].header['FOCRATIO'] = OPTICDIC[hdul[0].header['OPTIC']]["FOCALLEN"]/OPTICDIC[hdul[0].header['OPTIC']]["APATURE"]
                    print(f"The 'FOCRATIO' is set {hdul[0].header['FOCRATIO']}...")
                
                    print(hdul[0].header['OPTIC']+'_'+hdul[0].header['CCDNAME'])
                    # hdul[0].header['PIXSCALE'] = PIXSCALEDIC[hdul[0].header['OPTIC']\
                    #                                          +'_'+hdul[0].header['CCDNAME']]
                    #print(OPTICDIC['OPTIC']['FOCALLEN'], CCDDIC['CCDNAME']['PIXSIZE'])
                    hdul[0].header['PIXSCALE'] = calPixScale(OPTICDIC[hdul[0].header['OPTIC']]['FOCALLEN'],
                                                    CCDDIC[hdul[0].header['CCDNAME']]['PIXSIZE'])
                    print(f"The 'PIXSCALE' is set {hdul[0].header['PIXSCALE']}...")
                
                if (not 'TELESCOP' in hdul[0].header):
                    hdul[0].header['TELESCOP'] = "-"
                    print(f"The 'TELESCOP' is set {hdul[0].header['TELESCOP']}...")
                
                if (not 'XBINNING' in hdul[0].header)\
                    and (hdul[0].header["CCDNAME"] == "STX-16803") :
                    if hdul[0].header['NAXIS1'] == 4096 \
                        or  hdul[0].header['NAXIS2'] == 4096 :
                        hdul[0].header['XBINNING'] = 1
                        hdul[0].header['YBINNING'] = 1
                        hdul[0].header['TELESCOP'] = "-"   
                
                    elif hdul[0].header['NAXIS1'] == 2048 \
                        or  hdul[0].header['NAXIS2'] == 2048 :
                        hdul[0].header['XBINNING'] = 2
                        hdul[0].header['YBINNING'] = 2
                        hdul[0].header['TELESCOP'] = "-"
                
                    elif hdul[0].header['NAXIS1'] == 1024 \
                        or  hdul[0].header['NAXIS2'] == 1024 :
                        hdul[0].header['XBINNING'] = 3
                        hdul[0].header['YBINNING'] = 3
                        hdul[0].header['TELESCOP'] = "-"
                hdul[0].header['XBINNING'] = int(hdul[0].header['XBINNING'])
                hdul[0].header['YBINNING'] = int(hdul[0].header['YBINNING'])
                print(f"The 'XBINNING', 'YBINNING' are set {hdul[0].header['XBINNING']}, \
                        {hdul[0].header['YBINNING']},...")
                    
                if not "CCD-TEMP" in hdul[0].header :
                    hdul[0].header['CCD-TEMP'] = 'N'
                    print(f"The 'CCD-TEMP' is set {hdul[0].header['CCD-TEMP']}...")

                if "EXPOSURE" in hdul[0].header :
                    if not "EXPTIME" in hdul[0].header :
                        hdul[0].header["EXPTIME"] = hdul[0].header["EXPOSURE"]
                        print(f"The 'EXPTIME' is set {hdul[0].header['EXPOSURE']}...")
                elif "EXPTIME" in hdul[0].header :
                    hdul[0].header["EXPOSURE"] = hdul[0].header["EXPTIME"]
                    print(f"The 'EXPOSURE' is set {hdul[0].header['EXPTIME']}...")
                else :
                    hdul[0].header["EXPTIME"] = 'N'
                    hdul[0].header["EXPOSURE"] = 'N'
                    print(f"The 'EXPTIME' and 'EXPOSURE' are set 'N'...")

                #hdul[0].header['GAIN'] = GAINDIC[CCDNAME]
                hdul[0].header['GAIN'] = CCDDIC[hdul[0].header['CCDNAME']]['GAIN']
                print(f"The 'GAIN' is set {hdul[0].header['GAIN']}...")
                #hdul[0].header['EGAIN'] = GAINDIC[CCDNAME]
                hdul[0].header['EGAIN'] = CCDDIC[hdul[0].header['CCDNAME']]['GAIN']
                print(f"The 'EGAIN' is set {hdul[0].header['EGAIN']}...")
                #hdul[0].header['RDNOISE'] = RDNOISEDIC[CCDNAME]
                hdul[0].header['RDNOISE'] = CCDDIC[hdul[0].header['CCDNAME']]['RDNOISE']
                print(f"The 'RDNOISE' is set {hdul[0].header['RDNOISE']}...")
                
                hdul[0].header['FLIPSTAT'] = " "
                print(f"The 'FLIPSTAT' is set {hdul[0].header['FLIPSTAT']}...")
                
                for checkKEY in checkKEYs: 
                    print(f"{checkKEY}: ", hdul[0].header[checkKEY])

                hdul.flush()  # changes are written back to original.fits
                print('*'*30)
                print(f"The header of {fpath.name} is updated..")
            
    except:
        pass