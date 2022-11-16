"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

이 파일은 fits 파일의 헤더에 누락된 정보를 넣어주는 것입니다.
base_dir 폴더 안에 있는 모든 fit 파일에 대해서 바로 상위 디렉토리 명을 참조하여 
OPTIC, CCDNAME 등의 정보를 줍니다.
"""
#%%
from datetime import datetime
from astropy.io import fits
import os
import Python_utilities
import astro_utilities

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
base_dir = "../CCD_new_files/"
#base_dir = "../Rne_2022/"

fullnames = astro_utilities.getFullnameListOfallFiles(base_dir)

print ("fullnames: {}".format(fullnames))
#######################################################


#######################################################
# set gain and readout noise 
gain = 0
rdnoise = 0
binning = 1
gains = {"STF-8300M": 0.37, 
        "STX-16803": 1.27, 
        "STL-11000": 0.8, 
        "QSI683ws": 0.13 } 

rdnoises = {"STF-8300M": 9.3, 
            "STX-16803": 9.0, 
            "STL-11000": 9.6, 
            "QSI683ws": 8.0 } 

pixscale = {"STF-8300M": None, 
            "RiLA600_STX-16803_1bin": 0.512, 
            "RiLA600_STX-16803_2bin": 1.024, 
            "STL-11000": None, 
            "QSI683ws": None } 
#######################################################

#%%

checkKEYs = ["OPTIC", "OBJECT", "FLIPSTAT", "CCDNAME", 
            "GAIN", "EGAIN", "RDNOISE",
            "XBINNING", "YBINNING"]

changeKEYs = ["OPTIC", "OBJECT", "FLIPSTAT", "CCDNAME", 
            "GAIN", "EGAIN", "RDNOISE",
            "XBINNING", "YBINNING"]

n = 0    
for fullname in fullnames[:] :
#fullname = fullnames[0]
    #######################################################
    # print the processing status

    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames))*100, os.path.basename(__file__)))
    print ("Starting...   fullname: {}".format(fullname))
    ######################################################


    if fullname[-4:] == ".fit" or fullname[-4:] == ".new" :
        print('Starting......\n{0} ...'.format(fullname))
        fullname_el = fullname.split('/')
        foldername_el = fullname_el[-2].split('_')
        object_name = foldername_el[0]
        optic_name = foldername_el[5]
        ccd_name = foldername_el[6]
        flipstat = "        "

        try: 
            with fits.open('{0}'.format(fullname), mode="append") as hdul :
            #with fits.open('{0}'.format(fullname), mode="update") as hdul :
                for checkKEY in checkKEYs: 
                    if not checkKEY in hdul[0].header :
                        hdul[0].header.append(checkKEY, 
                                        '', 
                                        '{} value added'.format(checkKEY))                    
                hdul.flush()  # changes are written back to original.fits
                print('*'*30)
                astro_utilities.write_log(log_file, 
                    '{1} ::: fits header is append with {0} ...'\
                    .format(fullname, datetime.now()))
        except Exception as err :
            print("X"*60)
            Python_utilities.write_log(err_log_file,
                '{2} ::: \n{1} with {0} ...'\
                .format(fullname, err, datetime.now()))


        try :            
            # Change something in hdul.
            with fits.open('{0}'.format(fullname), mode="update") as hdul :
                hdul[0].header['OBJECT'] = object_name
                hdul[0].header.append('COMMENT', 
                                       'change HEADER OBJECT {0}'.format(object_name), 
                                       'change HEADER OBJECT {0}'.format(object_name))

                hdul[0].header['OPTIC'] = optic_name
                hdul[0].header.append('COMMENT',
                                      'change HEADER OPTIC {0}'.format(optic_name),
                                      'change HEADER OPTIC {0}'.format(optic_name))

                hdul[0].header['CCDNAME'] = ccd_name
                hdul[0].header.append('COMMENT',
                                      'change HEADER CCDNAME {0}'.format(ccd_name),
                                      'change HEADER CCDNAME {0}'.format(ccd_name))

                if "-8300" in hdul[0].header['INSTRUME'] :
                    hdul[0].header['GAIN'] = gains["STF-8300M"]
                    hdul[0].header.append('COMMENT',
                                      'change HEADER GAIN {0}'.format(gains["STF-8300M"]),
                                      'change HEADER GAIN {0}'.format(gains["STF-8300M"]))
                    hdul[0].header['EGAIN'] = gains["STF-8300M"]
                    hdul[0].header.append('COMMENT',
                                      'change HEADER EGAIN {0}'.format(gains["STF-8300M"]),
                                      'change HEADER EGAIN {0}'.format(gains["STF-8300M"]))
                    hdul[0].header['RDNOISE'] = rdnoises["STF-8300M"]
                    hdul[0].header.append('COMMENT',
                                      'change HEADER RDNOISE {0}'.format(rdnoises["STF-8300M"]),
                                      'change HEADER RDNOISE {0}'.format(rdnoises["STF-8300M"]))
                    hdul[0].header['XBINNING'] = int(hdul[0].header['XPIXSZ'] / 5.4)
                    hdul[0].header.append('COMMENT',
                                      'change HEADER XBINNING {0}'.format(int(hdul[0].header['XPIXSZ'] / 5.4)),
                                      'change HEADER XBINNING {0}'.format(int(hdul[0].header['XPIXSZ'] / 5.4)))
                    
                    hdul[0].header['YBINNING'] = int(hdul[0].header['YPIXSZ'] / 5.4)
                    hdul[0].header.append('COMMENT',
                                      'change HEADER YBINNING {0}'.format(int(hdul[0].header['YPIXSZ'] / 5.4)),
                                      'change HEADER YBINNING {0}'.format(int(hdul[0].header['YPIXSZ'] / 5.4)))


                elif "16803" in hdul[0].header['INSTRUME'] :
                    hdul[0].header['GAIN'] = gains["STX-16803"]
                    hdul[0].header.append('COMMENT',
                                      'change HEADER GAIN {0}'.format(gains["STX-16803"]),
                                      'change HEADER GAIN {0}'.format(gains["STX-16803"]))
                    hdul[0].header['EGAIN'] = gains["STX-16803"]
                    hdul[0].header.append('COMMENT',
                                      'change HEADER EGAIN {0}'.format(gains["STX-16803"]),
                                      'change HEADER EGAIN {0}'.format(gains["STX-16803"]))
                    
                    hdul[0].header['RDNOISE'] = rdnoises["STX-16803"]
                    hdul[0].header.append('COMMENT',
                                      'change HEADER RDNOISE {0}'.format(rdnoises["STX-16803"]),
                                      'change HEADER RDNOISE {0}'.format(rdnoises["STX-16803"]))

                    if 'XPIXSZ' in hdul[0].header['INSTRUME'] :
                        binning = int(hdul[0].header['XPIXSZ'] / 9)
                        
                    elif hdul[0].header['NAXIS1'] == 4096 or hdul[0].header['NAXIS2'] == 4096 :
                        binning = 1
                    elif hdul[0].header['NAXIS1'] == 2048 or hdul[0].header['NAXIS2'] == 2048 :
                        binning = 2
                    elif hdul[0].header['NAXIS1'] == 1024 or hdul[0].header['NAXIS2'] == 1024 :
                        binning = 3    

                    hdul[0].header['XBINNING'] = binning
                    hdul[0].header['YBINNING'] = binning
                    hdul[0].header.append('COMMENT',
                                        'change HEADER XBINNING {0}'.format(binning),
                                        'change HEADER XBINNING {0}'.format(binning))
                    hdul[0].header.append('COMMENT',
                                        'change HEADER YBINNING {0}'.format(binning),
                                        'change HEADER YBINNING {0}'.format(binning))

                elif "11000" in hdul[0].header['INSTRUME'] :
                    hdul[0].header['GAIN'] = gains["STL-11000"]
                    hdul[0].header.append('COMMENT',
                                      'change HEADER GAIN {0}'.format(gains["STL-11000"]),
                                      'change HEADER GAIN {0}'.format(gains["STL-11000"]))
                    hdul[0].header['EGAIN'] = gains["STL-11000"]
                    hdul[0].header.append('COMMENT',
                                      'change HEADER EGAIN {0}'.format(gains["STL-11000"]),
                                      'change HEADER EGAIN {0}'.format(gains["STL-11000"]))
                    
                    hdul[0].header['RDNOISE'] = rdnoises["STX-16803"]
                    hdul[0].header.append('COMMENT',
                                      'change HEADER RDNOISE {0}'.format(rdnoises["STL-11000"]),
                                      'change HEADER RDNOISE {0}'.format(rdnoises["STL-11000"]))
                    hdul[0].header['XBINNING'] = int(hdul[0].header['XPIXSZ'] / 9)
                    hdul[0].header.append('COMMENT',
                                      'change HEADER XBINNING {0}'.format(int(hdul[0].header['XPIXSZ'] / 9)),
                                      'change HEADER XBINNING {0}'.format(int(hdul[0].header['XPIXSZ'] / 9)))
        
                    hdul[0].header['YBINNING'] = int(hdul[0].header['YPIXSZ'] / 9)
                    hdul[0].header.append('COMMENT',
                                      'change HEADER YBINNING {0}'.format(int(hdul[0].header['YPIXSZ'] / 9)),
                                      'change HEADER YBINNING {0}'.format(int(hdul[0].header['YPIXSZ'] / 9)))
                
                elif "683" in hdul[0].header['INSTRUME'] :
                    hdul[0].header['GAIN'] = gains["QSI683ws"]
                    hdul[0].header.append('COMMENT',
                                      'change HEADER GAIN {0}'.format(gains["QSI683ws"]),
                                      'change HEADER GAIN {0}'.format(gains["QSI683ws"]))
                    hdul[0].header['EGAIN'] = gains["QSI683ws"]
                    hdul[0].header.append('COMMENT',
                                      'change HEADER EGAIN {0}'.format(gains["QSI683ws"]),
                                      'change HEADER EGAIN {0}'.format(gains["QSI683ws"]))
                    
                    hdul[0].header['RDNOISE'] = rdnoises["STX-16803"]
                    hdul[0].header.append('COMMENT',
                                      'change HEADER RDNOISE {0}'.format(rdnoises["QSI683ws"]),
                                      'change HEADER RDNOISE {0}'.format(rdnoises["QSI683ws"]))
                    hdul[0].header['XBINNING'] = int(hdul[0].header['XPIXSZ'] / 5.4)
                    hdul[0].header.append('COMMENT',
                                      'change HEADER XBINNING {0}'.format(int(hdul[0].header['XPIXSZ'] / 5.4)),
                                      'change HEADER XBINNING {0}'.format(int(hdul[0].header['XPIXSZ'] / 5.4)))
                    
                    hdul[0].header['YBINNING'] = int(hdul[0].header['YPIXSZ'] / 5.4)
                    hdul[0].header.append('COMMENT',
                                      'change HEADER YBINNING {0}'.format(int(hdul[0].header['YPIXSZ'] / 5.4)),
                                      'change HEADER YBINNING {0}'.format(int(hdul[0].header['YPIXSZ'] / 5.4)))
                
                hdul.flush()  # changes are written back to original.fits
                print('*'*30)
                print('*'*35)
                astro_utilities.write_log(log_file, 
                    '{1} ::: fits header is append with {0} ...'\
                    .format(fullname, datetime.now()))
    
        except Exception as err :
            print("X"*60)
            Python_utilities.write_log(err_log_file,
                '{2} ::: \n{1} with {0} ...'\
                .format(fullname, err, datetime.now()))