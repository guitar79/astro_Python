"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

conda env list
source activate astro_Python

conda install astropy

ModuleNotFoundError: No module named 'ccdproc'
conda install -c conda-forge ccdproc

2019.09.29  modify - missing 'IMAGETYP' on APT
"""

import os
from datetime import datetime
from astropy.io import fits
import shutil 
import Python_utilities
import astro_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))

master_file_dir_name = 'master_files_Python/'
processing_dir_name = 'processing_Python/'
integration_dir_name = 'integration_Python/'
alignment_dir_name = 'alignment_Python/'

base_dir = "../CCD_new_files/"
#base_dir = "../../../3TB1/CCD_obs"

wcs_one_dir_name = "../CCD_wcs_one/"
target_duplicate_files_dir = "../CCD_duplicate_files/"

if not os.path.exists('{0}'.format(target_duplicate_files_dir)):
    os.makedirs('{0}'.format(target_duplicate_files_dir))
if not os.path.exists('{0}'.format(wcs_one_dir_name)):
    os.makedirs('{0}'.format(wcs_one_dir_name))

                
fullnames = astro_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))
#fullname = fullnames[0]

n = 0   
for fullname in fullnames[:]:
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames))*100, os.path.basename(__file__)))
    print ("Starting...   fullname: {}".format(fullname))

    try :
        if fullname[-4:].lower() == ".txt" \
            or fullname[-4:].lower() == "xisf" \
            or fullname[-4:].lower() == ".zip" \
            or fullname[-4:].lower() == ".png" \
            or fullname[-4:].lower() == ".log" \
            or fullname[-4:].lower() == "seal" \
            or fullname[-4:].lower() == "tiff" \
            or fullname[-4:].lower() == "xosm" :
            os.remove("{}".format(fullname))
        
        elif (fullname[-4:].lower() == ".fit" or fullname[-4:].lower() == "fits" or fullname[-4:].lower() == ".new") \
            and (os.path.isfile('{}'.format(fullname))):

            fullname_el = fullname.split("/")
            new_filename = fullname_el[-1]
            new_foldername = wcs_one_dir_name
            print ("new_filename: {}".format(new_filename))
            print ("new_foldername: {}".format(new_foldername))
            
            if not os.path.exists('{0}'.format(new_foldername)):
                os.makedirs('{0}'.format(new_foldername))
                Python_utilities.write_log(log_file, \
                     '{1} ::: {0} is created'.format(new_foldername, datetime.now()))    
        
            if os.path.exists('{0}{1}'.format(new_foldername, new_filename)):
                Python_utilities.write_log(log_file, 
                     '{0}{1} is already exist...'.format(new_foldername, new_filename))
                shutil.move(r"{}".format(fullname), r"{}{}".format(target_duplicate_files_dir, new_filename))
                print ("move {}".format(fullname), "{}{}".format(target_duplicate_files_dir, new_filename))
                    
            else : 
                #os.rename(fullname, '{0}{1}'.format(new_foldername, new_filename))
                shutil.move(r"{}".format(fullname), r"{}{}".format(new_foldername, new_filename))
                Python_utilities.write_log(log_file, \
                         '{0} is moved to {1}{2}'.format(fullname, new_foldername, new_filename))
                fits.setval('{0}{1}'.format(new_foldername, new_filename), \
                        'NOTES', value='modified by guitar79@naver.com')
                #fits.setval('{0}{1}'.format(new_foldername, new_filename), \
                #        'observer', value='Kiehyun Park')
                
                hdul = fits.open("{0}{1}".format(new_foldername, new_filename))
                
                print("hdul[0].header.tostring: {}".format(hdul[0].header.tostring))
                fits_info1 = hdul[0].header.tostring()
                fits_info = fits_info1.replace("'", "'\'")
                print("fits_info: {}".format(fits_info))
                
                print("*"*60)
                Python_utilities.write_log(log_file, \
                     '{3} ::: {0} is moved to {1}{2} ...'\
                         .format(fullname, new_foldername, new_filename, datetime.now()))
    except Exception as err :
        print("X"*60)
        Python_utilities.write_log(err_log_file, \
                 '{2} ::: {0} with move {1} '.format(err, fullname, datetime.now()))
    
                
#############################################################################
#############################################################################
#############################################################################
import shutil 
for i in range(4) : 
    fullnames = astro_utilities.getFullnameListOfallsubDirs(base_dir)
    print ("fullnames: {}".format(fullnames))
    
    for fullname in fullnames[:] :
        fullname_el = fullname.split("/")
        if fullname_el[-1] == master_file_dir_name[:-1] : 
            #shutil.rmtree(r"{}".format(fullname))
            print ("rmtree {}\n".format(fullname))
    
        # Check is empty..
        if len(os.listdir(fullname)) == 0 :
            shutil.rmtree(r"{}".format(fullname)) # Delete..
            print ("rmtree {}\n".format(fullname))