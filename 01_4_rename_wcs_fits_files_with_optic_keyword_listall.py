"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

conda env list
source activate astro_Python

conda install astropy

ModuleNotFoundError: No module named 'ccdproc'
conda install -c conda-forge ccdproc
"""

import os
from datetime import datetime
from astropy.io import fits
import shutil 
import astro_utilities

add_log = True
if add_log == True :
    log_file = 'astro_Python.log'
    err_log_file = 'astro_Python_err.log'

master_file_dir_name = 'master_file_Python/'
processing_dir_name = 'processing_Python/'
integration_dir_name = 'integration_Python/'
alignment_dir_name = 'alignment_Python/'

destination_base_dir_name = "../CCD_obs_raw/"
#solved_base_dir_name = "../CCD_obs_solved/"
target_duplicate_files_dir = "../CCD_duplicate_files/"
base_dir = "../CCD_new_files/"
base_dir = "../CCD_wcs_one/"

if not os.path.exists('{0}'.format(target_duplicate_files_dir)):
    os.makedirs('{0}'.format(target_duplicate_files_dir))
#if not os.path.exists('{0}'.format(solved_base_dir_name)):
    #os.makedirs('{0}'.format(solved_base_dir_name))
if not os.path.exists('{0}'.format(destination_base_dir_name)):
    os.makedirs('{0}'.format(destination_base_dir_name))
                
fullnames = astro_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))
#fullname = fullnames[10]
   
for fullname in fullnames[:]:
    if fullname[-4:].lower() == ".txt" \
        or fullname[-4:].lower() == "xisf" \
        or fullname[-4:].lower() == ".zip" \
        or fullname[-4:].lower() == ".png" \
        or fullname[-4:].lower() == ".log" \
        or fullname[-4:].lower() == "seal" \
        or fullname[-4:].lower() == "tiff" \
        or fullname[-4:].lower() == "xosm" :
        os.remove("{}".format(fullname))
    
    elif (fullname[-4:].lower() == ".fit" or fullname[-4:].lower() == "fits")\
          and os.path.isfile('{}'.format(fullname)):
        
        try :
            print ("Starting...   fullname: {}".format(fullname))
            fits.setval('{}'.format(fullname), \
                            'NOTES', value='modified by guitar79@naver.com')
            hdul = fits.open("{}".format(fullname))
            print("hdul[0].header.tostring: {}".format(hdul[0].header.tostring))
            fits_info1 = hdul[0].header.tostring()
            fits_info = fits_info1.replace("'", "'\'")
            print("fits_info: {}".format(fits_info))
            print("*"*60)
            astro_utilities.write_log(log_file, \
                '{1} ::: {0} fits info modified ...'\
                .format(fullname, datetime.now()))                
            new_filename = astro_utilities.get_new_filename(fullname)
            #new_foldername = astro_utilities.get_new_foldername(new_filename)
            #if new_filename[-4:] ==".fits" :
            new_foldername = astro_utilities.get_new_foldername_from_filename(new_filename)
            
            print ("new_filename: {}".format(new_filename))
            new_foldername = "{}{}".format(destination_base_dir_name, new_foldername)
            print ("new_foldername: {}".format(new_foldername))
            
            if not os.path.exists('{0}'.format(new_foldername)):
                os.makedirs('{0}'.format(new_foldername))
                astro_utilities.write_log(log_file, \
                     '{1} ::: {0} is created'.format(new_foldername, datetime.now()))    
        
            if new_filename[-6:].lower() == "_-.fit" :
                if os.path.exists('{0}{1}_wcs.fit'.format(new_foldername, new_filename[:-6])):
                    astro_utilities.write_log(log_file, 
                         '{0}{1}_wcs.fit is already exist...'.format(new_foldername, new_filename))
                    shutil.move(r"{}".format(fullname), r"{}{}".format(target_duplicate_files_dir, new_filename))
                    print ("move {}".format(fullname), "{}{}".format(target_duplicate_files_dir, new_filename))
                else : 
                    os.rename(fullname, '{0}{1}'.format(new_foldername, new_filename))
                    astro_utilities.write_log(log_file, \
                             '{0} is moved to {1}{2}'.format(fullname, new_foldername, new_filename))
                    
            elif new_filename[-8:].lower() == "_wcs.fit" : 
                if os.path.exists('{0}{1}_-.fit'.format(new_foldername, new_filename[:-8])):
                    shutil.move(r'{0}{1}_-.fit'.format(new_foldername, new_filename[:-8]), \
                                r"{0}{1}_-.fit".format(target_duplicate_files_dir, new_filename[:-8]))
                if os.path.exists('{0}{1}_wcs.fit'.format(new_foldername, new_filename[:-8])):
                    shutil.move(r'{0}{1}_wcs.fit'.format(new_foldername, new_filename[:-8]), \
                                r"{0}{1}_wcs.fit".format(target_duplicate_files_dir, new_filename[:-8]))

                os.rename(fullname, '{0}{1}'.format(new_foldername, new_filename))
                astro_utilities.write_log(log_file, \
                    '{0} is moved to {1}{2}'.format(fullname, new_foldername, new_filename))
            
            elif fullname[-4:].lower() == ".fit" \
                or fullname[-4:].lower() == "fits" : 
                os.rename(fullname, '{0}{1}'.format(new_foldername, new_filename))
                astro_utilities.write_log(log_file, \
                    '{0} is moved to {1}{2}'.format(fullname, new_foldername, new_filename))
                                 
        except Exception as err :
            print("X"*60)
            astro_utilities.write_log(err_log_file, \
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
