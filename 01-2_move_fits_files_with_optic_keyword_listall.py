"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

이 파일은 base_dir 폴더 안에 있는 모든 fit 파일에 대해서 
fits header에 있는 정보를 바탕으로  
destination_base_dir_name 안에 규칙적으로 폴더를 만들어서 저장합니다. 

"""
#%%
import os
from datetime import datetime
from astropy.io import fits
import shutil 
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

destination_base_dir_name = "../CCD_obs_raw/"
target_duplicate_files_dir = "../CCD_duplicate_files/"


if not os.path.exists('{0}'.format(target_duplicate_files_dir)):
    os.makedirs('{0}'.format(target_duplicate_files_dir))

if not os.path.exists('{0}'.format(destination_base_dir_name)):
    os.makedirs('{0}'.format(destination_base_dir_name))
                
### make all file list...
fullnames = astro_utilities.getFullnameListOfallFiles(base_dir)
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))

#%%
n = 0   
for fullname in fullnames[:]:
    #fullname = fullnames[10]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames))*100, os.path.basename(__file__)))
    print ("Starting...   fullname: {}".format(fullname))

    try :
        if fullname[-4:].lower() in [".txt", "xisf", ".zip", ".png", ".log",
                                      "seal", "tiff", ".axy", "atch", "lved",
                                      "rdls", "xyls", "corr", "xosm", ".ini",
                                      ".wcs", ".csv"] \
                                and os.path.isfile('{}'.format(fullname)):
            os.remove("{}".format(fullname))
            print("{} is removed".format(fullname))

        elif fullname[-4:].lower() in [".fit", "fits", ".new", ".tmp"]\
                            and os.path.isfile('{}'.format(fullname)):

            fits.setval('{}'.format(fullname), \
                            'NOTES', value='modified by guitar79@naver.com')
            hdul = fits.open("{}".format(fullname))
            print("hdul[0].header.tostring: {}".format(hdul[0].header.tostring))
            fits_info1 = hdul[0].header.tostring()
            fits_info = fits_info1.replace("'", "'\'")
            print("fits_info: {}".format(fits_info))
            print("*"*60)
            
            Python_utilities.write_log(log_file, \
                        '{1} ::: {0} fits info modified ...'\
                        .format(fullname, datetime.now()))                
            
            new_filename = astro_utilities.get_new_filename(fullname)
            new_foldername = astro_utilities.get_new_foldername_from_filename(new_filename)
            print ("new_filename: {}".format(new_filename))
            new_foldername = "{}{}".format(destination_base_dir_name, new_foldername)
            print ("new_foldername: {}".format(new_foldername))
            
            if not os.path.exists('{0}'.format(new_foldername)):
                os.makedirs('{0}'.format(new_foldername))
                Python_utilities.write_log(log_file, \
                     '{1} ::: {0} is created'.format(new_foldername, datetime.now()))    
        
            if new_filename[-6:].lower() == "_-.fit" :
                if os.path.exists('{0}{1}_wcs.fit'.format(new_foldername, new_filename[:-6])):
                    Python_utilities.write_log(log_file, 
                         '{0}{1}_wcs.fit is already exist...'.format(new_foldername, new_filename))
                    os.rename(r"{}".format(fullname), 
                              r"{}{}".format(target_duplicate_files_dir, new_filename))
                    #shutil.move(r"{}".format(fullname), 
                    #            r"{}{}".format(target_duplicate_files_dir, new_filename))
                    print ("move {}".format(fullname), 
                           "{}{}".format(target_duplicate_files_dir, new_filename))
                else : 
                    os.rename(r'{0}'.format(fullname), 
                              r'{0}{1}'.format(new_foldername, new_filename))
                    #shutil.move(r'{0}'.format(fullname), 
                    #            r'{0}{1}'.format(new_foldername, new_filename))
                    Python_utilities.write_log(log_file, \
                             '{0} is moved to {1}{2}'.format(fullname, new_foldername, new_filename))
                    
            elif new_filename[-8:].lower() == "_wcs.fit" : 
                if os.path.exists('{0}{1}_-.fit'.format(new_foldername, new_filename[:-8])):
                    os.rename(r'{0}{1}_-.fit'.format(new_foldername, new_filename[:-8]), \
                                r"{0}{1}_-.fit".format(target_duplicate_files_dir, new_filename[:-8]))
                    #shutil.move(r'{0}{1}_-.fit'.format(new_foldername, new_filename[:-8]), \
                    #            r"{0}{1}_-.fit".format(target_duplicate_files_dir, new_filename[:-8]))
                if os.path.exists('{0}{1}_wcs.fit'.format(new_foldername, new_filename[:-8])):
                    os.rename(r'{0}{1}_wcs.fit'.format(new_foldername, new_filename[:-8]), \
                                r"{0}{1}_wcs.fit".format(target_duplicate_files_dir, new_filename[:-8]))
                    #shutil.move(r'{0}{1}_wcs.fit'.format(new_foldername, new_filename[:-8]), \
                    #            r"{0}{1}_wcs.fit".format(target_duplicate_files_dir, new_filename[:-8]))
                shutil.move(r'{0}'.format(fullname), r'{0}{1}'.format(new_foldername, new_filename))
                Python_utilities.write_log(log_file, \
                    '{0} is moved to {1}{2}'.format(fullname, new_foldername, new_filename))
            
            elif fullname[-4:].lower() == ".fit" \
                or fullname[-4:].lower() == "fits" : 
                os.rename(r"{}".format(fullname), r"{}{}".format(new_foldername, new_filename))
                #shutil.move(r"{}".format(fullname), r"{}{}".format(new_foldername, new_filename))
                Python_utilities.write_log(log_file, \
                    '{0} is moved to {1}{2}'.format(fullname, new_foldername, new_filename))
                                 
    except Exception as err :
        print("X"*60)
        Python_utilities.write_log(err_log_file, \
                '{2} ::: {0} with move {1} '.format(err, fullname, datetime.now()))

 #%%   
#############################################################################
#Check and delete empty folder....
#############################################################################

master_file_dir_name = 'master_file_Python/'

for i in range(4) : 
    fullnames = Python_utilities.getFullnameListOfallsubDirs(base_dir)
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
