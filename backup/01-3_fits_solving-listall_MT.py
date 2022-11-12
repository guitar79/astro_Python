# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

No module named 'ccdproc'
conda install -c conda-forge ccdproc
"""
#%%
import os
import numpy as np
import shutil
from datetime import datetime 
from concurrent.futures import ThreadPoolExecutor, as_completed
from astropy.io import fits
from time import time
import Python_utilities
import astro_utilities
import threading

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

#######################################################
# read all files in base directory for processing

base_dir = "../CCD_new_files/"
base_dir = "../CCD_obs_raw/STX-16803_2bin/"

destination_base_dir_name = "../CCD_obs_raw/"
target_duplicate_files_dir = "../CCD_duplicate_files/"

if not os.path.exists('{0}'.format(target_duplicate_files_dir)):
    os.makedirs('{0}'.format(target_duplicate_files_dir))

if not os.path.exists('{0}'.format(destination_base_dir_name)):
    os.makedirs('{0}'.format(destination_base_dir_name))

# make all fits file list...
fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))

fullnames_fit = [w for w in fullnames if (w.endswith(".fit") or w.endswith(".fits"))]
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames_fit): {}".format(len(fullnames_fit)))
#######################################################

#%%
#########################################
def threaded_process(items_chunk):
    """ Your main process which runs in thread for each chunk"""
    for fullname in items_chunk:                                               
        try:                                                                    
            astro_utilities.KevinSolver(fullname)
        except Exception:                                                       
            print('error with {}'.format(fullname))  


n_threads = 20
# Splitting the items into chunks equal to number of threads
array_chunk = np.array_split(fullnames_fit, n_threads)
print("len(array_chunk): {}".format(len(array_chunk)))
#print("array_chunk: {}".format(array_chunk))
print("len(array_chunk[0]): {}".format(len(array_chunk[0])))
print("array_chunk[0]: {}".format(array_chunk[0]))

#%%
thread_list = []
for thr in range(n_threads):
    for fullname in array_chunk[thr]:
        #thread = threading.Thread(target=threaded_process, args=(array_chunk[thr]),)
        thread = threading.Thread(target=astro_utilities.KevinSolver, args=(fullname,))
        thread.daemon = True
        thread_list.append(thread)
        thread_list[thr].start()
        #thread.start()

    for thread in thread_list:
        thread.join()

#%%
#############################################################################
#Check existence tmp file and rename ...
#############################################################################
fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))

fullnames_wcs = [w for w in fullnames if ((w.endswith(".tmp")) or (w.endswith(".new")))]

#print ("fullnames_wcs: {}".format(fullnames_wcs))
print ("len(fullnames_wcs): {}".format(len(fullnames_wcs)))

#%%
n = 0
for fullname in fullnames_wcs[:] :
#fullname = fullnames[5]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames_wcs), (n/len(fullnames_wcs))*100, os.path.basename(__file__)))
    print ("Starting...\nfullname: {}".format(fullname))

    try:
        if os.path.isfile('{}'.format(fullname)):
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
        
    except Exception as err:
        Python_utilities.write_log(err_log_file,
                    '{2} ::: {0} There is no {1} '.format(err, fullname, datetime.now()))    
