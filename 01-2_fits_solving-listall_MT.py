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
from time import time
import Python_utilities
import astro_utilities
import threading

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))

#########################################
def threaded_process(fullnames):
    """ Your main process which runs in thread for each chunk"""
    for fullname in fullnames: 
        try:                                                                    
            astro_utilities.KevinSolver(fullname)
        except Exception:                                                       
            print('error with {}}'.format(fullname))  
#########################################

base_dir = "../CCD_new_files/"

### make all fits file list...
fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
fullnames_fit = [w for w in fullnames if ".fit" in w]
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames_fit): {}".format(len(fullnames_fit)))

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
#print("array_chunk[0]: {}".format(array_chunk[0]))

thread_list = []
for thr in range(n_threads):
    for fullname in array_chunk[thr]:
        thread = threading.Thread(target=astro_utilities.KevinSolver, args=(fullname))
        #thread = threading.Thread(target=threaded_process, args=(array_chunk[thr]),)
        thread.daemon = True
        #thread_list.append(thread)
        #thread_list[thr].start()
        thread.start()

#for thread in thread_list:
#    thread.join()

#%%
#############################################################################
#Check existence tmp file and rename ...
#############################################################################
fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
fullnames_tmp = [w for w in fullnames if ".tmp" in w]
#print ("fullnames_tmp: {}".format(fullnames_tmp))

#%%
n = 0
for fullname in fullnames_tmp[:] :
#fullname = fullnames[5]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames_tmp))*100, os.path.basename(__file__)))
    print ("Starting...\nfullname: {}".format(fullname))

    try:
        shutil.move(r"{}".format(fullname), \
                        r"{}.fit".format(fullname[:-4]))

    except Exception as err:
        Python_utilities.write_log(err_log_file,
                    '{2} ::: {0} There is no {1} '.format(err, fullname, datetime.now()))      
