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

from datetime import datetime
import os
import shutil 
import Python_utilities
import astro_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))

base_dir = "../CCD_new_files/"
#base_dir = "../../../3TB1/CCD_obs"
base_dir = "../CCD_obs_raw/QSI683ws_1bin/Light_FSQ106ED-x.73"
base_dir = "../2016/"

target_dir = "../CCD_wcs_one/"

fullnames = astro_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))
    
n = 0
for fullname in fullnames[:2000] :
#fullname = fullnames[0]
    n += 1
    print('#'*40,
           "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames))*100, os.path.basename(__file__)))
    try :
        if fullname[-6:].lower() == "_-.fit" :
            fullname_el = fullname.split("/")
            shutil.copy2(r"{}".format(fullname), r"{}{}".format(target_dir, fullname_el[-1]))
            print ("copy {}".format(fullname), "{}{}".format(target_dir, fullname_el[-1]))
    
    except Exception as err :
        print("X"*60)
        python_utilities.write_log(err_log_file, 
            '{2} ::: \n{1} with {0} ...'\
            .format(fullname, err, datetime.now()))