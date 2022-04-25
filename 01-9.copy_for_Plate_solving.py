# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

No module named 'ccdproc'
conda install -c conda-forge ccdproc
"""
import os
import subprocess
import shutil
import Python_utilities
import astro_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))

base_dir = "../CCD_obs_raw/"
save_dir = "../CCD_wcs_one/"

if not os.path.exists(save_dir):
    os.makedirs(save_dir)
fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))

n = 0
for fullname in fullnames[:] :
#fullname = fullnames[5]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames))*100, os.path.basename(__file__)))
    print ("Starting...   fullname: {}".format(fullname))

    if fullname[-4:].lower() == ".fit" \
        and fullname[-7:].lower() != "wcs.fit":
        
        fullname_el = fullname.split("/")
        filename_el = fullname_el[-1].split("_")

        if filename_el[1].lower() == "light" :

            shutil.copy2(r"{}".format(fullname), r"{}{}".format(save_dir, fullname_el[-1]))
            print ("copy {}".format(fullname), "{}{}".format(save_dir, fullname_el[-1]))