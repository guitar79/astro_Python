# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com
"""

from glob import glob
import os
import subprocess
from datetime import datetime
from astropy.io import fits
import shutil 
import Python_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))

base_dir = "../CCD_new_files/"
#base_dir = "../CCD_obs_raw/"
#base_dir = "../CCD_wcs_one/"
#base_dir = "../CCD_obs_raw/STL-11000M_2bin/"
#base_dir = "../CCD_obs_raw/ST-8300M_1bin/"

### make all fits file list...
fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
fullnames_fit = [w for w in fullnames if ".fit" in w]
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))


n = 0
for fullname in fullnames_fit[:] :
#fullname = fullnames[5]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames_fit), (n/len(fullnames_fit))*100, os.path.basename(__file__)))
    print ("Starting...\nfullname: {}".format(fullname))

    fullname_el = fullname.split("/")
    filename_el = fullname_el[-1].split("_")

    if os.path.exists('{0}.wcs'.format(fullname[-4:])) \
        or os.path.exists('{0}.ini'.format(fullname[-4:])):
        print("{0}.wcs is already exist...".format(fullname[-4:]))

    else :  
        hdul = fits.open(fullname)
        print("fits file is opened".format(fullname_el[-1]))

        if "light" in hdul[0].header["IMAGETYP"].lower() :
            print("{} is light frame".format(fullname_el[-1]))

            with subprocess.Popen(['astap', 
                        '-f', 
                        '{0}'.format(fullname), 
                        '-analyse2',
                        '-update'],
                        stdout=subprocess.PIPE) as proc :
                print(proc.stdout.read())
            
            #if os.path.exists("{}.tmp".format(fullname[:-4])):
            #    shutil.copy(r"{}.tmp".format(fullname[:-4]), \
            #                r"{}.fit".format(fullname[:-4]))
            #    print(r"{}.fit is created...".format(fullname[:-4]))

        else :
            print("{} is not light frame".format(fullname_el[-1]))

#############################################################################
#Check existence tmp file and rename ...
#############################################################################

fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))

fullnames_tmp = [w for w in fullnames if ".tmp" in w]

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