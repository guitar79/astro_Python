# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com
"""
#%%
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
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))
    
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

#%%
n = 0
for fullname in fullnames_fit[:] :
#fullname = fullnames[5]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames_fit), (n/len(fullnames_fit))*100, os.path.basename(__file__)))
    print ("Starting...\nfullname: {}".format(fullname))

    fullname_el = fullname.split("/")
    filename_el = fullname_el[-1].split("_")

    hdul = fits.open(fullname)
    print("fits file is opened...".format(fullname_el[-1]))

    # check Light frame
    if "light" in hdul[0].header["IMAGETYP"].lower() :
        print("{} is light frame...".format(fullname_el[-1]))
        
        # check HFD data in fits header
        if "HFD" in hdul[0].header :
        #if "CD1_1" in hdul[0].header :
            print("{0} is already solved by ASTAP...".format(fullname_el[-1]))

        else : 
            print("{0} is being solved by ASTAP...".format(fullname_el[-1]))
            with subprocess.Popen(['astap', 
                        '-f', 
                        '{0}'.format(fullname), 
                        '-analyse2',
                        '-update'],
                        stdout=subprocess.PIPE) as proc :
                print(proc.stdout.read())
            

            if os.path.exists("{}".format(fullname)):
                hdul = fits.open(fullname)
                print("fits file is opened...".format(fullname_el[-1]))
                if "CD1_1" in hdul[0].header :
                    print("{0} is already solved...".format(fullname_el[-1]))
                
                else : 
                    print("{0} is being solved by local Astrometry...".format(fullname_el[-1]))
                    try :
                        with subprocess.Popen(['solve-field', 
                                                '-O', #--overwrite: overwrite output files if they already exist
                                                #'--scale-units', 'arcsecperpix', #pixel scale
                                                #'--scale-low', '0.1', '--scale-high', '0.40', #pixel scale
                                                '-g', #--guess-scale: try to guess the image scale from the FITS headers
                                                '--cpulimit', '30',  #will make it give up after 30 seconds.
                                                #'-p', # --no-plots: don't create any plots of the results
                                                #'-D', '{0}'.format(save_dir_name), 
                                                '{0}'.format(fullname)], 
                                                stdout=subprocess.PIPE) as proc :
                            print(proc.stdout.read())
                    
                    except Exception as err :
                        Python_utilities.write_log(log_file,
                            '{1} ::: {2} with {0} ...'\
                            .format(fullname, datetime.now(), err))

    else :
        print("{} is not light frame...".format(fullname_el[-1]))

#%%
#############################################################################
#Check existence tmp file and rename ...
#############################################################################
fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))

fullnames_tmp = [w for w in fullnames if ".tmp" in w]

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
