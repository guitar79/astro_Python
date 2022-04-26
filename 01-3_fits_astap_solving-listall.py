# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

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
#base_dir = "../CCD_wcs_one/"
base_dir = "../CCD_obs_raw/"

fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))

fullnames_fit = [w for w in fullnames if ".fit" in w]

n = 0
for fullname in fullnames_fit[:] :
#fullname = fullnames[5]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames_fit))*100, os.path.basename(__file__)))
    print ("Starting...\nfullname: {}".format(fullname))

    fullname_el = fullname.split("/")
    filename_el = fullname_el[-1].split("_")
    
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
            
    else :
        print("{} is not light frame".format(fullname_el[-1]))


    '''
    if fullname[-4:].lower() == ".fit" :
        try : 

            if (not os.path.isfile(r'{0}.wcs'.format(fullname[:-4]))) or False :
                print("Trying image sovser using ASTAP...")

            else : 
                
                hdul = fits.open(fullname)
                print("fits file is opened".format(fullname_el[-1]))

                if "light" in hdul[0].header["IMAGETYP"].lower() :
                    print("{} is light frame".format(fullname_el[-1]))

                    if not os.path.isfile(r'{0}.wcs'.format(fullname[:-4])):
                        with subprocess.Popen(['astap', 
                                    '-f', 
                                    '{0}'.format(fullname), 
                                    '-update'
                                    '-analyse2'
                                    ],
                                    stdout=subprocess.PIPE) as proc :
                            print(proc.stdout.read())
                
                    else : 
                        print("{} is already solved...".format(fullname_el[-1]))
                        
                else :
                    print("{} is not light frame".format(fullname_el[-1]))
        
            #astro_utilities.ASTAPSolver(fullname)

        except Exception as err:
            Python_utilities.write_log(err_log_file,
                    '{2} ::: {0} with solve {1} '.format(err, fullname, datetime.now()))
        '''
#####################################################################
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
                    '{2} ::: {0} with solve {1} '.format(err, fullname, datetime.now()))
        