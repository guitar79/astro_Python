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
import _Python_utilities
import pandas as pd
import _astro_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))

BASEDIR = "../A1_CCD_new_files/"
#BASEDIR = "../CCD_wcs_one/"
BASEDIR = "../A3_CCD_obs_raw/"

fullnames = _Python_utilities.getFullnameListOfallFiles(BASEDIR)
print ("fullnames: {}".format(fullnames))

fullnames_wcs = [w for w in fullnames if ".wcs" in w]

n = 0
for fullname in fullnames_wcs[:] :
#fullname = fullnames[5]
#fullname = fullnames_wcs[5]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames))*100, os.path.basename(__file__)))
    print ("Starting...\nfullname: {}".format(fullname))

    fullname_el = fullname.split("/")
    filename_el = fullname_el[-1].split("_")
    
    if not os.path.isfile(r'{0}.fit'.format(fullname[:-4])):
        print("There is no file:\n{0}.fit".format(fullname[:-4]))

    else : 
        f = open("{}".format(fullname), 'r')
        wcs_line = f.readline()

        #Index(['SIMPLE', '=', 'T', 'Unnamed: 3'], dtype='object')
        df_wcs = pd.read_fwf("{}".format(fullname), widths=(8,2,21,49))
        print("df_wcs:\n {}".format(df_wcs))

        hdul = fits.open('{0}.fit'.format(fullname[:-4]))
        print("fits file is opened".format(fullname_el[-1]))

        for idx, row in df_wcs.iterrows():
            print("index: {}".format(idx))
            print("df", df_wcs.loc[idx, 'SIMPLE'], df_wcs.loc[idx, "T"], type(df_wcs.loc[idx, "T"]))
            print("hdul", hdul[0].header[df_wcs.loc[idx, 'SIMPLE']], type(hdul[0].header[df_wcs.loc[idx, 'SIMPLE']]))
            print("value check", hdul[0].header[df_wcs.loc[idx, "SIMPLE"]] == df_wcs.loc[idx, "T"])
            print("==============")


        

    #_astro_utilities.ASTAPSolver(fullname)
'''
except Exception as err:
    _Python_utilities.write_log(err_log_file,
            '{2} ::: {0} with solve {1} '.format(err, fullname, datetime.now()))
            '''
