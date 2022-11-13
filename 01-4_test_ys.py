# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user


#first time
cd ~/Downloads/ && git clone https://github.com/ysBach/ysvisutilpy && cd ysvisutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/ysfitsutilpy && cd ysfitsutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/ysphotutilpy && cd ysphotutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/SNUO1Mpy && cd SNUO1Mpy && git pull && pip install -e . && cd ..

# second time...
cd ~/Downloads/ysvisutilpy && git pull && pip install -e . 
cd ~/Downloads/ysfitsutilpy && git pull && pip install -e . 
cd ~/Downloads/ysphouutilpy && git pull && pip install -e . 
cd ~/Downloads/SNUO1Mpy && git pull && pip install -e . 

"""
#%%
from glob import glob
import numpy as np
import os
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData
from astropy.io import fits
import Python_utilities
import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu


from pathlib import Path
from snuo1mpy import Preprocessor
import numpy as np

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
#%%
#######################################################
# read all files in base directory for processing
base_dir = "../Rne_2022/KLEOPATRA_Light_-_2022-10-27_-_RiLA600_STX-16803_-_2bin/"
#base_dir = "../RnE_2022/"

master_dir = "master_files/"

#%%
base_dirs = Python_utilities.getFullnameListOfsubDir(base_dir)
base_dirs = [w for w in base_dirs if not (w.endswith(master_dir) \
                or w.endswith(".fits"))]
print ("base_dirs: {}".format(base_dirs))

#%%
base_dir = "../Rne_2022/KLEOPATRA_Light_-_2022-10-27_-_RiLA600_STX-16803_-_2bin/"

fullnames = Python_utilities.getFullnameListOfallFiles("{}".format(base_dir))
fullnames = [w for w in fullnames \
            if ((w.endswith(".fit") or w.endswith(".fits")))
                and not ("bias" in w.lower()) or ("dark" in w.lower()) or ("flat" in w.lower())]
print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))

#%%
f_name1 = fullnames[10]
f_name2 = fullnames[48]
print(f_name1, "\n", f_name2)

#%%
hdul1 = fits.open("{}".format(f_name1))
hdul2 = fits.open("{}".format(f_name2))

print(hdul1[0].header.tostring())
print(hdul1[0].header.tostring())print(hdul2[0].header.tostring())
# %%
