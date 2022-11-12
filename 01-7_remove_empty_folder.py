"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

이 파일은 base_dir 폴더 안에 있는 모든 빈 디렉토리를 지워줍니다.

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
base_dir = "../Rne_2022/"

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
