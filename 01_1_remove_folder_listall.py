# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

"""
import os
import astro_filename_utility

add_log = True
if add_log == True :
    log_file = 'astro_Python.log'
    err_log_file = 'astro_Python_err.log'

master_file_dir_name = 'master_files_Python/'
processing_dir_name = 'processing_Python/'
integration_dir_name = 'integration_Python/'
alignment_dir_name = 'alignment_Python/'

c_method = 'median'    

add_log = True
if add_log == True :
    log_file = 'get_new_filename.log'
    err_log_file = 'get_new_filename_err.log'

base_dir = '../CCD_new_files'
base_dirs = ['../CCD_new_files']

fullnames = astro_filename_utility.getFullnameListOfallsubDirs(base_dir)
print ("fullnames: {}".format(fullnames))
    
import shutil 

for fullname in fullnames[:] :
    fullname_el = fullname.split("/")
    if fullname_el[-1] == master_file_dir_name[:-1] : 
        #shutil.rmtree(r"{}".format(fullname))
        print ("rmtree {}\n".format(fullname))

    # Check is empty..
    if len(os.listdir(fullname)) == 0 :
        shutil.rmtree(r"{}".format(fullname)) # Delete..
        print ("rmtree {}\n".format(fullname))