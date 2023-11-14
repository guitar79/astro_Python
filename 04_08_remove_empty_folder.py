"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

이 파일은 BASEDIR 폴더 안에 있는 모든 빈 디렉토리를 지워줍니다.

"""
#%%
import os
import _Python_utilities
import _astro_utilities

#%%
#######################################################
# for log file
#######################################################
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
#######################################################
BASEDIR = "../Rne_2022/"
BASEDIR = "/mnt/Rdata/NMSC/GK2A/AMI/PRIMARY/L1B/COMPLETE/EA/"

result = _Python_utilities.removeAllEmptyDirs(BASEDIR)


# %%
