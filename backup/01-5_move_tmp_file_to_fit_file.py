# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

"""
import os
import shutil
from datetime import datetime
import Python_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))


BASEDIR = "../CCD_obs_raw/"

#####################################################################
fullnames = Python_utilities.getFullnameListOfallFiles(BASEDIR)
print ("fullnames: {}".format(fullnames))

fullnames_tmp = [w for w in fullnames if ".tmp" in w]

n = 0
for fullname in fullnames_tmp[:] :
#fullname = fullnames[5]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames_tmp), (n/len(fullnames_tmp))*100, os.path.basename(__file__)))
    print ("Starting...\nfullname: {}".format(fullname))

    try:
        shutil.move(r"{}".format(fullname), \
                        r"{}.fit".format(fullname[:-4]))

    except Exception as err:
        Python_utilities.write_log(err_log_file,
                    '{2} ::: {0} with solve {1} '.format(err, fullname, datetime.now()))