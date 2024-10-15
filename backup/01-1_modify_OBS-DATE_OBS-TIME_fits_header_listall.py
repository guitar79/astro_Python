"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from datetime import datetime
from astropy.io import fits
import os
import _Python_utilities
import _astro_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))

BASEDIR = "../CCD_new_files/"
#BASEDIR = "../CCD_wcs_One/"

fullnames = sorted(_astro_utilities.getFullnameListOfallFiles(BASEDIR))

print ("fullnames: {}".format(fullnames))

#%%
n = 0    
for fullname in fullnames[:] :
#fullname = fullnames[0]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames))*100, os.path.basename(__file__)))
    print ("Starting...   fullname: {}".format(fullname))

    try :
        if fullname[-4:] == ".fit" or fullname[-4:] == ".new" :
            print('Starting......\n{0} ...'.format(fullname))
            fullname_el = fullname.split('/')
            foldername_el = fullname_el[-2].split('_')
            
            # Change something in hdul.
            with fits.open('{0}'.format(fullname), mode="update") as hdul :
                #
                #old_object_name = hdul[0].header['OBJECT']
                if 'DATE-OBS' in hdul[0].header and 'TIME-OBS' in hdul[0].header : 
                    if len(hdul[0].header['DATE-OBS']) == 10 :
                        hdul[0].header['DATE-OBS'] += 'T{}'.format(hdul[0].header['TIME-OBS'])

                hdul.flush()  # changes are written back to original.fits
                print('*'*60)
                _Python_utilities.write_log(log_file,
                    '{1} ::: fits header is update with {0} ...'\
                    .format(fullname, datetime.now()))
    
    except Exception as err :
        print("X"*60)
        _Python_utilities.write_log(err_log_file,
            '{2} ::: \n{1} with {0} ...'\
            .format(fullname, err, datetime.now()))