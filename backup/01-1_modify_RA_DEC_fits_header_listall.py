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
#
object_RA = 5.2
object_DEC = 34.2
#OBJCTRA_OBJCTDEC

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

            with fits.open('{0}'.format(fullname), mode="append") as hdul :
                #hdul =  fits.open('{0}'.format(fullname))
                if not 'RA' in hdul[0].header :
                    hdul[0].header.append('RA', 
                                       '{}'.format(object_RA), 
                                       'RA')
                    print('{0} is added at RA...'.format(object_RA))
                    _Python_utilities.write_log(log_file, 
                        '{1} ::: RA is appended at {0}...'\
                        .format(fullname, datetime.now()))
                if not 'DEC' in hdul[0].header :
                    hdul[0].header.append('DEC', 
                                       '{}'.format(object_DEC), 
                                       'DEC')
                    print('{0} is added at DEC...'.format(object_DEC))
                    _Python_utilities.write_log(log_file, 
                        '{1} ::: DEC is appended at {0}...'\
                        .format(fullname, datetime.now()))
                    
                hdul.flush()  # changes are written back to original.fits
                print('*'*30)
                _Python_utilities.write_log(log_file, 
                    '{1} ::: fits header is appended with {0} ...'\
                    .format(fullname, datetime.now()))
            
            # Change something in hdul.
            with fits.open('{0}'.format(fullname), mode="update") as hdul :
                
                if 'RA' in hdul[0].header and 'DEC' in hdul[0].header :
                    hdul[0].header['RA'] = object_RA
                    hdul[0].header['DEC'] = object_DEC
                    print('{0} is added at RA...'.format(object_RA))
                    print('{0} is added at DEC...'.format(object_DEC))
                    _Python_utilities.write_log(log_file, 
                        '{1} ::: RA, DEC is appended at {0}...'\
                        .format(fullname, datetime.now()))
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