"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from datetime import datetime
from astropy.io import fits
import os
import Python_utilities
import astro_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))

base_dir = "../CCD_new_files/"
#base_dir = "../CCD_wcs_One/"

fullnames = sorted(astro_utilities.getFullnameListOfallFiles(base_dir))

print ("fullnames: {}".format(fullnames))

#%%
focal_length = 910
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
            #with fits.open('{0}'.format(fullname), mode="update") as hdul :
                if not 'FOCALLEN' in hdul[0].header :
                    hdul[0].header.append('FOCALLEN', 
                                       '{0}'.format(focal_length), 
                                       'Focal length (mm)')
                    astro_utilities.write_log(log_file, 
                        '{1} ::: FOCALLEN is appended at {0}...'\
                        .format(fullname, datetime.now()))
                    
                hdul.flush()  # changes are written back to original.fits
                print('*'*30)
                astro_utilities.write_log(log_file, 
                    '{1} ::: fits header is append with {0} ...'\
                    .format(fullname, datetime.now()))
            
            # Change something in hdul.
            with fits.open('{0}'.format(fullname), mode="update") as hdul :
                #
                #old_object_name = hdul[0].header['OBJECT']
                if hdul[0].header['FOCALLEN'] != focal_length :
                    hdul[0].header['FOCALLEN'] = focal_length

                hdul.flush()  # changes are written back to original.fits
                print('*'*60)
                Python_utilities.write_log(log_file,
                    '{1} ::: fits header is update with {0} ...'\
                    .format(fullname, datetime.now()))
    
    except Exception as err :
        print("X"*60)
        Python_utilities.write_log(err_log_file,
            '{2} ::: \n{1} with {0} ...'\
            .format(fullname, err, datetime.now()))