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

base_dir = "../CCD_new_files/"
#base_dir = "../CCD_obs_raw/"
#base_dir = "../../../3TB1/CCD_obs"

fullnames = astro_utilities.getFullnameListOfallFiles(base_dir)
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
            object_name = foldername_el[0]
            optic_name = foldername_el[5]
            ccd_name = foldername_el[6]
            flipstat = "        "

            with fits.open('{0}'.format(fullname), mode="append") as hdul :
            #with fits.open('{0}'.format(fullname), mode="update") as hdul :
                if not 'OPTIC' in hdul[0].header :
                    hdul[0].header.append('OPTIC', 
                                       '{0}'.format(optic_name), 
                                       'OPTIC information')
                    astro_utilities.write_log(log_file, 
                        '{1} ::: OPTIC information is appended at {0}...'\
                        .format(fullname, datetime.now()))
                if not 'OBJECT' in hdul[0].header :
                    hdul[0].header.append('OBJECT', 
                                       '{0}'.format(object_name), 
                                       'OBJECT information')
                        
                if not 'FLIPSTAT' in hdul[0].header :
                    hdul[0].header.append('FLIPSTAT', 
                                       'FLIPSTAT {}'.format(flipstat), 
                                       'FLIPSTAT {}'.format(flipstat))
                    hdul[0].header.append('COMMENT', 
                                        'append FLIPSTAT {}'.format(flipstat), 
                                        'append FLIPSTAT {}'.format(flipstat))
                
                if not 'CCDNAME' in hdul[0].header :
                    hdul[0].header.append('CCDNAME', 
                                       'CCDNAME {}'.format(ccd_name), 
                                       'CCDNAME {}'.format(ccd_name))
                    hdul[0].header.append('CCDNAME', 
                                        'append CCDNAME {}'.format(ccd_name), 
                                        'append CCDNAME {}'.format(ccd_name))
                    
                hdul.flush()  # changes are written back to original.fits
                print('*'*30)
                astro_utilities.write_log(log_file, 
                    '{1} ::: fits header is append with {0} ...'\
                    .format(fullname, datetime.now()))
            
            # Change something in hdul.
            with fits.open('{0}'.format(fullname), mode="update") as hdul :
                #
                #old_object_name = hdul[0].header['OBJECT']
                hdul[0].header['OBJECT'] = object_name
                #hdul[0].header.update('OBJECT', 
                #                   '{0}'.format(object_name), 
                #                   'OBJECT information')
                
                hdul[0].header.append('COMMENT', 
                                       'change HEADER OBJECT {0}'.format(object_name), 
                                       'change HEADER OBJECT {0}'.format(object_name))

                # old_optic_name = hdul[0].header['OPTIC']
                hdul[0].header['OPTIC'] = optic_name
                # hdul[0].header.update('OPTIC',
                #                   '{0}'.format(optic_name),
                #                   'OPTIC information')

                hdul[0].header.append('COMMENT',
                                      'change HEADER OPTIC from {0}'.format(optic_name),
                                      'change HEADER OPTIC from {0}'.format(optic_name))

                # old_optic_name = hdul[0].header['INSTRUME']
                hdul[0].header['CCDNAME'] = ccd_name
                # hdul[0].header.update('INSTRUME',
                #                   '{0}'.format(instrument_name),
                #                   'INSTRUME information')

                hdul[0].header.append('COMMENT',
                                      'change HEADER CCDNAME from {0}'.format(ccd_name),
                                      'change HEADER CCDNAME from {0}'.format(ccd_name))

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