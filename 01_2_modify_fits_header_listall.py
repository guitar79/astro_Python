"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

conda env list
source activate astro_Python

conda install astropy

ModuleNotFoundError: No module named 'ccdproc'
conda install -c conda-forge ccdproc

"""

from datetime import datetime
from astropy.io import fits
import os
import astro_utilities

log_file = os.path.basename(__file__)[:-3]+".log"
err_log_file = os.path.basename(__file__)[:-3]+"_err.log"
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))


c_method = 'median'    

base_dir = "../CCD_new_files"
base_dir = "../../../3TB1/CCD_obs"

fullnames = astro_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))

n = 0    
for fullname in fullnames[:] :
#fullname = fullnames[0]
    n += 1
    if fullname[-4:] == ".fit" or fullname[-4:] == ".new" :
        print("\n{2:.01f}%  ({0}/{1})".format(n, len(fullnames), (n/len(fullnames))*100), '#'*40 )
        print('Starting......\n{0} ...'.format(fullname))
        fullname_el = fullname.split('/')
        foldername_el = fullname_el[-2].split('_')
        object_name = foldername_el[0]
        optic_name = foldername_el[5]

    
        try :
            #with fits.open('{0}'.format(fullname), mode="append") as hdul :
            with fits.open('{0}'.format(fullname), mode="update") as hdul :
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
                    astro_utilities.write_log(log_file, 
                        '{1} ::: OPTIC information is appended at {0}...'\
                        .format(fullname, datetime.now()))
            
                # Change something in hdul.
                if not hdul[0].header['OPTIC'] : 
                    hdul[0].header['OPTIC'] = '{0}'.format(optic_name)
                    astro_utilities.write_log(log_file, 
                        '{1} ::: OPTIC information is modified at {0}...'\
                        .format(fullname, datetime.now()))
                elif not '{0}'.format(optic_name).lower() in hdul[0].header['OPTIC'].lower() : 
                    hdul[0].header['OPTIC'] = '{0}'.format(optic_name)
                    astro_utilities.write_log(log_file, 
                        '{1} ::: OPTIC information is modified at {0}...'\
                        .format(fullname, datetime.now()))
                else : 
                    hdul[0].header['OPTIC'] = '{0}'.format(optic_name)
                    astro_utilities.write_log(log_file, 
                        '{1} ::: OPTIC information is modified at {0}...'\
                        .format(fullname, datetime.now()))
                hdul[0].header.append('COMMENT', 
                                       'add HEADER OPTIC {0}'.format(optic_name), 
                                       'add HEADER OPTIC {0}'.format(optic_name))
                
                # Change something in hdul.
                if not hdul[0].header['OBJECT'] : 
                    hdul[0].header['OBJECT'] = '{0}'.format(object_name)
                    astro_utilities.write_log(log_file, 
                        '{1} ::: OBJECT information is modified at {0}...'\
                        .format(fullname, datetime.now()))
                elif not '{0}'.format(object_name).lower() in hdul[0].header['OBJECT'].lower() : 
                    hdul[0].header['OBJECT'] = '{0}'.format(object_name)
                    astro_utilities.write_log(log_file, 
                        '{1} ::: OBJECT information is modified at {0}...'\
                        .format(fullname, datetime.now()))
                else : 
                    hdul[0].header['OBJECT'] = '{0}'.format(object_name)
                    astro_utilities.write_log(log_file, 
                        '{1} ::: OBJECT information is modified at {0}...'\
                        .format(fullname, datetime.now()))
                hdul[0].header.append('COMMENT', 
                                       'add HEADER OBJECT {0}'.format(object_name), 
                                       'add HEADER OBJECT {0}'.format(object_name))
                
                if not 'FLIPSTAT' in hdul[0].header :
                    print ("There is no 'FLIPSTAT' in the file")
                else : 
                    if hdul[0].header['FLIPSTAT'] != '        ' :
                        hdul[0].header.append('COMMENT', 
                                           'modified FLIPSTAT from {0}'.format(hdul[0].header['FLIPSTAT']), 
                                           'modified FLIPSTAT from {0}'.format(hdul[0].header['FLIPSTAT']))
                        hdul[0].header['FLIPSTAT'] = '        '
                        astro_utilities.write_log(log_file, 
                        '{1} ::: FLIPSTAT information is modified at {0}...'\
                        .format(fullname, datetime.now()))
            
                hdul.flush()  # changes are written back to original.fits
                print('#'*60)
                astro_utilities.write_log(log_file, 
                    '{1} ::: fits header is update with {0} ...'\
                    .format(fullname, datetime.now()))
    
        except Exception as err :
            print("X"*60)
            astro_utilities.write_log(err_log_file, 
                '{2} ::: \n{1} with {0} ...'\
                .format(fullname, err, datetime.now()))