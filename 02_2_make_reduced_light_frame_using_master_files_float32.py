# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user
"""
#%%
from glob import glob
import numpy as np
import os
import astropy.units as u
from ccdproc import CCDData, ccd_process
import Python_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))

base_dr = "../Post_processing/M35_Light_-_2018-10-31_-_TMB130ss_STF-8300M_-_1bin/"
base_dr = "../RnE_2022/KLEOPATRA_Light_-_2022-10-12_-_GSON300_STF-8300M_-_1bin/"

### make all fits file list...
fullnames = Python_utilities.getFullnameListOfallFiles("{}/input".format(base_dr))
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))

c_method = 'median'
master_dr = "master_files/"
reduced_dr = "readuced_files/"

if not os.path.exists('{0}'.format("{}{}".format(base_dr, reduced_dr))):
    os.makedirs("{}{}".format(base_dr, reduced_dr))

# Open master bias, dark, and flat
bias = CCDData.read('{}{}Master_Bias_{}_f32.fit'.\
            format(base_dr, master_dr, c_method),
            unit='adu')
dark0 = CCDData.read('{}{}Master_Dark0_{}_f32.fit'.\
            format(base_dr, master_dr, c_method),
            unit='adu')
dark = CCDData.read('{}{}Master_Dark_{}_f32.fit'.\
            format(base_dr, master_dr, c_method),
            unit='adu')


#%%
n = 0
fullnames_light = [w for w in fullnames if ("_bias_" not in w.lower()) and ("_dark_" not in w.lower()) and ("_flat_" not in w.lower())]
print ("len(fullnames_light): {}".format(len(fullnames_light)))

object_list = fullnames_light

for chl in['L', 'R', 'G', 'B', 'H', 'S', 'O'] :
    

    flat0 = CCDData.read('{}{}Master_Flat0_{}_{}_f32.fit'.\
            format(base_dr, master_dr, chl, c_method),
            unit='adu')
    flat = CCDData.read('{}{}Master_Flat0_{}_{}_f32.fit'.\
            format(base_dr, master_dr, chl, c_method),
            unit='adu')

    object_list = [w for w in fullnames_light if "_{}_".format(chl) in w.upper()]

    try : 

        # Reduce each object image separately.
        # Then save it with prefix 'p_' which means "preprocessed"
        for fullname in object_list:
            
            n += 1
            print('#'*40,
                "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames_light), (n/len(fullnames_light))*100, os.path.basename(__file__)))
            print ("Starting...\nfullname: {}".format(fullname))

            fullname_el = fullname.split("/")
            obj = CCDData.read(fullname, unit='adu')
            reduced = ccd_process(obj,                    # The input image
                                master_bias = bias,       # Master bias
                                dark_frame = dark0,        # dark
                                master_flat = flat0,      # non-calibrated flat
                                min_value = flat0.data.min(),        # flat.min() should be 30000
                                exposure_key = 'exptime', # exposure time keyword
                                exposure_unit = u.s,      # exposure time unit
                                dark_scale = True)        # whether to scale dark frame
            reduced.data = np.array(reduced.data, 
                    dtype=np.float32)
            reduced.write("{}{}{}_reduced_f32.fit".\
                        format(base_dr, reduced_dr, fullname_el[-1][:-4]), overwrite =True)
            print("{}{}{}_reduced_f32.fit is created...".\
                        format(base_dr, reduced_dr, fullname_el[-1][:-4]))

    except Exception as err: 
        print ('Error messgae .......')
        print (err)
            
