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
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData
import Python_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))

base_dir = "../Post_process/KLEOPATRA1_Light_-_2022-10-17_-_GSON300_STF-8300M_-_1bin/"

### make all fits file list...
fullnames = Python_utilities.getFullnameListOfallFiles("{}".format(base_dir))
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))

c_method = 'median'
master_dir = "master_files/"

if not os.path.exists('{0}'.format("{}{}".format(base_dir, master_dir))):
    os.makedirs("{}{}".format(base_dir, master_dir))

#%%
fullnames_bias = [w for w in fullnames if "_bias_" in w.lower()]
print ("len(fullnames_bias): {}".format(len(fullnames_bias)))

if os.path.exists('{}{}Master_Bias_{}_f32.fit'.\
        format(base_dir, master_dir, c_method)) :

    # Open master bias
    bias = CCDData.read('{}{}Master_Bias_{}_f32.fit'.\
            format(base_dir, master_dir, c_method),
            unit='adu')
    
    print('{0}{1}Master_Dark_{2}_f32.fit is already exist..'.\
        format(base_dir, master_dir, c_method))
else:
    bias_list = fullnames_bias

    bias = combine(bias_list,       # ccdproc does not accept numpy.ndarray, but only python list.
                method=c_method,         # default is average so I specified median.
                unit='adu')              # unit is required: it's ADU in our case.

    bias.data = np.array(bias.data, 
                    dtype=np.float32)
                    
    print('type(bias.data)', type(bias.data),
        'bias.data.mean()', bias.data.mean(), 
        'bias.data.std()', bias.data.std(), 
        'bias.data.max()', bias.data.max(), 
        'bias.data.min()', bias.data.min(), 
        'bias', bias)

    # Save file
    bias.write('{}{}Master_Bias_{}_f32.fit'.\
                format(base_dir, master_dir, c_method), 
                overwrite =True)
    print('{}{}Master_Bias_{}_f32.fit is created'.\
                format(base_dir, master_dir, c_method))

#%%
fullnames_dark = [w for w in fullnames if "_dark_" in w.lower()]
print ("len(fullnames_dark): {}".format(len(fullnames_dark)))

if os.path.exists('{}{}Master_Dark_{}_f32.fit'.\
        format(base_dir, master_dir, c_method))\
        and os.path.exists('{}{}Master_Dark0_{}_f32.fit'.\
        format(base_dir, master_dir, c_method)) :

    # Open master bias, dark, and flat
    dark0 = CCDData.read('{}{}Master_Dark0_{}_f32.fit'.\
                format(base_dir, master_dir, c_method),
                unit='adu')
    dark = CCDData.read('{}{}Master_Dark_{}_f32.fit'.\
                format(base_dir, master_dir, c_method),
                unit='adu')

    print('{0}{1}Master_Dark_{2}_f32.fit and {0}{1}Master_Dark0_{2}_f32.fit are already exist..'.\
        format(base_dir, master_dir, c_method))
else:
    dark_list = fullnames_dark

    dark0 = combine(dark_list,       # ccdproc does not accept numpy.ndarray, but only python list.
                method=c_method,         # default is average so I specified median.
                unit='adu')
    dark0.data = np.array(dark0.data, 
                    dtype=np.float32)
                    
    print('type(dark0.data)', type(dark0.data),
        'dark0.data.mean()', dark0.data.mean(), 
        'dark0.data.std()', dark0.data.std(), 
        'dark0.data.max()', dark0.data.max(), 
        'dark0.data.min()', dark0.data.min(), 
        'dark0', dark0)
    dark0.write('{}{}Master_Dark0_{}_f32.fit'.\
            format(base_dir, master_dir, c_method), 
            overwrite =True)
    print('{}{}Master_Dark0_{}_f32.fit is created'.\
            format(base_dir, master_dir, c_method))

    dark = ccd_process(dark0, master_bias=bias) 
    dark.data = np.array(dark.data, 
                    dtype=np.float32)
                    
    print('type(dark.data)', type(dark.data),
        'dark.data.mean()', dark.data.mean(), 
        'dark.data.std()', dark.data.std(), 
        'dark.data.max()', dark.data.max(), 
        'dark.data.min()', dark.data.min(), 
        'dark', dark)

    # Save file
    dark.write('{}{}Master_Dark_{}_f32.fit'.\
            format(base_dir, master_dir, c_method), 
            overwrite =True)
    print('{}{}Master_Dark_{}_f32.fit is created'.\
            format(base_dir, master_dir, c_method))



#%%
fullnames_flat = [w for w in fullnames if "_flat_" in w.lower()]
print ("len(fullnames_flat): {}".format(len(fullnames_flat)))

for chl in['L', 'R', 'G', 'B', 'H', 'S', 'O'] :
    
    if os.path.exists('{}{}Master_Flat_{}_{}_f32.fit'.\
                    format(base_dir, master_dir, chl, c_method)) :
        print('{}{}Master_Flat_{}_{}_f32.fit is already exist...'.\
                    format(base_dir, master_dir, chl, c_method))
    else:

        flat_list = [w for w in fullnames_flat if "_{}_".format(chl) in w.upper()]
        
        try : 
            flat0 = combine(flat_list,       # ccdproc does not accept numpy.ndarray, but only python list.
            method=c_method,         # default is average so I specified median.
            unit='adu')
            flat0.data = np.array(flat0.data, 
                    dtype=np.float32)
                    
            print('type(flat0.data)', type(flat0.data), 
                'flat0.data.mean()', flat0.data.mean(), 
                'flat0.data.std()', flat0.data.std(), 
                'flat0.data.max()', flat0.data.max(), 
                'flat0.data.min()', flat0.data.min(), 
                'flat0', flat0)
            flat0.write('{}{}Master_Flat0_{}_{}_f32.fit'.\
                    format(base_dir, master_dir, chl, c_method), 
                    overwrite =True)
            print('{}{}Master_Flat0_{}_{}_f32.fit is created'.\
                    format(base_dir, master_dir, chl, c_method))
            # This dark isn't bias subtracted yet, so let's subtract bias:               
            # Open master bias and dark

            # Subtract bias and dark
            flat = ccd_process(flat0,                  # The input image (median combined flat)
                            master_bias = bias,       # Master bias
                            dark_frame = dark0,        # dark
                            exposure_key='exptime', # the header keyword for exposure time
                            exposure_unit=u.s,      # the unit of exposure time
                            dark_scale=True)        # whether to scale dark frame

            print('flat', flat.data.dtype, flat.data.min(), flat.data.max(), flat)

            flat.data = np.array(flat.data, 
                            dtype=np.float32)
                    
            # Save file
            flat.write('{}{}Master_Flat_{}_{}_f32.fit'.\
                        format(base_dir, master_dir, chl, c_method), 
                        overwrite =True)
            print('{}{}Master_Flat_{}_{}_f32.fit is created'.\
                        format(base_dir, master_dir, chl, c_method))
        except Exception as err: 
            print ('Error messgae .......')
            print (err)
