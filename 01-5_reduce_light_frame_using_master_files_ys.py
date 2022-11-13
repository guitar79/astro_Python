# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user


#first time
cd ~/Downloads/ && git clone https://github.com/ysBach/ysvisutilpy && cd ysvisutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/ysfitsutilpy && cd ysfitsutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/ysphotutilpy && cd ysphotutilpy && git pull && pip install -e . && cd ..
cd ~/Downloads/ && git clone https://github.com/ysBach/SNUO1Mpy && cd SNUO1Mpy && git pull && pip install -e . && cd ..

# second time...
cd ~/Downloads/ysvisutilpy && git pull && pip install -e . 
cd ~/Downloads/ysfitsutilpy && git pull && pip install -e . 
cd ~/Downloads/ysphouutilpy && git pull && pip install -e . 
cd ~/Downloads/SNUO1Mpy && git pull && pip install -e . 

"""
#%%
from glob import glob
from pathlib import Path
import numpy as np
import os
import astropy.units as u
from ccdproc import CCDData, ccd_process
import Python_utilities

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

#######################################################
# for log file
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
base_dir = Path("../RnE_2022/KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin/")
base_dir = Path("../RnE_2022/KLEOPATRA_Light_-_2022-10-27_-_RiLA600_STX-16803_-_2bin/")

### make all fits file list...
fullnames = Python_utilities.getFullnameListOfallFiles("{}".format(base_dir))
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))

c_method = 'median'
master_dir = "master_files_ys/"
reduced_dir = "readuced_files/"

if not os.path.exists('{0}'.format("{}{}".format(base_dir, reduced_dir))):
    os.makedirs("{}{}".format(base_dir, reduced_dir))

# Open master bias, dark, and flat
#bias = CCDData.read('{}{}Master_Bias_{}_f32.fit'.\
#            format(base_dir, master_dir, c_method),
#            unit='adu')
#dark0 = CCDData.read('{}{}Master_Dark0_{}_f32.fit'.\
#            format(base_dir, master_dir, c_method),
#            unit='adu')
#dark = CCDData.read('{}{}Master_Dark_{}_f32.fit'.\
#            format(base_dir, master_dir, c_method),
#            unit='adu')

#%%
from pathlib import Path
summary = yfu.make_summary(base_dir/"*.fit")
df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
df_light = df_light.reset_index(drop=True)

# %%
for _, row in df_light.iterrows():
    fpath = Path(row["file"])
    ccd = yfu.load_ccd(fpath)
    filt = ccd.header["FILTER"]
    expt = ccd.header["EXPTIME"]
    red = yfu.ccdred(
        ccd,
        output=Path(f"{base_dir}/reduced/{fpath.stem}.fits"),
        mdarkpath=str(base_dir/"master_files_ys/master_dark_{:.0f}sec.fits".format(expt)),
        mflatpath=str(base_dir/"master_files_ys/master_flat_{}.fits".format(filt.upper())),
        # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
        overwrite=True
    )


# %%

# %%

# %%

# %%

#%%
n = 0
fullnames_light = [w for w in fullnames if not (("_bias_" in w.lower()) \
            or "_dark_" in w.lower() or "_flat_" in w.lower())]
print ("fullnames_light: {}".format(fullnames_light))
print ("len(fullnames_light): {}".format(len(fullnames_light)))

#%%


#%%
object_list = fullnames_light

for chl in['L', 'R', 'G', 'B', 'H', 'S', 'O', 'R', 'V'] :

    try : 
        #flat0 = CCDData.read('{}{}Master_Flat0_{}_{}_f32.fit'.\
        #    format(base_dir, master_dir, chl, c_method),
        #    unit='adu')
        flat = CCDData.read('{}{}master_flat_{}.fits'.\
            format(base_dir, master_dir, chl),
            unit='adu')

        object_list = [w for w in fullnames_light if "_{}_".format(chl) in w.upper()]

        # Reduce each object image separately.
        # Then save it with prefix 'p_' which means "preprocessed"
        for fullname in object_list:
            
            n += 1
            print('#'*40,
                "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames_light), (n/len(fullnames_light))*100, os.path.basename(__file__)))
            print ("Starting...\nfullname: {}".format(fullname))

            fullname_el = fullname.split("/")
            filename_el = fullname_el[-1].split("_")
            dark = CCDData.read('{}{}master_dark_{}.fits'.\
            format(base_dir, master_dir, filename_el[4]),
            unit='adu')
            
            obj = CCDData.read(fullname, unit='adu')
            reduced = ccd_process(obj,                    # The input image
                                #master_bias = bias,       # Master bias
                                dark_frame = dark,        # dark
                                master_flat = flat,      # non-calibrated flat
                                min_value = flat.data.min(),        # flat.min() should be 30000
                                exposure_key = 'EXPTIME', # exposure time keyword
                                #exposure_key = 'EXPOSURE', # exposure time keyword
                                exposure_unit = u.s,      # exposure time unit
                                dark_scale = True)        # whether to scale dark frame
            reduced.data = np.array(reduced.data, 
                    dtype=np.float32)
            print("reduced.data.shape:{}".format(reduced.data.shape))
            reduced.write("{}{}{}_reduced_f32.fit".\
                        format(base_dir, reduced_dir, fullname_el[-1][:-4]), overwrite =True)
            print("{}{}{}_reduced_f32.fit is created...".\
                        format(base_dir, reduced_dir, fullname_el[-1][:-4]))
            #break
    except Exception as err: 
        print ('Error messgae .......')
        print (err)
    
            
