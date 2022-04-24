# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.stats import sigma_clip, sigma_clipped_stats
from photutils import DAOStarFinder
#from photutils import CircularAperture as CircAp
#from photutils import CircularAnnulus as CircAn
#from photutils import aperture_photometry as APPHOT

#base_dir = '../CCD_obs_raw/QSI683ws_1bin/Light_FSQ106ED-x73/46P-WIRTANEN_Light_-_2018-12-23_-_FSQ106ED-x73_QSI683ws_-_1bin/'
#filename = '46P-WIRTANEN_Light_L_2018-12-23-14-01-36_080sec_FSQ106ED-x73_QSI683ws_-30C_1bin_wcs.fit'

from glob import glob
import os
from astropy.io import fits
import Python_utilities
import astro_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))

base_dir = "../CCD_obs_raw/"

fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))

n = 0
for fullname in fullnames[:] :
#fullname = fullnames[5]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames))*100, os.path.basename(__file__)))
    print ("Starting...   fullname: {}".format(fullname))
    fullname_el = fullname.split("/")
    filename_el = fullname_el[-1].split("_")

    if fullname[-4:].lower() in [".fit"] and filename_el[1].lower() == "light" \
        and not os.path.exists("{}{}_DAOStarFinder.csv".format(fullname[:-len(fullname_el[-1])], fullname_el[-1][:-4])): 
        hdu = fits.open("{}".format(fullname))
        img = np.array(hdu[0].data/65535.0, dtype=np.float32)

        # if image value < 10^(-6), replace the pixel as 10^(-6)
        img[img < 1.e-6] = 1.e-6

        FWHM   = 1.5
        sky_th = 100    # sky_th * sky_sigma will be used for detection lower limit : sky_th = 5

        # What if we do "sigma-clip" than MAD?
        sky_a, sky_m, sky_s_sc = sigma_clipped_stats(img) # default is 3-sigma, 5 iters
        thresh = sky_th * sky_s_sc
        print('3 sigma 5 iters clipped case:')
        print('{0:8.6f} x {1:4d} = {2:8.6f}\n'.format(sky_s_sc, sky_th, thresh))
        try : 
            DAOfind = DAOStarFinder(fwhm=FWHM, threshold=thresh, 
                            #sharplo=0.2, sharphi=3.0,  # default values: sharplo=0.2, sharphi=1.0,
                            #roundlo=-1.0, roundhi=1.0,  # default values: roundlo=-1.0, roundhi=1.0,
                            #sigma_radius=1.5,          # default values: sigma_radius=1.5,
                            ratio=0.9,                 # 1.0: circular gaussian:  ratio=1.0,
                            exclude_border=True)       # To exclude sources near edges : exclude_border=True

            # The DAOStarFinder object ("DAOfind") gets at least one input: the image.
            # Then it returns the astropy table which contains the aperture photometry results:
            DAOfound = DAOfind(img)

            if len(DAOfound) == 0 :
                print ('No star was founded using DAOStarFinder\n'*3)
            else : 
            
                # Use the object "found" for aperture photometry:
                print('{} stars were founded'.format(len(DAOfound)))
                #print('DAOfound \n', DAOfound)
                DAOfound.pprint(max_width=1800)
                # save XY coordinates:
                DAOfound.write("{}{}_DAOStarFinder.csv".format(fullname[:-len(fullname_el[-1])], fullname_el[-1][:-4]), overwrite=True, format='ascii.fast_csv')
        except Exception as err :
            print("X"*60)
            Python_utilities.write_log(err_log_file, \
                '{0}::: {1} '.format(fullname, err))