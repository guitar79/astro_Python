# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from glob import glob
from pathlib import Path
from datetime import datetime, timedelta
import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.stats import sigma_clip
from ccdproc import combine, ccd_process, CCDData
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

import _astro_utilities
import _Python_utilities

plt.rcParams.update({'figure.max_open_warning': 0})

import warnings
warnings.filterwarnings('ignore')

#%%
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
count_stars = False
verbose = True
tryagain = False
trynightsky = False
tryASTROMETRYNET = True
file_age = 365
downsample = 4
#######################################################
BASEDIR = Path("/mnt/Rdata/ASTRO_data")  

PROJECDIR = BASEDIR / "C1-Variable"
TODODIR = PROJECDIR / "-_-_-_2016-_-_RiLA600_STX-16803_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2017-01_-_RiLA600_STX-16803_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2017-03_-_RiLA600_STX-16803_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2017-05_-_RiLA600_STX-16803_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2017-06_-_RiLA600_STX-16803_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2021-10_-_RiLA600_STX-16803_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2022-01_-_RiLA600_STX-16803_-_2bin"

# PROJECDIR = BASEDIR / "C2-Asteroid"
# TODODIR = PROJECDIR / "-_-_-_2022-_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2022-_-_RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2022-_-_RiLA600_STX-16803_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2023-_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2023-_-_RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2023-_-_RiLA600_STX-16803_-_2bin"

# PROJECDIR = BASEDIR / "C3-EXO"
# TODODIR = PROJECDIR / "-_-_-_2024-05_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2024-05_-_RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2024-06_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2024-06_-_RiLA600_STX-16803_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2024-09_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2024-09_-_RiLA600_ASI6200MMPro_-_2bin"
# TODODIR = PROJECDIR / "-_-_-_2024-11_-_GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "-_-_-_2024-11_-_RiLA600_ASI6200MMPro_-_3bin"

# PROJECDIR = BASEDIR / "C4-Spectra"
# TODODIR = PROJECDIR / "-_-_-_2024-05_TEC140_ASI183MMPro_-_1bin"

# PROJECDIR = BASEDIR / "C5-Test"
# TODODIR = PROJECDIR / "-_-_-_-_GSON300_STF-8300M_-_1bin"

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(TODODIR))
if verbose == True :
    print ("DOINGDIRs: ", format(DOINGDIRs))
    print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

try : 
    BDFDIR = [x for x in DOINGDIRs if "CAL-BDF" in str(x)]
    print ("BDFDIR: ", format(BDFDIR))
    BDFDIR = Path(BDFDIR[0])    
except : 
    BDFDIR = TODODIR
    pass

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_" in str(x)])
# print ("DOINGDIRs: ", format(DOINGDIRs))
# print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

# filter_str = 'TT-ARI'
# DOINGDIRs = [x for x in DOINGDIRs if filter_str in str(x)]
# remove = 'BIAS'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'DARK'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'FLAT'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
if verbose == True :
    print ("DOINGDIRs: ", DOINGDIRs)
    print ("len(DOINGDIRs): ", len(DOINGDIRs))
#######################################################
#%%
#####################################################################
# Observed location
LOCATION = dict(lon=127.005, lat=37.308889, elevation=101)
Suwon =  EarthLocation(lon=127.005 * u.deg, 
                                 lat=37.308889 * u.deg, 
                                 height=101 * u.m)
observatory_code = "P64"

# Used for any `astropy.SkyCoord` object:
SKYC_KW = dict(unit=u.deg, frame='icrs')
#######################################################
# Initial guess of FWHM in pixel
FWHM_INIT = 4

# Photometry parameters
# R_AP = 1.5*FWHM_INIT # Aperture radius
# R_IN = 4*FWHM_INIT   # Inner radius of annulus
# R_OUT = 6*FWHM_INIT  # Outer radius of annulus

Mag_Low = 11.5
Mag_High = 15

Mag_target = 12.5
Mag_delta = 2
ERR_Max = 0.5

coord_deltas = np.arange(0.00001, 0.00050, 0.00001)
#######################################################

#%%
###################################################################
# 각 프로젝트 'PROJECDIR' 아래에 'TODODIR' 이 위치하며 
# 전처리를 수행할 파일은 단 하나의 '*CAL-DBF*' 폴더에 위치한다.
### combine_BDF 함수를 이용하여 BIAS, DARK, FLAT 파일을 그룹별로 합성하여 
### master 파일을 만든다.
###################################################################
try :
    _astro_utilities.combine_BDF(BDFDIR,  # 전처리를 수행할 파일이 들어있는 경로이다.
                tryagain = tryagain,  #bool 형태로 입력한다. 이미 파일이 존재하는 경우 다시 시도할지를 결정한다.
                file_age = file_age, #int 형태로 day 단위로 입력한다. 이보다 오래된 파일은 덮어 쓴다.
                verbose = verbose,  
                )
except Exception as err :
    print("X"*60)
    print(str(err))
    pass

###################################################################
# 각 관측대상/관측일 별로 폴더에 나누어 저장한 폴더별로 측광을 수행한다.
###################################################################

for DOINGDIR in DOINGDIRs[:1] :

    # find '/mnt/Rdata/ASTRO_data/C3-EXO/-_-_-_2024-11_-_GSON300_STF-8300M_-_1bin/XO-6b_LIGHT_-_2025-01-09_-_GSON300_STF-8300M_-_1bin/' -type f -name '*.fit' -exec astap -f '{}' -wcs -analyse2 -update \;
    # find '/mnt/Rdata/ASTRO_data/C3-EXO/-_-_-_2024-11_-_GSON300_STF-8300M_-_1bin/Qatar-9b_LIGHT_-_2025-01-09_-_GSON300_STF-8300M_-_1bin/' -type f -name '*.fit' -exec solve-field -O --downsample 4 --nsigma 30 --cpulimit 20 '{}' \;
    # print(f"astap -f {str(fpath)} -fov {hfov} -wcs -analyse2 -update")    
    # os.system(f"astap -f {str(fpath)} -fov {hfov} -wcs -analyse2 -update")
     
    try : 
        _astro_utilities.solving_fits_file(DOINGDIR,
                downsample = downsample,
                count_stars = count_stars,
                tryagain = tryagain,
                tryASTAP = True, # default False 
                tryASTROMETRYNET = tryASTROMETRYNET,
                verbose = verbose,
                )    
    except Exception as err :
        print("X"*60)
        _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
        pass

    try :
        _astro_utilities.ccd_Reduction(DOINGDIR,
                BDFDIR,
                tryagain = False,
                trynightsky = trynightsky,
                file_age = file_age,
                verbose = verbose,
                )
    except Exception as err :
        print("X"*60)
        _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
        pass

    try :
        _astro_utilities.solving_fits_file(DOINGDIR,
                SOLVINGDIR = _astro_utilities.reduced_dir,
                downsample = downsample,
                tryagain = tryagain,
                tryASTAP = True, # default False
                tryASTROMETRYNET = tryASTROMETRYNET,
                verbose = verbose,     
                )     
    except Exception as err :
        print("X"*60)
        _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
        pass

    try :
        _astro_utilities.diff_Photometry_PS1(DOINGDIR,
                tryagain = tryagain, # True, 
                LOCATION = LOCATION,
                SKYC_KW = SKYC_KW,
                FWHM_INIT = FWHM_INIT,
                Mag_Low = Mag_Low,
                Mag_High = Mag_High,
                Mag_target = Mag_target,
                Mag_delta = Mag_delta,
                ERR_Max = ERR_Max,
                READINGDIR =  _astro_utilities.reduced_dir,
                # READINGDIR =  _astro_utilities.reduced_nightsky_dir,
                file_age = file_age,
                verbose = verbose,   
                ) 
    except Exception as err :
        print("X"*60)
        _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
        pass

    if PROJECDIR.parts[-1] == "C2-Asteroid" :
        try :
            _astro_utilities.plot_light_curve_asteroids_using_csv(DOINGDIR,
                Mag_target = Mag_target,
                FWHM_INIT = FWHM_INIT,
                READINGDIR = _astro_utilities.reduced_dir,
                #  READINGDIR = _astro_utilities.reduced_nightsky_dir,
                verbose = verbose, 
                )
        except Exception as err :
            print("X"*60)
            _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
            pass
    
    else :
        try :
            _astro_utilities.plot_light_curve_variables_using_csv(DOINGDIR,
                Mag_target = Mag_target,
                FWHM_INIT = FWHM_INIT,
                READINGDIR = _astro_utilities.reduced_dir,
                #  READINGDIR = _astro_utilities.reduced_nightsky_dir,
                verbose = verbose, 
                )
        except Exception as err :
            print("X"*60)
            _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
            pass

    if trynightsky == True : 
        try :
            _astro_utilities.solving_fits_file(DOINGDIR,
                SOLVINGDIR = _astro_utilities.reduced_nightsky_dir,
                downsample = downsample,
                tryagain = tryagain,
                tryASTAP = True, # default False
                tryASTROMETRYNET = tryASTROMETRYNET,  
                verbose = verbose,    
                ) 
        except Exception as err :
            print("X"*60)
            _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
            pass

        try :
            _astro_utilities.diff_Photometry_PS1 (DOINGDIR,
                tryagain = tryagain, # True, 
                LOCATION = LOCATION,
                SKYC_KW = SKYC_KW,
                FWHM_INIT = FWHM_INIT,
                Mag_Low = Mag_Low,
                Mag_High = Mag_High,
                Mag_target = Mag_target,
                Mag_delta = Mag_delta,
                ERR_Max = ERR_Max,
                # READINGDIR =  _astro_utilities.reduced_dir,
                READINGDIR =  _astro_utilities.reduced_nightsky_dir,
                file_age = file_age,
                verbose = verbose, 
                )
        except Exception as err :
            print("X"*60)
            _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
            pass

        if PROJECDIR.parts[-1] == "C2-Asteroid" :
            try :
                _astro_utilities.plot_light_curve_asteroids_using_csv(DOINGDIR,
                    Mag_target = Mag_target,
                    FWHM_INIT = FWHM_INIT,
                    # READINGDIR = _astro_utilities.reduced_dir,
                    READINGDIR = _astro_utilities.reduced_nightsky_dir,
                    verbose = verbose,
                    )
            except Exception as err :
                print("X"*60)
                _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
                pass

        else :
            try : 
                _astro_utilities.plot_light_curve_variables_using_csv(DOINGDIR,
                    Mag_target = Mag_target,
                    FWHM_INIT=FWHM_INIT,
                    #  READINGDIR = _astro_utilities.reduced_dir,
                    READINGDIR = _astro_utilities.reduced_nightsky_dir,
                    verbose = verbose,
                    )
            except Exception as err :
                print("X"*60)
                _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
                pass    
    
    _Python_utilities.write_log(err_log_file, f"{str(DOINGDIR)} is finighed..", verbose=verbose)