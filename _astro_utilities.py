# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

"""
#%%
from pathlib import Path
from astropy.io import fits
import subprocess
from datetime import datetime, timedelta
import os
import numpy as np
import pandas as pd
from ccdproc import combine
import shutil
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip, sigma_clipped_stats
from photutils.centroids import centroid_com

from photutils.aperture import CircularAperture as CAp
from photutils.aperture import CircularAnnulus as CAn

import ysfitsutilpy as yfu

import _Python_utilities

#%%
#########################################
#directory variables
#########################################

c_method = "median"

CCD_obs_raw_dir = "CCD_obs_raw"
CCD_NEW_dir = "CCD_new_files"
CCD_NEWUP_dir = "CCD_newUpdated_files"
CCD_duplicate_dir = "CCD_duplicate_files"

master_dir = "master_files_ys"
reduced_dir = "reduced"
reduced_dir2 = "reduced2"
REDUC_nightsky_dir = "REDUC_nightsky"
solved_dir = "solved"
solved_dir2 = "solved2"
DAOfinder_result_dir = "DAOfinder_result"
IRAFfinder_result_dir = "IRAFfinder_result"
APh_result_dir = "APh_result"
Asteroid_result_dir = "Asteroid_result"
Asteroid_diff_Phot_dir = "Asteroid_diff_Phot"
Inst_Mag_dir = "Inst_Mag_result"
Diff_Phot_dir = "Diff_Phot_result"
Exoplanet_diff_Phot_dir = "Exoplanet_diff_Phot"

master_file_dir = 'master_file_Python/'
processing_dir = 'processing_Python/'
integration_dir = 'integration_Python/'
alignment_dir = 'alignment_Python/'

#######################################################
# OBS instruments information 
#######################################################

#CCDNAME, PIXSIZE, GAIN, RENOISE    
CCDDIC = {"ST-8300M": {"PIXSIZE":5.4, 
                        "GAIN":0.37,
                        "RDNOISE":9.3}, 
        "STF-8300M": {"PIXSIZE":5.4, 
                        "GAIN":0.37,
                        "RDNOISE":9.3}, 
        "QSI683ws": {"PIXSIZE":5.4, 
                        "GAIN":0.13,
                        "RDNOISE":8.0},
        "STL-11000M": {"PIXSIZE":9.0, 
                        "GAIN":0.8,
                        "RDNOISE":9.6},
        "STX-16803": {"PIXSIZE":9.0, 
                        "GAIN":1.27,
                        "RDNOISE":9.0},
        "QHY8": {"PIXSIZE":5.4, 
                "GAIN": "-",
                "RDNOISE":"-"},
        "ATR3CMOS26000KPA": {"PIXSIZE":3.76, 
                "GAIN": "-",
                "RDNOISE":"-"},
        "TT-2600CP": {"PIXSIZE":3.76, 
                "GAIN": "-",
                "RDNOISE":"-"},
        "ASI183MMPro": {"PIXSIZE":2.4, 
                "GAIN": "-",
                "RDNOISE":"-"},
        "ASI6200MMPro": {"PIXSIZE":3.76, 
                "GAIN": "-",
                "RDNOISE":"-"},
                        }

OPTICDIC = {"TMB130ss": {"APATURE" : 130, 
                         "FOCALLEN" : 910},
            "TMB130ss-x75": {"APATURE" : 130, 
                         "FOCALLEN" : 910*0.75}, 
            "RiLA600": {"APATURE" : 600, 
                        "FOCALLEN" : 3000}, 
            "RILA600": {"APATURE" : 600, 
                        "FOCALLEN" : 3000}, 
            "GSON300": {"APATURE" : 300, 
                        "FOCALLEN" : 1200}, 
            "OON300": {"APATURE" : 300, 
                        "FOCALLEN" : 1200}, 
            "FS60CB": {"APATURE" : 60, 
                       "FOCALLEN" : 355}, 
            "SVX80T": {"APATURE": 80,
                        "FOCALLEN": 480},
            "SVX80T-x80": {"APATURE":80,
                        "FOCALLEN": 480*0.8},
            "FSQ106ED": {"APATURE": 106,
                         "FOCALLEN": 530},
            "FSQ106ED-x73": {"APATURE": 106,
                         "FOCALLEN": 530*0.73},
            "FSQ106ED-x72": {"APATURE": 106,
                         "FOCALLEN": 530*0.72},
            "TEC140": {"APATURE": 140,
                       "FOCALLEN": 980},
            "TEC140-x75": {"APATURE": 140,
                       "FOCALLEN": 980*0.75},
            "TEC140-x72": {"APATURE": 140,
                       "FOCALLEN": 980*0.72},
                       }


#######################################################

#%%
#########################################
#calPixScale
#########################################
from warnings import warn
from astropy.visualization import (
    ImageNormalize,
    LinearStretch,
    ZScaleInterval,
    simple_norm,
)

__all__ = ["znorm", "zimshow", "norm_imshow"]

def znorm(image, **kwargs):
    return ImageNormalize(image, interval=ZScaleInterval(**kwargs))

def zimshow(
    ax,
    image,
    stretch=LinearStretch(),
    cmap=None,
    origin="lower",
    zscale_kw={},
    **kwargs
):
    im = ax.imshow(
        image,
        norm=znorm(image, stretch=stretch, **zscale_kw),
        origin=origin,
        cmap=cmap,
        **kwargs
    )
    return im

def norm_imshow(
    ax,
    data,
    origin="lower",
    stretch="linear",
    power=1.0,
    asinh_a=0.1,
    min_cut=None,
    max_cut=None,
    min_percent=None,
    max_percent=None,
    percent=None,
    clip=True,
    log_a=1000,
    invalid=-1.0,
    zscale=False,
    vmin=None,
    vmax=None,
    **kwargs
):
    """Do normalization and do imshow"""
    if vmin is not None and min_cut is not None:
        warn("vmin will override min_cut.")

    if vmax is not None and max_cut is not None:
        warn("vmax will override max_cut.")

    if zscale:
        zs = ImageNormalize(data, interval=ZScaleInterval())
        min_cut = vmin = zs.vmin
        max_cut = vmax = zs.vmax

    if vmin is not None or vmax is not None:
        im = ax.imshow(data, origin=origin, vmin=vmin, vmax=vmax, **kwargs)
    else:
        im = ax.imshow(
            data,
            origin=origin,
            norm=simple_norm(
                data=data,
                stretch=stretch,
                power=power,
                asinh_a=asinh_a,
                min_cut=min_cut,
                max_cut=max_cut,
                min_percent=min_percent,
                max_percent=max_percent,
                percent=percent,
                clip=clip,
                log_a=log_a,
                invalid=invalid
            ),
            **kwargs)
    return im

def sky_fit(all_sky, method='mode', sky_nsigma=3, sky_iter=5, \
            mode_option='sex', med_factor=2.5, mean_factor=1.5):
    '''
    Estimate sky from given sky values.
    Parameters
    ----------
    all_sky : ~numpy.ndarray
        The sky values as numpy ndarray format. It MUST be 1-d for proper use.
    method : {"mean", "median", "mode"}, optional
        The method to estimate sky value. You can give options to "mode"
        case; see mode_option.
        "mode" is analogous to Mode Estimator Background of photutils.
    sky_nsigma : float, optinal
        The input parameter for sky sigma clipping.
    sky_iter : float, optinal
        The input parameter for sky sigma clipping.
    mode_option : {"sex", "IRAF", "MMM"}, optional.
        sex  == (med_factor, mean_factor) = (2.5, 1.5)
        IRAF == (med_factor, mean_factor) = (3, 2)
        MMM  == (med_factor, mean_factor) = (3, 2)
    Returns
    -------
    sky : float
        The estimated sky value within the all_sky data, after sigma clipping.
    std : float
        The sample standard deviation of sky value within the all_sky data,
        after sigma clipping.
    nsky : int
        The number of pixels which were used for sky estimation after the
        sigma clipping.
    nrej : int
        The number of pixels which are rejected after sigma clipping.
    -------
    '''
    sky = all_sky.copy()
    if method == 'mean':
        return np.mean(sky), np.std(sky, ddof=1)

    elif method == 'median':
        return np.median(sky), np.std(sky, ddof=1)

    elif method == 'mode':
        sky_clip   = sigma_clip(sky, sigma=sky_nsigma,
                                maxiters=sky_iter, #iters=sky_iter,
                                )
        sky_clipped= sky[np.invert(sky_clip.mask)]
        nsky       = np.count_nonzero(sky_clipped)
        mean       = np.mean(sky_clipped)
        med        = np.median(sky_clipped)
        std        = np.std(sky_clipped, ddof=1)
        nrej       = len(all_sky) - len(sky_clipped)

        if nrej < 0:
            raise ValueError('nrej < 0: check the code')

        if nrej > nsky: # rejected > survived
            raise Warning('More than half of the pixels rejected.')

        if mode_option == 'IRAF':
            if (mean < med):
                sky = mean
            else:
                sky = 3 * med - 2 * mean

        elif mode_option == 'MMM':
            sky = 3 * med - 2 * mean

        elif mode_option == 'sex':
            if (mean - med) / std > 0.3:
                sky = med
            else:
                sky = (2.5 * med) - (1.5 * mean)
        else:
            raise ValueError('mode_option not understood')

        return sky, std, nsky, nrej
#%%
#########################################
#calPixScale
#########################################

def calPixScale (
    F_length,
    Pix_size,
    binn) :
    '''
        Parameters
        ----------
        F_length : float or int
            Focal Length of Telescope with out accesery (mm)
        
        Pix_Size : float
            pixel size of detector (um), 
        
        binn : int
            binning number, 
    
        Pixel scale : Pix_Size  /   Telescope Focal Length   )   X 206.265  
            (arcsec / pixel)        
    '''

    PIXScale = Pix_size * binn / (F_length ) *  206.265
    return PIXScale

#%%
#########################################
#KvinFitsMover
#########################################
def KevinFitsNewFname(
    fpath,
    fnameKEYs = ["OBJECT", "IMAGETYP", "FILTER", "DATE-OBS", 
            "EXPOSURE", "OPTIC", "CCDNAME", "CCD-TEMP", "XBINNING"],
    ):
    '''
        Parameters
        ----------
        fpath : string
            The fullname of input file...
        
        fnameKEYs : list
            KEY of fits file header for update
    '''
    
    fpath = Path(fpath)
    hdul = fits.open(str(fpath))
    print("fpath: ", fpath)
    for fnameKEY in fnameKEYs: 
        print(f"{fnameKEY}: ", hdul[0].header[fnameKEY])
        try :
            ccdtemp = str(int(hdul[0].header["CCD-TEMP"]))
        except : 
            ccdtemp = "N"
        print("ccdtemp: ", ccdtemp)

        new_fname = hdul[0].header["OBJECT"]+"_"+hdul[0].header["IMAGETYP"]+"_"+hdul[0].header["FILTER"]+"_"
        new_fname += hdul[0].header["DATE-OBS"][:19].replace("T","-").replace(":","-")+"_"
        new_fname += str(int(hdul[0].header["EXPOSURE"]))+"sec_"
        new_fname += hdul[0].header["OPTIC"]+"_"+hdul[0].header["CCDNAME"]+"_"       
        new_fname += ccdtemp+"c_"+str(int(hdul[0].header["XBINNING"]))+"bin.fit"
        #new_fname += fpath.ext
        print("new_fname: ", new_fname)
        hdul.close()
    return new_fname


#%%
#########################################
#KvinFitsUpdater
#########################################
def KevinFitsUpdater(
    fpath,
    checkKEYs = ["OBJECT", "TELESCOP", "OPTIC", "CCDNAME", 'FILTER',
            #"GAIN", "EGAIN", "RDNOISE", 
            "PIXSCALE", "FOCALLEN", "APATURE", "CCD-TEMP",
            'XPIXSZ', 'YPIXSZ',
            "XBINNING", "YBINNING", "FLIPSTAT", "EXPTIME", "EXPOSURE"],
    imgtype_update=False,
    fil_update=True,
    ):
    '''
        Parameters
        ----------
        fpath : string
            The fullname of input file...
        checkKEYs : dictionary
            KEY of fits file header for update
    '''
    
    fpath = Path(fpath)

    foldername_el = fpath.parts[-2].split('_')
    fname_el = fpath.parts[-1].split('_')
    print("foldername_el", foldername_el)
    print("fname_el", fname_el)
    object_name = foldername_el[0].replace(" ","")
    print("object_name", object_name)
    image_type = foldername_el[1]
    #filter_name = fname_el[2]
    #print("filter_name", filter_name)
    optic_name = foldername_el[5]
    print("optic_name", optic_name)
    ccd_name = foldername_el[6]
    print("ccd_name", ccd_name)
    with fits.open(str(fpath), mode="append") as hdul :
        for checkKEY in checkKEYs: 
            if not checkKEY in hdul[0].header :
                hdul[0].header.append(checkKEY, 
                                '', 
                                f"The keyword '{checkKEY}' is added.") 
            print(f"{checkKEY}: ", hdul[0].header[checkKEY])

        hdul.flush()  # changes are written back to original.fits

    # Change something in hdul.
    with fits.open(str(fpath), mode="update") as hdul :
        
        ###########################
        #### "OBJECT"
        hdul[0].header["OBJECT"] = object_name   #delete upper()
        print(f"The 'OBJECT' is set {object_name}")

        ###########################
        #### "FILTER"
        # hdul[0].header["FILTER"] = object_name.upper()
        # print(f"The 'OBJECT' is set {object_name.upper()}")

        ###########################
        #### "CCDNAME"
        if 'INSTRUME' in hdul[0].header :
            if 'qsi' in hdul[0].header['INSTRUME'].lower() :     
                CCDNAME = 'QSI683ws'
            elif 'st-8300' in hdul[0].header['INSTRUME'].lower() : 
                CCDNAME = 'ST-8300M'
            elif 'qhy8' in hdul[0].header['INSTRUME'].lower() : 
                CCDNAME = 'QHY8'
            elif 'stf-8300' in hdul[0].header['INSTRUME'].lower() : 
                CCDNAME = 'STF-8300M'
            elif '11000' in hdul[0].header['INSTRUME'] : 
                CCDNAME = 'STL-11000M'
            elif '16803' in hdul[0].header['INSTRUME'] : 
                CCDNAME = 'STX-16803'
            elif 'ASI183MM Pro' in hdul[0].header['INSTRUME'] : 
                CCDNAME = 'ASI183MMPro'
            elif 'ASI6200MM Pro' in hdul[0].header['INSTRUME'] : 
                CCDNAME = 'ASI6200MMPro'
            elif "SBIG" in hdul[0].header['INSTRUME'] :
                if hdul[0].header['XPIXSZ'] == 5.4 \
                        or hdul[0].header['XPIXSZ'] == 10.8 :
                    CCDNAME = 'STF-8300M'
                elif hdul[0].header['XPIXSZ'] == 9.0 \
                        or hdul[0].header['XPIXSZ'] == 18.0 :
                    if hdul[0].header['NAXIS1'] == 2048 \
                            or  hdul[0].header['NAXIS1'] == 4096 :
                        CCDNAME = 'STX-16803'
                    elif hdul[0].header['NAXIS1'] == 4008 \
                            or  hdul[0].header['NAXIS1'] == 2672 \
                            or  hdul[0].header['NAXIS1'] == 2004 \
                            or  hdul[0].header['NAXIS1'] == 1336 :
                        CCDNAME = 'STL-11000M'
                    else:
                        CCDNAME = ccd_name
            elif "CMOS26000" in hdul[0].header['INSTRUME'] or \
                "ToupTek" in hdul[0].header['INSTRUME'] :
                CCDNAME = "TT-2600CP"
                hdul[0].header["FILTER"] = "-"

            else :
                #CDNAME = hdul[0].header['INSTRUME']
                CCDNAME = ccd_name
        else :
            CCDNAME = ccd_name
        print("CCDNAME", CCDNAME)

        hdul[0].header["CCDNAME"] = CCDNAME
        print(f"The 'CCDNAME' is set {hdul[0].header['CCDNAME']}...")

        ###########################
        #### 'DATE-OBS'
        if len(hdul[0].header['DATE-OBS']) == 10 \
            and 'TIME-OBS' in hdul[0].header : 
            hdul[0].header['DATE-OBS'] += 'T' + hdul[0].header['TIME-OBS']
            print(f"The 'DATE-OBS' is set {hdul[0].header['DATE-OBS']}")
        
        ###########################
        #### 'IMAGETYP'
        if imgtype_update == True :
            hdul[0].header["IMAGETYP"] = image_type.upper()
            print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")
        if not "IMAGETYP" in hdul[0].header :
            hdul[0].header["IMAGETYP"] = image_type  
        elif "ze" in hdul[0].header["IMAGETYP"].lower() \
                or "bi" in hdul[0].header["IMAGETYP"].lower() :
            hdul[0].header["IMAGETYP"] = "BIAS"
            print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")
        elif "da" in hdul[0].header["IMAGETYP"].lower() :
            hdul[0].header["IMAGETYP"] = "DARK"
            print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")
        elif "fl" in hdul[0].header["IMAGETYP"].lower() :
            hdul[0].header["IMAGETYP"] = "FLAT"
            print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")
        elif "obj" in hdul[0].header["IMAGETYP"].lower() \
                or "lig" in hdul[0].header["IMAGETYP"].lower() :
            hdul[0].header["IMAGETYP"] = "LIGHT"
            print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")

        if "BIAS" in hdul[0].header["IMAGETYP"] \
            or "DARK" in hdul[0].header["IMAGETYP"] :
            for _KEY in ['FILTER', 'OPTIC', 'FOCALLEN', 'APATURE', 'PIXSCALE',] :
                hdul[0].header[_KEY] = "-"
                print(f"The '{_KEY}' is set {hdul[0].header[_KEY]}")

        if "FLAT" in hdul[0].header["IMAGETYP"] \
            or "LIGHT" in hdul[0].header["IMAGETYP"] :
            filter_name = fname_el[2]
            
            if not "FILTER" in hdul[0].header :
                hdul[0].header["FILTER"] = filter_name   #delete uper()
            if hdul[0].header["FILTER"] != filter_name \
                and fil_update==True :
                hdul[0].header["FILTER"] = filter_name
            print(f"FILTER is set {hdul[0].header['FILTER']}")
            if not "OPTIC" in hdul[0].header :
                hdul[0].header["OPTIC"] = optic_name
                print(f"The 'OPTIC' is set {hdul[0].header['OPTIC']}")
            elif  hdul[0].header["OPTIC"] != optic_name :
                hdul[0].header["OPTIC"] = optic_name
                print(f"The 'OPTIC' is set {hdul[0].header['OPTIC']}")

            hdul[0].header['FOCALLEN'] = OPTICDIC[hdul[0].header['OPTIC']]['FOCALLEN']
            print(f"The 'FOCALLEN' is set {hdul[0].header['FOCALLEN']}...")
            hdul[0].header['FOCRATIO'] = OPTICDIC[hdul[0].header['OPTIC']]["FOCALLEN"]/OPTICDIC[hdul[0].header['OPTIC']]["APATURE"]
            print(f"The 'FOCRATIO' is set {hdul[0].header['FOCRATIO']}...")
        
        ###########################
        #### 
        if (not 'TELESCOP' in hdul[0].header):
            hdul[0].header['TELESCOP'] = "-"
            print(f"The 'TELESCOP' is set {hdul[0].header['TELESCOP']}...")
        
        ###########################
        #### 
        if (not 'XBINNING' in hdul[0].header)\
            and (hdul[0].header["CCDNAME"] == "STX-16803") :
            if hdul[0].header['NAXIS1'] == 4096 \
                or  hdul[0].header['NAXIS2'] == 4096 :
                hdul[0].header['XBINNING'] = 1
                hdul[0].header['YBINNING'] = 1
                hdul[0].header['TELESCOP'] = "-"   
        
            elif hdul[0].header['NAXIS1'] == 2048 \
                or  hdul[0].header['NAXIS2'] == 2048 :
                hdul[0].header['XBINNING'] = 2
                hdul[0].header['YBINNING'] = 2
                hdul[0].header['TELESCOP'] = "-"
        
            elif hdul[0].header['NAXIS1'] == 1024 \
                or  hdul[0].header['NAXIS2'] == 1024 :
                hdul[0].header['XBINNING'] = 3
                hdul[0].header['YBINNING'] = 3
                hdul[0].header['TELESCOP'] = "-"
        hdul[0].header['XBINNING'] = int(hdul[0].header['XBINNING'])
        hdul[0].header['YBINNING'] = int(hdul[0].header['YBINNING'])
        print(f"The 'XBINNING', 'YBINNING' are set {hdul[0].header['XBINNING']}, \
                {hdul[0].header['YBINNING']},...")

        ###########################
        ####
        if (not 'XPIXSZ' in hdul[0].header) \
                and CCDNAME == 'STX-16803' :
            hdul[0].header['XPIXSZ'] = 9 * hdul[0].header['XBINNING']
            hdul[0].header['YPIXSZ'] = 9 * hdul[0].header['YBINNING']
            print(f"The 'XPIXSZ' and 'YPIXSZ' are set {9 * hdul[0].header['XBINNING']} \
                and {9 * hdul[0].header['YBINNING']}...")
        #hdul[0].header['GAIN'] = GAINDIC[CCDNAME]
        #hdul[0].header['GAIN'] = CCDDIC[hdul[0].header['CCDNAME']]['GAIN']
        #print(f"The 'GAIN' is set {hdul[0].header['GAIN']}...")
        #hdul[0].header['EGAIN'] = GAINDIC[CCDNAME]
        #hdul[0].header['EGAIN'] = CCDDIC[hdul[0].header['CCDNAME']]['GAIN']
        #print(f"The 'EGAIN' is set {hdul[0].header['EGAIN']}...")
        #hdul[0].header['RDNOISE'] = RDNOISEDIC[CCDNAME]
        #hdul[0].header['RDNOISE'] = CCDDIC[hdul[0].header['CCDNAME']]['RDNOISE']
        #print(f"The 'RDNOISE' is set {hdul[0].header['RDNOISE']}...")
        
        ###########################
        ####     
        if not "CCD-TEMP" in hdul[0].header :
            hdul[0].header['CCD-TEMP'] = 'N'
            print(f"The 'CCD-TEMP' is set {hdul[0].header['CCD-TEMP']}...")

        ###########################
        #### 
        if "EXPOSURE" in hdul[0].header :
            if not "EXPTIME" in hdul[0].header :
                hdul[0].header["EXPTIME"] = hdul[0].header["EXPOSURE"]
                print(f"The 'EXPTIME' is set {hdul[0].header['EXPOSURE']}...")
        elif "EXPTIME" in hdul[0].header :
            hdul[0].header["EXPOSURE"] = hdul[0].header["EXPTIME"]
            print(f"The 'EXPOSURE' is set {hdul[0].header['EXPTIME']}...")
        else :
            hdul[0].header["EXPTIME"] = 'N'
            hdul[0].header["EXPOSURE"] = 'N'
            print(f"The 'EXPTIME' and 'EXPOSURE' are set 'N'...")

        ###########################
        #### 
        if not "OPTIC" in hdul[0].header :
            hdul[0].header["OPTIC"] = optic_name
        if not "FOCALLEN" in hdul[0].header :
            hdul[0].header["FOCALLEN"] = OPTICDIC[hdul[0].header['OPTIC']]['FOCALLEN']
        if not "APATURE" in hdul[0].header :
            hdul[0].header["APATURE"] = OPTICDIC[hdul[0].header['OPTIC']]["APATURE"]

        print(hdul[0].header['OPTIC']+'_'+hdul[0].header['CCDNAME'])
        if "FLAT" in hdul[0].header["IMAGETYP"] or \
            "LIGHT" in hdul[0].header["IMAGETYP"] :
            if not "PIXSCALE" in hdul[0].header :
                hdul[0].header["PIXSCALE"] = calPixScale(hdul[0].header['FOCALLEN'], 
                                                            hdul[0].header['XPIXSZ'],
                                                            hdul[0].header['XBINNING'],)
            hdul[0].header["PIXSCALE"] = calPixScale(hdul[0].header['FOCALLEN'], 
                                                            hdul[0].header['XPIXSZ'],
                                                            hdul[0].header['XBINNING'],)
        hdul[0].header['FLIPSTAT'] = " "
        print(f"The 'FLIPSTAT' is set {hdul[0].header['FLIPSTAT']}...")
        
        for checkKEY in checkKEYs: 
            print(f"{checkKEY}: ", hdul[0].header[checkKEY])

        hdul.flush()  # changes are written back to original.fits
        print('*'*30)
        print(f"The header of {fpath.name} is updated..")

    return hdul

#%%
class KevinFitsHeader():
    def __init__(self, fpath):
        self.fpath = Path(fpath)
        self.checkKEYs = ["OBJECT", "TELESCOP", "OPTIC", "CCDNAME", 'FILTER',
            "GAIN", "EGAIN", "RDNOISE", "FOCALLEN", "PIXSCALE",
            'XPIXSZ', 'YPIXSZ', "XBINNING", "YBINNING", "FLIPSTAT"]
        '''
        Parameters
        ----------
        fpath : string
            The fullname of input file...

        '''
    def append_header(self):
        with fits.open(str(self.fpath), mode="append") as self.hdul :
            for self.checkKEY in self.checkKEYs: 
                if not self.checkKEY in self.hdul[0].header :
                    self.hdul[0].header.append(self.checkKEY, 
                                    '', 
                                    f"The keyword '{self.checkKEY}' is added.") 
                print(f"{self.checkKEY}: ", self.hdul[0].header[self.checkKEY])

            self.hdul.flush()  # changes are written back to original.fits
        return self.hdul

    def update_header(self): 
        self.foldername_el = self.fpath.parts[-2].split('_')
        self.fname_el = self.fpath.parts[-1].split('_')
        print("foldername_el", self.foldername_el)
        print("fname_el", self.fname_el)
        self.object_name = self.foldername_el[0]
        self.filter_name = self.fname_el[2]
        self.optic_name = self.foldername_el[5]
        self.ccd_name = self.foldername_el[6]
        print("object_name", self.object_name)
        print("filter_name", self.filter_name)
        print("optic_name", self.optic_name)
        print("ccd_name", self.ccd_name)
   
        # Change something in hdul.
        with fits.open(str(self.fpath), mode="update") as self.hdul :
            
            #if object_name != hdul[0].header["OBJECT"] : 
            self.hdul[0].header["OBJECT"] = self.object_name.upper()
            print(f"The 'OBJECT' is set {self.object_name.upper()}")
        
            if len(self.hdul[0].header['DATE-OBS']) == 10 \
                and 'TIME-OBS' in self.hdul[0].header : 
                self.hdul[0].header['DATE-OBS'] += 'T' + self.hdul[0].header['TIME-OBS']
                print(f"The 'DATE-OBS' is set {self.hdul[0].header['DATE-OBS']}")

            if "ze" in self.hdul[0].header["IMAGETYP"].lower() \
                    or "bi" in self.hdul[0].header["IMAGETYP"].lower() :
                self.hdul[0].header["IMAGETYP"] = "BIAS"
                print(f"The 'IMAGETYP' is set {self.hdul[0].header['IMAGETYP']}")
            elif "da" in self.hdul[0].header["IMAGETYP"].lower() :
                self.hdul[0].header["IMAGETYP"] = "DARK"
                print(f"The 'IMAGETYP' is set {self.hdul[0].header['IMAGETYP']}")
            elif "fl" in self.hdul[0].header["IMAGETYP"].lower() :
                self.hdul[0].header["IMAGETYP"] = "FLAT"
                print(f"The 'IMAGETYP' is set {self.hdul[0].header['IMAGETYP']}")
            elif "da" in self.hdul[0].header["IMAGETYP"].lower() \
                    or "lig" in self.hdul[0].header["IMAGETYP"].lower() :
                self.hdul[0].header["IMAGETYP"] = "LIGHT"
                print(f"The 'IMAGETYP' is set {self.hdul[0].header['IMAGETYP']}")
            
            if "BIAS" in self.hdul[0].header["IMAGETYP"] \
                or "DARK" in self.hdul[0].header["IMAGETYP"] :
                self.hdul[0].header["FILTER"] = "-"
                print(f"The 'FILTER' is set {self.hdul[0].header['FILTER']}")
                self.hdul[0].header['OPTIC'] = "-"
                print(f"The 'OPTIC' is set {self.hdul[0].header['OPTIC']}")

            if "FLAT" in self.hdul[0].header["IMAGETYP"] \
                or "LIGHT" in self.hdul[0].header["IMAGETYP"] :
                self.hdul[0].header["FILTER"] = self.filter_name.upper()
                print(f"The 'FILTER' is set {self.hdul[0].header['FILTER']}")
                self.hdul[0].header["OPTIC"] = self.optic_name
                print(f"The 'OPTIC' is set {self.hdul[0].header['OPTIC']}")

            try : 
                if 'qsi' in self.hdul[0].header['INSTRUME'] : 
                    self.CCDNAME = 'QSI683ws'
                elif 'st-8300' in self.hdul[0].header['INSTRUME'] : 
                    self.CCDNAME = 'ST-8300M'
                elif 'stf-8300' in self.hdul[0].header['INSTRUME'] : 
                    self.CCDNAME = 'STF-8300M'
                elif '11000' in self.hdul[0].header['INSTRUME'] : 
                    self.CCDNAME = 'STL-11000M'
                elif '16803' in self.hdul[0].header['INSTRUME'] : 
                    self.CCDNAME = 'STX-16803'
                elif "SBIG" in self.hdul[0].header['INSTRUME'] :
                    if self.hdul[0].header['XPIXSZ'] == 5.4 \
                            or self.hdul[0].header['XPIXSZ'] == 10.8 :
                        self.CCDNAME = 'STF-8300M'
                    elif self.hdul[0].header['XPIXSZ'] == 9.0 \
                            or self.hdul[0].header['XPIXSZ'] == 18.0 :
                        if self.hdul[0].header['NAXIS1'] == 2048 \
                                or self.hdul[0].header['NAXIS1'] == 4096 :
                            self.CCDNAME = 'STX-16803'
                        elif self.hdul[0].header['NAXIS1'] == 4008 \
                                or  self.hdul[0].header['NAXIS1'] == 2672 \
                                or  self.hdul[0].header['NAXIS1'] == 2004 \
                                or  self.hdul[0].header['NAXIS1'] == 1336 :
                            self.CCDNAME = 'STL-11000M'
                else :
                    self.CCDNAME = self.ccd_name    
            except :
                self.CCDNAME = self.ccd_name
            print("CCDNAME", self.CCDNAME)

            #self.hdul[0].header["CCDNAME"] = self.CCDNAME
            #print(f"The 'CCDNAME' is set {self.CCDNAME}...")

            if (not 'XPIXSZ' in self.hdul[0].header) \
                and self.CCDNAME == 'STX-16803' :
                self.hdul[0].header['XPIXSZ'] = 9 * self.hdul[0].header['XBINNING']
                self.hdul[0].header['YPIXSZ'] = 9 * self.hdul[0].header['YBINNING']
                print(f"The 'XPIXSZ' and 'YPIXSZ' are set {9 * self.hdul[0].header['XBINNING']} and {9 * self.hdul[0].header['YBINNING']}...")

            if not "CCD-TEMP" in self.hdul[0].header :
                self.hdul[0].header['CCD-TEMP'] = 'N'
                print(f"The 'CCD-TEMP' is set {self.hdul[0].header['CCD-TEMP']}...")

            if not "EXPOSURE" in self.hdul[0].header :
                self.hdul[0].header["EXPOSURE"] = self.hdul[0].header["EXPTIME"]
                print(f"The 'EXPOSURE' is set {self.hdul[0].header['EXPTIME']}...")

            self.hdul[0].header['GAIN'] = CCDDIC[self.CCDNAME]["GAIN"]
            print(f"The 'GAIN' is set {self.hdul[0].header['GAIN']}...")
            self.hdul[0].header['EGAIN'] = CCDDIC[self.CCDNAME]["EGAIN"]
            print(f"The 'EGAIN' is set {self.hdul[0].header['EGAIN']}...")
            self.hdul[0].header['RDNOISE'] = CCDDIC[self.CCDNAME]["RDNOISE"]
            print(f"The 'RDNOISE' is set {self.hdul[0].header['RDNOISE']}...")
            self.hdul.flush()  # changes are written back to original.fits
            print('*'*30)
            print(f"The header of {self.fpath.name} is updated..")
        return self.hdul



#%%
def fits_newpath(
        fpath,
        rename_by,
        mkdir_by=None,
        header=None,
        delimiter='_',
        fillnan="",
        fileext='.fits'
):
    ''' Gives the new path of the FITS file from header.
    Parameters
    ----------
    fpath : path-like
        The path to the original FITS file.
    rename_by : list of str, optional
        The keywords of the FITS header to rename by.
    mkdir_by : list of str, optional
        The keys which will be used to make subdirectories to classify files.
        If given, subdirectories will be made with the header value of the
        keys.
    header : Header object, optional
        The header to extract `rename_by` and `mkdir_by`. If `None`, the
        function will do ``header = fits.getheader(fpath)``.
    delimiter : str, optional
        The delimiter for the renaming.
    fillnan : str, optional
        The string that will be inserted if the keyword is not found from the
        header.
    fileext : str, optional
        The extension of the file name to be returned. Normally it should be
        ``'.fits'`` since this function is `fits_newname`, but you may prefer,
        e.g., ``'.fit'`` for some reason. If `fileext` does not start with a
        period (``"."``), it is automatically added to the final file name in
        front of the ``fileext``.
    Returns
    -------
    newpath : path
        The new path.
    '''

    if header is None:
        hdr = fits.getheader(fpath)
    else:
        hdr = header.copy()

    # First make file name without parent path
    hdrvals = []
    for k in rename_by:
        try:
            hdrvals.append(str(hdr[k]))
        except KeyError:
            hdrvals.append(fillnan)

    if not fileext.startswith('.'):
        fileext = f".{fileext}"

    newname = delimiter.join(list(hdrvals))  # just in case, re-listify...
    newname = newname + fileext
    newpath = Path(fpath.parent)

    if mkdir_by is not None:
        for k in mkdir_by:
            newpath = newpath / hdr[k]

    newpath = newpath / newname

    return newpath

#%%
#########################################
#KevinPSolver
#########################################
def KevinSolver(fpath, 
                    #solved_dir,
                    **kwargs
                    #downsample,
                    #pixscale,
                    ):
    """
    Parameters
    ----------
    fpath : path-like
        The path to the original FITS file.

    solved dir: string
        The directory where the output file

    pixscale : int

    """
    if not 'downsample' in kwargs :
        downsample = 1
    else: 
        downsample = kwargs['downsample']
    print("downsample: ", downsample)

    if not 'pixscale' in kwargs :
        pixscale = 1.0
    else : 
        pixscale = kwargs['pixscale']
    print(f"pixscale: {pixscale:.03f}, L: {pixscale*0.97:.03f}, U: {pixscale*1.03:.03f}")

    fpath = Path(fpath)
    
    try :
        if not (fpath.parent/f'{fpath.stem}.fits').exists() :
            # solve command.
            # solve-field fullname.fit -O --cpulimit 120 --nsigma 15 -u app -L 1.2 -U 1.3 -N new_filename.fits -p --no-plots -D output_directory {0}
            with subprocess.Popen(['solve-field', 
                                '-O', #--overwrite: overwrite output files if they already exist
                                #'-g', #--guess-scale: try to guess the image scale from the FITS headers
                                '--cpulimit', '10',  #will make it give up after 30 seconds.
                                #'--nsigma', '15',
                                '--downsample', f'{str(downsample)}',
                                '-u', 'app', #'--scale-units', 'arcsecperpix', #pixel scale
                                '-L', f'{pixscale*0.99:.03f}', 
                                '-U', f'{pixscale*1.01:.03f}',   
                                '-N', f'{fpath.parent/ fpath.stem}.fits', #--new-fits <filename>: output filename of the new FITS file containingthe WCS header; "none" to not create this file
                                #-p', 
                                '--no-plots',#: don't create any plots of the results
                                #'-D', str(SOLVEDDIR),
                                str(fpath)
                                ], 
                                stdout=subprocess.PIPE) as proc :
                print(proc.stdout.read())
    
    except Exception as err :
        print('{1} ::: {2} with {0} ...'\
            .format(fpath, datetime.now(), err))  

    try :
        if not (fpath.parent/f'{fpath.stem}.fits').exists():
            #https://www.hnsky.org/astap.htm#astap_command_line
            with subprocess.Popen(['astap', 
                                '-f', str(fpath), 
                                # '-o', 
                                #'-fov',
                                '-z', f'{str(downsample)}',
                                '-wcs',
                                '-analyse2',
                                '-update',],
                                stdout=subprocess.PIPE) as proc :
                print(proc.stdout.read())
            
            if (fpath.parent/f'{fpath.stem}.wcs').exists() \
                and (fpath.parent/f'{fpath.stem}.ini').exists():
                shutil.copy(str(fpath), f"{fpath.parent / fpath.stem}.fits")
                print(f"{fpath.name} is solved using ASTAP...")
        
    except Exception as err :
        print('{1} ::: {2} with {0} ...'\
            .format(fpath, datetime.now(), err))

    if fpath.exists() and (fpath.parent/f'{fpath.stem}.fits').exists():
        os.remove(str(fpath))
        print(f"{str(fpath)} is removed...")


#%%
#########################################
# LOCALPSolver
#########################################
def LOCALPSolver(fpath, 
                    #solved_dir,
                    **kwargs
                    #downsample,
                    #pixscale,
                    ):
    """
    Parameters
    ----------
    fpath : path-like
        The path to the original FITS file.

    solved dir: string
        The directory where the output file

    pixscale : int

    """
    fpath = Path(fpath)
    if not 'downsample' in kwargs :
        downsample = 1
    else: 
        downsample = kwargs['downsample']
    print("downsample: ", downsample)

    if not 'pixscale' in kwargs :
        hdul = fits.open(fpath)
        if 'PIXSCALE' in hdul[0].header:
            pixscale = hdul[0].header['PIXSCALE']
        else : 
            pixscale = calPixScale(hdul[0].header['FOCALLEN'], hdul[0].header['XPIXSZ'])
        hdul.close()
    else : 
        pixscale = kwargs['pixscale']
    print(f"pixscale: {pixscale:.03f}, L: {pixscale*0.97:.03f}, U: {pixscale*1.03:.03f}")
    
    # try :
        # solve command.
        # solve-field fullname.fit -O --cpulimit 120 --nsigma 15 -u app -L 1.2 -U 1.3 -N new_filename.fits -p --no-plots -D output_directory {0}
    with subprocess.Popen(['solve-field', 
                                '-O', #--overwrite: overwrite output files if they already exist
                                #'-g', #--guess-scale: try to guess the image scale from the FITS headers
                                '--cpulimit', '10',  #will make it give up after 30 seconds.
                                #'--nsigma', '15',
                                '--downsample', f'{str(downsample)}',
                                '-u', 'app', #'--scale-units', 'arcsecperpix', #pixel scale
                                '-L', f'{pixscale*0.99:.03f}', 
                                '-U', f'{pixscale*1.01:.03f}',   
                                '-N', f'{fpath.parent/ fpath.stem}.new', #--new-fits <filename>: output filename of the new FITS file containingthe WCS header; "none" to not create this file
                                #-p', 
                                '--no-plots',#: don't create any plots of the results
                                #'-D', str(SOLVEDDIR),
                                str(fpath)
                                ], 
                                stdout=subprocess.PIPE) as proc :
        print(proc.stdout.read())
    
    # except Exception as err :
    #     print('{1} ::: {2} with {0} ...'\
    #         .format(fpath, datetime.now(), err))  

    if fpath.exists() and (fpath.parent/f'{fpath.stem}.new').exists():
        #os.remove(str(fpath))
        #print(str(fpath))
        #shutil.move(str(fpath.parent/f'{fpath.stem}.new'), str(fpath.parent/f'{fpath.stem}.fits'))
        shutil.move(str(fpath.parent/f'{fpath.stem}.new'), str(fpath))
        print(str(fpath.parent/f'{fpath.stem}.new'), str(fpath))
        #print(f"{str(fpath)} is removed...")

#%%
#########################################
#ASTAPPSolver
#########################################
def ASTAPSolver(fpath, 
                    #solved_dir,
                    **kwargs
                    #downsample,
                    #pixscale,
                    ):
    """
    Parameters
    ----------
    fpath : path-like
        The path to the original FITS file.

    solved dir: string
        The directory where the output file

    pixscale : int

    """
    fpath = Path(fpath)
    if not 'downsample' in kwargs :
        downsample = 1
    else: 
        downsample = kwargs['downsample']
    print("downsample: ", downsample)
    
    if not 'pixscale' in kwargs :
        hdul = fits.open(fpath)
        if 'PIXSCALE' in hdul[0].header:
            pixscale = hdul[0].header['PIXSCALE']
        else : 
            pixscale = calPixScale(hdul[0].header['FOCALLEN'], hdul[0].header['XPIXSZ'])
        hdul.close()
    else : 
        pixscale = kwargs['pixscale']
    print(f"pixscale: {pixscale:.03f}, L: {pixscale*0.97:.03f}, U: {pixscale*1.03:.03f}")
    
    #try :
    #https://www.hnsky.org/astap.htm#astap_command_line
    with subprocess.Popen(['astap', 
                        '-f', str(fpath), 
                        #'-o', 
                        #'-fov',
                        '-z', f'{str(downsample)}',
                        '-wcs',
                        '-analyse2',
                        '-update',],
                        stdout=subprocess.PIPE) as proc :
        print(proc.stdout.read())
    
    if (fpath.parent/f'{fpath.stem}.wcs').exists() \
        and (fpath.parent/f'{fpath.stem}.ini').exists() \
        and (fpath.parent/f'{fpath.stem}.tmp').exists() :
        #shutil.move(str(fpath), f"{fpath.parent / fpath.stem}.fits")
        print(f"{fpath.parent / fpath.stem}.tmp")
        print(f"{fpath.parent / fpath.stem}.fits")
        #shutil.copy(f"{fpath.parent / fpath.stem}.tmp", f"{fpath.parent / fpath.stem}.fits")
        #shutil.move(f"{fpath.parent / fpath.stem}.tmp", f"{fpath.parent / fpath.stem}.fits")
        shutil.move(f"{fpath.parent / fpath.stem}.tmp", f"{fpath}")
        print(f"{fpath.name} is solved using ASTAP...")
    
    #except Exception as err :
        # print('{1} ::: {2} with {0} ...'\
        #     .format(fpath, datetime.now(), err))

    # if fpath.exists() and (fpath.parent/f'{fpath.stem}.fits').exists():
    #     os.remove(str(fpath))
    #     print(f"{str(fpath)} is removed...")



#%%
#########################################
# checkPSolve
#########################################
def checkPSolve(fpath, 
                    #solved_dir,
                    **kwargs,
                    #downsample,
                    #pixscale,
                    ):
    """
    Parameters
    ----------
    fpath : path-like
        The path to the original FITS file.

    solved dir: string
        The directory where the output file

    pixscale : int

    """

    fpath = Path(fpath)
    hdul = fits.open(fpath)
    PSKeys = ["CD1_1", "CD1_2", "CD2_1", "CD2_2", ]
    
    chk = 0
    SOLVE = False
    ASTAP = False
    LOCAL = False

    for PSKey in PSKeys :
        if PSKey in hdul[0].header : 
            chk += 1
    if chk == 4 : 
        SOLVE = True
        LOCAL = False
        ASTAP = False
        try : 
            for comment in hdul[0].header["COMMENT"]:
                if "scale:" in comment :
                    LOCAL = True
        except :
            LOCAL = False
                        
        if "PLTSOLVD" in hdul[0].header:
            try : 
                ASTAP = hdul[0].header["PLTSOLVD"]
            except : 
                ASTAP = False
    else : 
        SOLVE = False
        ASTAP = False
        LOCAL = False

    hdul.close()
    return SOLVE, ASTAP, LOCAL

#%%
#########################################
# makingAstrometrySH
#########################################
def makingAstrometrySH(fpath, 
                    solved_dir,
                    downsample,
                    pixscale,
                    ): 
    """
    Parameters
    ----------
    fpath : path-like
        The path to the original FITS file.

    solved dir: string
        The directory where the output file

    pixscale : int

    """
    if downsample is None :
        pixscale = 1
    print("downsample: ", downsample)

    if pixscale is None :
        pixscale = 1.0
    print("pixscale: ", pixscale)
    

    fpath = Path(fpath)
    SOLVEDDIR = Path(solved_dir)
    
    result = f'solve-field -O --cpulimit 10 --nsigma 15 -p {str(fpath)} '
    if downsample is not None:
        result += f'--downsample {str(downsample)} '
    if pixscale is not None:
        result += f'-u app -L {pixscale*0.9:.02f} -U {pixscale*1.1:.02f} '
    if solved_dir is not None:
        result += f'-D {str(SOLVEDDIR)} '
    result += '\n'
    
    print("result:", result)
    return result

#%%        
# =============================================================================
def get_RADEC_offset(hdul):
    """
    Parameters
    ----------
    hdul : fits file
    """
    Suwon = location = EarthLocation(lon=127.005 * u.deg, lat=37.308889 * u.deg, height=101 * u.m)
    #print('Starting get_new_foldername ...\n{0}'.format(filename))   
    w = WCS(hdul[0].header)
    t_obs = Time(hdul[0].header["DATE-OBS"]) + hdul[0].header["EXPOSURE"] * u.s / 2  # middle of observation time
    cent_coord = yfu.center_radec(ccd_or_header=hdul[0].header, center_of_image=True)
    print(f"calcualted (RA, DEC): ({cent_coord.ra}, {cent_coord.dec})")
    print(f"in header (RA, DEC): ({hdul[0].header['RA']*u.deg}, {hdul[0].header['DEC']*u.deg})")
    offset_RA = (cent_coord.ra - hdul[0].header['RA']*u.deg).to(u.arcmin)
    offset_DEC = (cent_coord.dec - hdul[0].header['DEC']*u.deg).to(u.arcmin)
    
    altaz = AltAz(obstime=t_obs, location=Suwon)

    cent_aa = cent_coord.transform_to(altaz)
    print(f"calculated (Az, Alt): ({cent_aa.az}, {cent_aa.alt})")
    print(f"in header (Az, Alt): ({hdul[0].header['CENTAZ ']*u.deg}, {hdul[0].header['CENTALT']*u.deg})")
    offset_AZ = (cent_aa.az - hdul[0].header['CENTAZ']*u.deg).to(u.arcmin)
    offset_ALT = (cent_aa.alt - hdul[0].header['CENTALT']*u.deg).to(u.arcmin)
    
    return offset_RA, offset_DEC, offset_AZ, offset_ALT

#%%        
# =============================================================================
def get_new_foldername_from_filename(filename):
    """
    Parameters
    ----------
    filename : str, path-like
        The path to the original FITS file.
    """
    #print('Starting get_new_foldername ...\n{0}'.format(filename))   
    filename_stem = filename.split(".")
    filename_el = filename_stem[-2].split("_")
    #print("filename_el: ", filename_el)
    timez = 9

    if int(filename_el[3][17:19])>=60 :
        obs_UT = datetime.strptime("{}59".format(filename_el[3][:17]), '%Y-%m-%d-%H-%M-%S')
    else:
         obs_UT = datetime.strptime(filename_el[3], '%Y-%m-%d-%H-%M-%S')
    obs_LST = obs_UT + timedelta(hours = timez)
    if obs_LST.hour < 12 :
        obs_LST = obs_LST - timedelta(days = 1)
    filename_el[3] = obs_LST.strftime('%Y-%m-%d-%H-%M-%S')
    if filename_el[1].upper() == 'BIAS':
        new_foldername = '{6}_{8}/Cal/-_{1}_-_{3}_-_-_{6}_-_{8}/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    elif filename_el[1].upper() == 'DARK' :
        new_foldername = '{6}_{8}/Cal/-_{1}_-_{3}_{4}_-_{6}_-_{8}/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    elif filename_el[1].upper() == 'FLAT' :
        new_foldername = '{6}_{8}/Cal_{5}/-_{1}_-_{3}_-_{5}_{6}_-_{8}/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    else : 
        new_foldername = '{6}_{8}/LIGHT_{5}/{0}_{1}_-_{3}_-_{5}_{6}_-_{8}/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    #write_log(log_file, 
    #            '{1} ::: \nNew foldername is {0} ...'\
    #            .format(new_foldername, datetime.now()))    
    return new_foldername


def get_new_foldername(filename):
    #log_file = 'get_new_foldername.log'
    print('Starting get_new_foldername ...\n{0}'.format(filename))
    
    filename_el1 = filename.split("bin")
    filename_el = filename_el1[0].split("_")
    
    if filename_el[1].lower() == 'BIAS':
        new_foldername = '{6}_{8}bin/Cal/-_{3}_-_{1}_-_{4}_-_{6}_-_{8}bin/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    elif filename_el[1].lower() == 'DARK' :
        new_foldername = '{6}_{8}bin/Cal/-_{3}_-_{1}_-_{4}_-_{6}_-_{8}bin/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    elif filename_el[1].lower() == 'FLAT' :
        new_foldername = '{6}_{8}bin/Cal_{5}/-_{3}_-_{1}_-_{5}_{6}_-_{8}bin/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    else : 
        new_foldername = '{6}_{8}bin/LIGHT_{5}/{0}_{1}_-_{3}_-_{5}_{6}_-_{8}bin/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    #write_log(log_file, 
    #            '{1} ::: \nNew foldername is {0} ...'\
    #            .format(new_foldername, datetime.now()))    
    return new_foldername


#%%                                
def connectMariaDB():
    #import pymysql
    import pymysql.cursors
    #conda install pymysql
    
    #mariaDB info
    db_host = 'parksparks.iptime.org'
    db_user = 'root'
    db_pass = 'rlgusl01'
    db_name = 'CCD_obs'
    db_port = 3307
        
    conn = pymysql.connect(host = db_host,
                          port = db_port,
                          user = db_user, password = db_pass,
                          db = db_name, charset = 'utf8mb4',
                          cursorclass = pymysql.cursors.DictCursor)
    
    return conn


#%%
def subp_solve_field(fullname, save_dir_name, sub_start_time): 
    import subprocess
    print('-'*60)
    print(fullname)
    with subprocess.Popen(['solve-field', 
                           '-O', #--overwrite: overwrite output files if they already exist
                           #'--scale-units', 'arcsecperpix', #pixel scale
                           #'--scale-low', '0.1', '--scale-high', '0.40', #pixel scale
                           '-g', #--guess-scale: try to guess the image scale from the FITS headers
                           #'-p', # --no-plots: don't create any plots of the results
                           '-D', '{0}'.format(save_dir_name), 
                           '{0}'.format(fullname)], 
                          stdout=subprocess.PIPE) as proc :
        print(proc.stdout.read())

        '''
        solve-field -O fullname
       '''
    return 0