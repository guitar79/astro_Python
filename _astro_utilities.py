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
import matplotlib.pyplot as plt
import shutil
from ccdproc import combine

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip, sigma_clipped_stats

from photutils.aperture import CircularAperture as CAp
from photutils.aperture import CircularAnnulus as CAn
from photutils import DAOStarFinder, detect_threshold
from photutils.centroids import centroid_com

from astroquery.jplhorizons import Horizons
from astroquery.imcce import Skybot
from astroquery.astrometry_net import AstrometryNet

import ysfitsutilpy as yfu
import ysphotutilpy as ypu

import _Python_utilities

#%%
#########################################
#directory variables
#########################################

c_method = "median"

CCD_obs_raw_dir = "A3_CCD_obs_raw"
CCD_NEW_dir = "A1_CCD_new_files"
CCD_NEWUP_dir = "A2_CCD_newUpdated_files"
CCD_duplicate_dir = "A4_CCD_duplicate_files"

bad_fits_dir = "Bad_fits"
master_dir = "master_files_ys"
reduced_dir = "reduced"
reduced_dir2 = "reduced2"
reduced_nightsky_dir = "reduced_nightsky"
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

GSHS = EarthLocation(lon=127.005 * u.deg, lat=37.308889 * u.deg, height=101 * u.m),
#######################################################
# OBS instruments information 
#######################################################

#CCDNAME, PIXSIZE, GAIN, RENOISE    
CCDDIC = {  "ST-8300M"    : {"PIXSIZE":5.4,  "GAIN":0.37,    "RDNOISE":9.3}, 
            "STF-8300M"   : {"PIXSIZE":5.4,  "GAIN":0.37,    "RDNOISE":9.3,  "RESOLUTION":(3352, 2532)}, 
            "QSI683ws"    : {"PIXSIZE":5.4,  "GAIN":0.13,    "RDNOISE":8.0,  "RESOLUTION":(3326, 2504)},
            "STL-11000M"  : {"PIXSIZE":9.0,  "GAIN":0.8,     "RDNOISE":9.6,  "RESOLUTION":(4008, 2672)},
            "STX-16803"   : {"PIXSIZE":9.0,  "GAIN":1.27,    "RDNOISE":9.0,  "RESOLUTION":(4096, 4096)},
            "QHY8"        : {"PIXSIZE":5.4,  "GAIN": "-",    "RDNOISE":"-",  "RESOLUTION":(3032, 2016)},
        "ATR3CMOS26000KPA": {"PIXSIZE":3.76, 
                "GAIN": "-",
                "RDNOISE":"-",
                "RESOLUTION":(6224, 4168)},
        "TT-2600CP": {"PIXSIZE":3.76, 
                "GAIN": "-",
                "RDNOISE":"-",
                "RESOLUTION":(6224, 4168)},
        "ASI183MMPro": {"PIXSIZE":2.4, 
                "GAIN": "-",
                "RDNOISE":"-",
                "RESOLUTION":(5496, 3672)},
        "ASI6200MMPro": {"PIXSIZE":3.76, 
                "GAIN": "-",
                "RDNOISE":"-",
                "RESOLUTION":(9576, 6388)},
        "ASI2600MC": {"PIXSIZE":3.76, 
                "GAIN": "-",
                "RDNOISE":"-",
                "RESOLUTION":(6248, 4176)},
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
            "WORC300": {"APATURE": 300,
                       "FOCALLEN": 2400},
            "WORC300x80": {"APATURE": 300,
                       "FOCALLEN": 2400*0.8},
                       }

fnameKEYs = ["OBJECT", "IMAGETYP", "FILTER", "DATE-OBS", 
            "EXPOSURE", "OPTIC", "CCDNAME", "CCD-TEMP", "XBINNING"]
#######################################################

#%%
#########################################
# Visiulization
#########################################
from warnings import warn
from astropy.visualization import (

    ImageNormalize,
    LinearStretch,
    ZScaleInterval,
    simple_norm,
)

def znorm(image, stretch=LinearStretch(), **kwargs):
    return ImageNormalize(image, interval=ZScaleInterval(**kwargs), stretch=stretch)

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
        binn,
        verbose = False,
        **kwargs
    ) :
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
    if verbose == True :
        print("PIXScale :",PIXScale)
    return PIXScale
 
#%%
#########################################
#KvinFitsMover
#########################################
def KevinFitsNewFname(
    fpath,
    fnameKEYs = ["OBJECT", "IMAGETYP", "FILTER", "DATE-OBS", 
            "EXPOSURE", "OPTIC", "CCDNAME", "CCD-TEMP", "XBINNING"],
    verbose = False,
    **kwargs
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
    if verbose == True :
        print("fpath: ", fpath)
    hdul = fits.open(str(fpath))
    for fnameKEY in fnameKEYs: 
        if verbose == True :
            print(f"{fnameKEY}: ", hdul[0].header[fnameKEY])
        try :
            ccdtemp = str(int(hdul[0].header["CCD-TEMP"]))
        except : 
            ccdtemp = "N"
        if verbose == True :
            print("ccdtemp: ", ccdtemp)
        new_fname = hdul[0].header["OBJECT"]+"_"+hdul[0].header["IMAGETYP"]+"_"+hdul[0].header["FILTER"]+"_"
        new_fname += hdul[0].header["DATE-OBS"][:19].replace("T","-").replace(":","-")+"_"
        new_fname += str(int(hdul[0].header["EXPOSURE"]))+"sec_"
        new_fname += hdul[0].header["OPTIC"]+"_"+hdul[0].header["CCDNAME"]+"_"       
        new_fname += ccdtemp+"c_"+str(int(hdul[0].header["XBINNING"]))+"bin.fit"
        #new_fname += fpath.ext
        if verbose == True :
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
                "GAIN", "EGAIN", "RDNOISE", 
                "PIXSCALE", "FOCALLEN", "APATURE", "CCD-TEMP",
                'XPIXSZ', 'YPIXSZ',
                "XBINNING", "YBINNING", "FLIPSTAT", "EXPTIME", "EXPOSURE"],
    imgtype_update = False,
    fil_update = False,
    verbose = False,
    **kwargs
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
    if verbose == True :
        print("foldername_el", foldername_el)
        print("fname_el", fname_el)
    object_name = foldername_el[0].replace(" ","")
    if verbose == True :
        print("object_name", object_name)
    image_type = foldername_el[1]
    if verbose == True :
        print("image_type", image_type)
    filter_name = fname_el[2].upper()
    if verbose == True :
        print("filter_name", filter_name)
    optic_name = foldername_el[5]
    if verbose == True :
        print("optic_name", optic_name)
    ccd_name = foldername_el[6]
    if verbose == True :
        print("ccd_name", ccd_name)
    with fits.open(str(fpath), mode="append") as hdul :
        for checkKEY in checkKEYs: 
            if not checkKEY in hdul[0].header :
                hdul[0].header.append(checkKEY, 
                                '', 
                                f"The keyword '{checkKEY}' is added.") 
            if verbose == True :
                print(f"{checkKEY}: ", hdul[0].header[checkKEY])

        hdul.flush()  # changes are written back to original.fits

    # Change something in hdul.
    with fits.open(str(fpath), mode="update") as hdul :
        
        ###########################
        #### "OBJECT"
        hdul[0].header["OBJECT"] = object_name   #delete upper()
        if verbose == True :
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
            elif 'ASI6200MM' in hdul[0].header['INSTRUME'] : 
                CCDNAME = 'ASI6200MMPro'
            elif 'ASI2600MC' in hdul[0].header['INSTRUME'] : 
                CCDNAME = 'ASI2600MC'   
                hdul[0].header["FILTER"] = "-" 
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
        if verbose == True :
            print("CCDNAME", CCDNAME)

        hdul[0].header["CCDNAME"] = CCDNAME
        if verbose == True :
            print(f"The 'CCDNAME' is set {hdul[0].header['CCDNAME']}...")

        ###########################
        #### 'DATE-OBS'
        if len(hdul[0].header['DATE-OBS']) == 10 \
            and 'TIME-OBS' in hdul[0].header : 
            hdul[0].header['DATE-OBS'] += 'T' + hdul[0].header['TIME-OBS']
            if verbose == True :
                print(f"The 'DATE-OBS' is set {hdul[0].header['DATE-OBS']}")
        
        ###########################
        #### 'IMAGETYP'
        if imgtype_update == True :
            hdul[0].header["IMAGETYP"] = image_type.upper()
            if verbose == True :
                print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")
        if not "IMAGETYP" in hdul[0].header :
            hdul[0].header["IMAGETYP"] = image_type  
        elif "ze" in hdul[0].header["IMAGETYP"].lower() \
                or "bi" in hdul[0].header["IMAGETYP"].lower() :
            hdul[0].header["IMAGETYP"] = "BIAS"
            if verbose == True :
                print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")
            hdul[0].header["OBJECT"] = "-"
            if verbose == True :
                print(f"The 'OBJECT' is set {hdul[0].header['OBJECT']}")
        elif "da" in hdul[0].header["IMAGETYP"].lower() :
            hdul[0].header["IMAGETYP"] = "DARK"
            if verbose == True :
                print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")
            hdul[0].header["OBJECT"] = "-"
            if verbose == True :
                print(f"The 'OBJECT' is set {hdul[0].header['OBJECT']}")
        elif "fl" in hdul[0].header["IMAGETYP"].lower() :
            hdul[0].header["IMAGETYP"] = "FLAT"
            if verbose == True :
                print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")
            hdul[0].header["OBJECT"] = "-"
            if verbose == True :
                print(f"The 'OBJECT' is set {hdul[0].header['OBJECT']}")
        elif "obj" in hdul[0].header["IMAGETYP"].lower() \
                or "lig" in hdul[0].header["IMAGETYP"].lower() :
            hdul[0].header["IMAGETYP"] = "LIGHT"
            if verbose == True :
                print(f"The 'IMAGETYP' is set {hdul[0].header['IMAGETYP']}")

        if "BIAS" in hdul[0].header["IMAGETYP"] \
            or "DARK" in hdul[0].header["IMAGETYP"] :
            for _KEY in ['FILTER', 'OPTIC', 'FOCALLEN', 'APATURE', 'PIXSCALE',] :
                hdul[0].header[_KEY] = "-"
                if verbose == True :
                    print(f"The '{_KEY}' is set {hdul[0].header[_KEY]}")

        if "FLAT" in hdul[0].header["IMAGETYP"] \
            or "LIGHT" in hdul[0].header["IMAGETYP"] :
            if not "OPTIC" in hdul[0].header :
                hdul[0].header["OPTIC"] = optic_name
                if verbose == True :
                    print(f"The 'OPTIC' is set {hdul[0].header['OPTIC']}")
            elif  hdul[0].header["OPTIC"] != optic_name :
                hdul[0].header["OPTIC"] = optic_name
                if verbose == True :
                    print(f"The 'OPTIC' is set {hdul[0].header['OPTIC']}")
            if not "FOCALLEN" in hdul[0].header :
                hdul[0].header["FOCALLEN"] = OPTICDIC[hdul[0].header['OPTIC']]['FOCALLEN']
            if not "APATURE" in hdul[0].header :
                hdul[0].header["APATURE"] = OPTICDIC[hdul[0].header['OPTIC']]["APATURE"]

            if not "FILTER" in hdul[0].header :
                hdul[0].header["FILTER"] = foldername_el[2].replace(" ","")
                if verbose == True :
                    print(f"FILTER is set {hdul[0].header['FILTER']}")
            elif hdul[0].header["FILTER"] != filter_name and fil_update==True :
                hdul[0].header["FILTER"] = filter_name
                if verbose == True :
                    print(f"FILTER is updateed to {hdul[0].header['FILTER']}")

            hdul[0].header['FOCALLEN'] = OPTICDIC[hdul[0].header['OPTIC']]['FOCALLEN']
            if verbose == True :
                print(f"The 'FOCALLEN' is set {hdul[0].header['FOCALLEN']}...")
            hdul[0].header['FOCRATIO'] = OPTICDIC[hdul[0].header['OPTIC']]["FOCALLEN"]/OPTICDIC[hdul[0].header['OPTIC']]["APATURE"]
            if verbose == True :
                print(f"The 'FOCRATIO' is set {hdul[0].header['FOCRATIO']}...")

            if not "PIXSCALE" in hdul[0].header :
                hdul[0].header["PIXSCALE"] = calPixScale(hdul[0].header['FOCALLEN'], 
                                                            hdul[0].header['XPIXSZ'],
                                                            hdul[0].header['XBINNING'],)
            # hdul[0].header["PIXSCALE"] = calPixScale(hdul[0].header['FOCALLEN'], 
            #                                                 hdul[0].header['XPIXSZ'],
            #                                                 hdul[0].header['XBINNING'],)
        
        ##########################
        if (not 'TELESCOP' in hdul[0].header):
            hdul[0].header['TELESCOP'] = "-"
            if verbose == True :
                print(f"The 'TELESCOP' is set {hdul[0].header['TELESCOP']}...")
        ##########################
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
        if verbose == True :
            print(f"The 'XBINNING', 'YBINNING' are set {hdul[0].header['XBINNING']}, \
                {hdul[0].header['YBINNING']},...")

        ###########################
        ####
        if (not 'XPIXSZ' in hdul[0].header) \
                and CCDNAME == 'STX-16803' :
            hdul[0].header['XPIXSZ'] = 9 * hdul[0].header['XBINNING']
            hdul[0].header['YPIXSZ'] = 9 * hdul[0].header['YBINNING']
            if verbose == True :
                print(f"The 'XPIXSZ' and 'YPIXSZ' are set {9 * hdul[0].header['XBINNING']} \
                    and {9 * hdul[0].header['YBINNING']}...")
        if not 'GAIN' in hdul[0].header :
            hdul[0].header['GAIN'] = CCDDIC[hdul[0].header['CCDNAME']]['GAIN']
            if verbose == True :
                print(f"The 'GAIN' is set {hdul[0].header['GAIN']}...")
        #hdul[0].header['EGAIN'] = GAINDIC[CCDNAME]
        #hdul[0].header['EGAIN'] = CCDDIC[hdul[0].header['CCDNAME']]['GAIN']
        #print(f"The 'EGAIN' is set {hdul[0].header['EGAIN']}...")
        if not 'RDNOISE' in hdul[0].header :
            hdul[0].header['RDNOISE'] = CCDDIC[hdul[0].header['CCDNAME']]['RDNOISE']
            if verbose == True :
                print(f"The 'RDNOISE' is set {hdul[0].header['RDNOISE']}...")
        
        ###########################
        ####     
        if not "CCD-TEMP" in hdul[0].header :
            hdul[0].header['CCD-TEMP'] = 'N'
            if verbose == True :
                print(f"The 'CCD-TEMP' is set {hdul[0].header['CCD-TEMP']}...")

        ###########################
        #### 
        if "EXPOSURE" in hdul[0].header :
            if not "EXPTIME" in hdul[0].header :
                hdul[0].header["EXPTIME"] = hdul[0].header["EXPOSURE"]
                if verbose == True :
                    print(f"The 'EXPTIME' is set {hdul[0].header['EXPOSURE']}...")
        elif "EXPTIME" in hdul[0].header :
            hdul[0].header["EXPOSURE"] = hdul[0].header["EXPTIME"]
            if verbose == True :
                print(f"The 'EXPOSURE' is set {hdul[0].header['EXPTIME']}...")
        else :
            hdul[0].header["EXPTIME"] = 'N'
            hdul[0].header["EXPOSURE"] = 'N'
            if verbose == True :
                print(f"The 'EXPTIME' and 'EXPOSURE' are set 'N'...")

        ###########################
        #### 
        
        if verbose == True :
            print(hdul[0].header['OPTIC']+'_'+hdul[0].header['CCDNAME'])
               
        hdul[0].header['FLIPSTAT'] = " "
        if verbose == True :
            print(f"The 'FLIPSTAT' is set {hdul[0].header['FLIPSTAT']}...")
        
        for checkKEY in checkKEYs: 
            if verbose == True :
                print(f"{checkKEY}: ", hdul[0].header[checkKEY])

        hdul.flush()  # changes are written back to original.fits
        if verbose == True :
            print('*'*30)
            print(f"The header of {fpath.name} is updated..")

    return hdul


#%%
def fits_newpath(
        fpath,
        rename_by,
        mkdir_by=None,
        header=None,
        delimiter='_',
        fillnan="",
        fileext='.fit',
        verbose = False,
        **kwargs
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
                    solved_dir = None,
                    downsample = 4,
                    pixscale = None,
                    cpulimit = 30,
                    SOLVE = False, 
                    nsigma = 15,
                    tryASTAP = True, 
                    tryLOCAL = True,
                    tryASTROMETRYNET = False,
                    makeLOCALsh = False,
                    verbose = False,
                    **kwargs
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
    from astropy.coordinates import SkyCoord
    fpath = Path(fpath)

    if fpath.exists() and (fpath.parent/f'{fpath.stem}.new').exists():
        #print(str(fpath))
        shutil.move(str(fpath.parent/f'{fpath.stem}.new'), str(fpath))
        if verbose == True :
            print(str(fpath.parent/f'{fpath.stem}.new'), str(fpath))
            #print(f"{str(fpath)} is removed...")

    SOLVE, ASTAP, LOCAL = checkPSolve(fpath)
    if verbose == True :
        print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)
    if not SOLVE :
        hdul = fits.open(fpath)
        
        if pixscale is None :
            if 'PIXSCALE' in hdul[0].header:
                pixscale = hdul[0].header['PIXSCALE']
            else : 
                pixscale = calPixScale(hdul[0].header['FOCALLEN'], 
                                            hdul[0].header['XPIXSZ'],
                                            hdul[0].header['XBINNING'])

        if verbose == True :
            print(f"pixscale: {pixscale:.03f}, L: {pixscale*0.97:.03f}, U: {pixscale*1.03:.03f}")
        vfov = pixscale * hdul[0].header['NAXIS1']/3600    
        hfov = pixscale * hdul[0].header['NAXIS2']/3600
        if verbose == True :
            print("vfov , hfov  :", vfov , hfov )
        
        if "RA" in hdul[0].header and "DEC" in hdul[0].header :
            try : 
                spd = 180 + hdul[0].header['DEC']
                ra = hdul[0].header['RA']
                dec = hdul[0].header['DEC']
            except :
                ra_h, ra_m, ra_s = hdul[0].header['RA'].split(':')
                dec_d, dec_m, dec_s = hdul[0].header['DEC'].split(':')
                ra_new = f"{ra_h}h{ra_m}m{ra_s}s"
                dec_new = f"{dec_d}d{dec_m}m{dec_s}s"
                c = SkyCoord(ra=ra_new, dec=dec_new)
                spd = 180 + c.dec.deg
                ra = c.ra.deg
                dec = c.dec.deg
        elif "OBJECTRA" in hdul[0].header and "OBJECTDEC" in hdul[0].header :
            try : 
                spd = 180 + hdul[0].header['OBJECTDEC']
                ra = hdul[0].header['OBJECTRA']
                dec = hdul[0].header['OBJECTDEC']
            except :
                ra_h, ra_m, ra_s = hdul[0].header['OBJECTRA'].split(':')
                dec_d, dec_m, dec_s = hdul[0].header['OBJECTDEC'].split(':')
                ra_new = f"{ra_h}h{ra_m}m{ra_s}s"
                dec_new = f"{dec_d}d{dec_m}m{dec_s}s"
                c = SkyCoord(ra=ra_new, dec=dec_new)
                spd = 180 + c.dec.deg
                ra = c.ra.deg
                dec = c.dec.deg
        elif "OBJCTRA" in hdul[0].header and "OBJCTDEC" in hdul[0].header :
            try : 
                spd = 180 + hdul[0].header['OBJCTDEC']
                ra = hdul[0].header['OBJCTRA']
                dec = hdul[0].header['OBJCTDEC']
            except :
                ra_h, ra_m, ra_s = hdul[0].header['OBJCTRA'].split(' ')
                dec_d, dec_m, dec_s = hdul[0].header['OBJCTDEC'].split(' ')
                ra_new = f"{ra_h}h{ra_m}m{ra_s}s"
                dec_new = f"{dec_d}d{dec_m}m{dec_s}s"
                c = SkyCoord(ra=ra_new, dec=dec_new)
                spd = 180 + c.dec.deg
                ra = c.ra.deg
                dec = c.dec.deg
        else :
            spd = 360
            ra = 0
            dec = 0
        hdul.close()

        # trying plate solving using ASTAP twice...

        if tryASTAP == True : 
            if verbose == True :
                print(f"Trying to solve using ASTAP:\n   {fpath.parent/fpath} ") 
            ASTAP_solve_cmd = f"astap -f {str(fpath)} "
            # ASTAP_solve_cmd += f"-z {str(downsample)} "
            ASTAP_solve_cmd += f"-wcs -analyse2 -update "
            if verbose == True :
                print(ASTAP_solve_cmd)
            os.system(ASTAP_solve_cmd)

            SOLVE, ASTAP, LOCAL = checkPSolve(fpath)
            if verbose == True :
                print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)
            if not SOLVE :
                if verbose == True :
                    print(f"Trying again to solve using ASTAP (different option) :\n   {fpath.parent/fpath} ")  
                ASTAP_solve_cmd = f"astap -f {str(fpath)} "
                ASTAP_solve_cmd += f"-z {str(downsample)} "
                ASTAP_solve_cmd += f"-wcs -analyse2 -update "
                if verbose == True :
                    print(ASTAP_solve_cmd)
                os.system(ASTAP_solve_cmd)    
            else : 
                return 0
    else : 
        return 0
            
    ######################################
    ###### LOCAL SOLVER
    ######################################
    SOLVE, ASTAP, LOCAL = checkPSolve(fpath)
    if verbose == True :
        print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)
    if not SOLVE :
        if tryLOCAL == True : 
            if verbose == True :
                print(f"Trying to solve using LOCAL:\n   {fpath.parent/fpath} ")
            LOCAL_solve_cmd = f"solve-field -O --cpulimit {cpulimit} "
            # LOCAL_solve_cmd += f"-g --nsigma {nsigma} "
            LOCAL_solve_cmd += f"--downsample {str(downsample)} "
            # LOCAL_solve_cmd += f"--scale-units app -L {pixscale*0.95:.03f} -H {pixscale*1.05:.03f} "
            # LOCAL_solve_cmd += f"--ra {ra} --dec {dec} "
            LOCAL_solve_cmd += f"--no-plots {str(fpath)} "
            if verbose == True :
                print(LOCAL_solve_cmd)
            os.system(LOCAL_solve_cmd)

            if fpath.exists() and (fpath.parent/f'{fpath.stem}.new').exists():
                #print(str(fpath))
                shutil.move(str(fpath.parent/f'{fpath.stem}.new'), str(fpath))
                if verbose == True :
                    print(str(fpath.parent/f'{fpath.stem}.new'), str(fpath))
                
            SOLVE, ASTAP, LOCAL = checkPSolve(fpath)
            if verbose == True :
                print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)
            if not SOLVE :
                if tryLOCAL == True : 
                    if verbose == True :
                        print(f"Trying again to solve using LOCAL (different option) :\n   {fpath.parent/fpath} ")
                    LOCAL_solve_cmd = f"solve-field -O --cpulimit {cpulimit} "
                    LOCAL_solve_cmd += f"-g --nsigma {nsigma} "
                    LOCAL_solve_cmd += f"--downsample {str(downsample)} "
                    LOCAL_solve_cmd += f"--scale-units app -L {pixscale*0.95:.03f} -H {pixscale*1.05:.03f} "
                    LOCAL_solve_cmd += f"--ra {ra} --dec {dec} "
                    LOCAL_solve_cmd += f"--no-plots {str(fpath)} "
                    if verbose == True :
                        print(LOCAL_solve_cmd)
                    os.system(LOCAL_solve_cmd)

                    if fpath.exists() and (fpath.parent/f'{fpath.stem}.new').exists():
                        #print(str(fpath))
                        shutil.move(str(fpath.parent/f'{fpath.stem}.new'), str(fpath))
                        if verbose == True :
                            print(str(fpath.parent/f'{fpath.stem}.new'), str(fpath))
            else : 
                return 0 
    else : 
        return 0 
     
    SOLVE, ASTAP, LOCAL = checkPSolve(fpath)
    if verbose == True :
        print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)

    if not SOLVE :
        if tryASTROMETRYNET == True : 
            if verbose == True :
                print(f"Trying to solve using astrometry:\n  {fpath} ")
            from astroquery.astrometry_net import AstrometryNet
            ast = AstrometryNet()
            # ger from nova.astrometry.net
            ast.api_key = 'bldvwzzuvktnwfph' #must changed...
            
            submission_id = None
            solve_timeout = 600
            try_again = True                

            try :                
                while try_again:
                    try:
                        if not submission_id:
                            wcs_header = ast.solve_from_image(str(fpath),
                                                force_image_upload=True,
                                                solve_timeout = solve_timeout,
                                                submission_id=submission_id)
                        else:
                            wcs_header = ast.monitor_submission(submission_id,
                                                                solve_timeout = solve_timeout)
                    except TimeoutError as e:
                        submission_id = e.args[1]
                    else:
                        # got a result, so terminate
                        try_again = False

                if not wcs_header:
                    # Code to execute when solve fails
                    if verbose == True :
                        print("fits file solving failure...")

                else:
                    # Code to execute when solve succeeds
                    if verbose == True :
                        print("fits file solved successfully...")

                    with fits.open(str(fpath), mode='update') as hdul:
                        for card in wcs_header :
                            try: 
                                if verbose == True :
                                    print(card, wcs_header[card], wcs_header.comments[card])
                                hdul[0].header.set(card, wcs_header[card], wcs_header.comments[card])
                            except : 
                                if verbose == True :
                                    print(card)
                        hdul.flush

                    if verbose == True :
                        print(str(fpath)+" is created...")
            
            except Exception as err: 
                if verbose == True :
                    print("Err :", err)
                    pass
    else :
        return 0

    SOLVE, ASTAP, LOCAL = checkPSolve(fpath)
    if verbose == True :
        print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)

    if not SOLVE :
        if makeLOCALsh == True :
            #solve-field -O -g --cpulimit 15 --nsigma 15 --downsample 4 -u app -L 0.6 -U 0.63 --no-plots
            solve_sh += f"solve-field -O -g --cpulimit {cpulimit} --nsigma {nsigma} --downsample {downsample} -u app -L {pixscale*0.95:.03f} -U {pixscale*1.05:.03f} --no-plots {str(fpath)}\n"
            solve_sh += f"mv {fpath.parent/fpath.stem}.new {str(fpath)}\n"
            # solve_sh += f"rm {fpath.parent/fpath.stem}-indx.xyls\n"
            # solve_sh += f"rm {fpath.parent/fpath.stem}-indx.axy\n"
            # solve_sh += f"rm {fpath.parent/fpath.stem}.rdls\n"
            # solve_sh += f"rm {fpath.parent/fpath.stem}.corr\n"
            # solve_sh += f"rm {fpath.parent/fpath.stem}.solved\n"
            # solve_sh += f"rm {fpath.parent/fpath.stem}.match\n"
            # solve_sh += f"rm {fpath.parent/fpath.stem}.axy\n"
            # solve_sh += f"rm {fpath.parent/fpath.stem}.wcs\n"
            # print("result:", solve_sh)

            with open(f"__{datetime.now().strftime('%Y%m%d')}_todo_astrometry_solve.sh", 'a') as f:
                f.write(solve_sh) 
    else : 
        return 0

    return 0
        

#%%
#########################################
# checkPSolve
#########################################
def checkPSolve(fpath, 
                verbose = False,
                **kwargs,
                ):
    """
    Parameters
    ----------
    fpath : path-like
        The path to the original FITS file.

    return
    ----------
    SOLVE, ASTAP, LOCAL :  bool, bool, bool

    """

    fpath = Path(fpath)
    if verbose == True :
        print(str(fpath))
    if fpath.exists() and (fpath.parent/f'{fpath.stem}.new').exists():
        shutil.move(str(fpath.parent/f'{fpath.stem}.new'), str(fpath))
        if verbose == True :
            print(str(fpath.parent/f'{fpath.stem}.new'), str(fpath))
            #print(f"{str(fpath)} is removed...")

    hdul = fits.open(fpath)
    PSKeys = ["CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2",
              "CD1_1", "CD1_2", "CD2_1", "CD2_2", 
              "A_0_0", "A_0_1", "A_1_0","A_1_1",
              "PC1_1", "PC1_2", "PC2_1", "PC2_2", 
              ]
    
    chk = 0
    SOLVE = False
    ASTAP = False
    LOCAL = False

    for PSKey in PSKeys :
        if PSKey in hdul[0].header : 
            chk += 1
    if chk > 5 : 
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
    remove_ext  = [".ini", ".axy", ".corr", ".match", ".rdls", ".tmp", 
                   ".solved", "-indx.xyls", ".solved", "-objs.png", "-ngc.png", "-indx.png", 
                   ".wcs"]
    for ext in remove_ext : 
        if (fpath.parent / f"{fpath.stem}{ext}").exists() :
            os.remove(fpath.parent / f"{fpath.stem}{ext}")
            if verbose == True :
                print(f"{fpath.parent}/{fpath.stem}{ext} is removed...")
    
    if fpath.exists() and (fpath.parent/f'{fpath.stem}.new').exists():
        if verbose == True :
            print(str(fpath.parent/f'{fpath.stem}.new'), str(fpath))
            #print(f"{str(fpath)} is removed...")

        shutil.move(str(fpath.parent/f'{fpath.stem}.new'), str(fpath))
        
    return SOLVE, ASTAP, LOCAL


#%%
#########################################
#get_RADEC_offset
#########################################
def get_RADEC_offset(fpath,
                    verbose = False,
                    **kwarg
                    ):
    """
    Parameters
    ----------
    hdul : fits file
    """
    fpath = Path(fpath)
    hdul = fits.open(str(fpath))
    
    SOLVE, ASTAP, LOCAL = checkPSolve(fpath)
    if verbose == True :
        print(SOLVE, ASTAP, LOCAL)
    
    if not SOLVE :
        if verbose == True :
            print(f"The file is not solved...\n{fpath}")
        return None, None
    else : 
        #print('Starting get_new_foldername ...\n{0}'.format(filename))   
        w = WCS(hdul[0].header)
        t_obs = Time(hdul[0].header["DATE-OBS"]) + hdul[0].header["EXPOSURE"] * u.s / 2  # middle of observation time
        cent_coord = yfu.center_radec(ccd_or_header=hdul[0].header, center_of_image=True)
        if verbose == True :
            print(f"calcualted (RA, DEC): ({cent_coord.ra}, {cent_coord.dec})")
            print(f"in header (RA, DEC): ({hdul[0].header['RA']*u.deg}, {hdul[0].header['DEC']*u.deg})")
        offset_RA = (cent_coord.ra - hdul[0].header['RA']*u.deg).to(u.arcmin)
        offset_DEC = (cent_coord.dec - hdul[0].header['DEC']*u.deg).to(u.arcmin)
        
        if verbose == True :
            print(f"(offset_RA, offset_DEC) = ({offset_RA}, {offset_DEC})")
    return offset_RA, offset_DEC

#%%
#########################################
#get_AZALT_offset
#########################################
def get_AZALT_offset(fpath,
                    location = EarthLocation(lon=127.005 * u.deg, lat=37.308889 * u.deg, height=101 * u.m),
                    verbose = False,
                    **kwarg
                    ):
    """
    Parameters
    ----------
    hdul : fits file
    """
    fpath = Path(fpath)
    hdul = fits.open(str(fpath))
    
    SOLVE, ASTAP, LOCAL = checkPSolve(fpath)
    if verbose == True :
        print(SOLVE, ASTAP, LOCAL)
    
    if not SOLVE :
        if verbose == True :
            print(f"The file is not solved...\n{fpath}")
        return None, None, None, None
    else : 
        #print('Starting get_new_foldername ...\n{0}'.format(filename))   
        w = WCS(hdul[0].header)
        t_obs = Time(hdul[0].header["DATE-OBS"]) + hdul[0].header["EXPOSURE"] * u.s / 2  # middle of observation time
        cent_coord = yfu.center_radec(ccd_or_header=hdul[0].header, center_of_image=True)
        if verbose == True :
            print(f"calcualted (RA, DEC): ({cent_coord.ra}, {cent_coord.dec})")
            print(f"in header (RA, DEC): ({hdul[0].header['RA']*u.deg}, {hdul[0].header['DEC']*u.deg})")
        offset_RA = (cent_coord.ra - hdul[0].header['RA']*u.deg).to(u.arcmin)
        offset_DEC = (cent_coord.dec - hdul[0].header['DEC']*u.deg).to(u.arcmin)
        
        altaz = AltAz(obstime=t_obs, location=location)

        cent_aa = cent_coord.transform_to(altaz)
        if verbose == True :
            print(f"calculated (Az, Alt): ({cent_aa.az}, {cent_aa.alt})")
            print(f"in header (Az, Alt): ({hdul[0].header['CENTAZ ']*u.deg}, {hdul[0].header['CENTALT']*u.deg})")
        offset_AZ = (cent_aa.az - hdul[0].header['CENTAZ']*u.deg).to(u.arcmin)
        offset_ALT = (cent_aa.alt - hdul[0].header['CENTALT']*u.deg).to(u.arcmin)
        if verbose == True :
            print(f"(offset_RA, offset_DEC) = ({offset_RA}, {offset_DEC})")
            print(f"(offset_AZ, offset_ALT) = ({offset_AZ}, {offset_ALT})")
    return offset_RA, offset_DEC, offset_AZ, offset_ALT

#%% 
#########################################
# get_new_foldername
#########################################
def get_new_foldername_from_filename1(filename, 
                                     timez=9,
                                     verbose = False,
                                     **kwarg
                                     ):
    """Generates a new folder name based on a FITS filename.

    Parameters
    ----------
    filename : str, path-like
        Path to the FITS file.
    timez : int, optional
        Time zone offset, by default 9.

    Returns
    -------
    str
        New folder name.
    """

    fpath = Path(filename)
    filename_el = fpath.stem.split("_")

    try:
        obs_ut = datetime.strptime(filename_el[3], '%Y-%m-%d-%H-%M-%S')
    except ValueError:  # Handle seconds > 60
        obs_ut = datetime.strptime(filename_el[3][:17] + "59", '%Y-%m-%d-%H-%M-%S')

    obs_lst = obs_ut + timedelta(hours=timez)
    if obs_lst.hour < 12:
        obs_lst = obs_lst - timedelta(days=1)
    
    date_obs = obs_lst.strftime('%Y-%m-%d')

    imgtype = filename_el[1].upper()
    if imgtype == 'BIAS' or imgtype == "DARK":
        new_foldername = f"{filename_el[6]}_{filename_el[8]}/Cal/-_{imgtype}_-_{date_obs}_-_-_{filename_el[6]}_-_{filename_el[8]}/"
    elif imgtype == 'FLAT':
        new_foldername = f"{filename_el[6]}_{filename_el[8]}/Cal_{filename_el[5]}/-_{imgtype}_-_{date_obs}_-_{filename_el[5]}_{filename_el[6]}_-_{filename_el[8]}/"
    else:  # LIGHT
        new_foldername = f"{filename_el[6]}_{filename_el[8]}/LIGHT_{filename_el[5]}/{filename_el[0]}_{imgtype}_-_{date_obs}_-_{filename_el[5]}_{filename_el[6]}_-_{filename_el[8]}/"

    return new_foldername



def get_new_foldername_from_filename(filename,
                                        verbose = False,
                                        timez = 9,
                                        **kwarg
                                        ):
    """
    Parameters
    ----------
    filename : str, path-like
        The path to the original FITS file.
    """
    filename_stem = filename.split(".")
    filename_el = filename_stem[-2].split("_")
    if verbose == True :
        print('Starting get_new_foldername ...\n{0}'.format(filename))   
        print("filename_el: ", filename_el)

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

#%%
#########################################
# get_new_foldername
#########################################
from pathlib import Path

def get_new_foldername1(filename,
                        verbose = False,
                        **kwargs
                        ):
    '''
    Generates a new folder name based on a FITS filename.

    Parameters
    ----------
    filename : str
        FITS filename.

    Returns
    -------
    str
        New folder name.
    '''
    parts = filename.split("_")
    date_obs = parts[3][:10]  # Extract date part
    img_type = parts[1].upper()  # Image type
    binning = parts[-1].split("bin")[0]  # Binning
    filter_name = parts[2]
    
    if img_type in ('BIAS', 'DARK'):
        new_foldername = f"{parts[6]}_{binning}bin/Cal/-_{date_obs}_-_{img_type}_-_{parts[4]}_-_{parts[6]}_-_{binning}bin/"
    elif img_type == 'FLAT':
        new_foldername = f"{parts[6]}_{binning}bin/Cal_{filter_name}/-_{date_obs}_-_{img_type}_-_{filter_name}_{parts[6]}_-_{binning}bin/"
    else:  # Assumed LIGHT
        new_foldername = f"{parts[6]}_{binning}bin/LIGHT_{filter_name}/{parts[0]}_{img_type}_-_{date_obs}_-_{filter_name}_{parts[6]}_-_{binning}bin/"

    return new_foldername


def get_new_foldername(filename,
                        verbose = False,
                        **kwargs):
    #log_file = 'get_new_foldername.log'
    if verbose == True :
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
    if verbose == True :
        print("new_foldername :", new_foldername)    
    return new_foldername




#%%
#########################################
#combine_BDF
#########################################
def combine_BDF(BDFDIR,
                tryagain = False,
                verbose = False,
                **kwargs
                ):
    MASTERDIR = BDFDIR / master_dir
    if verbose == True :
        print(f"Starting: {str(BDFDIR.parts[-1])}")
        print ("MASTERDIR: ", format(MASTERDIR))
    
    if not MASTERDIR.exists():
        os.makedirs("{}".format(str(MASTERDIR)))
        if verbose == True :
            print("{} is created...".format(str(MASTERDIR)))

    summary = yfu.make_summary(BDFDIR/"*.fit*", 
                                verify_fix=True,
                                ignore_missing_simple=True,
                                verbose = verbose,
                                )
    
    if summary is not None :
        if verbose == True :
            #print(summary)
            print("len(summary):", len(summary))
            print("summary:", summary)
            #print(summary["file"][0])

        if (MASTERDIR / "master_bias.fits").exists() and tryagain == False :
            if verbose == True :
                print("bias file is already exist....")
        else :
            summary_bias = summary.loc[summary["IMAGETYP"] == "BIAS"].copy()
            summary_bias.reset_index(inplace=True)
            if verbose == True :
                print("summary_bias", summary_bias)

            bias_fits = summary_bias["file"]
            if verbose == True :
                # print("type(bias_fits)", type(bias_fits))
                print("len(bias_fits)", len(bias_fits))
                # print("bias_fits", bias_fits)

            bias_comb = yfu.group_combine(
                            bias_fits.tolist(),
                            type_key = ["IMAGETYP"],
                            type_val = ["BIAS"],
                            group_key = ["EXPTIME"],
                            fmt = "master_bias.fits",  # output file name format
                            outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                            combine = "med",
                            memlimit = 2.e+10,
                            verbose = False,
                        )

        summary_dark = summary.loc[summary["IMAGETYP"] == "DARK"].copy()
        summary_dark.reset_index(inplace=True)
        if verbose == True :
            print("summary_dark", summary_dark)

        if 'EXPTIME' in summary_dark :
            check_exptimes = summary_dark['EXPTIME'].drop_duplicates()
            check_exptimes = check_exptimes.reset_index(drop=True)
            if verbose == True :
                print("check_exptimes", check_exptimes)

            for exptime in check_exptimes :
                if (MASTERDIR / f"master_dark_{exptime:.0f}sec.fits" ).exists() and tryagain == False :
                    if verbose == True :
                        print(f"master_dark_{exptime:.0f}sec.fits already exist....")
                else :
                    summary_dark_each = summary_dark.loc[summary_dark['EXPTIME'] == exptime]
                    dark_fits = summary_dark_each['file']
                    if verbose == True :
                        # print("type(dark_fits)", type(dark_fits))
                        print("len(dark_fits)", len(dark_fits))
                        # print("dark_fits", dark_fits)

                    dark_comb = yfu.group_combine(
                                dark_fits.tolist(),
                                type_key = ["IMAGETYP"],
                                type_val = ["DARK"],
                                group_key = ["EXPTIME"],
                                fmt = "master_dark_{:.0f}sec.fits",  # output file name format
                                outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                                combine = "med",
                                memlimit = 2.e+10,
                                verbose = False,
                            )
        dark_fpaths = sorted(list((MASTERDIR).glob('*dark*.fit*')))
        if verbose == True : 
            print(f"dark_fpaths: {dark_fpaths}")
            print(f"len(dark_fpaths): {len(dark_fpaths)}")
                    
        summary_flat = summary.loc[summary["IMAGETYP"] == "FLAT"].copy()
        summary_flat.reset_index(inplace=True)

        if 'FILTER' in summary_flat :
            check_filters = summary_flat['FILTER'].drop_duplicates()
            check_filters = check_filters.reset_index(drop=True)
            if verbose == True :
                print("check_filters", check_filters)

            for filter in check_filters :
                if (MASTERDIR / f"master_flat_{filter:s}.fits" ).exists() and tryagain == False :
                    if verbose == True :
                        print(f"master_flat_{filter:s}.fits already exist....")
                else : 
                    summary_flat_each = summary_flat.loc[summary_flat['FILTER'] == filter]
                    flat_fits = summary_flat_each['file']
                    if verbose == True :
                        # print("type(flat_fits)", type(flat_fits))
                        print("len(flat_fits)", len(flat_fits))
                        # print("flat_fits", flat_fits)

                    try : 
                        flat_comb= yfu.group_combine(
                                        flat_fits.tolist(),
                                        type_key = ["IMAGETYP"],
                                        type_val = ["FLAT"],
                                        group_key = ["FILTER"],
                                        fmt = "master_flat_{:s}_norm.fits",  # output file name format
                                        scale="med_sc", #norm
                                        scale_to_0th=False, #norm
                                        outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                                        combine = "med",
                                        memlimit = 2.e+10,
                                        verbose=verbose,
                                    )
                    except : 
                        flat_comb = yfu.group_combine(
                                    flat_fits.tolist(),
                                    type_key = ["IMAGETYP"],
                                    type_val = ["FLAT"],
                                    group_key = ["FILTER"],
                                    fmt = "master_flat_{:s}.fits",  # output file name format
                                    #scale="med_sc", #norm
                                    #scale_to_0th=False, #norm
                                    outdir = MASTERDIR,  # output directory (will automatically be made if not exist)
                                    combine = "med",
                                    memlimit = 2.e+10,
                                    verbose=verbose,
                                )
        flat_fpaths = sorted(list((MASTERDIR).glob('*flat*.fit*')))
        if verbose == True : 
            print(f"flat_fpaths: {dark_fpaths}")
            print(f"len(flat_fpaths): {len(flat_fpaths)}")

    return 0


#%%
#########################################
# count_Num_stars
#########################################
def count_Num_stars(fpath, 
                    FWHM = 6,
                    verbose = False,
                    **kwargs
                    ):
    """
    Parameters
    ----------
    fpath : path-like
        The path to the original FITS file.

    """
    fpath = Path(fpath)
    hdul = fits.open(fpath)

    ## thres
    # thresh = detect_threshold(data=hdul[0].data, nsigma=3)
    # thresh = thresh[0][0]
    # print('detect_threshold', thresh)

    avg, med, std = sigma_clipped_stats(hdul[0].data)  # by default, 3-sigma 5-iteration.
    thresh = 5. * std
    # print('detect_threshold', thresh)

    DAOfind = DAOStarFinder(
                            fwhm = FWHM, 
                            threshold = thresh, 
                            # sharplo = 0.2, sharphi = 1.0,  # default values: sharplo=0.2, sharphi=1.0,
                            # roundlo = -1.0, roundhi = 1.0,  # default values -1 and +1
                            # sigma_radius = 1.5,           # default values 1.5
                            # ratio = 1.0,                  # 1.0: circular gaussian
                            exclude_border = True         # To exclude sources near edges
                            )
    # The DAOStarFinder object ("DAOfind") gets at least one input: the image.
    # Then it returns the astropy table which contains the aperture photometry results:
    DAOfound = DAOfind(hdul[0].data)
    if DAOfound is None :
        return 0
    else : 
        return len(DAOfound)


#%%
#########################################
# move_bad_fits
#########################################
def move_bad_fits(fpath, 
                    FWHM = 6,
                    verbose = False,
                    **kwargs
                    ):
    """
    Parameters
    ----------
    fpath : path-like
        The path to the original FITS file.

    """

    fpath = Path(fpath)
    hdul = fits.open(fpath)

    ## thres
    # thresh = detect_threshold(data=hdul[0].data, nsigma=3)
    # thresh = thresh[0][0]
    # print('detect_threshold', thresh)

    avg, med, std = sigma_clipped_stats(hdul[0].data)  # by default, 3-sigma 5-iteration.
    thresh = 5. * std
    # print('detect_threshold', thresh)

    DAOfind = DAOStarFinder(
                            fwhm = FWHM, 
                            threshold = thresh, 
                            # sharplo = 0.2, sharphi = 1.0,  # default values: sharplo=0.2, sharphi=1.0,
                            # roundlo = -1.0, roundhi = 1.0,  # default values -1 and +1
                            # sigma_radius = 1.5,           # default values 1.5
                            # ratio = 1.0,                  # 1.0: circular gaussian
                            exclude_border = True         # To exclude sources near edges
                            )
    # The DAOStarFinder object ("DAOfind") gets at least one input: the image.
    # Then it returns the astropy table which contains the aperture photometry results:
    DAOfound = DAOfind(hdul[0].data)
    if DAOfound is None :
        return 0
    else : 
        return len(DAOfound)

#%%
#########################################
#move_bad_fits
#########################################
def move_bad_fits(fpath,
                    verbose = False,
                    **kwarg,
                    ):
    fpath = Path(fpath)
    BADFITSDIR = fpath.parent / "Bad_fits"
    Num_stars = count_Num_stars(fpath, 
                            FWHM = 6,
                            )
    if verbose == True :
        print("Num_stars :", Num_stars)
    if Num_stars < 3 :
        if not BADFITSDIR.exists():
            os.makedirs(str(BADFITSDIR))
            if verbose == True :
                print(f"{BADFITSDIR} is created...")

        shutil.move(str(fpath), str(BADFITSDIR / fpath.name))
        if verbose == True :
            print(f"{str(fpath.name)} is moved to Bad_fits...")
    
    return 0


#%%
#########################################
#solving_fits_fils
#########################################
def solving_fits_file(DOINGDIR,
                SOLVINGDIR = None,
                downsample = 4,
                # count_stars= False,
                tryagain = False,  
                tryASTAP = True, 
                tryLOCAL = True,
                tryASTROMETRYNET = False,  
                verbose = False,
                **kwarg,
                ):
    DOINGDIR = Path(DOINGDIR)
    if verbose == True :
        print("DOINGDIR", DOINGDIR)

    if SOLVINGDIR is None :
        SOLVINGDIR = DOINGDIR
    else : 
        SOLVINGDIR = DOINGDIR / SOLVINGDIR

    # BADFITSDIR = DOINGDIR / bad_fits_dir
    # if not BADFITSDIR.exists() :
    #     os.mkdir(str(BADFITSDIR))
    #     if verbose == True :
    #         print(f"{str(BADFITSDIR)} is created...")

    summary = yfu.make_summary(SOLVINGDIR/"*.fit*",
                                    verify_fix=True,
                                    )
    if summary is not None :
        if verbose == True :
            print("len(summary):", len(summary))
            print("summary:", summary)
            #print(summary["file"][0])  
        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)
        if verbose == True :
            print("df_light:\n{}".format(df_light))

        for _, row  in df_light.iterrows():
            fpath = Path(row["file"])
            if verbose == True :
                print("fpath :" ,fpath)
            ####
            if 'reduced' in fpath.parts[-2] and \
                not (fpath.parents[1] / fpath.name).exists() :
                os.remove(fpath)
            else : 
                # hdul = fits.open(fpath)
            
                try : 
                    KevinSolver(fpath, 
                                # solved_dir = None,
                                downsample = downsample,
                                # pixscale = None ,
                                # cpulimit = 15,
                                tryASTAP = tryASTAP, 
                                tryLOCAL = tryLOCAL,
                                tryASTROMETRYNET = tryASTROMETRYNET, 
                                makeLOCALsh = False,
                                verbose = verbose,
                                )
                except Exception as err :
                    if verbose == True :
                        print("X"*60)
                        print(err)
                
    return 0

#%%
#########################################
#ccd_Reducuction
#########################################

def ccd_Reduction(DOINGDIR,
                    BDFDIR,
                    tryagain = False,
                    trynightsky = True,
                    file_age = 365,
                    verbose = False,
                    **kwarg) :
    DOINGDIR = Path(DOINGDIR)
    if verbose == True :
        print(f"Starting: {str(DOINGDIR.parts[-1])}")
    
    MASTERDIR = BDFDIR / master_dir
    sMASTERDIR = DOINGDIR / master_dir
    REDUCEDDIR = DOINGDIR / reduced_dir

    if not sMASTERDIR.exists():
        os.makedirs(str(sMASTERDIR))
        if verbose == True :
            print("{} is created...".format(str(sMASTERDIR)))

    if not REDUCEDDIR.exists():
        os.makedirs(str(REDUCEDDIR))
        if verbose == True :
            print("{} is created...".format(str(REDUCEDDIR)))

    summary = yfu.make_summary(DOINGDIR/"*.fit*", 
                            verify_fix=True,
                            ignore_missing_simple=True,
                            )
    if summary is not None :
        if verbose == True :
            #print(summary)
            print("len(summary):", len(summary))
            print("summary:", summary)
            #print(summary["file"][0])

        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)

        summary_master = yfu.make_summary(MASTERDIR/"*.fit*", 
                           verbose = False,
                           )
        if verbose == True :
            print("summary_master", summary_master)

        summary_master_dark = summary_master.loc[summary_master["IMAGETYP"] == "DARK"].copy()
        summary_master_dark.reset_index(inplace=True)
        if verbose == True :
            print("summary_master_dark", summary_master_dark)

        if 'EXPTIME' in summary_master_dark :
            check_exptimes = summary_master_dark['EXPTIME'].drop_duplicates()
            check_exptimes = check_exptimes.reset_index(drop=True)
            if verbose == True :
                print("check_exptimes", check_exptimes)

        for _, row in df_light.iterrows():

            fpath = Path(row["file"])     
            ccd = yfu.load_ccd(fpath)
            filt = ccd.header["FILTER"]
            expt = ccd.header["EXPTIME"]
            
            idx = abs(summary_master_dark['EXPTIME'] - expt).idxmin()
            if verbose == True :
                print(idx)

            if (REDUCEDDIR / fpath.name).exists() :
                if verbose == True :
                    print(f"The reduction file already exists...\n{fpath.name}")
                fpath_age = _Python_utilities.get_file_age(str(REDUCEDDIR / fpath.name))
                if tryagain == False and fpath_age.days < file_age :
                    if verbose == True :
                        print("*"*10)
                        print(f"The reduction file is younger than {file_age} days...\n{fpath.name}")
                    pass
                
                else :
                    try : 
                        if not (MASTERDIR / f"master_flat_{filt.upper()}.fits").exists() :
                            if verbose == True :
                                print(f"{MASTERDIR}/master_flat_{filt.upper()}.fits is not exists...")
                        else :
                            if (MASTERDIR / f"master_dark_{expt:.0f}sec.fits").exists() :
                                if verbose == True :
                                    print(f"Trying Reduction with master_dark_{expt:.0f}sec.fits ...")

                                red = yfu.ccdred(
                                    ccd,
                                    output=Path(f"{REDUCEDDIR/ fpath.name}"),
                                    mdarkpath=str(MASTERDIR / "master_dark_{:.0f}sec.fits".format(expt)),
                                    mflatpath=str(MASTERDIR / "master_flat_{}.fits".format(filt.upper())),
                                    # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                                    overwrite=True,
                                    )
                            else : 
                                if verbose == True :
                                    print(f"Trying Reduction with master_dark_{summary_master_dark['EXPTIME'][idx]:.0f}sec.fits is not exists...")
                                red = yfu.ccdred(
                                    ccd,
                                    output=Path(f"{REDUCEDDIR/ fpath.name}"),
                                    mdarkpath=str(MASTERDIR / f"master_dark_{summary_master_dark['EXPTIME'][idx]:.0f}sec.fits"),
                                    mflatpath=str(MASTERDIR / f"master_flat_{filt.upper()}.fits"),
                                    dark_scale = True,
                                    exptime_dark = summary_master_dark['EXPTIME'][idx],
                                    # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                                    overwrite=True,
                                    )
                            if verbose == True :
                                print (f"Reduce Reduce {fpath.name} +++...")

                    except Exception as err: 
                        if verbose == True :
                            print("Err :", err)
                            pass
            else :
       
                if (MASTERDIR / f"master_dark_{expt:.0f}sec.fits").exists() :
                    if verbose == True :
                        print(f"Reduce with master_dark_{expt:.0f}sec.fits ...")
                    try :
                        red = yfu.ccdred(
                            ccd,
                            output=Path(f"{REDUCEDDIR/ fpath.name}"),
                            mdarkpath=str(MASTERDIR / "master_dark_{:.0f}sec.fits".format(expt)),
                            mflatpath=str(MASTERDIR / "master_flat_{}_norm.fits".format(filt.upper())),
                            # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                            overwrite=True,
                            )
                    except : 
                        red = yfu.ccdred(
                            ccd,
                            output=Path(f"{REDUCEDDIR/ fpath.name}"),
                            mdarkpath=str(MASTERDIR / "master_dark_{:.0f}sec.fits".format(expt)),
                            mflatpath=str(MASTERDIR / "master_flat_{}.fits".format(filt.upper())),
                            # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                            overwrite=True,
                            )
                else : 
                    if verbose == True :
                        print(f"Reduce with master_dark_{summary_master_dark['EXPTIME'][idx]:.0f}sec.fits is not exists...")
                    try : 
                        red = yfu.ccdred(
                            ccd,
                            output=Path(f"{REDUCEDDIR/ fpath.name}"),
                            mdarkpath=str(MASTERDIR / f"master_dark_{summary_master_dark['EXPTIME'][idx]:.0f}sec.fits"),
                            mflatpath=str(MASTERDIR / f"master_flat_{filt.upper()}_norm.fits"),
                            dark_scale = True,
                            exptime_dark = summary_master_dark['EXPTIME'][idx],
                            # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                            overwrite=True,
                            )
                    except :
                        red = yfu.ccdred(
                            ccd,
                            output=Path(f"{REDUCEDDIR/ fpath.name}"),
                            mdarkpath=str(MASTERDIR / f"master_dark_{summary_master_dark['EXPTIME'][idx]:.0f}sec.fits"),
                            mflatpath=str(MASTERDIR / f"master_flat_{filt.upper()}.fits"),
                            dark_scale = True,
                            exptime_dark = summary_master_dark['EXPTIME'][idx],
                            # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                            overwrite=True,
                            )
                        
                if verbose == True :
                    print (f"Reduce Reduce {fpath.name} +++...")


    if trynightsky == True : 
        REDUCNSKYDIR = DOINGDIR / reduced_nightsky_dir
        if not REDUCNSKYDIR.exists():
            os.makedirs("{}".format(str(REDUCNSKYDIR)))
            if verbose == True :
                print("{} is created...".format(str(REDUCNSKYDIR)))
    
        summary = yfu.make_summary(REDUCEDDIR /"*.fit*")
        if summary is not None :
            if verbose == True :
                print("len(summary):", len(summary))
                print("summary:", summary)
                #print(summary["file"][0])   

            df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
            df_light = df_light.reset_index(drop=True)

            if 'FILTER' in df_light :
                check_filters = df_light['FILTER'].drop_duplicates()
                check_filters = check_filters.reset_index(drop=True)
                if verbose == True :
                    print("check_filters", check_filters)

                for filt in check_filters:
                #for filt in ["V"]:
                    df_light_filt = df_light.loc[df_light["FILTER"] == filt].copy()
                    
                    if df_light_filt.empty:
                        if verbose == True :
                            print(f"The dataframe(df_light_filt) {filt} is empty")
                        pass
                    else:
                        if verbose == True :
                            print("len(df_light_filt):", len(df_light_filt))
                            print("df_light_filt:", df_light_filt)
                    
                        if (sMASTERDIR / f"nightskyflat-{filt}.fits").exists() :
                            if verbose == True :
                                print(f"{sMASTERDIR}/nightskyflat-{filt}.fits is already exists")
                            fpath_age = _Python_utilities.get_file_age(sMASTERDIR / f"nightskyflat-{filt}.fits")
                            if tryagain == False and fpath_age.days < file_age :
                                if verbose == True :
                                    print("*"*10)
                                    print(f"The file is younger than {file_age} days...\n{fpath.name}")
                                pass

                            else :
                                if verbose == True :
                                    print("*"*10)
                                    print(f"The file is older than {file_age} days...\nTry gagin {fpath.name}")      
                                File_Num = 80
                                if len(df_light_filt["file"]) > File_Num :
                                    combine_lst = df_light_filt["file"].tolist()[:File_Num]
                                else : 
                                    combine_lst = df_light_filt["file"].tolist()
                                try : 
                                    ccd = yfu.imcombine(
                                                        combine_lst, 
                                                        combine="med",
                                                        scale="avg", 
                                                        scale_to_0th=False, #norm
                                                        reject="sc", 
                                                        sigma=2.5,
                                                        verbose=verbose,
                                                        memlimit = 2.e+11,
                                                        )
                                except :
                                    ccd = yfu.imcombine(
                                                        combine_lst, 
                                                        combine="med",
                                                        scale="avg", 
                                                        scale_to_0th=False, #norm
                                                        reject="sc", 
                                                        # sigma=2.5,
                                                        verbose=verbose,
                                                        memlimit = 2.e+11,
                                                        )
                                ccd.write(sMASTERDIR / f"nightskyflat-{filt}.fits", overwrite=True)
                                print (f"nightskyflat-{filt}.fits is created +++...")
                        else : 
                            File_Num = 80
                            if len(df_light_filt["file"]) > File_Num :
                                combine_lst = df_light_filt["file"].tolist()[:File_Num]
                            else : 
                                combine_lst = df_light_filt["file"].tolist()
                            try : 
                                ccd = yfu.imcombine(
                                                combine_lst, 
                                                combine="med",
                                                scale="avg", 
                                                scale_to_0th=False, #norm
                                                reject="sc", 
                                                sigma=2.5,
                                                verbose=verbose,
                                                memlimit = 2.e+11,
                                                )
                            except :
                                ccd = yfu.imcombine(
                                                combine_lst, 
                                                combine="med",
                                                scale="avg", 
                                                scale_to_0th=False, #norm
                                                reject="sc", 
                                                # sigma=2.5,
                                                verbose=verbose,
                                                memlimit = 2.e+11,
                                                )
                        ccd.write(sMASTERDIR / f"nightskyflat-{filt}.fits", overwrite=True)
                        if verbose == True :
                            print (f"nightskyflat-{filt}.fits is created +++...")

            for _, row in df_light.iterrows():
                fpath = Path(row["file"])
                ccd = yfu.load_ccd(REDUCEDDIR / fpath.name)
                filt = row["FILTER"]
                if (REDUCNSKYDIR / fpath.name).exists() and tryagain == False:
                    if verbose == True :
                        print(f"Nightsky reduction file already exists...\n{fpath}")
                    fpath_age = _Python_utilities.get_file_age(REDUCNSKYDIR / fpath.name)
                    if tryagain == False or fpath_age.days < file_age:
                        if verbose == True :
                            print("*"*10)
                            print(f"The file is younger than {file_age} days...\n{fpath}")
                    else :
                        if verbose == True :
                            print("*"*10)
                            print(f"The file is older than {file_age} days...\nTry gagin {fpath.name}")   
                        try:    
                            ccd = yfu.ccdred(
                                            ccd, 
                                            mflatpath = sMASTERDIR / f"nightskyflat-{filt}.fits",
                                            output = REDUCNSKYDIR / fpath.name
                                            # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                                        )
                        except Exception as err: 
                            if verbose == True :
                                print("Err :", err)
                            pass 
                else :
                    try:    
                        ccd = yfu.ccdred(
                                        ccd, 
                                        mflatpath = sMASTERDIR / f"nightskyflat-{filt}.fits",
                                        output = REDUCNSKYDIR / fpath.name
                                        # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                                    )
                    except Exception as err: 
                        if verbose == True :
                            print("Err :", err)
                        pass
    return 0


#%%
#%%    
# #########################################
# mag_inst
# #########################################
def mag_inst(flux, ferr):
    m_inst = -2.5 * np.log10(flux)
    merr   = 2.5/ np.log(10) * ferr / flux
    return m_inst, merr

def linf(x, a, b):
    return a + b*x


#%%
#########################################
#diff_Photometry_PS1
#########################################
def diff_Photometry_PS1 (DOINGDIR,
                        tryagain = False,
                        LOCATION = None,
                        SKYC_KW = None,
                        FWHM_INIT = 6,
                        Mag_Low = 11.5,
                        Mag_High = 15,
                        Mag_target = 12.5,
                        Mag_delta = 2,
                        ERR_Minimum = 0.5,
                        READINGDIR = reduced_dir,
                        # READINGDIR =  reduced_nightsky_dir,
                        file_age = 365,
                        verbose = False,
                        **kwarg
                        ) :
    from astropy.wcs import WCS
    from astropy.time import Time
    from astropy.nddata import Cutout2D
    from photutils.detection import DAOStarFinder
    from astropy.stats import sigma_clipped_stats

    from photutils.aperture import CircularAperture as CAp
    from photutils.aperture import CircularAnnulus as CAn
    from photutils.aperture import aperture_photometry as apphot
    import seaborn as sns
    from scipy.optimize import curve_fit

    FWHM = FWHM_INIT
    R_AP = 1.5 * FWHM_INIT # Aperture radius
    R_IN = 4 * FWHM_INIT   # Inner radius of annulus
    R_OUT = 6 * FWHM_INIT 
    
    DOINGDIR = Path(DOINGDIR)
    if verbose == True :
        print("DOINGDIR", DOINGDIR)

    READINGDIR = DOINGDIR / READINGDIR

    DIFFPRESULTDIR = DOINGDIR / f"{READINGDIR.parts[-1]}_DPhot_Mag{Mag_target}_fw{FWHM_INIT}"
    if not DIFFPRESULTDIR.exists():
        os.makedirs("{}".format(str(DIFFPRESULTDIR)))
        if verbose == True :
            print("{} is created...".format(str(DIFFPRESULTDIR)))

    summary = yfu.make_summary(READINGDIR/"*.fit*",
                                verify_fix=True,
                                ignore_missing_simple=True,
                                verbose = verbose,
                                )
    if summary is not None : 
        if verbose == True :
            print("len(summary):", len(summary))
            #print("summary:", summary)
            #print(summary["file"][0])
        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)
        if verbose == True :
            print("df_light:\n{}".format(df_light))

        for _, row  in df_light.iterrows():
            try:
                fpath = Path(row["file"])

                if (DIFFPRESULTDIR/f"{fpath.stem}_result_photometry.csv").exists() :
                    if verbose == True :
                        print(f"{fpath.stem}_result_photometry.csv is already exist...")
                    fpath_age = _Python_utilities.get_file_age(DIFFPRESULTDIR/f"{fpath.stem}_result_photometry.csv")
                    if verbose == True :
                        print("fpath_age :", fpath_age)
                    if (tryagain == False or fpath_age.days < file_age):
                        if verbose == True :
                            print("*"*10)
                            print(f"The file is younger than {file_age}...\n{fpath}")
                        return 0
                else : 
                    if verbose == True :
                        print(f"{fpath.stem}_result_photometry.csv is being reprocess...")
                        print("*"*20)
                        # print(f"Starting {fpath.name}...")
                        print(f"Starting {fpath}...")
                    hdul = fits.open(fpath)
                    ccd = yfu.load_ccd(fpath)
                    flt = hdul[0].header["filter"]

                    SOLVE, ASTAP, LOCAL = checkPSolve(fpath)
                    if verbose == True :
                        print(SOLVE, ASTAP, LOCAL)
                    
                    if SOLVE :
                        wcs = WCS(hdul[0].header)
                        # It is used as a rough estimate, so no need to be accurate:
                        #PIX2ARCSEC = 0.62*u.arcsec
                            
                        if hdul[0].header['CCDNAME'] == 'ASI6200MMPro' :
                            val_figsize = (12, 9)
                            val_fraction = 0.035
                            hdul[0].header["RDNOISE"] = 1.2
                        
                        if hdul[0].header['CCDNAME'] == 'STF-8300M' :
                            val_figsize = (12, 9)
                            val_fraction = 0.035
                            hdul[0].header["GAIN"] = 0.37,
                            hdul[0].header["RDNOISE"] = 9.3

                        if hdul[0].header['CCDNAME'] == 'STX-16803' :
                            val_figsize=(10, 9)
                            val_fraction = 0.0455
                            hdul[0].header["GAIN"] = 1.27
                            hdul[0].header["RDNOISE"] = 9.0    

                        # It is used as a rough estimate, so no need to be accurate:
                        if 'PIXSCALE' in hdul[0].header:
                            PIX2ARCSEC = hdul[0].header['PIXSCALE']
                        else : 
                            PIX2ARCSEC = calPixScale(hdul[0].header['FOCALLEN'], 
                                                            hdul[0].header['XPIXSZ'],
                                                            hdul[0].header['XBINNING'])
                        rdnoise = hdul[0].header["RDNOISE"]
                        if "EGAIN" in hdul[0].header :
                            gain    = hdul[0].header["EGAIN"]
                        elif "GAIN" in hdul[0].header :
                            gain    = hdul[0].header["GAIN"]
                        if verbose == True :
                            print(f"rdnoise : {rdnoise}, gain : {gain}, PIX2ARCSEC : {PIX2ARCSEC}")
                        
                        # D.2. Find the observation time and exposure time to set the obs time
                        t_start = Time(hdul[0].header['DATE-OBS'], format='isot')
                        t_expos = hdul[0].header['EXPTIME'] * u.s
                        t_middle = t_start + t_expos / 2 # start time + 0.5 * exposure time
                        if verbose == True :
                            print(f"t_start: {t_start}, t_expos: {t_expos}, t_middle: {t_middle}")
                        
                        # Get the radius of the smallest circle which encloses all the pixels
                        rad = yfu.fov_radius(header=hdul[0].header,
                                            unit=u.deg)
                        if verbose == True :
                            print("rad: {}".format(rad))  # 시야각(FOV)으로 구한 반지름

                        cent_coord = yfu.center_radec(ccd_or_header=hdul[0].header, 
                                                            center_of_image=True)
                        if verbose == True :
                            print("cent_coord: {}".format(cent_coord))

                        pos_sky = SkyCoord(cent_coord, unit='deg')
                        pos_pix = pos_sky.to_pixel(wcs=wcs)
                        if verbose == True :
                            print("pos_sky: {}".format(pos_sky))
                            print("pos_pix: {}".format(pos_pix))

                        #%%
                        ps1 = ypu.PanSTARRS1(cent_coord.ra, cent_coord.dec, radius=rad,
                        column_filters={"rmag":f"{Mag_target-Mag_delta}..{Mag_target+Mag_delta}",
                                                "e_rmag":"<0.10", "nr":">5"}
                                                )
                        PS1_stars_all = ps1.query()
                        if verbose == True :
                            print("type(PS1_stars_all) :", type(PS1_stars_all))
                            print("len(PS1_stars_all) :", len(PS1_stars_all))

                        df_stars_all = PS1_stars_all.to_pandas()

                        isnear = ypu.organize_ps1_and_isnear(
                                            ps1, 
                                            # header=ccd.header+ccd.wcs.to_header(), 
                                            ccd.header+ccd.wcs.to_header(), 
                                            # bezel=5*FWHM_INIT*PIX2ARCSEC.value,
                                            # nearby_obj_minsep=5*FWHM_INIT*PIX2ARCSEC.value,
                                            bezel=5*FWHM_INIT*PIX2ARCSEC,
                                            nearby_obj_minsep=5*FWHM_INIT*PIX2ARCSEC,
                                            group_crit_separation=6*FWHM_INIT
                                        )
                        df_stars = ps1.queried.to_pandas()
                        if verbose == True :
                            print("len(df_stars):", len(df_stars))
                        df_stars = df_stars.dropna(subset=["gmag", "rmag"])
                        if verbose == True :
                            print("len(df_stars):", len(df_stars))

                        pos_stars_all = np.array([df_stars_all["RAJ2000"].array, df_stars_all["DEJ2000"].array]).T
                        pos_stars_all = SkyCoord(pos_stars_all, **SKYC_KW).to_pixel(wcs)
                        pos_stars_all = np.transpose(pos_stars_all)
                        # pos_stars_all   # PS1 query 모든 별

                        pos_stars = np.array([df_stars["RAJ2000"].array, df_stars["DEJ2000"].array]).T
                        pos_stars = SkyCoord(pos_stars, **SKYC_KW).to_pixel(wcs)
                        pos_stars = np.transpose(pos_stars)
                        
                        # pos_stars     # PS1 query 중 비교 측광에 사용될 별

                        ap_stars = CAp(positions=pos_stars, r=R_IN)
                        ap_stars_all = CAp(positions=pos_stars_all, r=R_IN)
                        #apert
                        an_stars = CAn(positions=pos_stars, r_in=R_IN, r_out=R_OUT)
                        an_stars_all = CAn(positions=pos_stars_all, r_in=R_IN, r_out=R_OUT)
                        
                        #%%
                        fig, axs = plt.subplots(1, 1, figsize=val_figsize,
                        subplot_kw={'projection': wcs},
                        sharex=False, sharey=False, gridspec_kw=None)

                        im = zimshow(axs, hdul[0].data, )
                        axs.coords.grid(True, color='white', ls=':')
                        axs.coords['ra'].set_axislabel('Right Ascension (J2000)', minpad=0.5, fontsize=8)
                        axs.coords['ra'].set_ticklabel_position('bl')
                        axs.coords['dec'].set_axislabel('Declination (J2000)', minpad=0.4, fontsize=8)
                        axs.coords['dec'].set_ticklabel_position('bl')
                        axs.coords['ra'].set_major_formatter('hh:mm')
                        axs.coords['dec'].set_major_formatter('dd:mm')
                        axs.coords['ra'].display_minor_ticks(True)
                        axs.coords['dec'].display_minor_ticks(True)
                        axs.coords['ra'].set_minor_frequency(2)
                        axs.coords['dec'].set_minor_frequency(2)
                        axs.tick_params(labelsize=8)

                        for i in range(len(pos_stars)):
                            axs.text(pos_stars[i][0], pos_stars[i][1], f"Star #{str(i)}", fontsize=6, color='w')

                        ap_stars_all.plot(axs, color='w', lw=1)
                        ap_stars.plot(axs, color='r', lw=1)

                        axs.set_title(f"fname: {fpath.name}\n Comparison Stars of PS1 (red tag, Magnitude : {Mag_target}+-{Mag_delta})", fontsize=10,)

                        cbar = plt.colorbar(im, ax = axs, fraction=0.035, pad=0.04, )
                        cbar.ax.tick_params(labelsize=8)

                        axs.annotate(f'Number of star(s): {len(pos_stars)}', fontsize=8,
                            xy=(0, 0), xytext=(-10, -50), va='top', ha='left',
                            xycoords='axes fraction', textcoords='offset points')

                        plt.tight_layout()
                        plt.savefig(f"{DIFFPRESULTDIR/fpath.stem}_PS1_comparison.png")

                        # plt.show()
                        plt.close()

                        #%%
                        fig, axs = plt.subplots(1, 1, figsize=val_figsize,
                                                subplot_kw={'projection': wcs},
                                                sharex=False, sharey=False, gridspec_kw=None)

                        im = zimshow(axs, hdul[0].data, )

                        _phot_stars = []

                        for i, row in df_stars.iterrows():
                            pos_star = SkyCoord(row["RAJ2000"], row["DEJ2000"],
                                                **SKYC_KW).to_pixel(wcs)
                            ap = CAp([pos_star[0], pos_star[1]],
                                    r=R_AP)
                            an = CAn([pos_star[0], pos_star[1]],
                                    r_in=R_IN, r_out=R_OUT)
                            _phot_star = ypu.apphot_annulus(hdul[0].data,
                                                            ap, an,
                                                            error=yfu.errormap(hdul[0].data))
                            _phot_star[f"{flt}mag"] = row[f"{flt}mag"]
                            _phot_star[f"e_{flt}mag"] = row[f"e_{flt}mag"]
                            _phot_star["gmag"] = row["gmag"]
                            _phot_star["e_gmag"] = row["e_gmag"]
                            _phot_star["rmag"] = row["rmag"]
                            _phot_star["e_rmag"] = row["e_rmag"]
                            _phot_star["grcolor"] = row["grcolor"]
                            _phot_star["e_grcolor"] = row["e_grcolor"]
                            _phot_star["id"] = i
                            _phot_star["objID"] = int(row["objID"])
                            _phot_stars.append(_phot_star)
                            axs.text(pos_star[0]+10, pos_star[1]+10, f"star {i}:{row[f'{flt}mag']:.01f}",
                                    fontsize=8, color="w")
                            ap.plot(axs, color="orange")
                            # an.plot(axs, color="w")

                        axs.coords.grid(True, color='white', ls=':')
                        axs.coords['ra'].set_axislabel('Right Ascension (J2000)', minpad=0.5, fontsize=8)
                        axs.coords['ra'].set_ticklabel_position('bl')
                        axs.coords['dec'].set_axislabel('Declination (J2000)', minpad=0.4, fontsize=8)
                        axs.coords['dec'].set_ticklabel_position('bl')
                        axs.coords['ra'].set_major_formatter('hh:mm')
                        axs.coords['dec'].set_major_formatter('dd:mm')
                        axs.coords['ra'].display_minor_ticks(True)
                        axs.coords['dec'].display_minor_ticks(True)
                        axs.coords['ra'].set_minor_frequency(2)
                        axs.coords['dec'].set_minor_frequency(2)
                        axs.tick_params(labelsize=8)

                        cbar = plt.colorbar(im, ax = axs, fraction=0.035, pad=0.04, )

                        axs.set_title(f"fname: {fpath.name}\n {flt} magnitude of PS1 comparison stars (Magnitude : {Mag_target}+-{Mag_delta})", fontsize=10,)
                        axs.annotate(f'Number of star(s): {len(pos_stars)}', fontsize=8,
                                xy=(0, 0), xytext=(-10, -50), va='top', ha='left',
                                xycoords='axes fraction', textcoords='offset points')

                        plt.tight_layout()
                        plt.savefig(f"{DIFFPRESULTDIR}/{fpath.stem}_PS1_magnitude.png")

                        # plt.show()
                        plt.close()

                        #%%
                        df_phot_stars = pd.concat(_phot_stars)
                        df_phot_stars_na = df_phot_stars.dropna()
                        print(len(df_phot_stars_na))

                        df_phot_stars_na = df_phot_stars[df_phot_stars["merr"] < ERR_Minimum]
                        # phot_stars_na = phot_stars_na.set_index('id', drop=True)
                        df_phot_stars_na = df_phot_stars_na.reset_index(drop=True)
                        print(len(df_phot_stars_na))
                        # print(df_phot_stars_na)
                        df_phot_stars_na

                        #%%
                        merr_total1 = np.sqrt((df_phot_stars_na["merr"])**2 + (df_phot_stars_na[f"e_{flt}mag"])**2)

                        # === Calculate zero point and errors
                        _xx = np.linspace(Mag_Low, Mag_High)
                        zeropt_med = np.median(df_phot_stars_na["mag"] - df_phot_stars_na[f"{flt}mag"])
                        zeropt_avg = np.average(df_phot_stars_na["mag"] - df_phot_stars_na[f"{flt}mag"],
                                                weights=1/merr_total1**2)
                        dzeropt = np.max([1/np.sqrt(np.sum(1/(merr_total1)**2)),
                                        np.std((df_phot_stars_na[f"e_{flt}mag"] - df_phot_stars_na["merr"]), ddof=1)/np.sqrt(len(df_phot_stars_na[f"{flt}mag"]))])
                        merr_total2 = np.sqrt(np.sqrt(merr_total1**2 + dzeropt**2))

                        # === Find fitting lines
                        # Search for the usage of scipy.optimize.curve_fit.
                        poptm, _ = curve_fit(linf, df_phot_stars_na[f"{flt}mag"],
                                            df_phot_stars_na["mag"],
                                            sigma= df_phot_stars_na["merr"], absolute_sigma=True)
                        poptc, _ = curve_fit(linf, df_phot_stars_na["grcolor"],
                                            df_phot_stars_na["mag"] - df_phot_stars_na[f"{flt}mag"],
                                            sigma=merr_total2, absolute_sigma=True)

                        #%%
                        fig, axs = plt.subplots(2, 3, figsize=(15, 6), sharex=False, sharey=False,
                                        gridspec_kw={'height_ratios': [1, 3]})
                        
                        errkw = dict(marker="", ls="", ecolor="gray", elinewidth=0.5)

                        def plot_common(ax, x, y, xerr, yerr, title="", xlabel="", ylabel="", ylim=None):
                            ax.plot(x, y, '+')
                            ax.errorbar(x, y, xerr=xerr, yerr=yerr, **errkw)
                            ax.axhline(zeropt_med, color="r", lw=1, label=f"$Z = {{{zeropt_med:.3f}}} ± {{{dzeropt:.3f}}}$\n(median value)")
                            ax.axhline(zeropt_avg, color="b", lw=1, label=f"$Z = {{{zeropt_avg:.3f}}} ± {{{dzeropt:.3f}}}$\n(average value)")
                            ax.hlines([zeropt_med + dzeropt, zeropt_med - dzeropt, zeropt_avg + dzeropt, zeropt_avg - dzeropt],
                                    *ax.get_xlim(), color=["r","r","b","b"], lw=1, ls=":")
                            ax.set(title=title, xlabel=xlabel, ylabel=ylabel, ylim=ylim)
                            # ax.legend(fontsize=8, loc='best')

                        # 상단 행
                        plot_common(axs[0, 0], df_phot_stars_na[f"{flt}mag"], df_phot_stars_na["mag"] - df_phot_stars_na[f"{flt}mag"],
                                    df_phot_stars_na[f"e_{flt}mag"], df_phot_stars_na["merr"],
                                    ylabel=f"${{{flt}}}_{{inst}} - {{{flt}}}_{{PS1}}$",
                                    ylim=(zeropt_med-0.8, zeropt_med+0.8))

                        plot_common(axs[0, 1], df_phot_stars_na["grcolor"], df_phot_stars_na["mag"] - df_phot_stars_na[f"{flt}mag"],
                                    df_phot_stars_na[f"e_grcolor"], merr_total2,
                                    title=f"${{{flt}}}_{{inst}} - {{{flt}}}_{{PS1}} = (z + k'X) + (k''X + k)C$",
                                    ylabel=f"${{{flt}}}_{{inst}} - {{{flt}}}_{{PS1}}$",
                                    ylim=(zeropt_med-0.8, zeropt_med+0.8))
                        axs[0, 1].plot(axs[0, 1].get_xlim(), linf(np.array(axs[0, 1].get_xlim()), *poptc),
                                    "g-", lw=1, label=f"$y = {{{poptc[1]:+.3f}}}x {{{poptc[0]:+.3f}}}$\n(curve_fit)")
                        # axs[0, 1].legend(fontsize=8, loc='best')

                        data = df_phot_stars_na[["mag", "merr", "grcolor", "e_grcolor"]]
                        sns.heatmap(data.corr(), annot=True, cmap='coolwarm', vmin=-1, vmax=1, center=0, ax=axs[0, 2])
                        axs[0, 2].set(title='Correlation Heatmap')

                        # 하단 행
                        plot_common(axs[1, 0], df_phot_stars_na[f"{flt}mag"], df_phot_stars_na["mag"] - df_phot_stars_na[f"{flt}mag"],
                                    df_phot_stars_na[f"e_{flt}mag"], df_phot_stars_na["merr"],
                                    xlabel=f"${{{flt}}}_{{PS1}}$ (PS1 to {flt} filter by Tonry+2012)",
                                    ylabel=f"${{{flt}}}_{{inst}} - {{{flt}}}_{{PS1}}$")

                        plot_common(axs[1, 1], df_phot_stars_na["grcolor"], df_phot_stars_na["mag"] - df_phot_stars_na[f"{flt}mag"],
                                    df_phot_stars_na[f"e_grcolor"], merr_total2,
                                    xlabel="$g - r$ (PS1)",
                                    ylabel=f"${{{flt}}}_{{inst}} - {{{flt}}}_{{PS1}}$")
                        axs[1, 0].legend(fontsize=8, loc='best')

                        axs[1, 1].plot(axs[1, 1].get_xlim(), linf(np.array(axs[1, 1].get_xlim()), *poptc),
                                    "g-", lw=1, label=f"$y = {{{poptc[1]:+.3f}}}x {{{poptc[0]:+.3f}}}$\n(curve_fit)")
                        axs[1, 1].legend(fontsize=8, loc='best')

                        axs[1, 2].plot(_xx, _xx + zeropt_med,
                                    label=f"${{{flt}}}_{{inst}} = {{{flt}}}_{{PS1}}$ {zeropt_med:+.03f}\n(median vlaue)",
                                    color="r", lw=1, ls="-")
                        axs[1, 2].plot(axs[1, 2].get_xlim(), linf(np.array(axs[1, 2].get_xlim()), *poptm),
                                    "g-", lw=1, label=f"$y = {{{poptm[1]:.3f}}}x {{{poptm[0]:+.3f}}}$\n(curve_fit)")
                        axs[1, 2].plot(df_phot_stars_na[f"{flt}mag"], df_phot_stars_na["mag"], '+')
                        axs[1, 2].errorbar(df_phot_stars_na[f"{flt}mag"],
                                    df_phot_stars_na["mag"],
                                    xerr=df_phot_stars_na[f"e_{flt}mag"],
                                    yerr=df_phot_stars_na["merr"],
                                    **errkw)
                        axs[1, 2].set(
                                    xlabel=f"${{{flt}}}_{{PS1}}$ (PS1 to {flt} filter by Tonry+2012)",
                                    ylabel =f"${{{flt}}}_{{inst}}$",
                                )
                        axs[1, 2].legend(fontsize=8, loc='best')
                        axs[1, 2].axis('square')

                        # ID 텍스트 추가
                        for _, row in df_phot_stars_na.iterrows():
                            for i in range(2):
                                for j in range(2):
                                    axs[i, j].text(row[f"{flt}mag" if j == 0 else "grcolor"],
                                                row["mag"] - row[f"{flt}mag"], int(row["id"]), fontsize=8, clip_on=True)
                            axs[1, 2].text(row[f"{flt}mag"], row["mag"], int(row["id"]), fontsize=8, clip_on=True)

                        # x축 레이블 숨기기 (상단 행)
                        for ax in axs[0, :2]:
                            ax.tick_params(labelbottom=False)

                        plt.suptitle(f"fname: {fpath.name}\n PS1 check for differential photometry (Magnitude : {Mag_target}+-{Mag_delta})", fontsize=10)
                        plt.tight_layout()
                        plt.savefig(f"{DIFFPRESULTDIR}/{fpath.stem}_standardization_extended.png")

                        # plt.show()
                        plt.close()


                        #%%
                        FWHM = FWHM_INIT
                        avg, med, std = sigma_clipped_stats(hdul[0].data)  # by default, 3-sigma 5-iteration.
                        thresh = 5. * std

                        DAOfind = DAOStarFinder(
                                                fwhm = FWHM,
                                                threshold=thresh,   # In reality, FWHM must be measured a priori using, e.g., ``ginga``
                                                # sharplo=0.2, sharphi=1.0,   # default values 0.2 and 1.0
                                                # roundlo=-1.0, roundhi=1.0,  # default values -1 and +1
                                                # sigma_radius=1.5,           # default values 1.5
                                                # ratio=1.0,                  # 1.0: circular gaussian
                                                exclude_border=True         # To exclude sources near edges
                                                )

                        DAOfound = DAOfind(hdul[0].data)
                        if len(DAOfound) > 2000 :
                            from photutils import detect_threshold
                            thresh_snr = detect_threshold(data=hdul[0].data, nsigma=3,)
                            print('type(thresh_snr) :', type(thresh_snr))
                            print('thresh_snr.shape :', thresh_snr.shape)
                            print('detect_threshold', thresh_snr)
                            thresh = thresh_snr[0][0]

                            DAOfind = DAOStarFinder(
                                                fwhm = FWHM,
                                                threshold=thresh,   # In reality, FWHM must be measured a priori using, e.g., ``ginga``
                                                # sharplo=0.2, sharphi=1.0,   # default values 0.2 and 1.0
                                                # roundlo=-1.0, roundhi=1.0,  # default values -1 and +1
                                                # sigma_radius=1.5,           # default values 1.5
                                                # ratio=1.0,                  # 1.0: circular gaussian
                                                exclude_border=True         # To exclude sources near edges
                                                )
                            DAOfound = DAOfind(hdul[0].data)

                        print("len(DAOfound) :",len(DAOfound))
                        print(DAOfound.colnames)

                        # DAOfound.write(f"{DIFFPRESULTDIR/fpath.stem}_DAOStarfinder_fwhm_{FWHM}.csv",
                        #                             overwrite = True,
                        #                             format='ascii.fast_csv')
                        df_DAO = DAOfound.to_pandas()
                        print(type(df_DAO))
                        df_DAO

                        pos = np.transpose((DAOfound['xcentroid'], DAOfound['ycentroid']))
                        apert = CAp(pos, r=R_AP)
                        annul = CAn(positions=pos, r_in= R_IN, r_out=R_OUT)
                        #%%
                        fig, axs = plt.subplots(1, 1, figsize=val_figsize,
                        subplot_kw={'projection': wcs},
                        sharex=False, sharey=False, gridspec_kw=None)

                        im = zimshow(axs, hdul[0].data, )
                        axs.set_title('World coordinate system', fontsize=9)
                        axs.coords.grid(True, color='white', ls=':')
                        axs.coords['ra'].set_axislabel('Right Ascension (J2000)', minpad=0.5, fontsize=8)
                        axs.coords['ra'].set_ticklabel_position('bl')
                        axs.coords['dec'].set_axislabel('Declination (J2000)', minpad=0.4, fontsize=8)
                        axs.coords['dec'].set_ticklabel_position('bl')
                        axs.coords['ra'].set_major_formatter('hh:mm')
                        axs.coords['dec'].set_major_formatter('dd:mm')
                        axs.coords['ra'].display_minor_ticks(True)
                        axs.coords['dec'].display_minor_ticks(True)
                        axs.coords['ra'].set_minor_frequency(2)
                        axs.coords['dec'].set_minor_frequency(2)
                        axs.tick_params(labelsize=8)

                        annul.plot(axs, color="r")
                        for i in range(len(pos)):
                            axs.text(pos[i][0], pos[i][1], f"Star #{str(i)}", fontsize=6, color='w')

                        annul.plot(axs, color="r")

                        cbar = plt.colorbar(im, ax = axs, fraction=0.035, pad=0.04, )
                        cbar.ax.tick_params(labelsize=8)

                        axs.set_title(f"fname: {fpath.name}\n Result of DAOFinder", fontsize=10,)

                        axs.annotate(f'FWHM: {FWHM}', fontsize=8,
                            xy=(0, 0), xytext=(-10, -30), va='top', ha='left',
                            xycoords='axes fraction', textcoords='offset points')

                        axs.annotate(f'Sky threshold: {thresh:.02f}', fontsize=8,
                            xy=(0, 0), xytext=(-10, -40), va='top', ha='left',
                            xycoords='axes fraction', textcoords='offset points')

                        axs.annotate(f'Number of star(s): {len(DAOfound)}', fontsize=8,
                            xy=(0, 0), xytext=(-10, -50), va='top', ha='left',
                            xycoords='axes fraction', textcoords='offset points')

                        plt.tight_layout()
                        plt.savefig(f"{DIFFPRESULTDIR}/{fpath.stem}_DAOStarfinder_fwhm_{FWHM}.png")

                        # plt.show()
                        plt.close()

                        #%%
                        apphot_result = apphot(hdul[0].data, apert, method='center')
                        print(type(apphot_result))
                        # df_apphot = pd.DataFrame()
                        # apphot_result
                        df_apphot = apphot_result.to_pandas()
                        print(type(df_apphot))
                        df_apphot

                        ap_area  = apert.area
                        ap_area

                        # since our `annul` has many elements,
                        mask_apert = (apert.to_mask(method='center'))
                        mask_annul = (annul.to_mask(method='center'))

                        mag_ann  = np.zeros(len(apphot_result))
                        merr_ann = np.zeros(len(apphot_result))

                        #%%
                        for i in range(len(apphot_result)):
                            annul_weighted = mask_annul[i].multiply(hdul[0].data)
                            sky_non0   = np.nonzero(annul_weighted)
                            sky_pixel  = annul_weighted[sky_non0]

                            msky, sky_std, nsky, nrej = sky_fit(sky_pixel, method='mode',
                                                                                mode_option='sex')


                            flux_star = apphot_result['aperture_sum'][i] - msky * ap_area  # total - sky

                            flux_err  = np.sqrt(apphot_result['aperture_sum'][i] * gain    # Poissonian (star + sky)
                                                + ap_area * rdnoise**2 # Gaussian
                                                + (ap_area * (gain * sky_std))**2 / nsky )

                            mag_ann[i], merr_ann[i] = mag_inst(flux_star, flux_err)
                            df_apphot.at[i, 'msky'] = msky
                            df_apphot.at[i, 'sky_std'] = sky_std
                            df_apphot.at[i, 'nsky'] = nsky
                            df_apphot.at[i, 'nrej'] = nrej
                            df_apphot.at[i, 'flux_star'] = flux_star
                            df_apphot.at[i, 'flux_err'] = flux_err
                            df_apphot.at[i, 'mag_ann'] = mag_ann[i]
                            df_apphot.at[i, 'merr_ann'] = merr_ann[i]

                        df_apphot['filename'] = fpath.stem
                        df_apphot['t_start'] = t_start
                        df_apphot['t_expos'] = t_expos
                        df_apphot['t_middle'] = t_middle
                        df_apphot['filter'] = flt
                        df_apphot["zeropt_med"] = zeropt_med
                        df_apphot["zeropt_avg"] = zeropt_avg
                        df_apphot["e_zeropt"] = dzeropt

                        df_apphot[f"{flt}_magnitude"] = df_apphot["mag_ann"] - df_apphot["zeropt_med"]

                        df_apphot['filename'] = fpath.stem
                        df_apphot['t_start'] = t_start
                        df_apphot['t_expos'] = t_expos
                        df_apphot['t_middle'] = t_middle
                        df_apphot['filter'] = flt
                        df_apphot["zeropt_med"] = zeropt_med
                        df_apphot["zeropt_avg"] = zeropt_avg
                        df_apphot["e_zeropt"] = dzeropt

                        df_apphot[f"{flt}_magnitude"] = df_apphot["mag_ann"] - df_apphot["zeropt_med"]

                        sky_coord = wcs.pixel_to_world(df_apphot['xcenter'], df_apphot['ycenter'])
                        sky_coord
                        print(type(sky_coord))

                        # df_apphot["RA2000"] = sky_coord.ra
                        # df_apphot["RA2000"]
                        df_RADEC = pd.DataFrame({"RAJ2000": sky_coord.ra.degree, "DEJ2000": sky_coord.dec.degree})
                        # df_RADEC
                        #type(df_RADEC["RA2000"][0])
                        df_apphot = pd.concat([df_apphot, df_RADEC], axis=1,)

                        df_apphot.to_csv(f"{DIFFPRESULTDIR}/{fpath.stem}_result_photometry.csv")
                        
                        df_apphot_sub = df_apphot.dropna()
                        print(len(df_apphot_sub))
                        df_apphot_sub = df_apphot_sub.loc[(df_apphot_sub["merr_ann"] < ERR_Minimum)]
                        df_apphot_sub

                        #%%
                        fig, axs = plt.subplots(2, 2, figsize=(10, 8),
                                                sharex=False, sharey=False, gridspec_kw=None)

                        for idx, row in df_apphot_sub.iterrows():
                            im0 = axs[0, 0].errorbar(df_apphot_sub["id"],
                                        df_apphot_sub[f"{flt}_magnitude"], yerr=df_apphot_sub["merr_ann"],
                                        marker='x',
                                        ls='none',
                                        #ms=10,
                                        capsize=3)

                        axs[0, 0].invert_yaxis()
                        axs[0, 0].set(
                            xlabel='Star ID',
                            ylabel=f"${{{flt}}}_{{obs}}$"
                            )

                        style = {'edgecolor': 'white', 'linewidth': 3}
                        im1 = axs[0, 1].hist(df_apphot_sub[f"{flt}_magnitude"],
                                    **style)
                        axs[0, 1].set(
                            xlabel=f"${{{flt}}}_{{obs}}$",
                            ylabel="number of stars"
                            )

                        # 상관관계 계산
                        data =  df_apphot_sub[[f"{flt}_magnitude", "merr_ann"]]
                        corr = data.corr()

                        # 히트맵 그리기
                        im2 = sns.heatmap(corr, annot=True, cmap='coolwarm',
                                            vmin=-1, vmax=1, center=0, ax = axs[1, 0])
                        axs[1, 0].set(
                            title = 'Correlation Heatmap',
                            )

                        axs[1, 1].scatter(df_apphot_sub[f"{flt}_magnitude"], df_apphot_sub["merr_ann"], marker='x',)
                        axs[1, 1].errorbar(x=df_apphot_sub[f"{flt}_magnitude"], y=df_apphot_sub["merr_ann"],
                                    yerr=None, xerr=df_apphot_sub["merr_ann"], fmt="o", color="gray", capsize=3, alpha=0.5)
                        axs[1, 1].set(
                            title = "Correlation between Magnitude and Error",
                            xlabel=f"${{{flt}}}_{{obs}}$",
                            ylabel="Error",
                            )

                        plt.suptitle(f"fname: {fpath.name}\n Result of differential photometry (Magnitude : {Mag_target}+-{Mag_delta})", fontsize=10,)

                        plt.tight_layout()
                        plt.savefig(f"{DIFFPRESULTDIR}/{fpath.stem}_Result_of_differential_photometry.png")

                        # plt.show()
                        plt.close()

            except Exception as err: 
                if verbose == True :
                    print("Err :", err)
                    pass 
    return 0

#%%
#########################################
#plot_light_curve_variables_using_csv
#########################################
def plot_light_curve_variables_using_csv(DOINGDIR,
                                        READINGDIR = reduced_dir,
                                        FWHM_INIT = 6,
                                        Mag_target = 12.5,
                                        ERR_Minimum = 0.5,
                                        coord_deltas = np.arange(0.00001, 0.00050, 0.00001),
                                        #  READINGDIR = _astro_utilities.reduced__nightsky_dir,
                                        verbose = False,
                                        **kwarg
                                        ) :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)

    READINGDIR = DOINGDIR / READINGDIR
    # READINGDIR = DOINGDIR / _astro_utilities.reduced_nightsky_dir

    DIFFPRESULTDIR = DOINGDIR / f"{READINGDIR.parts[-1]}_DPhot_Mag{Mag_target}"
    DIFFPRESULTDIR = DOINGDIR / f"{READINGDIR.parts[-1]}_DPhot_Mag{Mag_target}_fw{FWHM_INIT}"
    LIGHTCUEVEDIR = DOINGDIR / "LightCurve"
    if not LIGHTCUEVEDIR .exists():
        os.makedirs("{}".format(str(LIGHTCUEVEDIR )))
        print("{} is created...".format(str(LIGHTCUEVEDIR )))

    csv_in_dir = sorted(list(DIFFPRESULTDIR.glob('*result_photometry.csv')))
    if len(csv_in_dir) == 0 : 
        print("len(csv_in_dir):", len(csv_in_dir))
    else : 
        try: 
            df = pd.DataFrame()
            for fpath in csv_in_dir[:]:
                fpath = Path(fpath)
                print(f"starting... {fpath}")
                df_csv = pd.read_csv(fpath)
                df = pd.concat([df, df_csv], axis=0)

            print("len(df):", len(df))
            print("df.columns:", df.columns)

            df['t_middle_dt'] = pd.to_datetime(df['t_middle'])
            df = df.drop(columns=['Unnamed: 0'], axis=0)
            df = df.reset_index(drop=True)
            targ_name = DOINGDIR.parts[-1].split("_")[0]
            targ_name = targ_name.replace("-"," ")
            print("targ_name :", targ_name)

            from astroquery.simbad import Simbad
            result_table = Simbad.query_object(targ_name)
            
            if not result_table :
                print("there is no result...")
            else : 
                # print(result_table.columns)
                print("result_table :", result_table)

                # type(result_table['RA'][0])
                # result_table['RA'][0].split(" ")
                targ_sky = SkyCoord(ra=result_table['RA'][0],
                    dec=result_table['DEC'][0], 
                    unit=(u.hourangle, u.degree),
                    frame='icrs')
                for coord_delta in coord_deltas :
                    df_targ = df.loc[(df["RAJ2000"] > targ_sky.ra.value*(1-coord_delta)) \
                                    & (df["RAJ2000"] < targ_sky.ra.value*(1+coord_delta)) \
                                    & (df["DEJ2000"] > targ_sky.dec.value*(1-coord_delta))\
                                    & (df["DEJ2000"] < targ_sky.dec.value*(1+coord_delta))\
                                    & (df["merr_ann"] < ERR_Minimum)]

                    if df_targ.empty :
                        print("df_targ is empty")
                    else : 
                        df_targ.to_csv(f"{LIGHTCUEVEDIR}/{READINGDIR.parts[-2]}_{READINGDIR.parts[-1]}_DPhot_Mag{Mag_target}_fw{FWHM_INIT}_light_curve_{coord_delta:.05f}.csv")
                        print(f"{READINGDIR.parts[-2]}_{READINGDIR.parts[-1]}_DPhot_Mag{Mag_target}_fw{FWHM_INIT}_light_curve_{coord_delta:.05f}.csv is created...")

                        check_filters = df_targ['filter'].drop_duplicates()
                        check_filters.reset_index(drop=True)
                        # print("check_filters)", check_filters)

                        ttime = Time(df_targ["t_middle_dt"])

                        fig, axs = plt.subplots(2, 1, figsize=(12, 8), 
                                sharex=False, sharey=False, gridspec_kw=None)

                        for chl in check_filters :
                            # if f'{chl}_magnitude' in df_targ:    
                            df_targ_chl = df_targ.loc[df_targ["filter"] == chl].copy()
                            if not df_targ_chl.empty :
                                # print(df_targ_chl)
                                # ttime = Time(df_targ_chl["t_middle_dt"])
                                im0 = axs[0].errorbar(Time(df_targ_chl["t_middle_dt"]).mjd, 
                                        df_targ_chl[f'{chl}_magnitude'], yerr=abs(df_targ_chl["merr_ann"]),
                                        marker='x',
                                        ls='none',
                                        #ms=10,
                                        capsize=3,
                                        label=f'{chl}_magnitude')
                                axs[1].errorbar(df_targ_chl["t_middle_dt"], 
                                        df_targ_chl['flux_star'], yerr=abs(df_targ_chl["flux_err"]),
                                        marker='x',
                                        ls='none',
                                        #ms=10,
                                        capsize=3,
                                        label=f'flux_star at {chl}')

                        axs[0].invert_yaxis()

                        axs[0].set(
                            xlabel='Time (MJD)',
                            ylabel="Magnitude",
                            # ylim=(10.8+1, 10.8-1),
                            # ylim=(11.25+1.2, 11.25-1.2),   
                            # ylim=(10.75+.6, 10.75-.6),   
                            # ylim=(10.8+.9, 10.8-.9), 
                        )
                        axs[0].legend()
                        axs[0].grid(linestyle=':')

                        axs[0].set_title(f"light curve of {targ_name}", fontsize=12,)
                        axs[0].annotate(f'Coord: {targ_sky} +-{coord_delta}', fontsize=8,
                                    xy=(0, 0), xytext=(-10, -30), va='top', ha='left',
                                    xycoords='axes fraction', textcoords='offset points')

                        axs[1].set(
                            xlabel='Time (date)',
                            ylabel="flux",
                        ) 
                        axs[1].legend()
                        axs[1].grid(linestyle=':')

                        axs[1].set_title(f"light curve of {targ_name}", fontsize=12,)
                        axs[1].annotate(f'Coord: {targ_sky} +-{coord_delta}', fontsize=8,
                                    xy=(0, 0), xytext=(-10, -30), va='top', ha='left',
                                    xycoords='axes fraction', textcoords='offset points')

                        plt.tight_layout()
                        plt.savefig(f"{LIGHTCUEVEDIR}/{READINGDIR.parts[-2]}_{READINGDIR.parts[-1]}_DPhot_Mag{Mag_target}_fw{FWHM_INIT}_light_curve_{coord_delta:.05f}.png")

                        # plt.show()
                        # plt.close()
        except Exception as err: 
            if verbose == True :
                print("Err :", err)
                pass 
        
    return 0

#%%
#########################################
#plot_light_curve_Asteroids_using_csv
#########################################
def plot_light_curve_asteroids_using_csv(DOINGDIR,
                                        READINGDIR = reduced_dir,
                                        FWHM_INIT = 6,
                                        Mag_target = 12.5,
                                        ERR_Minimum = 0.5,
                                        coord_deltas = np.arange(0.00001, 0.00050, 0.00001),
                                        # READINGDIR = reduced__nightsky_dir,
                                        LOCATION = dict(lon=127.005, lat=37.308889, elevation=101),
                                        verbose = False,
                                        **kwarg
                                        ) :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)

    READINGDIR = DOINGDIR / READINGDIR
    # READINGDIR = DOINGDIR / _astro_utilities.reduced_nightsky_dir

    # DIFFPRESULTDIR = DOINGDIR / f"{READINGDIR.parts[-1]}_DPhot_Mag{Mag_target}"
    DIFFPRESULTDIR = DOINGDIR / f"{READINGDIR.parts[-1]}_DPhot_Mag{Mag_target}_fw{FWHM_INIT}"
    LIGHTCUEVEDIR = DOINGDIR / "LightCurve"
    if not LIGHTCUEVEDIR .exists():
        os.makedirs("{}".format(str(LIGHTCUEVEDIR )))
        print("{} is created...".format(str(LIGHTCUEVEDIR )))

    csv_in_dir = sorted(list(DIFFPRESULTDIR.glob('*result_photometry.csv')))
    if len(csv_in_dir) == 0 : 
        print("len(csv_in_dir):", len(csv_in_dir))
    else : 
        df = pd.DataFrame()
        for fpath in csv_in_dir[:]:
            fpath = Path(fpath)
            print(f"starting... {fpath}")
            df_csv = pd.read_csv(fpath)
            df = pd.concat([df, df_csv], axis=0)

        print("len(df):", len(df))
        print("df.columns:", df.columns)

        df['t_middle_dt'] = pd.to_datetime(df['t_middle'])
        df = df.drop(columns=['Unnamed: 0'], axis=0)
        df = df.reset_index(drop=True)

        targ_name = DOINGDIR.parts[-1].split("_")[0]
        targ_name = targ_name.replace("-"," ")
        targ_name = ''.join([i for i in targ_name  if not i.isdigit()])
        print("targ_name :", targ_name)

        check_ttimes = df[['t_middle_dt']].drop_duplicates()
        check_ttimes = check_ttimes.reset_index(drop=True)
        check_ttimes

        try : 
            df_targ = pd.DataFrame()
            for idx, row in check_ttimes.iterrows() :
                print(idx, row)
                targ_ttime = Time(row['t_middle_dt'])

                obj = Horizons(id=targ_name, location=LOCATION, epochs=targ_ttime.jd)
                result_table = obj.ephemerides()
                print("result_table : {}".format(result_table ))

                pos_sky = SkyCoord(result_table ["RA"][0], result_table ["DEC"][0], unit='deg')
                print("pos_sky: {}".format(pos_sky))

                if not result_table :
                    print("there is no result...")
                else : 
                    # print(result_table.columns)
                    print("result_table :", result_table)

                    # targ_sky = SkyCoord(ra=result_table['RA'][0],
                    #                     dec=result_table['DEC'][0], 
                    #                     unit=(u.hourangle, u.degree),
                    #                     frame='icrs')
                    # print("targ_sky :", targ_sky)       
                    targ_sky = pos_sky
                    print("targ_sky :", targ_sky)

                    for coord_delta in coord_deltas :
                        df_one = df.loc[(df["RAJ2000"] > targ_sky.ra.value*(1-coord_delta)) \
                                        & (df["RAJ2000"] < targ_sky.ra.value*(1+coord_delta)) \
                                        & (df["DEJ2000"] > targ_sky.dec.value*(1-coord_delta))\
                                        & (df["DEJ2000"] < targ_sky.dec.value*(1+coord_delta))\
                                        & (df['t_middle_dt'] == row['t_middle_dt'])]
                        print("df_one :", df_one)
                        df_targ = pd.concat([df_targ, df_one], axis=0)

                        if df_targ.empty :
                            print("df_targ is empty")
                        else : 
                            df_targ.to_csv(f"{LIGHTCUEVEDIR}/{READINGDIR.parts[-2]}_{READINGDIR.parts[-1]}_light_curve_{coord_delta}.csv")
                            print(f"{READINGDIR.parts[-2]}_{READINGDIR.parts[-1]}_light_curve_{coord_delta}.csv is created...")
                            
                            check_filters = df_targ['filter'].drop_duplicates()
                            check_filters.reset_index(drop=True)
                            # print("check_filters)", check_filters)

                            ttime = Time(df_targ["t_middle_dt"])

                            fig, axs = plt.subplots(2, 1, figsize=(12, 8), 
                                    sharex=False, sharey=False, gridspec_kw=None)

                            for chl in check_filters :
                                # if f'{chl}_magnitude' in df_targ:    
                                df_targ_chl = df_targ.loc[df_targ["filter"] == chl].copy()
                                if not df_targ_chl.empty :
                                    # print(df_targ_chl)
                                    # ttime = Time(df_targ_chl["t_middle_dt"])
                                    im0 = axs[0].errorbar(Time(df_targ_chl["t_middle_dt"]).mjd, 
                                            df_targ_chl[f'{chl}_magnitude'], yerr=abs(df_targ_chl["merr_ann"]),
                                            marker='x',
                                            ls='none',
                                            #ms=10,
                                            capsize=3,
                                            label=f'{chl}_magnitude')
                                    axs[1].errorbar(df_targ_chl["t_middle_dt"], 
                                            df_targ_chl['flux_star'], yerr=abs(df_targ_chl["flux_err"]),
                                            marker='x',
                                            ls='none',
                                            #ms=10,
                                            capsize=3,
                                            label=f'flux_star at {chl}')

                            axs[0].invert_yaxis()

                            axs[0].set(
                                xlabel='Time (MJD)',
                                ylabel="Magnitude",
                                # ylim=(10.8+1, 10.8-1),
                                # ylim=(11.25+1.2, 11.25-1.2),   
                                # ylim=(10.75+.6, 10.75-.6),   
                                # ylim=(10.8+.9, 10.8-.9), 
                            )
                            axs[0].legend()
                            axs[0].grid(linestyle=':')

                            axs[0].set_title(f"light curve of {targ_name}", fontsize=12,)
                            axs[0].annotate(f'Coord: {targ_sky} +-{coord_delta}', fontsize=8,
                                        xy=(0, 0), xytext=(-10, -30), va='top', ha='left',
                                        xycoords='axes fraction', textcoords='offset points')

                            axs[1].set(
                                xlabel='Time (date)',
                                ylabel="flux",
                            ) 
                            axs[1].legend()
                            axs[1].grid(linestyle=':')

                            axs[1].set_title(f"light curve of {targ_name}", fontsize=12,)
                            axs[1].annotate(f'Coord: {targ_sky} +-{coord_delta}', fontsize=8,
                                        xy=(0, 0), xytext=(-10, -30), va='top', ha='left',
                                        xycoords='axes fraction', textcoords='offset points')

                            plt.tight_layout()
                            plt.savefig(f"{LIGHTCUEVEDIR}/{READINGDIR.parts[-2]}_{READINGDIR.parts[-1]}_DPhot_Mag{Mag_target}_fw{FWHM_INIT}_light_curve_{coord_delta:.05f}.png")

                            # plt.show()
                            # plt.close()
        except Exception as err: 
            if verbose == True :
                print("Err :", err)
                pass
            
#%%
#########################################
#checkAsteroids
#########################################
def checkAsteroids(DOINGDIR,
        summary,
        ):
    ''' Gives the new path of the FITS file from header.
    Parameters
    ----------
    DOINGDIR: pathlike
        The path to the original .
    summary : dataframe
        
    Returns
    -------
    
    '''
    #####################################################################
    # Our object (will be queried to JPL HORIZONS)
    #OBJID = '2159' # 

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
    FWHM_INIT = 6

    # Photometry parameters
    R_AP = 1.5*FWHM_INIT # Aperture radius
    R_IN = 4*FWHM_INIT   # Inner radius of annulus
    R_OUT = 6*FWHM_INIT  # Outer radius of annulus

    Mag_UP = 17
    #######################################################

    ASTRESULTDIR = DOINGDIR / Asteroid_result_dir
    if not ASTRESULTDIR.exists():
        os.makedirs("{}".format(str(ASTRESULTDIR)))
        print("{} is created...".format(str(ASTRESULTDIR)))

    DOINGDIR = DOINGDIR / reduced_nightsky_dir

    summary_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
    summary_light = summary_light.reset_index(drop=True) 

    df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
    df_light = df_light.reset_index(drop=True)
    print("df_light:\n{}".format(df_light))

    for _, row  in df_light.iterrows():
        fpath = Path(row["file"])
        fpath_age = _Python_utilities.get_file_age(fpath)
        hdul = fits.open(fpath)

        SOLVE, ASTAP, LOCAL = checkPSolve(fpath)
        print(SOLVE, ASTAP, LOCAL)
        
        if SOLVE :
            wcs = WCS(hdul[0].header)
            # It is used as a rough estimate, so no need to be accurate:
            #PIX2ARCSEC = 0.62*u.arcsec
            if 'PIXSCALE' in hdul[0].header:
                PIX2ARCSEC = hdul[0].header['PIXSCALE']
            else : 
                PIX2ARCSEC = calPixScale(hdul[0].header['FOCALLEN'], 
                                                hdul[0].header['XPIXSZ'],
                                                hdul[0].header['XBINNING'])

            # D.2. Find the observation time and exposure time to set the obs time
            t_start = Time(hdul[0].header['DATE-OBS'], format='isot')
            t_expos = hdul[0].header['EXPTIME'] * u.s
            t_middle = t_start + t_expos / 2 # start time + 0.5 * exposure time
            #print(f"t_start: {t_start}, t_expos: {t_expos}, t_middle: {t_middle}")
            
            cent_coord = yfu.center_radec(ccd_or_header=hdul[0].header, 
                                        center_of_image=True)
            results_ast = Skybot.cone_search(cent_coord, 
                                            50*u.arcmin, 
                                            t_middle)
            #print(results_ast.pprint(max_width=80) )

            offset_RA = (cent_coord.ra.to(u.deg) - hdul[0].header['RA']*u.deg).to(u.arcmin)
            offset_DEC = (cent_coord.dec.to(u.deg) - hdul[0].header['DEC']*u.deg).to(u.arcmin) 
            altaz = AltAz(obstime=t_middle, location=Suwon)   
            cent_aa = cent_coord.transform_to(altaz)
            offset_AZ = (cent_aa.az.to(u.deg) - hdul[0].header['CENTAZ']*u.deg).to(u.arcmin)
            offset_ALT = (cent_aa.alt.to(u.deg) - hdul[0].header['CENTALT']*u.deg).to(u.arcmin)
            
            df_ast = results_ast.to_pandas()
            df_ast

            df_targ = df_ast[df_ast['V'] < Mag_UP]
            df_targ = df_targ.sort_values(by=['V'])
            df_targ = df_targ.reset_index(drop=True)
            df_targ

            if df_targ.empty:
                pass
            else:
                df_targ_eph = pd.DataFrame()

                for i, row in df_targ.iterrows() :
                    try : 
                        #print("type(row)", type(row))
                        #Query the ephemerides of this target! 
                        obj = Horizons(id=row['Number'], 
                                    location=observatory_code, 
                                    epochs=t_middle.jd)
                        obj_ephem = obj.ephemerides()
                        #print(obj_ephem)
                        df_eph = obj_ephem.to_pandas()
                        df_targ_eph = pd.concat([df_targ_eph, df_eph], axis = 0)
                    except : 
                        continue

                #print(df_targ_eph)
                df_targ_eph = df_targ_eph.reset_index(drop=True)
                df_targ_eph = pd.concat([df_targ, df_targ_eph], axis = 1)
                print("df_targ_eph :", df_targ_eph)

                duplicated_columns_list = []
                list_of_all_columns = list(df_targ_eph.columns)
                for column in list_of_all_columns:
                    if list_of_all_columns.count(column) > 1 and not column in duplicated_columns_list:
                        duplicated_columns_list.append(column)
                duplicated_columns_list

                for column in duplicated_columns_list:
                    list_of_all_columns[list_of_all_columns.index(column)] = column
                    list_of_all_columns[list_of_all_columns.index(column)] = column + '_1'

                df_targ_eph.columns = list_of_all_columns
                #print(df_targ_eph.columns)
                df_targ_eph.to_csv(f"{ASTRESULTDIR/fpath.stem}_AST_Mag{Mag_UP}.csv")
                df_targ_eph.dropna(subset = ['RA', 'DEC', 'V', 'RA_1', 'DEC_1', 'V_1'], inplace=True)
                df_targ_eph[['RA', 'DEC', 'V', 'RA_1', 'DEC_1', 'V_1']]
                
                # RADEC_targ = np.array([df_targ_eph['RA'], df_targ_eph["DEC"]]).T
                # RADEC_targ
                # pos_targ_init = SkyCoord(RADEC_targ, 
                #         **SKYC_KW).to_pixel(wcs, origin=1, mode='wcs')
                # print("pos_targ_init:", pos_targ_init)

                pos_targ_init = SkyCoord(df_targ_eph['RA']*u.deg, df_targ_eph["DEC"]*u.deg, 
                                        **SKYC_KW).to_pixel(wcs, origin=0, mode='wcs')
                pos_targ_init = np.array(pos_targ_init).T
                print("pos_targ_init:", pos_targ_init)

            if hdul[0].header['CCDNAME'] == 'STF-8300M' :
                val_figsize = (13, 5.2)
                val_fraction = 0.035
            if hdul[0].header['CCDNAME'] == 'STX-16803' :
                val_figsize=(12, 6.2)
                val_fraction = 0.0455

            fig_set = plt.figure(figsize=val_figsize)
            ax1 = plt.subplot2grid((1,2), (0,0),
                                fig=fig_set)
            im1 = zimshow(ax1, hdul[0].data, )
            ax1.set_title('Pixel coordinate system', fontsize=9)
            ax1.tick_params(labelsize=8)
            plt.colorbar(im1, ax = ax1, fraction=val_fraction, pad=0.04)

            ax2 = plt.subplot2grid((1,2), (0,1),
                                projection=wcs,
                                fig=fig_set)
            im2 = zimshow(ax2, hdul[0].data, )
            ax2.set_title('World coordinate system', fontsize=9)
            ax2.coords.grid(True, color='white', ls=':')
            ax2.coords['ra'].set_axislabel('Right Ascension (J2000)', minpad=0.5, fontsize=8)
            ax2.coords['ra'].set_ticklabel_position('bl')
            ax2.coords['dec'].set_axislabel('Declination (J2000)', minpad=0.4, fontsize=8)
            ax2.coords['dec'].set_ticklabel_position('bl')
            ax2.coords['ra'].set_major_formatter('hh:mm')
            ax2.coords['dec'].set_major_formatter('dd:mm')
            ax2.coords['ra'].display_minor_ticks(True)
            ax2.coords['dec'].display_minor_ticks(True)
            ax2.coords['ra'].set_minor_frequency(1)
            ax2.coords['dec'].set_minor_frequency(1)
            ax2.tick_params(labelsize=8)

            if df_targ.empty:
                pass
            else:
                targ_ap = CAp(pos_targ_init,
                        r=R_AP, 
                        )
                targ_an = CAn(pos_targ_init,
                        r_in=R_IN,
                        r_out=R_OUT)
                
                #targ_ap.plot(ax1, color="r")
                targ_an.plot(ax1, color="r")
                #targ_ap.plot(ax2, color="r")
                targ_an.plot(ax2, color="r")

                ax1.annotate(f"{pos_targ_init}",
                        xy=(0, 0), xytext=(0, -0.1),
                        xycoords='axes fraction',
                        va='top', ha='left',
                        fontsize=7)

                ax2.annotate(f"{df_targ_eph[['Number', 'RA', 'DEC', 'V']]}",
                        xy=(0, 0), xytext=(0, -0.1),
                        xycoords='axes fraction',
                        va='top', ha='left',
                        fontsize = 6)
            plt.colorbar(im2, ax = ax2, fraction=val_fraction, pad=0.04)
            plt.suptitle(f"fname: {fpath.name}")

            ax2.annotate(f"image center (RA, DEC): ({cent_coord.ra:.03f}, {cent_coord.dec:.03f})\ntelescope center (RA, DEC): ({hdul[0].header['RA']*u.deg:.03f}, {hdul[0].header['RA']*u.deg:.03f})\noffset (RA, DEC): ({offset_RA:.03f}, {offset_DEC:.03f})\noffset (AZ, ALT): ({offset_AZ:.03f}, {offset_ALT:.03f})",
                        xy=(0, 0), xytext=(0.6, -0.1),
                        xycoords='axes fraction',
                        va='top', ha='left',
                        fontsize = 6)
            plt.tight_layout()
            plt.savefig(f"{ASTRESULTDIR/fpath.stem}_AST_Mag{Mag_UP}.png")

            if df_targ.empty:
                pass
            else:
                cutsizes = 49
                for i, row in df_targ_eph.iterrows():
                
                    #1. cut asteroia area
                    #print(i)
                    cut_hdu = Cutout2D(
                                data = hdul[0].data,
                                position = (pos_targ_init[i]),
                                size=(cutsizes, cutsizes) #cut ccd
                                )
                    avg, med, std = sigma_clipped_stats(cut_hdu.data)  # by default, 3-sigma 5-iteration.

                    fig_set = plt.figure(figsize=(8, 5.5))
                    
                    ax11 = plt.subplot2grid((2, 2), (0,0),
                                fig=fig_set)
                    im11 = zimshow(ax11, cut_hdu.data)
                    ax11.plot(round(cutsizes/2), round(cutsizes/2), 'rx')
                    ax11.set_ylabel('pixels')
                    ax11.grid(ls=':')
                    ax11.set_title(f'Asteroid area image', fontsize=9)
                    ax11.annotate(   f"mean: {np.mean(cut_hdu.data):.01f}, std: {np.std(cut_hdu.data):.01f} \nmax: {np.max(cut_hdu.data):.01f}, min: {np.min(cut_hdu.data):.01f} \nNumber of Pixel: {np.shape(cut_hdu.data)[0]:.0f}x{np.shape(cut_hdu.data)[1]:.0f}",
                        xy=(0, 0), xytext=(0.1, -0.20),
                        xycoords='axes fraction',
                        va='top', ha='left',
                        fontsize=8)
                    plt.colorbar(im11,
                                ax=ax11,
                                label="ADU",
                                fraction=0.0455, pad=0.04)
                    #print("Image size is: ", cut_hdu.data.shape)

                    #2. Get center dx, dy
                    thresh_3sig = med + 3 * std
                    mask_3sig = (cut_hdu.data < thresh_3sig)
                    center = centroid_com(
                                data = cut_hdu.data,
                                mask = mask_3sig
                                )

                    centerdx = center[0] - ((cutsizes+1)/2)
                    centerdy = center[1] - ((cutsizes+1)/2)
                    # print("type(center):", type(center))
                    # print("center:", center)
                    # print("center dx, dy:", centerdx, centerdy)

                    ax12 = plt.subplot2grid((2,2), (0,1),
                                fig=fig_set)
                    ax12.grid(ls=':')
                    ax12.set_title(f'The new center of asteroid', fontsize=9)
                    im12 = ax12.imshow(mask_3sig.astype(int),
                        origin="lower")
                    im12 = ax12.imshow(cut_hdu.data,
                            alpha=0.4,
                            origin="lower")
                    ax12.plot(*center, 'rx')
                    ax12.annotate(f"center: {center[0]:.02f}, {center[1]:.02f}\ncenter dx, dy: {centerdx:.02f}, {centerdy:.02f}",
                            xy=(0, 0), xytext=(0.01, -0.20),
                            xycoords='axes fraction',
                            va='top', ha='left',
                            fontsize=8)
                    
                    ax11.annotate(f"asteroid No.{i}: {row['targetname']}, \n{row['datetime_str']}, V_mag {row['V']}",
                            xy=(1, 0), xytext=(-0.1, 1.33),
                            xycoords='axes fraction',
                            va='top', ha='left',
                            fontsize=8)
                    
                    plt.colorbar(im12,
                                ax=ax12,
                                label="ADU",
                                fraction=0.0455, pad=0.04)
                    plt.suptitle(f"{fpath.name}", 
                                fontsize=9)
                    
                    plt.tight_layout()
                    plt.savefig(f"{ASTRESULTDIR/fpath.stem}_AST_Mag{Mag_UP}_{i:02d}.png")
                    #plt.show()
                    plt.close()

    return 0


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