# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

2019.09.29  modify - missing 'IMAGETYP' on APT

ModuleNotFoundError: No module named 'ccdproc'
conda install -c condaforge ccdproc
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
from astropy.visualization import ZScaleInterval, ImageNormalize

import ysfitsutilpy as yfu
import ysfitsutilpy as yfu

from glob import glob
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ccdproc import combine, ccd_process, CCDData

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time



from astroquery.jplhorizons import Horizons
from astroquery.imcce import Skybot


import _astro_utilities
import _Python_utilities

from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip, sigma_clipped_stats
from photutils.centroids import centroid_com

from photutils.aperture import CircularAperture as CAp
from photutils.aperture import CircularAnnulus as CAn

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

master_file_dir = 'master_file_Python/'
processing_dir = 'processing_Python/'
integration_dir = 'integration_Python/'
alignment_dir = 'alignment_Python/'



# %%
def znorm(image, **kwargs):
    return ImageNormalize(image, interval=ZScaleInterval(**kwargs))

def zimshow(ax, image, **kwargs):
    return ax.imshow(image, norm=znorm(image, **kwargs), origin='lower')

#####################################################################

#%%
#########################################
#reduceLightFrame
#########################################
def reduceLightFrame(
        DOINGDIR,
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
    MASTERDIR = DOINGDIR / _astro_utilities.master_dir
    REDUCEDDIR = DOINGDIR / _astro_utilities.reduced_dir
    
    df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
    df_light = df_light.reset_index(drop=True)

    for _, row in df_light.iterrows():

        fpath = Path(row["file"])
        ccd = yfu.load_ccd(fpath)
        filt = ccd.header["FILTER"]
        expt = ccd.header["EXPTIME"]
        red = yfu.ccdred(
                ccd,
                output=Path(f"{REDUCEDDIR/ fpath.name}"),
                mdarkpath=str(MASTERDIR / "master_dark_{:.0f}sec.fits".format(expt)),
                mflatpath=str(MASTERDIR / "master_flat_{}_norm.fits".format(filt.upper())),
                # flat_norm_value=1,  # 1 = skip normalization, None = normalize by mean
                overwrite=True,
                )

    return 0

#%%
#########################################
#makeNightskyflatReduceLightFrame
#########################################
def makeNightskyflatReduceLightFrame(DOINGDIR,
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
    MASTERDIR = DOINGDIR / _astro_utilities.master_dir
    REDUC_nightsky = DOINGDIR / _astro_utilities.REDUC_nightsky_dir

    summary_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
    summary_light = summary_light.reset_index(drop=True) 

    for filt in ["b", "v", "r", "L", "R", "G", "B"]:
    #for filt in ["V"]:
        summary_light_filt = summary_light.loc[summary_light["FILTER"] == filt].copy()
        
        if summary_light_filt.empty:
            print("The dataframe(summary_light_filt) is empty")
            pass
        else:
        
            print("len(summary_light_filt):", len(summary_light_filt))
            print("summary_light_filt:", summary_light_filt)

            ccd = yfu.imcombine(
                summary_light_filt["file"].tolist(), 
                combine="med",
                scale="avg", 
                scale_to_0th=False, 
                reject="sc", 
                sigma=2.5,
                verbose=True,
                memlimit = 2.e+11,
                )
            ccd.write(MASTERDIR / f"nightskyflat-{filt}.fits", overwrite=True)

    for _, row in summary_light.iterrows():
        fpath = Path(row["file"])
        filt = row["FILTER"]
        ccd = yfu.ccdred(
            fpath, 
            mflatpath=str(MASTERDIR / f"nightskyflat-{filt}.fits"),
            output=REDUC_nightsky/fpath.name
        )
    
    return 0
    
#%%
#########################################
#solvingLightFrame
#########################################
def solvingLightFrame(DOINGDIR,
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
    
    DOINGDIR = DOINGDIR / _astro_utilities.REDUC_nightsky_dir

    summary_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
    summary_light = summary_light.reset_index(drop=True) 

    df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
    df_light = df_light.reset_index(drop=True)
    print("df_light:\n{}".format(df_light))

    for _, row  in df_light.iterrows():

        fpath = Path(row["file"])
        hdul = fits.open(fpath)
        
        if 'PIXSCALE' in hdul[0].header:
            PIXc = hdul[0].header['PIXSCALE']
        else : 
            PIXc = _astro_utilities.calPixScale(hdul[0].header['FOCALLEN'], 
                                                hdul[0].header['XPIXSZ'],
                                                hdul[0].header['XBINNING'])
        print("PIXc : ", PIXc)
        hdul.close()

        SOLVE, ASTAP, LOCAL = _astro_utilities.checkPSolve(fpath)
        print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)

        if ASTAP :
            print(f"{fpath.name} is solved by ASTAP")
        else : 
            print(f"{fpath.name} is solving now by ASTAP")
            solved = _astro_utilities.ASTAPSolver(fpath, 
                                                    #str(SOLVEDDIR), 
                                                    downsample = 2,
                                                    pixscale = PIXc,
                                                            )

        SOLVE, ASTAP, LOCAL = _astro_utilities.checkPSolve(fpath)
        print("SOLVE:", SOLVE, "ASTAP:", ASTAP, "LOCAL:", LOCAL)

        if LOCAL :
            print(f"{fpath.name} is solved by LOCAL")
        else : 
            print(f"{fpath.name} is solving now by LOCAL")
            if 'PIXSCALE' in hdul[0].header:
                PIXc = hdul[0].header['PIXSCALE']
            else : 
                PIXc = _astro_utilities.calPixScale(hdul[0].header['FOCALLEN'], 
                                                hdul[0].header['XPIXSZ'],
                                                hdul[0].header['XBINNING'])
            print("PIXc : ", PIXc)

            solved = _astro_utilities.LOCALPSolver(fpath, 
                                                    #str(SOLVEDDIR), 
                                                    downsample = 2,
                                                    pixscale = PIXc,
                                                            )
    return 0
    
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
    Suwon = location = EarthLocation(lon=127.005 * u.deg, 
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

    ASTRESULTDIR = DOINGDIR / _astro_utilities.Asteroid_result_dir
    if not ASTRESULTDIR.exists():
        os.makedirs("{}".format(str(ASTRESULTDIR)))
        print("{} is created...".format(str(ASTRESULTDIR)))

    DOINGDIR = DOINGDIR / _astro_utilities.REDUC_nightsky_dir

    summary_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
    summary_light = summary_light.reset_index(drop=True) 

    df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
    df_light = df_light.reset_index(drop=True)
    print("df_light:\n{}".format(df_light))

    for _, row  in df_light.iterrows():
        fpath = Path(row["file"])
        hdul = fits.open(fpath)

        SOLVE, ASTAP, LOCAL = _astro_utilities.checkPSolve(fpath)
        print(SOLVE, ASTAP, LOCAL)
        
        if SOLVE :
            wcs = WCS(hdul[0].header)
            # It is used as a rough estimate, so no need to be accurate:
            #PIX2ARCSEC = 0.62*u.arcsec
            if 'PIXSCALE' in hdul[0].header:
                PIX2ARCSEC = hdul[0].header['PIXSCALE']
            else : 
                PIX2ARCSEC = _astro_utilities.calPixScale(hdul[0].header['FOCALLEN'], 
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
