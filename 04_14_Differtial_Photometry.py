# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user


"""
#%%
import os
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clip, sigma_clipped_stats

from astropy.nddata import Cutout2D

from astropy.time import Time
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord

from photutils.aperture import CircularAperture as CAp
from photutils.aperture import CircularAnnulus as CAn
from photutils.centroids import centroid_com
from photutils import aperture_photometry as apphot

from ccdproc import CCDData, ccd_process

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

import _Python_utilities
import _astro_utilities

plt.rcParams.update({'figure.max_open_warning': 0})

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
#######################################################
BASEDIR = _astro_utilities.base_dir

BASEDIRs = sorted(_Python_utilities.getFullnameListOfsubDir(BASEDIR))
print ("BASEDIRs: {}".format(BASEDIRs))
print ("len(BASEDIRs): {}".format(len(BASEDIRs)))
#%%
#####################################################################
# Our object (will be queried to JPL HORIZONS)
OBJID = '216' # Kleopatra

# Observed location
LOCATION = dict(lon = 127.0, lat = 37.3, elevation = 130)

# It is used as a rough estimate, so no need to be accurate:
PIX2ARCSEC = 1.24 * u.arcsec

# Used for any `astropy.SkyCoord` object:
SKYC_KW = dict(unit = u.deg, frame = 'icrs')

# Initial guess of FWHM in pixel
FWHM_INIT = 6

# Photometry parameters
R_AP = 1.5 * FWHM_INIT # Aperture radius
R_IN = 4 * FWHM_INIT   # Inner radius of annulus
R_OUT = 6 * FWHM_INIT  # Outer radius of annulus
#####################################################################

#%%

for BASEDIR in BASEDIRs[0:8]:
    print ("Starting...\n{}".format(BASEDIR))

    BASEDIR = Path(BASEDIR)

    SOLVEDDIR = BASEDIR / _astro_utilities.solved_dir2
    AsteroidRESULTDIR = BASEDIR / _astro_utilities.Asteroid_result_dir

    if not AsteroidRESULTDIR.exists():
        os.makedirs("{}".format(str(AsteroidRESULTDIR)))
        print("{} is created...".format(str(AsteroidRESULTDIR)))

    #%%
    #summary = yfu.make_summary(BASEDIR/"*.fit*")
    summary = yfu.make_summary(SOLVEDDIR/"*.fit*")

    if summary.empty:
        print("The dataframe(summary) is empty")
        pass
    else:
        print("len(summary):", len(summary))
        print("summary:", summary)

    for filt in ["r","v", "b"]:
        summary_filt = summary.loc[summary["FILTER"] == filt].copy()
        
        if summary_filt.empty:
            print("The dataframe(summary_filt) is empty")
            pass
        else:
            print("summary_filt:", summary_filt)
            print("len(summary_filt):", len(summary_filt))

            df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
            df_light = df_light.reset_index(drop=True)
            print("df_light:\n{}".format(df_light))

            #for filt in ["r", "v", "b"]:
            for filt in ["r","v", "b"]:
                df_light_filt = df_light.loc[df_light["FILTER"] == filt].copy()
            
                if df_light_filt.empty:
                    print("The dataframe(df_light_filt) is empty")
                    pass
                else:
                    print("df_light_filt:", df_light_filt)
                    print("len(df_light_filt):", len(df_light_filt))
                    
                for _, row  in df_light_filt.iterrows():
                    fpath = Path(row["file"])
                    print("type(fpath)", type(fpath))
                    print("fpath", fpath)
            
                    #%%
                    # load as ccd
                    ccd = yfu.load_ccd(fpath, 
                                    unit="adu")
                    
                    ## Get the Ephemeris
                    _, eph, _ = ypu.horizons_query(
                                                OBJID, 
                                                epochs=Time(ccd.header["DATE-OBS"]).jd, 
                                                location=LOCATION)
                    print("eph:", eph)
                    eph.write("{}_eph-{}.csv"\
                            .format(str(AsteroidRESULTDIR / fpath.stem), OBJID),
                            overwrite = True)

                    # Initial Photometry of the Target
                    pos_targ_init = SkyCoord(eph["RA"], 
                                            eph["DEC"], 
                                            **SKYC_KW).to_pixel(ccd.wcs)
                    targ_ap = CAp([pos_targ_init[0][0], 
                            pos_targ_init[1][0]], 
                            r=R_AP)
                    targ_an = CAn([pos_targ_init[0][0], 
                            pos_targ_init[1][0]], 
                            r_in=R_IN, 
                            r_out=R_OUT)
                    print("pos_targ_init:", pos_targ_init)

                    phot_targ = ypu.apphot_annulus(ccd, 
                                                    targ_ap, 
                                                    targ_an, 
                                                    error=yfu.errormap(ccd))
                    print("phot_targ:", phot_targ)

                    #%%
                    # 별의 목록을 가져옴.
                    r_fov = yfu.fov_radius(ccd.header+ccd.wcs.to_header())
                    print("r_fov: ", r_fov)
                    
                    ps1 = ypu.PanSTARRS1(
                                        ccd.wcs.wcs.crval[0]*u.deg, 
                                        ccd.wcs.wcs.crval[1]*u.deg, 
                                        radius=r_fov,
                                        column_filters = {"{}mag".format(filt.lower()):"12.0..14.5", 
                                                        "e_{}mag".format(filt.lower()):"<0.10", 
                                                        "nr":">5"}
                                        )
                    print("ps1:", ps1)
                    #%%
                    # 가까이 붙어 있는 별은 지우자.
                    isnear = ypu.organize_ps1_and_isnear(
                                        ps1, 
                                        header = ccd.header+ccd.wcs.to_header(), 
                                        bezel = 5*FWHM_INIT * PIX2ARCSEC.value,
                                        nearby_obj_minsep= 5 * FWHM_INIT*PIX2ARCSEC.value,
                                        group_crit_separation = 6 * FWHM_INIT
                                    )
                    print("isnear:", isnear)

                    # 별의 목록
                    df_stars = ps1.queried.to_pandas()
                    
                    if df_stars.empty:
                        print("The dataframe(df_stars) is empty")
                        pass

                    else:
                        print(df_stars)
                        print("len(df_stars):", len(df_stars))

                        df_stars.to_csv("{}_stars.csv".format(str(AsteroidRESULTDIR / fpath.stem)))
                        print("df_stars:", df_stars)

                        #%%
                        fig, axs = plt.subplots(1, 1, 
                                                figsize=(8, 8), 
                                                sharex=False, 
                                                sharey=False, 
                                                gridspec_kw=None
                                                )

                        axs = plt.subplot(projection=ccd.wcs, label='overlay')

                        im = yvu.norm_imshow(axs, 
                                        ccd, zscale=True)

                        plt.colorbar(im, 
                                    ax=axs,
                                    fraction=0.042, #0.0455
                                    #aspect=10,
                                    pad=0.12)
                        overlay = axs.get_coords_overlay('fk5')
                                                    #overlay = ax.get_coords_overlay('icrs')
                        overlay.grid(True, color='white', ls=':', alpha=0.7)
                        overlay[0].set_axislabel('Right Ascension (J2000)')
                        overlay[1].set_axislabel('Declination (J2000)')

                        # sat tick label
                        lon, lat = axs.coords
                        lon.set_ticks(color='red')
                        lon.set_ticks_position('lbtr')
                        lon.set_ticklabel_position('lbtr')
                        lat.set_ticks(color='blue')
                        lat.set_ticks_position('lbtr')
                        lat.set_ticklabel_position('lbtr')

                        plt.title(f"Marking {eph['targetname'][0]} and {len(df_stars)} stars",
                                fontsize = 14, pad=50)

                        plt.annotate(f"{str(fpath.name)}",
                                fontsize=10, xy=(0, 0), xytext=(3, -50), va='top', ha='left',
                                xycoords='axes fraction', textcoords='offset points')
                        #
                        targ_ap.plot(axs, color="r")
                        targ_an.plot(axs, color="b")

                        # 각 별의 측광을 수행
                        _phot_stars = []
                        for idx, row in df_stars.iterrows():
                            print("Starting photometry star {}:".format(idx))

                            #별의 적경, 적위를 이미지 안에서의 픽셀 값으로 
                            pos_star = SkyCoord(row["RAJ2000"], 
                                                row["DEJ2000"], 
                                                **SKYC_KW).to_pixel(ccd.wcs)
                            ap = CAp([pos_star[0], 
                                    pos_star[1]], 
                                    r=R_AP)
                            an = CAn([pos_star[0], 
                                    pos_star[1]], 
                                    r_in=R_IN, 
                                    r_out=R_OUT)
                            _phot_star = ypu.apphot_annulus(ccd, ap, an, 
                                                            error=yfu.errormap(ccd))
                            _phot_star["{}mag".format(filt.upper())] = row["{}mag".format(filt.upper())]
                            _phot_star["e_{}mag".format(filt.upper())] = row["e_{}mag".format(filt.upper())]
                            _phot_star["grcolor"] = row["grcolor"]
                            _phot_star["e_grcolor"] = row["e_grcolor"]
                            _phot_star["id"] = idx
                            _phot_star["objID"] = int(row["objID"])
                            _phot_stars.append(_phot_star)
                            axs.text(pos_star[0]+10, 
                                    pos_star[1]+10, 
                                    f"star {idx}", 
                                    fontsize=8)
                            ap.plot(axs, color="orange")
                            an.plot(axs, color="w")
                        plt.title("Marking {} and {} stars".format(eph["targetname"][0], 
                                len(df_stars)), fontsize = 14)
                        plt.tight_layout()
                        plt.savefig("{}_stars.png".format(str(AsteroidRESULTDIR / fpath.stem)))
                        #plt.show()
                        #%%
                        phot_stars = pd.concat(_phot_stars)
                        # phot_stars = phot_stars.loc[phot_stars["objID"] != 110823405221754720].copy()  # star 15
                        # SEE THE LAST CELL IN THIS FILE FOR DESCRIPTION
                        print("phot_stars: ", phot_stars)
                        phot_stars.to_csv("{}_phot_stars.csv".format(str(AsteroidRESULTDIR / fpath.stem)))

                        # %%
                        # Centroid and Re-photometry
                        _phot_stars = []
                        for idx, row in df_stars.iterrows():
                            print("Starting RE-photometry star {}:".format(idx))

                            #1. 별의 적경, 적위를 이미지 안에서의 픽셀 값으로 
                            pos_star = SkyCoord(row["RAJ2000"], 
                                                row["DEJ2000"], 
                                                **SKYC_KW).to_pixel(ccd.wcs)
                            
                            #2. Loading and Cut Data
                            cutsizes = 32
                            cut_hdu = Cutout2D(
                                        data = ccd, 
                                        position = ([pos_star[0], pos_star[1]]), 
                                        size=(cutsizes, cutsizes) #cut ccd
                                        )
                            avg, med, std = sigma_clipped_stats(cut_hdu.data)  # by default, 3-sigma 5-iteration.
                            thresh_3sig = med + 3 * std
                            mask_3sig = (cut_hdu.data < thresh_3sig)
                            center = centroid_com(
                                        data = cut_hdu.data, 
                                        mask = mask_3sig
                                        )
                            
                            centerdx = int(cutsizes/2-center[0])
                            centerdy = int(cutsizes/2-center[1])

                            print("type(center):", type(center))
                            print("center:", center)
                            print("center dx, dy:", centerdx, centerdy)
                            print("center dx, dy:", centerdx, centerdy)
                            print("center dx, dy:", centerdx, centerdy)

                            #3. Loading and RE-Cut Data with New center
                            bigcutsizes = 100
                            bigcut_hdu = Cutout2D(
                                        data = ccd, 
                                        position = ([pos_star[0], pos_star[1]]), 
                                        size=(bigcutsizes, bigcutsizes) #cut ccd
                                        )
                            avg, med, std = sigma_clipped_stats(bigcut_hdu.data)  # by default, 3-sigma 5-iteration.

                            #4. re center Aperture and Annulus
                            bigcenter = [bigcutsizes/2 - centerdx, bigcutsizes/2 - centerdy]
                            
                            
                            ap = CAp(positions = bigcenter, 
                                    r=R_AP)
                            an = CAn(positions = bigcenter, 
                                    r_in=R_IN, 
                                    r_out=R_OUT)
                            
                            print("ap", ap)
                            print("type(ap)", type(ap))
                            print("an", an)
                            print("type(an)", type(an))

                            _phot_star = ypu.apphot_annulus(ccd, ap, an, 
                                                            error=yfu.errormap(ccd))

                            _phot_star["{}mag".format(filt.upper())] = row["{}mag".format(filt.upper())]
                            _phot_star["e_{}mag".format(filt.upper())] = row["e_{}mag".format(filt.upper())]
                            _phot_star["grcolor"] = row["grcolor"]
                            _phot_star["e_grcolor"] = row["e_grcolor"]
                            _phot_star["id"] = idx
                            _phot_star["objID"] = int(row["objID"])
                            _phot_stars.append(_phot_star)
                            axs.text(pos_star[0]+10, 
                                    pos_star[1]+10, 
                                    f"star {idx}", 
                                    fontsize=8)
                            ap.plot(axs, color="orange")
                            an.plot(axs, color="w")
                        plt.title("Marking {} and Re-centering {} stars".format(eph["targetname"][0], 
                                len(df_stars)), fontsize = 14)
                        plt.tight_layout()
                        plt.savefig("{}_stars_Re.png".format(str(AsteroidRESULTDIR / fpath.stem)))
                        #plt.show()
                        #%%
                        phot_stars = pd.concat(_phot_stars)
                        # phot_stars = phot_stars.loc[phot_stars["objID"] != 110823405221754720].copy()  # star 15
                        # SEE THE LAST CELL IN THIS FILE FOR DESCRIPTION
                        print("phot_stars: ", phot_stars)
                        phot_stars.to_csv("{}_phot_stars_Re.csv".format(str(AsteroidRESULTDIR / fpath.stem)))

                   