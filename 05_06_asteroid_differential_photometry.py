"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
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
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

from astroquery.jplhorizons import Horizons
from astroquery.imcce import Skybot

import ysfitsutilpy as yfu
import ysphotutilpy as ypu

import _astro_utilities
import _Python_utilities

from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip, sigma_clipped_stats
from photutils.centroids import centroid_com

from photutils.aperture import CircularAperture as CAp
from photutils.aperture import CircularAnnulus as CAn

#%%
#######################################################
# for log file
#######################################################
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
BASEDIR = Path("/mnt/Rdata/OBS_data") 
#BASEDIR = Path("/Volumes/OBS_data") 

DOINGDIR = Path(BASEDIR/ "asteroid" / "RiLA600_STX-16803_-_1bin")
DOINGDIR = Path(BASEDIR/ "asteroid" / "GSON300_STF-8300M_-_1bin")
# DOINGDIR = Path(BASEDIR/ "asteroid" / "good")
# DOINGDIR = Path(BASEDIR/"2023OA/asteroids_teacher")

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(DOINGDIR))

# filter_str = '678FREDEGUNDIS_LIGHT_-_2023-11-11_-_RiLA600_STX-16803_-_1bin'
# DOINGDIRs = [x for x in DOINGDIRs if filter_str in x]
# remove = 'BIAS'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'DARK'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'FLAT'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]

print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################

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

# %%
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)

    ASTRESULTDIR = DOINGDIR / _astro_utilities.Asteroid_result_dir

    if str(DOINGDIR.parts[-2]) == "RiLA600_STX-16803_-_1bin" :
        DOINGDIR = DOINGDIR / _astro_utilities.reduced_nightsky_dir
    if str(DOINGDIR.parts[-2]) == "GSON300_STF-8300M_-_1bin" :
        DOINGDIR = DOINGDIR / _astro_utilities.reduced_dir

    if not ASTRESULTDIR.exists():
        os.makedirs("{}".format(str(ASTRESULTDIR)))
        print("{} is created...".format(str(ASTRESULTDIR)))
    
    fits_in_dir = sorted(list(DOINGDIR.glob('*.fit*')))
    #print("fits_in_dir", fits_in_dir)
    print("len(fits_in_dir)", len(fits_in_dir))

    if len(fits_in_dir) == 0 :
        print(f"There is no fits fils in {DOINGDIR}")
        pass
    else : 
        print(f"Starting: {str(DOINGDIR.parts[-1])}")

        summary = yfu.make_summary(DOINGDIR/"*.fit*")
        print("len(summary):", len(summary))
        #print("summary:", summary)
        #print(summary["file"][0])
        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)
        #print("df_light:\n{}".format(df_light))

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
                    
                if hdul[0].header['CCDNAME'] == 'STF-8300M' :
                    val_figsize = (13, 5.2)
                    val_fraction = 0.035
                    hdul[0].header["GAIN"] = 0.37,
                    hdul[0].header["RDNOISE"] = 9.3

                if hdul[0].header['CCDNAME'] == 'STX-16803' :
                    val_figsize=(12, 6.2)
                    val_fraction = 0.0455
                    hdul[0].header["GAIN"] = 1.27
                    hdul[0].header["RDNOISE"] = 9.0    

                # It is used as a rough estimate, so no need to be accurate:
                PIX2ARCSEC = hdul[0].header["PIXSCALE"]
                rdnoise = hdul[0].header["RDNOISE"]
                gain    = hdul[0].header["GAIN"]

                print(rdnoise, gain, PIX2ARCSEC)
                

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
                
                # Get the radius of the smallest circle which encloses all the pixels
                rad = yfu.fov_radius(header=hdul[0].header, unit=u.deg)
                print("rad: {}".format(rad))

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
                        # try : 
                            #print("type(row)", type(row))
                            #Query the ephemerides of this target! 
                        obj = Horizons(id=row['Number'], 
                                    location=observatory_code, 
                                    epochs=t_middle.jd)
                        obj_ephem = obj.ephemerides()
                        #print(obj_ephem)
                        df_eph = obj_ephem.to_pandas()
                        df_targ_eph = pd.concat([df_targ_eph, df_eph], axis = 0)
                        # except : 
                        #     continue

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

                    fig_set = plt.figure(figsize=val_figsize)
                    ax1 = plt.subplot2grid((1,2), (0,0),
                                        fig=fig_set)
                    im1 = _astro_utilities.zimshow(ax1, hdul[0].data, )
                    ax1.set_title('Pixel coordinate system', fontsize=9)
                    ax1.tick_params(labelsize=8)
                    plt.colorbar(im1, ax = ax1, fraction=val_fraction, pad=0.04)

                    ax2 = plt.subplot2grid((1,2), (0,1),
                                        projection=wcs,
                                        fig=fig_set)
                    im2 = _astro_utilities.zimshow(ax2, hdul[0].data, )
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

                    if df_targ_eph.empty:
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

                    ax2.annotate(f"image center (RA, DEC): ({cent_coord.ra:.03f}, {cent_coord.dec:.03f})\ntelescope center (RA, DEC): ({hdul[0].header['RA']*u.deg:.03f}, {hdul[0].header['DEC']*u.deg:.03f})\noffset (RA, DEC): ({offset_RA:.03f}, {offset_DEC:.03f})\noffset (AZ, ALT): ({offset_AZ:.03f}, {offset_ALT:.03f})",
                                xy=(0, 0), xytext=(0.6, -0.1),
                                xycoords='axes fraction',
                                va='top', ha='left',
                                fontsize = 6)
                    plt.tight_layout()
                    plt.savefig(f"{ASTRESULTDIR/fpath.stem}_AST_Mag{Mag_UP}.png")
                    plt.close()

                    from photutils.detection import DAOStarFinder
                    from astropy.stats import sigma_clipped_stats

                    FWHM = FWHM_INIT
                    avg, med, std = sigma_clipped_stats(hdul[0].data)  # by default, 3-sigma 5-iteration.
                    thresh = 5. * std

                    DAOfind = DAOStarFinder(
                                            fwhm = FWHM,
                                            threshold = thresh,
                                            # sharplo = 0.2, sharphi = 1.0,  # default values: sharplo=0.2, sharphi=1.0,
                                            # roundlo = 0, roundhi = 1.0,  # default values -1 and +1
                                            # sigma_radius = 3,           # default values 1.5
                                            # ratio = 1.0,                  # 1.0: circular gaussian
                                            # exclude_border = True         # To exclude sources near edges
                                            )

                    DAOfound = DAOfind(hdul[0].data)
                    print("len(DAOfound) :",len(DAOfound))
                    print(DAOfound.colnames)
                    DAOfound
                    DAOfound.write(f"{ASTRESULTDIR/fpath.stem}_DAOStarfinder_fwhm_{FWHM}.csv",
                                                overwrite = True,
                                                format='ascii.fast_csv')
                    df_DAO = DAOfound.to_pandas()
                    print(type(df_DAO))
                    df_DAO

                    import numpy as np
                    pos = np.transpose((DAOfound['xcentroid'], DAOfound['ycentroid']))

                    from photutils.aperture import CircularAperture as CAp
                    from photutils.aperture import CircularAnnulus as CAn
                    apert = CAp(pos, r=R_AP)
                    #apert
                    annul = CAn(positions=pos, r_in= R_IN, r_out=R_OUT)
                    #annul

                    from astropy.wcs import WCS

                    wcs = WCS(hdul[0].header)
                    #print("wcs :", wcs)
                    #print("type(wcs) :", type(wcs))
                    #print("dir(wcs) :", dir(wcs))

                    wcs.pixel_n_dim

                    fig_set = plt.figure(figsize=val_figsize)

                    ax1 = plt.subplot2grid((1,2), (0,0),
                                        fig=fig_set)
                    im1 = _astro_utilities.zimshow(ax1, hdul[0].data, )
                    ax1.set_title('Pixel coordinate system', fontsize=9)
                    ax1.tick_params(labelsize=8)

                    ax2 = plt.subplot2grid((1,2), (0,1),
                                        projection=wcs,
                                        fig=fig_set)
                    im2 = _astro_utilities.zimshow(ax2, hdul[0].data, )
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
                    ax2.coords['ra'].set_minor_frequency(2)
                    ax2.coords['dec'].set_minor_frequency(2)
                    ax2.tick_params(labelsize=8)

                    annul.plot(ax1, color="r")
                    annul.plot(ax2, color="r")

                    cbar1 = plt.colorbar(im1, ax = ax1, fraction=val_fraction, pad=0.04)
                    cbar2 = plt.colorbar(im2, ax = ax2, fraction=val_fraction, pad=0.04, )
                    cbar1.ax.tick_params(labelsize=8)
                    cbar2.ax.tick_params(labelsize=8)

                    plt.suptitle(f"fname: {fpath.name}\n Result of DAOFinder", fontsize=10,)

                    ax1.annotate(f'FWHM: {FWHM}', fontsize=8,
                        xy=(0, 0), xytext=(-10, -30), va='top', ha='left',
                        xycoords='axes fraction', textcoords='offset points')

                    ax1.annotate(f'Sky threshold: {thresh:.02f}', fontsize=8,
                        xy=(0, 0), xytext=(-10, -40), va='top', ha='left',
                        xycoords='axes fraction', textcoords='offset points')

                    ax1.annotate(f'Number of star(s): {len(DAOfound)}', fontsize=8,
                        xy=(0, 0), xytext=(-10, -50), va='top', ha='left',
                        xycoords='axes fraction', textcoords='offset points')

                    plt.tight_layout()
                    plt.savefig(f"{ASTRESULTDIR/fpath.stem}_DAOStarfinder_fwhm_{FWHM}.png")

                    # plt.show()
                    # plt.close()

                    if df_targ_eph.empty:
                        pass
                    else:
                        
                        # Initialize PanSTARRS1 class
                        q = ypu.PanSTARRS1(
                            ra=cent_coord.ra.value,
                            dec=cent_coord.dec.value,
                            radius=rad,
                            column_filters={"gmag":"11.0..15.0", "e_gmag":"<0.10"},
                        )

                        # Query to the website (VizieR)
                        # This is where the most of the time is spent.
                        q.query()

                        # # Only select the stars within 50-pixel bezel in the FOV.
                        q.select_xyinFOV(hdul[0].header,
                                        #bezel=50
                                        bezel=5*FWHM_INIT*PIX2ARCSEC
                                        )

                        # # Remove objects not suitable for differential photometry (see description below)
                        q.drop_for_diff_phot(drop_by_Kron=True)

                        q_stars_orig = q.queried.copy()
                        pos_stars_orig = np.array([q_stars_orig["x"], q_stars_orig["y"]]).T
                        q_stars_orig

                        q.drop_star_groups(crit_separation = 6 * FWHM_INIT)
                        q_stars_diropped = q.queried.copy()

                        q_stars = q.queried.copy()
                        #pos_sky_targ_init, pos_pix_targ_init

                        pos_stars = np.array([q_stars["x"], q_stars["y"]]).T

                        df_stars = q_stars_diropped.to_pandas()

                        for idx, row in df_targ_eph.iterrows():
                            #try: 
                            ASTNAME = row['targetname']
                            cutsizes = 49
                            #1. cut asteroia area
                            #print(i)
                            try :
                                cut_hdu = Cutout2D(
                                            data = hdul[0].data,
                                            position = (pos_targ_init[idx]),
                                            size=(cutsizes, cutsizes) #cut ccd
                                            )
                                avg, med, std = sigma_clipped_stats(cut_hdu.data)  # by default, 3-sigma 5-iteration.

                                fig_set = plt.figure(figsize=(8, 5.5))
                                
                                ax11 = plt.subplot2grid((2, 2), (0,0),
                                            fig=fig_set)
                                im11 = _astro_utilities.zimshow(ax11, cut_hdu.data)
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
                                plt.savefig(f"{ASTRESULTDIR/fpath.stem}_AST_Mag{Mag_UP}_{ASTNAME}.png")
                                #plt.show()
                                # plt.close()
                            except : 
                                continue


                            from photutils.aperture import aperture_photometry as apphot
                            apphot_result = apphot(hdul[0].data, apert, method='center')
                            print(type(apphot_result))
                            apphot_result
                            df_apphot = apphot_result.to_pandas()
                            df_apphot
                            ap_area   = apert.area

                            # since our `annul` has many elements,
                            mask_apert = (apert.to_mask(method='center'))
                            mask_annul = (annul.to_mask(method='center'))

                            mag_ann  = np.zeros(len(apphot_result))
                            merr_ann = np.zeros(len(apphot_result))

                            #Returns magnitude from flux.
                            def mag_inst(flux, ferr):
                                m_inst = -2.5 * np.log10(flux)
                                merr   = 2.5/ np.log(10) * ferr / flux
                                return m_inst, merr

                            fig, ax = plt.subplots()

                            for i in range(len(apphot_result)):
                                annul_weighted = mask_annul[i].multiply(hdul[0].data)
                                sky_non0   = np.nonzero(annul_weighted)
                                sky_pixel  = annul_weighted[sky_non0]
                                msky, sky_std, nsky, nrej = _astro_utilities.sky_fit(sky_pixel, method='mode', 
                                                                                    mode_option='sex')

                                flux_star = apphot_result['aperture_sum'][i] - msky * ap_area  # total - sky

                                flux_err  = np.sqrt(apphot_result['aperture_sum'][i] * gain    # Poissonian (star + sky)
                                                    + ap_area * rdnoise**2 # Gaussian
                                                    + (ap_area * (gain * sky_std))**2 / nsky )

                                mag_ann[i], merr_ann[i] = mag_inst(flux_star, flux_err)
                                # print(i, msky, sky_std, nsky, nrej)
                                # print(i, mag_ann[idx], merr_ann[idx])
                                try:
                                    ax.errorbar(i, mag_ann[i], yerr=merr_ann[i],
                                                marker='x',
                                                #ms=10,
                                                capsize=3)
                                except:
                                    continue

                            ax.invert_yaxis()
                            ax.set_ylim(ymin=-20, ymax=0)

                            plt.xlabel('Star ID')
                            plt.ylabel('Instrumental mag')
                            plt.grid(ls=':')
                            
                            plt.savefig(f"{ASTRESULTDIR}/{fpath.name}_stars_IM.png")
                            # plt.show()
                            # plt.close()

                            for i, row in df_apphot.iterrows():
                                annul_weighted = mask_annul[i].multiply(hdul[0].data)
                                sky_non0   = np.nonzero(annul_weighted)
                                sky_pixel  = annul_weighted[sky_non0]
                                msky, sky_std, nsky, nrej = _astro_utilities.sky_fit(sky_pixel, method='mode', mode_option='sex')

                                flux_star = apphot_result['aperture_sum'][i] - msky * ap_area  # total - sky

                                flux_err  = np.sqrt(apphot_result['aperture_sum'][i] * gain    # Poissonian (star + sky)
                                                    + ap_area * rdnoise**2 # Gaussian
                                                    + (ap_area * (gain * sky_std))**2 / nsky )

                                mag_ann, merr_ann = mag_inst(flux_star, flux_err)
                                df_apphot.loc[i, 'msky'] = msky
                                df_apphot.loc[i, 'sky_std'] = sky_std
                                df_apphot.loc[i, 'nsky'] = nsky
                                df_apphot.loc[i, 'nrej'] = nrej
                                df_apphot.loc[i, 'flux_star'] = flux_star
                                df_apphot.loc[i, 'flux_err'] = flux_err
                                df_apphot.loc[i, 'mag_inst'] = mag_ann
                                df_apphot.loc[i, 'merr_inst'] = merr_ann

                            df_apphot.to_csv(f"{ASTRESULTDIR}/{fpath.stem}_m_inst.csv")
                            #df_apphot


                            # fig, axs = plt.subplots(1, 1, figsize=(8, 5), sharex=False, sharey=False, gridspec_kw=None)

                            # im = _astro_utilities.zimshow(axs, hdul[0].data)

                            # targ_an.plot(axs, color="r")
                            # targ_ap.plot(axs, color="b")
                            # axs.text(pos_targ_init[idx][0]+10, pos_targ_init[idx][1]+10, f"ASTNAME", fontsize=8)

                            # phot_targ = ypu.apphot_annulus(hdul[0].data, targ_ap, targ_an, error=yfu.errormap(hdul[0].data))
                            # phot_targ

                            # _phot_stars = []

                            # for i, row in df_stars.iterrows():
                            #     pos_star = SkyCoord(row["RAJ2000"], row["DEJ2000"],
                            #                         **SKYC_KW).to_pixel(wcs, origin=0, mode='wcs')
                            #     ap = CAp([pos_star[0], pos_star[1]], r=R_AP)
                            #     an = CAn([pos_star[0], pos_star[1]], r_in=R_IN, r_out=R_OUT)
                            #     _phot_star = ypu.apphot_annulus(hdul[0].data, ap, an, error=yfu.errormap(hdul[0].data))
                            #     _phot_star["Rmag"] = row["rmag"]
                            #     _phot_star["e_Rmag"] = row["e_rmag"]
                            #     _phot_star["grcolor"] = row["gmag"] - row["rmag"] #row["grcolor"]
                            #     _phot_star["e_grcolor"] = row["e_gmag"] + row["e_rmag"]
                            #     _phot_star["id"] = i
                            #     _phot_star["objID"] = int(row["objID"])
                            #     _phot_stars.append(_phot_star)
                            #     axs.text(pos_star[0]+10, pos_star[1]+10, f"star {i}", fontsize=8)
                            #     ap.plot(axs, color="orange")
                            #     an.plot(axs, color="w")
                            # plt.title(f"{fpath.name}")

                            # plt.tight_layout()
                            # plt.savefig(f"{ASTRESULTDIR}/{fpath.stem}_asteroid_PS1_query.png")
                            # # plt.show()
                            # # plt.close()

                            fig_set = plt.figure(figsize=val_figsize)

                            ax1 = plt.subplot2grid((1,2), (0,0),
                                                fig=fig_set)
                            im1 = _astro_utilities.zimshow(ax1, hdul[0].data, )
                            ax1.set_title('Pixel coordinate system', fontsize=9)
                            ax1.tick_params(labelsize=8)

                            ax2 = plt.subplot2grid((1,2), (0,1),
                                                projection=wcs,
                                                fig=fig_set)
                            im2 = _astro_utilities.zimshow(ax2, hdul[0].data, )

                            phot_targ = ypu.apphot_annulus(hdul[0].data, targ_ap, targ_an, error=yfu.errormap(hdul[0].data))
                            print("len(phot_targ) :", len(phot_targ))

                            _phot_stars = []

                            for i, row in df_stars.iterrows():
                                pos_star = SkyCoord(row["RAJ2000"], row["DEJ2000"],
                                                    **SKYC_KW).to_pixel(wcs, origin=0, mode='wcs')
                                ap = CAp([pos_star[0], pos_star[1]], r=R_AP)
                                an = CAn([pos_star[0], pos_star[1]], r_in=R_IN, r_out=R_OUT)
                                _phot_star = ypu.apphot_annulus(hdul[0].data, ap, an, error=yfu.errormap(hdul[0].data))
                                _phot_star["Rmag"] = row["rmag"]
                                _phot_star["e_Rmag"] = row["e_rmag"]
                                _phot_star["grcolor"] = row["gmag"] - row["rmag"] #row["grcolor"]
                                _phot_star["e_grcolor"] = row["e_gmag"] + row["e_rmag"]
                                _phot_star["id"] = i
                                _phot_star["objID"] = int(row["objID"])
                                _phot_stars.append(_phot_star)
                                ax1.text(pos_star[0]+10, pos_star[1]+10, f"star {i}", fontsize=8)
                                ap.plot(ax1, color="orange")
                                an.plot(ax1, color="w")
                                ax2.text(pos_star[0]+10, pos_star[1]+10, f"star {i}", fontsize=8)
                                ap.plot(ax2, color="orange")
                                an.plot(ax2, color="w")

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
                            ax2.coords['ra'].set_minor_frequency(2)
                            ax2.coords['dec'].set_minor_frequency(2)
                            ax2.tick_params(labelsize=8)

                            targ_an.plot(ax1, color="r")
                            targ_ap.plot(ax1, color="b")
                            targ_an.plot(ax2, color="r")
                            targ_ap.plot(ax2, color="b")
                            ax1.text(pos_targ_init[idx][0]+10, pos_targ_init[idx][1]+10, f"{ASTNAME}", fontsize=8)
                            ax2.text(pos_targ_init[idx][0]+10, pos_targ_init[idx][1]+10, f"{ASTNAME}", fontsize=8)

                            cbar1 = plt.colorbar(im1, ax = ax1, fraction=val_fraction, pad=0.04)
                            cbar2 = plt.colorbar(im2, ax = ax2, fraction=val_fraction, pad=0.04, )
                            cbar1.ax.tick_params(labelsize=8)
                            cbar2.ax.tick_params(labelsize=8)

                            plt.suptitle(f"fname: {fpath.name}\n Result of PS1 query for differential photometry", fontsize=10,)

                            ax1.annotate(f'FWHM: {FWHM}', fontsize=8,
                                xy=(0, 0), xytext=(-10, -30), va='top', ha='left',
                                xycoords='axes fraction', textcoords='offset points')

                            ax1.annotate(f'Sky threshold: {thresh:.02f}', fontsize=8,
                                xy=(0, 0), xytext=(-10, -40), va='top', ha='left',
                                xycoords='axes fraction', textcoords='offset points')

                            ax1.annotate(f'Number of star(s): {len(DAOfound)}', fontsize=8,
                                xy=(0, 0), xytext=(-10, -50), va='top', ha='left',
                                xycoords='axes fraction', textcoords='offset points')

                            plt.tight_layout()
                            plt.savefig(f"{ASTRESULTDIR}/{fpath.stem}_asteroid_PS1_query_{ASTNAME}.png")

                            # plt.show()
                            # plt.close()
                            

                            df_phot_stars = pd.concat(_phot_stars)
                            df_phot_stars = df_phot_stars.dropna()
                            #phot_stars = phot_stars.loc[phot_stars["objID"] != 125240299358089744].copy()  # star 15
                            # SEE THE LAST CELL IN THIS FILE FOR DESCRIPTION
                            #df_phot_stars = df_phot_stars.set_index('id', drop=True)
                            df_phot_stars
                            #print(phot_targ)
                            print(df_phot_stars)
                            #try: 
                            fig, axs = plt.subplots(1, 1, figsize=(7, 7),
                                sharex=False, sharey=False, gridspec_kw=None)

                            _xx = np.linspace(11, 17)
                            axs.plot(df_phot_stars["Rmag"], df_phot_stars["mag"], '+')
                            axs.axhline(phot_targ["mag"][idx],
                                        label=f"Asteroid #{ASTNAME}, instrumental mag :{float(phot_targ['mag'][idx]):.02f}")
                            axs.plot(_xx, _xx + np.median(df_phot_stars['mag'] - df_phot_stars['Rmag']),
                                    label=f"R_inst = R magnitude + ({np.median(df_phot_stars['mag'] - df_phot_stars['Rmag']):.02f})")

                            for _, row in df_phot_stars.iterrows():
                                axs.text(row["Rmag"], row["mag"], int(row["id"]), fontsize=8)

                            axs.set(
                                xlabel="R magnitude (PS1 to R_C filter by Tonry+2012)",
                                ylabel="R_inst"
                                    )
                            axs.axis('equal')
                            axs.grid(which = "major", alpha=0.5)
                            axs.grid(':', which = "minor", alpha=0.3)
                            axs.invert_yaxis()
                            axs.minorticks_on()
                            axs.tick_params(which='both', width=2)
                            axs.tick_params(which='major', length=7)
                            axs.tick_params(which='minor', length=4, color='r')
                            axs.legend(fontsize="8", 
                                    #loc ="ower right"
        )

                            axs.annotate(f"target name: {df_eph['targetname'][0]}\nR_mag: {np.median(df_phot_stars['Rmag'] - df_phot_stars['mag'])+ phot_targ['mag'][0]:.02f}",
                                        fontsize=8,
                                xy=(0, 0), xytext=(0, -0.15), va='top', ha='left',
                                xycoords='axes fraction')
                            
                            plt.tight_layout()
                            plt.savefig(f"{ASTRESULTDIR}/{fpath.stem}_standardization_{ASTNAME}.png")
                            # plt.show()
                            # plt.close()

                            # except Exception as err: 
                            #     print("Err :", err)
                            #     continue


                    