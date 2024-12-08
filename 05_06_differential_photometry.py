# -*- coding: utf-8 -*-
"""
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
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

import astropy.units as u

import ysfitsutilpy as yfu
import ysphotutilpy as ypu

import _astro_utilities
import _Python_utilities
import _tool_visualization

from astropy.nddata import Cutout2D
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clip, sigma_clipped_stats
from photutils.centroids import centroid_com

from photutils.aperture import CircularAperture as CAp
from photutils.aperture import CircularAnnulus as CAn
from photutils.aperture import aperture_photometry as apphot

from astroquery.simbad import Simbad
from urllib.parse import urlencode

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
# read all files in base directory for processing
BASEDIR = Path("/mnt/Rdata/ASTRO_data") 
DOINGDIR = Path(BASEDIR/ "asteroid" / "RiLA600_STX-16803_-_1bin")
DOINGDIR = Path(BASEDIR/ "asteroid" / "GSON300_STF-8300M_-_1bin")
# DOINGDIR = Path(BASEDIR/ "asteroid" / "good")
# DOINGDIR = Path(BASEDIR/"2023OA/asteroids_teacher")

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(DOINGDIR))
DOINGDIRs = sorted([x for x in DOINGDIR.iterdir() if x.is_dir()])
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

MASTERDIR = [x for x in DOINGDIRs if "CAL-BDF" in str(x)]
MASTERDIR = MASTERDIR[0]/ _astro_utilities.master_dir
print ("MASTERDIR: ", format(MASTERDIR))

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_" in str(x)])

print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

# filter_str = '2023-12'
# DOINGDIRs = [x for x in DOINGDIRs if filter_str in str(x)]
# remove = 'BIAS'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'DARK'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'FLAT'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# print ("DOINGDIRs: ", DOINGDIRs)
# print ("len(DOINGDIRs): ", len(DOINGDIRs))
#######################################################
#%%
#Returns magnitude from flux.
def mag_inst(flux, ferr):
    m_inst = -2.5 * np.log10(flux)
    merr   = 2.5/ np.log(10) * ferr / flux
    return m_inst, merr

import matplotlib
matplotlib.use('Agg')
#####################################################################
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

Mag_UP = 16
#######################################################

#%%
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    print("DOINGDIR", DOINGDIR)

    DIFFPRESULTDIR = DOINGDIR / _astro_utilities.Diff_Phot_dir

    if "RiLA600_STX-16803_" in str(DOINGDIR.parts[-2]) :
        DOINGDIR = DOINGDIR / _astro_utilities.reduced_nightsky_dir
    if str(DOINGDIR.parts[-2]) == "GSON300_STF-8300M_-_1bin" :
        DOINGDIR = DOINGDIR / _astro_utilities.reduced_dir

    if not DIFFPRESULTDIR.exists():
        os.makedirs("{}".format(str(DIFFPRESULTDIR)))
        print("{} is created...".format(str(DIFFPRESULTDIR)))
    
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
        # df_light = summary.loc[summary["FILTER"] == flt].copy()
        df_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        df_light = df_light.reset_index(drop=True)
        #print("df_light:\n{}".format(df_light))

        try:

            for _, row  in df_light.iterrows():

                fpath = Path(row["file"])
                hdul = fits.open(fpath)
                ccd = yfu.load_ccd(fpath)
                flt = hdul[0].header["filter"]

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

                    # print(rdnoise, gain, PIX2ARCSEC)
                    
                    # D.2. Find the observation time and exposure time to set the obs time
                    t_start = Time(hdul[0].header['DATE-OBS'], format='isot')
                    t_expos = hdul[0].header['EXPTIME'] * u.s
                    t_middle = t_start + t_expos / 2 # start time + 0.5 * exposure time
                    #print(f"t_start: {t_start}, t_expos: {t_expos}, t_middle: {t_middle}")
                    
                    cent_coord = yfu.center_radec(ccd_or_header=hdul[0].header, 
                                                center_of_image=True)

                    offset_RA = (cent_coord.ra.to(u.deg) - hdul[0].header['RA']*u.deg).to(u.arcmin)
                    offset_DEC = (cent_coord.dec.to(u.deg) - hdul[0].header['DEC']*u.deg).to(u.arcmin) 
                    altaz = AltAz(obstime=t_middle, location=Suwon)   
                    cent_aa = cent_coord.transform_to(altaz)
                    offset_AZ = (cent_aa.az.to(u.deg) - hdul[0].header['CENTAZ']*u.deg).to(u.arcmin)
                    offset_ALT = (cent_aa.alt.to(u.deg) - hdul[0].header['CENTALT']*u.deg).to(u.arcmin)

                    # Get the radius of the smallest circle which encloses all the pixels
                    rad = yfu.fov_radius(header=hdul[0].header, unit=u.deg)
                    print("rad: {}".format(rad))
                    

                    cent_coord = yfu.center_radec(ccd_or_header=hdul[0].header, 
                                        center_of_image=True)
                    pos_sky = SkyCoord(cent_coord, unit='deg')
                    pos_pix = pos_sky.to_pixel(wcs=wcs)

                    print("pos_sky: {}".format(pos_sky))
                    print("pos_pix: {}".format(pos_pix))


                    #%%
                    # Position of the telescope FOV center 
                    # (RA/DEC of the pixel at the center)
                    cent_coord = yfu.center_radec(ccd_or_header=hdul[0].header, 
                                                center_of_image=True)
                    print("cent_coord: {}".format(cent_coord))

                    # Get the radius of the smallest circle which encloses all the pixels
                    rad = yfu.fov_radius(header=hdul[0].header, 
                                        unit=u.deg)
                    print("rad: {}".format(rad))



                    r_fov = yfu.fov_radius(ccd.header+ccd.wcs.to_header())
                    print(r_fov)
                    ps1 = ypu.PanSTARRS1(ccd.wcs.wcs.crval[0]*u.deg, ccd.wcs.wcs.crval[1]*u.deg, radius=r_fov,
                                        column_filters={"rmag":"10.0..14.5", "e_rmag":"<0.10", "nr":">5"})

                    df_stars_all = ps1.query()
                    df_stars_all = df_stars_all.to_pandas()

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
                    df_stars = df_stars.dropna(subset=["gmag", "rmag"])
                    print(df_stars.columns)

                    pos_stars = np.array([df_stars["RAJ2000"].array, df_stars["DEJ2000"].array]).T
                    pos_stars = SkyCoord(pos_stars, **SKYC_KW).to_pixel(wcs)
                    pos_stars = np.transpose(pos_stars)
                    pos_stars

                    pos_stars_all = np.array([df_stars_all["RAJ2000"].array, df_stars_all["DEJ2000"].array]).T
                    pos_stars_all = SkyCoord(pos_stars_all, **SKYC_KW).to_pixel(wcs)
                    pos_stars_all = np.transpose(pos_stars_all)
                    pos_stars_all



                    # ap_stars = CAp(positions=pos_stars, r=20)
                    # ap_stars_all = CAp(positions=pos_stars_all, r=20)

                    # fig_set = plt.figure(figsize=val_figsize)
                    # ax1 = plt.subplot2grid((1,2), (0,0),
                    #                     fig=fig_set)
                    # im1 = _astro_utilities.norm_imshow(ax1, hdul[0].data, 
                    #                                     zscale=True)
                    # ax1.set_title('Pixel coordinate system', fontsize=9)
                    # ax1.tick_params(labelsize=8)
                    # ap_stars_all.plot(ax1, color='w', lw=2)
                    # ap_stars.plot(ax1, color='r', lw=2)
                    # plt.colorbar(im1, ax = ax1, fraction=val_fraction, pad=0.04)

                    # ax2 = plt.subplot2grid((1,2), (0,1),
                    #                     projection=wcs,
                    #                     fig=fig_set)
                    # im2 = _astro_utilities.norm_imshow(ax2, hdul[0].data, 
                    #                                     zscale=True)
                    # ax2.set_title('World coordinate system', fontsize=9)
                    # ap_stars_all.plot(ax2, color='w', lw=2)
                    # ap_stars.plot(ax2, color='r', lw=2)

                    # ax2.coords.grid(True, color='white', ls=':')
                    # ax2.coords['ra'].set_axislabel('Right Ascension (J2000)', minpad=0.5, fontsize=8)
                    # ax2.coords['ra'].set_ticklabel_position('bl')
                    # ax2.coords['dec'].set_axislabel('Declination (J2000)', minpad=0.4, fontsize=8)
                    # ax2.coords['dec'].set_ticklabel_position('bl')
                    # ax2.coords['ra'].set_major_formatter('hh:mm')
                    # ax2.coords['dec'].set_major_formatter('dd:mm')
                    # ax2.coords['ra'].display_minor_ticks(True)
                    # ax2.coords['dec'].display_minor_ticks(True)
                    # ax2.coords['ra'].set_minor_frequency(1)
                    # ax2.coords['dec'].set_minor_frequency(1)
                    # ax2.tick_params(labelsize=8)

                    # ax2.annotate(f"image center (RA, DEC): ({cent_coord.ra:.03f}, {cent_coord.dec:.03f})\ntelescope center (RA, DEC): ({hdul[0].header['RA']*u.deg:.03f}, {hdul[0].header['DEC']*u.deg:.03f})\noffset (RA, DEC): ({offset_RA:.03f}, {offset_DEC:.03f})\noffset (AZ, ALT): ({offset_AZ:.03f}, {offset_ALT:.03f})",
                    #             xy=(0, 0), xytext=(0.6, -0.1),
                    #             xycoords='axes fraction',
                    #             va='top', ha='left',
                    #             fontsize = 6)

                    # plt.colorbar(im2, ax = ax2, fraction=val_fraction, pad=0.04)
                    # plt.suptitle(f"fname: {fpath.name}\nPS1 Query: Dropping Nearby Stars")

                    # plt.tight_layout()
                    # plt.savefig(f"{DIFFPRESULTDIR/fpath.stem}_PS1_nearby.png")
                    # # plt.show()
                    # plt.close()

                    

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

                    _phot_stars = []

                    for i, row in df_stars.iterrows():
                        pos_star = SkyCoord(row["RAJ2000"], row["DEJ2000"], **SKYC_KW).to_pixel(wcs)
                        ap = CAp([pos_star[0], pos_star[1]], r=R_AP)
                        an = CAn([pos_star[0], pos_star[1]], r_in=R_IN, r_out=R_OUT)
                        _phot_star = ypu.apphot_annulus(hdul[0].data, ap, an, error=yfu.errormap(hdul[0].data))
                        _phot_star[f"{flt}mag"] = row[f"{flt}mag"]
                        _phot_star[f"e_{flt}mag"] = row[f"e_{flt}mag"]
                        _phot_star["grcolor"] = row["grcolor"]
                        _phot_star["e_grcolor"] = row["e_grcolor"]
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

                    cbar1 = plt.colorbar(im1, ax = ax1, fraction=val_fraction, pad=0.04)
                    cbar2 = plt.colorbar(im2, ax = ax2, fraction=val_fraction, pad=0.04, )
                    cbar1.ax.tick_params(labelsize=8)
                    cbar2.ax.tick_params(labelsize=8)

                    plt.suptitle(f"fname: {fpath.name}\n Result of PS1 query for differential photometry", fontsize=10,)

                    plt.tight_layout()
                    plt.savefig(f"{DIFFPRESULTDIR}/{fpath.stem}_PS1_query.png")

                    # plt.show()
                    plt.close()



                    phot_stars = pd.concat(_phot_stars)
                    # phot_stars = phot_stars.loc[phot_stars["objID"] != 110823405221754720].copy()  # star 15
                    # SEE THE LAST CELL IN THIS FILE FOR DESCRIPTION
                    print(len(phot_stars))
                    print(phot_stars)

                    phot_stars_na = phot_stars.dropna()
                    print(len(phot_stars_na))
                    print(phot_stars_na)

                    fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharex=False, sharey=False, gridspec_kw=None)

                    axs[0].plot(phot_stars_na[f"{flt}mag"], phot_stars_na["mag"] - phot_stars_na[f"{flt}mag"], '+')
                    axs[1].plot(phot_stars_na["grcolor"], phot_stars_na["mag"] - phot_stars_na[f"{flt}mag"], '+')
                    _xx = np.linspace(10, 15)
                    axs[2].plot(phot_stars_na[f"{flt}mag"], phot_stars_na["mag"], '+')
                    y0 = np.median(phot_stars_na["mag"] - phot_stars_na[f"{flt}mag"])
                    axs[2].plot(_xx, _xx + y0,
                            label=f"{flt}_inst = {flt} magnitude + ({y0:.01f})"
                            )

                    for _, row in phot_stars_na.iterrows():
                        axs[0].text(row[f"{flt}mag"], row["mag"] - row[f"{flt}mag"], int(row["id"]), fontsize=8)
                        axs[1].text(row["grcolor"], row["mag"] - row[f"{flt}mag"], int(row["id"]), fontsize=8)
                        axs[2].text(row[f"{flt}mag"], row["mag"], int(row["id"]), fontsize=8)
                        
                    axs[0].set(
                        xlabel=f"{flt} magnitude (PS1 to {flt} filter by Tonry+2012)",
                        ylabel=f"{flt}_inst - {flt}"
                    )
                    axs[1].set(
                        xlabel="g - r (PS1)",
                        ylabel=f"{flt}_inst - {flt}"
                    )
                    axs[2].set(
                        xlabel =f"{flt} magnitude (PS1 to {flt} filter by Tonry+2012)",
                        ylabel =f"{flt}_inst"
                    )
                    axs[2].legend()

                    plt.suptitle(f"fname: {fpath.name}\n PS1 check for differential photometry", fontsize=10,)

                    plt.tight_layout()
                    plt.savefig(f"{DIFFPRESULTDIR}/{fpath.stem}_standardization.png")

                    # plt.show()
                    plt.close()



                    FWHM = FWHM_INIT
                    avg, med, std = sigma_clipped_stats(hdul[0].data)  # by default, 3-sigma 5-iteration.
                    thresh = 5. * std

                    DAOfind = DAOStarFinder(
                                            fwhm = FWHM,
                                            threshold=thresh,   # In reality, FWHM must be measured a priori using, e.g., ``ginga``
                                            sharplo=0.2, sharphi=1.0,   # default values 0.2 and 1.0
                                            roundlo=-1.0, roundhi=1.0,  # default values -1 and +1
                                            sigma_radius=1.5,           # default values 1.5
                                            ratio=1.0,                  # 1.0: circular gaussian
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
                    plt.savefig(f"{DIFFPRESULTDIR/fpath.stem}_DAOStarfinder_fwhm_{FWHM}.png")

                    # plt.show()
                    plt.close()


                    apphot_result = apphot(hdul[0].data, apert, method='center')
                    print(type(apphot_result))
                    # df_apphot = pd.DataFrame()
                    # apphot_result
                    df_apphot = apphot_result.to_pandas()
                    print(type(df_apphot))
                    df_apphot



                    cutsizes = 49
                    for idx, row in df_apphot.iterrows():
                        #1. cut asteroia area
                        #print(idx)
                        try :

                            cut_hdu = Cutout2D(
                                        data = hdul[0].data,
                                        position = ((row['xcenter'],row['ycenter'])),
                                        size=(cutsizes, cutsizes) #cut ccd
                                        )
                            avg, med, std = sigma_clipped_stats(cut_hdu.data)  # by default, 3-sigma 5-iteration.

                            fig_set = plt.figure(figsize=(8, 5.5))
                            
                            ax11 = plt.subplot2grid((2, 2), (0,0),
                                        fig=fig_set)
                            #im11 = _astro_utilities.zimshow(ax11, cut_hdu.data)
                            im11 = ax11.imshow(cut_hdu.data,
                                            origin="lower")

                            ax11.plot(round(cutsizes/2), round(cutsizes/2), 'rx')
                            ax11.set_ylabel('pixels')
                            ax11.grid(ls=':')
                            ax11.set_title(f'Star area image', fontsize=9)
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
                            ax12.set_title(f'The new center of Star', fontsize=9)
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
                            
                            ax11.annotate(f"star ID {idx}: ",
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
                            # plt.savefig(f"{DIFFPRESULTDIR/fpath.stem}_Star_{idx:03d}.png")
                            # plt.show()
                            plt.close()
                        except : 
                            continue

                    

                    ap_area  = apert.area
                    ap_area

                    # since our `annul` has many elements,
                    mask_apert = (apert.to_mask(method='center'))
                    mask_annul = (annul.to_mask(method='center'))

                    mag_ann  = np.zeros(len(apphot_result))
                    merr_ann = np.zeros(len(apphot_result))



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
                        df_apphot.at[i, 'msky'] = msky
                        df_apphot.at[i, 'sky_std'] = sky_std
                        df_apphot.at[i, 'nsky'] = nsky
                        df_apphot.at[i, 'nrej'] = nrej
                        df_apphot.at[i, 'flux_star'] = flux_star
                        df_apphot.at[i, 'flux_err'] = flux_err
                        df_apphot.at[i, 'mag_ann'] = mag_ann[i]
                        df_apphot.at[i, 'merr_ann'] = merr_ann[i]

                    sky = wcs.pixel_to_world(df_apphot['xcenter'], df_apphot['ycenter'])
                    sky
                    # df_apphot.to_csv(f"{INSTRESULTDIR}/{fpath.stem}_m_inst.csv")
                    df_apphot


                    sky_coord = wcs.pixel_to_world(df_apphot['xcenter'], df_apphot['ycenter'])
                    sky_coord
                    print(type(sky_coord))
                    #sky_coord[0]

                    dir(sky_coord)
                    len(sky_coord.ra)
                    # df_apphot["RA2000"] = sky_coord.ra
                    # df_apphot["RA2000"]
                    df_RADEC = pd.DataFrame({"RAJ2000": sky_coord.ra.degree, "DEJ2000": sky_coord.dec.degree})
                    df_RADEC
                    #type(df_RADEC["RA2000"][0])

                    df_apphot = pd.concat([df_apphot, df_RADEC], axis=1,)
                    df_apphot['filename'] = fpath.stem
                    df_apphot['t_start'] = t_start
                    df_apphot['t_expos'] = t_expos
                    df_apphot['t_middle'] = t_middle
                    df_apphot["filter"] = f"{flt}"

                    df_apphot["y0"] = y0
                    df_apphot

                    
                    df_apphot[f"{flt}_magnitude"] = df_apphot["mag_ann"] - df_apphot["y0"]

                    df_apphot.to_csv(f"{DIFFPRESULTDIR}/{fpath.stem}_result_photometry.csv")
                    df_apphot[f"{flt}_magnitude"]

                    df_apphot_sub = df_apphot.dropna()
                    df_apphot_sub


                    fig, axs = plt.subplots(1, 2, figsize=(10, 5), sharex=False, sharey=False, gridspec_kw=None)

                    for idx, row in df_apphot_sub.iterrows():

                        axs[0].errorbar(df_apphot_sub["id"], 
                                    df_apphot_sub[f"{flt}_magnitude"], yerr=df_apphot_sub["merr_ann"],
                                    marker='x',
                                    ls='none',
                                    #ms=10,
                                    capsize=3)

                    axs[0].invert_yaxis()

                    style = {'edgecolor': 'white', 'linewidth': 3}
                    axs[1].hist(df_apphot_sub[f"{flt}_magnitude"], 
                                **style)
                        
                    axs[0].set(
                        xlabel='Star ID',
                        ylabel=f"{flt}mag"
                    )
                    axs[1].set(
                        xlabel=f"{flt}mag",
                        ylabel="number of stars"
                    )

                    plt.suptitle(f"fname: {fpath.name}\n Result of differential photometry", fontsize=10,)

                    plt.tight_layout()
                    plt.savefig(f"{DIFFPRESULTDIR}/{fpath.stem}_result_photometry_chart.png")

                    # plt.show()
                    plt.close()
               
        except Exception as err: 
            print("Err :", err)
            continue