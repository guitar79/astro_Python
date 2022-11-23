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
cd ~/Downloads/ysphotutilpy && git pull && pip install -e . 
cd ~/Downloads/SNUO1Mpy && git pull && pip install -e . 

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

import Python_utilities
import astro_utilities

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
# 
BASEDIR = "../RnE_2022/"
BASEDIR = "../RnE_2022/RiLA600_STX-16803_2bin/"

BASEDIR = astro_utilities.base_dir
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
BASEDIRs = sorted(Python_utilities.getFullnameListOfsubDir(BASEDIR))
print ("BASEDIRs: {}".format(BASEDIRs))

for BASEDIR in BASEDIRs[4:5]:
    print ("Starting...\n{}".format(BASEDIR))

    BASEDIR = Path(BASEDIR)

    SOLVEDDIR = BASEDIR / astro_utilities.solved_dir2
    AsteroidRESULTDIR = BASEDIR / astro_utilities.Asteroid_result_dir

    if not AsteroidRESULTDIR.exists():
        os.makedirs("{}".format(str(AsteroidRESULTDIR)))
        print("{} is created...".format(str(AsteroidRESULTDIR)))

    #%%
    summary = yfu.make_summary(SOLVEDDIR/"*.fits")
    #print(summary)
    #print("len(summary):", len(summary))
    #print(summary["file"][0])

    for filt in ["r", "v", "b"]:
        summary_filt = summary.loc[summary["FILTER"] == filt].copy()
        
        if summary_filt.empty:
            print("The dataframe(summary_filt) is empty")
            pass
        else:
            print("len(summary_filt):", len(summary_filt))
            print("summary_filt:", summary_filt)

            for fname in summary_filt["file"][:]:
                #fpath = summary["file"][1]
                print ("Starting...\nfname: {}".format(fname))
                fpath = Path(fname)
        
                #%%
                # load as ccd
                ccd = yfu.load_ccd(fname, 
                                unit="adu")
                
                #%%
                ## Get the Ephemeris
                _, eph, _ = ypu.horizons_query(
                                            OBJID, 
                                            epochs=Time(ccd.header["DATE-OBS"]).jd, 
                                            location=LOCATION)
                print("eph:", eph)
                #%%
                eph.write("{}_eph-{}.csv"\
                        .format(str(AsteroidRESULTDIR / fpath.stem), OBJID),
                        overwrite = True)

                #%%
                # Initial Photometry of the Target
                pos_targ_init = SkyCoord(eph["RA"], 
                                        eph["DEC"], 
                                        **SKYC_KW).to_pixel(ccd.wcs)
                ap = CAp([pos_targ_init[0][0], 
                        pos_targ_init[1][0]], 
                        r=R_AP)
                an = CAn([pos_targ_init[0][0], 
                        pos_targ_init[1][0]], 
                        r_in=R_IN, 
                        r_out=R_OUT)
                print("pos_targ_init:", pos_targ_init)

                phot_targ = ypu.apphot_annulus(ccd, 
                                                ap, 
                                                an, 
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
                #%%
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
                                            figsize=(10, 10), 
                                            sharex=False, 
                                            sharey=False, 
                                            gridspec_kw=None
                                            )

                    yvu.norm_imshow(axs, ccd, zscale=True)
                    # ap = CircularAperture([pos_targ_init[0][0], 
                    #                         pos_targ_init[1][0]], 
                    #                         r=R_AP)

                    # an = CircularAnnulus([pos_targ_init[0][0], 
                    #                         pos_targ_init[1][0]], 
                    #                         r_in=R_IN, 
                    #                         r_out=R_OUT)
                    ap.plot(axs, color="r")
                    an.plot(axs, color="b")

                    # 각 별의 측광을 수행
                    _phot_stars = []
                    for idx, row in df_stars.iterrows():
                        print("Starting photometry star {}:".format(idx))
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
                    plt.title("{} and {} stars".format(eph["targetname"][0], 
                            len(df_stars)), fontsize = 16)
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

                        pos_star = SkyCoord(row["RAJ2000"], 
                                            row["DEJ2000"], 
                                            **SKYC_KW).to_pixel(ccd.wcs)
                        ap = CAp([pos_star[0]-centerdx, 
                                pos_star[1]-centerdy], 
                                r = R_AP)
                        an = CAn([pos_star[0]-centerdx, 
                                pos_star[1]-centerdy], 
                                r_in = R_IN, 
                                r_out = R_OUT)
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
                    plt.title("{} and {} stars".format(eph["targetname"][0], 
                            len(df_stars)), fontsize = 16)
                    plt.tight_layout()
                    plt.savefig("{}_stars_Re.png".format(str(AsteroidRESULTDIR / fpath.stem)))
                    #plt.show()
                    #%%
                    phot_stars = pd.concat(_phot_stars)
                    # phot_stars = phot_stars.loc[phot_stars["objID"] != 110823405221754720].copy()  # star 15
                    # SEE THE LAST CELL IN THIS FILE FOR DESCRIPTION
                    print("phot_stars: ", phot_stars)
                    phot_stars.to_csv("{}_phot_stars_Re.csv".format(str(AsteroidRESULTDIR / fpath.stem)))
                        

                    # %%
                    # Standardization Plots
                    #import seaborn as sns
                    fig, axs = plt.subplots(1, 1, 
                                        figsize=(8, 8), 
                                        sharex=False, 
                                        sharey=False, 
                                        gridspec_kw=None)

                    _xx = np.linspace(12, 14.5)
                    axs.plot(phot_stars["Rmag"], phot_stars["mag"], '+')
                    axs.axhline(phot_targ["mag"].values, label="Kleopatra, instrumental mag")
                    axs.plot(_xx, _xx + np.median(phot_stars["mag"] - phot_stars["Rmag"]))
                    #axs.plot(x, x + np.median(y-x))
                    # IDEA: y = 1*x + c  -> c = mean(y - x)



                    # axs.plot(phot_stars["Rmag"], phot_stars["mag"], '+')
                    # sns.scatterplot(
                    #                     x = "Rmag",
                    #                     y = "mag", 
                    #                     marker = '+',
                    #                     ax = axs,
                    #                     data = phot_stars
                    #                     )
                    # p = sns.regplot(
                    #                 x = "Rmag",
                    #                 y = "mag", 
                    #                 ax = axs,
                    #                 data = phot_stars
                    #                 )
                    
                    
                    

                    for _, row in phot_stars.iterrows():
                        axs.text(row["Rmag"], row["mag"], int(row["id"]), fontsize=8)

                    axs.set(
                        xlabel="R magnitude (PS1 to R_C filter by Tonry+2012)",
                        ylabel="R_inst",
                        title = "R_mag of {} and {} stars".format(eph["targetname"][0], len(df_stars))
                    )
                    axs.grid(ls=':')
                    axs.legend()

                    plt.tight_layout()
                    #plt.show()
                    plt.savefig("{}_R_mag.png".format(str(AsteroidRESULTDIR / fpath.stem)))
                    
                    # %%
                    fig, axs = plt.subplots(1, 2, 
                                            figsize=(16, 8), 
                                            sharex=False, 
                                            sharey=False, 
                                            gridspec_kw=None)

                    axs[0].plot(phot_stars["Rmag"], phot_stars["mag"] - phot_stars["Rmag"], '+')
                    axs[1].plot(phot_stars["grcolor"], phot_stars["mag"] - phot_stars["Rmag"], '+')
                    for _, row in phot_stars.iterrows():
                        axs[0].text(row["Rmag"], row["mag"] - row["Rmag"], int(row["id"]), fontsize=8)
                        axs[1].text(row["grcolor"], row["mag"] - row["Rmag"], int(row["id"]), fontsize=8)
                        
                    axs[0].set(
                        xlabel="R magnitude (PS1 to R_C filter by Tonry+2012)",
                        ylabel="R_inst - R"
                    )
                    axs[0].grid(ls=':')
                    
                    axs[1].set(
                        xlabel="g - r (PS1)",
                        ylabel="R_inst - R"
                    )
                    axs[1].grid(ls=':')
                    plt.tight_layout()
                    #plt.show()
                    plt.savefig("{}_result.png".format(str(AsteroidRESULTDIR / fpath.stem)))

                    
            # %%
