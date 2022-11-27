# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 23:15:44 2018

@author: user

https://github.com/ysBach/SNU_AOclass/blob/master/Notebooks/04-Aperture_Phot_01.ipynb

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
from matplotlib.colors import LogNorm
from matplotlib import rcParams

import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clip, sigma_clipped_stats

from astropy.nddata import Cutout2D

from photutils.aperture import CircularAperture as CAp
from photutils.aperture import CircularAnnulus as CAn
from photutils.centroids import centroid_com
from photutils import aperture_photometry as apphot

import warnings

from ccdproc import CCDData, ccd_process

from astropy.time import Time
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

import Python_utilities
import astro_utilities

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
BASEDIR = "../RnE_2022/"
BASEDIR = "../RnE_2022/RiLA600_STX-16803_2bin/"

BASEDIR = astro_utilities.base_dir
#%%
#####################################################################
# Our object (will be queried to JPL HORIZONS)
#OBJID = '216' # Kleopatra

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
#######################################################

#%%
BASEDIRs = sorted(Python_utilities.getFullnameListOfsubDir(BASEDIR))
print ("BASEDIRs: {}".format(BASEDIRs))

for BASEDIR in BASEDIRs[:]:
    print ("Starting...\n{}".format(BASEDIR))

    BASEDIR = Path(BASEDIR)
    
    SOLVEDDIR = BASEDIR / astro_utilities.solved_dir2
    APhRESULTDIR = BASEDIR / astro_utilities.APh_result_dir

    if not APhRESULTDIR.exists():
        os.makedirs("{}".format(str(APhRESULTDIR)))
        print("{} is created...".format(str(APhRESULTDIR)))

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
            print("summary_filt:", summary_filt)
            print("len(summary_filt):", len(summary_filt))

            for fname in summary_filt["file"][:]:
                #fpath = summary["file"][1]
                print ("Starting...\nfname: {}".format(fname))
                fpath = Path(fname)

                #%%
                # load as ccd
                ccd = yfu.load_ccd(fpath, 
                                unit="adu")

                #%%
                # 별의 목록을 가져옴.
                r_fov = yfu.fov_radius(ccd.header+ccd.wcs.to_header())
                print("r_fov:", r_fov)
                
                ps1 = ypu.PanSTARRS1(ccd.wcs.wcs.crval[0]*u.deg, 
                                    ccd.wcs.wcs.crval[1]*u.deg, 
                                    radius = r_fov,
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

                    df_stars.to_csv("{}_stars.csv".format(str(APhRESULTDIR / fpath.stem)))
                    print("df_stars:", df_stars)

                    #%%
                    fig, axs = plt.subplots(1, 1, 
                                figsize=(10, 10), 
                                sharex=False, 
                                sharey=False, 
                                gridspec_kw=None
                                )
                    yvu.norm_imshow(axs, ccd, zscale=True)

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
                        
                        axs.text(pos_star[0]+10, pos_star[1]+10, 
                                f"star {idx}", fontsize=8)
                        ap.plot(axs, color="orange")
                        an.plot(axs, color="w")
                    plt.title("Marking Stars usigs star catalogue",
                                fontsize = 14)
                    plt.tight_layout()
                    plt.savefig("{}_stars.png".format(str(APhRESULTDIR / fpath.stem)))
                    #plt.show()
                    #%%
                    phot_stars = pd.concat(_phot_stars)
                    # phot_stars = phot_stars.loc[phot_stars["objID"] != 110823405221754720].copy()  # star 15
                    # SEE THE LAST CELL IN THIS FILE FOR DESCRIPTION
                    print("phot_stars: ", phot_stars)
                    phot_stars.to_csv("{}_phot_stars.csv".format(str(APhRESULTDIR / fpath.stem)))

                    #%%
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

                        #4. Putting Aperture and Annulus
                        bigcenter = [bigcutsizes/2 - centerdx, bigcutsizes/2 - centerdy]
                        #fwhm = 4
                        fwhm = FWHM_INIT
                        r_ap = 2 * fwhm
                        r_in = 4 * fwhm
                        r_out = 6 * fwhm
                        ap = CAp(positions = bigcenter, 
                                r = r_ap)
                        an = CAn(positions = bigcenter, 
                                r_in = r_in, r_out = r_out)
                        print("ap", ap)
                        print("type(ap)", type(ap))
                        print("an", an)
                        print("type(an)", type(an))

                        # 5. Estimating Sky
                        sky_mask = an.to_mask(method = 'center')

                        try:  # prior to photutils 0.7
                            sky_vals = sky_mask[0].multiply(bigcut_hdu.data)
                        except TypeError:
                            sky_vals = sky_mask.multiply(bigcut_hdu.data)

                        sky_vals = sky_vals[sky_vals > 0]
                        avg, med, std = sigma_clipped_stats(
                                                    sky_vals, 
                                                    sigma=3, 
                                                    maxiters=10, 
                                                    std_ddof=1
                                                    )

                        if med - avg < 0.3 * std:
                            msky = med
                        else:
                            msky = 2.5 * med - 1.5 * avg

                        print(f"Sky estimation: {msky:.3f} +- {std:.3f}")

                        phot = apphot(
                                    data = bigcut_hdu.data, 
                                    apertures = ap
                                    )
                        phot["sky"] = msky
                        
                        try:  # prior to photutils 0.7
                            phot["source_sum"] = phot["aperture_sum"] - ap.area() * phot["sky"]
                        
                        except TypeError:
                            phot["source_sum"] = phot["aperture_sum"] - ap.area * phot["sky"]
                            
                        phot["inst_mag"] = -2.5 * np.log10(phot["source_sum"] / ccd.header["EXPTIME"])
                        print("phot:", phot)

                        #%%
                        fig = plt.figure(figsize=(12, 12))
                        ax1 = plt.subplot(2, 2, 1)
                        ax1.imshow(cut_hdu.data, 
                                    origin='lower'
                                    )
                        ax1.set_ylabel('pixels')
                        ax1.grid(ls=':')
                        ax1.set_title('Star area image')
                        ax1.text( 0, -2, 
                                'sum: {0:.01f}, mean: {1:.01f}, std: {2:.01f} \nmax: {3:.01f}, min: {4:.01f} \nNumber of Pixel: {5:.0f}x{6:.0f}'\
                                .format(np.sum(cut_hdu.data), 
                                        np.mean(cut_hdu.data), 
                                        np.std(cut_hdu.data), 
                                        np.max(cut_hdu.data), 
                                        np.min(cut_hdu.data),
                                        np.shape(cut_hdu.data)[0],
                                        np.shape(cut_hdu.data)[1]),
                                va = 'top')

                        ax2 = plt.subplot(2, 2, 2)
                        ax2.grid(ls=':')
                        ax2.set_title('The center of Star')
                        ax2.imshow(mask_3sig.astype(int), 
                                    origin="lower")
                        ax2.imshow(cut_hdu.data, alpha=0.4, 
                                    origin="lower")
                        ax2.plot(*center, 'rx')
                        ax2.text(0, -2, 
                                'center: {0:.01f}, {1:.1f}'\
                                .format(center[0], center[1]),
                                va = 'top'
                                )
                        ax2.text(0, -3, 
                                "center dx, dy: {}, {}"\
                                .format(centerdx, centerdy),
                                va = 'top'
                                )
                                    
                        ax3 = plt.subplot(2, 2, 3)
                        ax3.grid(ls=':')
                        ax3.set_title('The result of photometry star {}'.format(idx))
                        ax3.imshow(bigcut_hdu.data,
                                    origin='lower' )
                        ap.plot(ax3, color='r', lw=2)
                        an.plot(ax3, color='w', lw=2)
                        ax3.plot(*bigcenter, 'rx')
                        ax3.text(0, -7, 
                                'center: {0:.01f}, {1:.1f}'\
                                .format(bigcenter[0], bigcenter[1]),
                                va = 'top'
                                )
                        ax3.text(0, -10, 
                                'aperture sum: {0:.01f}, sky: {1:.01f}, source sum: {2:.01f}\n initrument magnitude: {3:.01f}'\
                                .format(float(phot["aperture_sum"]),
                                        float(phot["sky"]),
                                        float(phot["source_sum"]),
                                        float(phot["inst_mag"])),
                                        va = 'top'
                                )
                        #ax3.colorbar()

                        ax4 = plt.subplot(2, 2, 4)
                        ax4.set_title('The histrogram of sky value')
                        ax4.hist(sky_vals, 50, histtype='step')
                        ax4.axvline(msky, ls=':', color='r')


                        plt.tight_layout()
                        plt.savefig("{}_star_{:02d}.png".format(str(APhRESULTDIR / fpath.stem), idx))
