# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

이 파일은 fits file의 header에 있는 plate solving 관련 keyword를 삭제해 준다.
필요할때만 사용하면 된다.
"""
#%%
from glob import glob
from pathlib import Path
import os
import numpy as np
import astropy.units as u

from astropy.io import fits
from astropy.modeling.polynomial import Polynomial1D
from astropy.modeling.models import Gaussian1D, Linear1D
from astropy.modeling.fitting import LinearLSQFitter, LevMarLSQFitter
# astroquery provides an interface to the NIST atomic line database
from astroquery.nist import Nist

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import ysfitsutilpy as yfu
#import ysphotutilpy as ypu
#import ysvisutilpy as yvu

import _astro_utilities
import _Python_utilities

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
#%%
#######################################################
BASEDIR = Path("/mnt/Rdata/ASTRO_data")  

# PROJECDIR = BASEDIR / "C1-Variable"
# TODODIR = PROJECDIR / "-_-_-_2016-_-_RiLA600_STX-16803_-_2bin"
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

PROJECDIR = BASEDIR / "C4-Spectra"
TODODIR = PROJECDIR / "-_-_-_2024-05_TEC140_ASI183MMPro_-_1bin"

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(TODODIR))
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

try : 
    BDFDIR = [x for x in DOINGDIRs if "CAL-BDF" in str(x)]
    print ("BDFDIR: ", format(BDFDIR))
    MASTERDIR = Path(BDFDIR[0]) / _astro_utilities.master_dir
    if not MASTERDIR.exists():
        os.makedirs("{}".format(str(MASTERDIR)))
        print("{} is created...".format(str(MASTERDIR)))
    print ("MASTERDIR: ", format(MASTERDIR))
except : 
    pass

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_" in str(x)])
# print ("DOINGDIRs: ", format(DOINGDIRs))
# print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

filter_str = '2024-05-01'
DOINGDIRs = [x for x in DOINGDIRs if filter_str in x]
# remove = 'BIAS'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'DARK'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
# remove = 'FLAT'
# DOINGDIRs = [x for x in DOINGDIRs if remove not in x]
print ("DOINGDIRs: ", DOINGDIRs)
print ("len(DOINGDIRs): ", len(DOINGDIRs))
#######################################################
#%%
#######################################################
ylow, yhigh = 1775, 1910 
wei_ylow, wei_yhigh = 1800, 1880
xhigh = 4200

# start by taking +/- 15 pixels
npixels_to_cut = 15
npixels = 20

guessed_wavelengths = [667.728, 640.2, 585.2]  # Ar, Ne, Ne
guessed_xvals = [600, 970, 1550]

# guessed_wavelengths = [763.511, 750.387, 727.294, 706.722, 696.543,
#                        667.728, 640.225, 614.306, 585.249, 
#                        476.487, 451.073, 420.067, 394.898
#                        ]  
guessed_wavelengths = [763.511, 750.387, 727.294, 
                       706.722, 696.543,
                       614.306, 585.249,
                       420.067, 394.989]
guessed_xvals = [1018, 1127, 1318, 
                 1482, 1560,
                 2200, 2420,
                 3770, 3965]

HBalmer = [656.279,	486.135, 434.0472, 410.1734, 397.0075] * u.nm
guessed_xvals = [1883, 3180, 3570, 3747, 3845]

calfile_idx = 0
#######################################################
#%%
# BDFDIR = Path(BDFDIR[0])
# print("BDFDIR", BDFDIR)

# summary_cal = yfu.make_summary(BDFDIR / "*.fit*")
# if summary_cal is not None : 
#     #print(summary)
#     print("len(summary_cal):", len(summary_cal))
#     print("summary_cal:", summary_cal)
#     #print(summary_cal["file"][0])
    
# fpath_cal = Path(summary_cal["file"][calfile_idx])
#%%
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    print(f"Starting {DOINGDIR}...")

    READINGDIR = DOINGDIR / "reduced"
    TRACERESULTDIR = DOINGDIR / f"{READINGDIR.parts[-1]}_trace_HBalmer"
    if not TRACERESULTDIR.exists() :
        os.mkdir(str(TRACERESULTDIR))
        print(f"{str(TRACERESULTDIR)} is created...")

    summary = yfu.make_summary(READINGDIR / "*.fit*")
    if summary is not None : 
        #print(summary)
        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])

        for _, row  in summary.iterrows():
            try:
                ######################################################################################################
                fpath = Path(row["file"])
                hdul = fits.open(fpath)
                print("*"*20)
                print(f"Starting {fpath.name}")

                ##############################################################################
                ##############################################################################
                ##############################################
                # ylow, yhigh = 1775, 1910 
                ##############################################
                result = f"ylow:{ylow}\nyhigh:{yhigh}\n" 
                result += f"wei_ylow:{wei_ylow}\nwei_yhigh:{wei_yhigh}\n"
                result += f"xhigh:{xhigh}\n"
                result += f"npixels_to_cut:{npixels_to_cut}\n"
                result += f"npixels:{npixels}\n"
                # result += f"calibration file:{fpath_cal}\n"
                result += f"HBalmer:{HBalmer}\n"
                result += f"guessed_xvals:{guessed_xvals}\n"

                yvals = np.argmax(hdul[0].data, axis=0)
                xvals = np.arange(hdul[0].data.shape[1])
                bad_pixels = (yvals < ylow) | (yvals > yhigh)

                fig, axs = plt.subplots(13, 1, figsize=(16, 36), 
                                        sharex=False, sharey=False, 
                                        gridspec_kw={'height_ratios': [14, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4]})

                im0 = _astro_utilities.zimshow(axs[0], 
                                hdul[0].data,
                                origin="lower",
                                cmap = "viridis"
                                )
                axs[0].set_title('The whole mage (zimshow)')
                axs[0].annotate(f"Min value: {hdul[0].data.min()}, Mean value: {hdul[0].data.mean():.02f}, Max value: {hdul[0].data.max()}", 
                        xy=(15, 15), xycoords='axes pixels')
                ##############################################################################
                ##############################################################################
                im1 = _astro_utilities.zimshow(axs[1], 
                                                hdul[0].data[ylow:yhigh, :],
                                                origin="lower",
                                                cmap = "viridis",
                                                )
                axs[1].set_title(f'The cut image with ylow ~ yhigh (zimshow): {ylow} ~ {yhigh} rows')

                ##############################################################################
                ##############################################################################
                im2 = axs[2].plot(hdul[0].data[int((ylow+yhigh)/2),:])
                axs[2].set_xlabel('pixel')
                axs[2].set_ylabel('ADU')
                axs[2].set_xlim((0,hdul[0].data.shape[1]))
                axs[2].set_title(f'Intensity of single row: {int((ylow+yhigh)/2)} row')

                ##############################################################################
                ##############################################################################
                # we use a cutout around the traced line, so the Y-values are from that cutout
                # the `repeat` command here is used to extend our Y-axis position values, which are 425, 426, ... 475
                # along the X-direction.  The indexing with [:, None] adds a "dummy" axis along the second (x) dimension,
                # then `repeat` copies our Y-axis values.  The resulting array has the same shape as our weight array,
                # which is hdul[0].data[425:475, :] minus the median
                yaxis = np.repeat(np.arange(ylow, yhigh)[:,None],
                                hdul[0].data.shape[1], axis=1)
                background = np.median(hdul[0].data)
                # moment 1 is the data-weighted average of the Y-axis coordinates
                weighted_yaxis_values = np.average(yaxis, axis=0,
                                                weights=hdul[0].data[ylow:yhigh,:] - background)
                # print("weighted_yaxis_values ", weighted_yaxis_values)
                # print("weighted_yaxis_values.shape :", weighted_yaxis_values.shape)

                ##############################################################################
                ##############################################################################
                ##############################################
                # wei_ylow, wei_yhigh = 1800, 1880
                ##############################################

                # we need to use the 'extent' keyword to have the axes correctly labeled
                im3 = _astro_utilities.zimshow(axs[3], hdul[0].data[ylow:yhigh,:],
                        extent=[0,hdul[0].data.shape[1],ylow,yhigh],
                        origin="lower",
                        cmap="viridis")
                # plt.gca().set_aspect(10) # we stretch the image out by 10x in the y-direction
                axs[3].plot(xvals[~bad_pixels], yvals[~bad_pixels], 'm+', label="Argmax", markersize=4, alpha=0.2)
                axs[3].plot(xvals, weighted_yaxis_values, 'cx', label="Weighted", markersize=3, alpha=0.2)
                axs[3].set_xlim((0,hdul[0].data.shape[1]))
                axs[3].set_ylim((ylow,yhigh))
                axs[3].set_title(f'Argmax and Weighted value on the image: {ylow} ~ {yhigh} rows')
                axs[3].legend(loc='best');

                ##############################################################################
                ##############################################################################
                bad_moments = (weighted_yaxis_values > wei_yhigh) | (weighted_yaxis_values < wei_ylow)

                im4 = axs[4].plot(xvals[~bad_pixels], yvals[~bad_pixels], 'm+', label="Argmax", markersize=4, alpha=0.2)
                axs[4].plot(xvals[~bad_moments], weighted_yaxis_values[~bad_moments], 'cx', label="Weighted", markersize=3, alpha=0.2)
                axs[4].set_xlim((0,hdul[0].data.shape[1]))
                axs[4].set_ylim((ylow,yhigh))
                axs[4].set_title(f'Argmax and Weighted value with wei_ylow ~ wei_yhigh: {wei_ylow} ~ {wei_yhigh} rows')
                axs[4].legend(loc='best');

                ##############################################################################
                ##############################################################################
                polymodel = Polynomial1D(degree=3)
                linfitter = LinearLSQFitter()

                fitted_polymodel = linfitter(polymodel, xvals[~bad_moments], weighted_yaxis_values[~bad_moments])
                fitted_polymodel

                im5 = _astro_utilities.zimshow(axs[5], hdul[0].data[ylow:yhigh,:], 
                        extent=[0,hdul[0].data.shape[1],ylow,yhigh],
                        origin="lower",
                        cmap="viridis")
                axs[5].plot(xvals[~bad_moments], weighted_yaxis_values[~bad_moments], 'x', label="Weighted", markersize=3, alpha=0.2)
                axs[5].plot(xvals, fitted_polymodel(xvals), color='r', label="fitted_polymodel");
                axs[5].set_title(f'fitted_polymodel wity ylow ~ yhigh on the imnage: {ylow} ~ {yhigh} rows')
                axs[5].annotate(f'fitted_polymodel parameters: {fitted_polymodel.parameters}', fontsize=8,
                                        xy=(0, 0), xytext=(-10, -30), va='top', ha='left',
                                        xycoords='axes fraction', textcoords='offset points')
                axs[5].legend(loc='best')

                ##############################################################################
                ##############################################################################
                #######################################
                # xhigh = 4200
                #######################################

                fitted_polymodel = linfitter(polymodel, xvals[(~bad_moments) & (xvals < xhigh)],
                                                weighted_yaxis_values[(~bad_moments) & (xvals < xhigh)])
                fitted_polymodel

                im6 = _astro_utilities.zimshow(axs[6], hdul[0].data[ylow:yhigh,:], extent=[0,hdul[0].data.shape[1],ylow,yhigh],
                                origin="lower",
                                cmap="viridis")
                axs[6].plot(xvals[~bad_moments], weighted_yaxis_values[~bad_moments], 'x', markersize=3, alpha=0.25)
                axs[6].plot(xvals, fitted_polymodel(xvals), color='r', label="fitted polymodel with xhigh");
                axs[6].set_title(f'Fitted polymodel with ylow ~ yhigh with xhigh on the image : {ylow} ~ {yhigh} rows, xhigh: {xhigh}')
                axs[6].axis((0,hdul[0].data.shape[1],ylow,yhigh));
                axs[6].annotate(f'fitted_polymodel parameters with xhigh : {fitted_polymodel.parameters}, xhigh ({xhigh})', fontsize=8,
                                        xy=(0, 0), xytext=(-10, -30), va='top', ha='left',
                                        xycoords='axes fraction', textcoords='offset points')
                axs[6].legend(loc='best')

                ##############################################################################
                ##############################################################################
                im7 = axs[7].plot(xvals[~bad_moments & (xvals < xhigh)],
                        weighted_yaxis_values[~bad_moments & (xvals < xhigh)] - fitted_polymodel(xvals[~bad_moments & (xvals < xhigh)]), 'x', markersize=3 )
                axs[7].plot(xvals[~bad_moments & (xvals > xhigh)],
                        weighted_yaxis_values[~bad_moments & (xvals > xhigh)] - fitted_polymodel(xvals[~bad_moments & (xvals > xhigh)]), 'r+', alpha=0.5)
                axs[7].set_ylim(-55, 55)
                axs[7].set_xlim((0,hdul[0].data.shape[1]))
                axs[7].set_ylabel("Residual (data-model)");
                axs[7].set_title(f'Residual (data-model)')
                axs[7].legend(loc='best')

                ##############################################################################
                ##############################################################################
                # we need to use the 'extent' keyword to have the axes correctly labeled
                #######################################
                # start by taking +/- 15 pixels
                # npixels_to_cut = 15
                #######################################
                im8 = _astro_utilities.zimshow(axs[8], hdul[0].data[ylow:yhigh,:], extent=[0,hdul[0].data.shape[1],ylow,yhigh],
                        origin="lower")
                axs[8].fill_between(xvals, fitted_polymodel(xvals)-npixels_to_cut,
                                fitted_polymodel(xvals)+npixels_to_cut,
                                color='orange', alpha=0.5, 
                                label=f'fitted polmodel with xhigh +/- {npixels_to_cut} pixels')
                axs[8].set_title(f'fitted polimodel with xhigh +/- {npixels_to_cut} pixels')
                axs[8].set_xlim((0,hdul[0].data.shape[1]));
                axs[8].set_ylim((ylow,yhigh));
                axs[8].legend(loc='best')

                ##############################################################################
                ##############################################################################
                trace_center = fitted_polymodel(xvals)
                cutouts = np.array([hdul[0].data[int(yval)-npixels_to_cut:int(yval)+npixels_to_cut, ii]
                                        for yval, ii in zip(trace_center, xvals)])

                im9 = _astro_utilities.zimshow(axs[9], hdul[0].data[ylow:yhigh,:], extent=[0,hdul[0].data.shape[1],ylow,yhigh],
                        origin="lower", cmap="viridis")
                axs[9].set_title("We go from this...")

                ##############################################################################
                ##############################################################################

                im10 = _astro_utilities.zimshow(axs[10], cutouts.T, 
                        origin="lower", cmap="viridis")
                axs[10].set_xlim((0,hdul[0].data.shape[1]))
                axs[10].set_title("...to this")

                ##############################################################################
                ##############################################################################
                mean_trace_profile = (cutouts - background).mean(axis=0)
                trace_profile_xaxis = np.arange(-npixels_to_cut, npixels_to_cut)

                lmfitter = LevMarLSQFitter()
                guess = Gaussian1D(amplitude=mean_trace_profile.max(), mean=0, stddev=5)
                fitted_trace_profile = lmfitter(model=guess, x=trace_profile_xaxis, y=mean_trace_profile)
                model_trace_profile = fitted_trace_profile(trace_profile_xaxis)

                im11 = axs[11].plot(trace_profile_xaxis, mean_trace_profile, label='data')
                axs[11].plot(trace_profile_xaxis, model_trace_profile, label='model')
                axs[11].legend(loc='best')
                axs[11].set_xlabel("Distance from center")
                axs[11].set_ylabel("Average source profile");
                axs[11].set_title(f'Average source profile')

                ##############################################################################
                ##############################################################################
                average_spectrum = (cutouts - background).mean(axis=1)
                trace_avg_spectrum = np.array([np.average(
                        hdul[0].data[int(yval)-npixels_to_cut:int(yval)+npixels_to_cut, ii] - background,
                        weights=mean_trace_profile)
                                                for yval, ii in zip(trace_center, xvals)])
                gaussian_trace_avg_spectrum = np.array([np.average(
                        hdul[0].data[int(yval)-npixels_to_cut:int(yval)+npixels_to_cut, ii] - background,
                        weights=model_trace_profile)
                                                for yval, ii in zip(trace_center, xvals)])

                im12 = axs[12].plot(average_spectrum, label="Direct Average")
                axs[12].plot(trace_avg_spectrum, label='Trace-weighted average')
                axs[12].plot(gaussian_trace_avg_spectrum, label='Gaussian-model-Trace-weighted average', alpha=0.5, linewidth=0.5, color='r')
                axs[12].set_xlabel("pixel")
                axs[12].set_xlim((0,hdul[0].data.shape[1]))
                axs[12].set_title(f'Spectrum (raw)')

                axs[1].set_aspect(2) # we stretch the image out by 10x in the y-direction 
                axs[3].set_aspect(2) # we stretch the image out by 10x in the y-direction 
                axs[5].set_aspect(2) # we stretch the image out by 10x in the y-direction 
                axs[6].set_aspect(2) # we stretch the image out by 10x in the y-direction 
                axs[7].set_aspect(2) # we stretch the image out by 10x in the y-direction 
                axs[8].set_aspect(2) # we stretch the image out by 10x in the y-direction 
                axs[9].set_aspect(2) # we stretch the image out by 10x in the y-direction 
                axs[10].set_aspect(4) # we stretch the image out by 10x in the y-direction 

                plt.suptitle(f'filename : {fpath.name}')

                plt.colorbar(im0, fraction=0.025, pad=0.03, location='bottom')
                plt.tight_layout(pad=2.0)
                plt.savefig(f"{TRACERESULTDIR}/{fpath.stem}_spctrum_location.png")
                                
                # plt.show()


                ######################################################################################################
                trace_model = fitted_polymodel
                xaxis = np.arange(hdul[0].data.shape[1])
                trace_center = trace_model(xaxis)
                yaxis = np.arange(-npixels_to_cut, npixels_to_cut)

                cutouts = np.array([hdul[0].data[int(yval)-npixels_to_cut:int(yval)+npixels_to_cut, ii]
                                        for yval, ii in zip(trace_center, xvals)])


                HBalmer = [656.279,	486.135, 434.0472, 410.1734, 397.0075] * u.nm
                guessed_xvals = [1883, 3180, 3570, 3747, 3845]

                improved_xval_guesses = [np.average(xaxis[g-npixels:g+npixels],
                                                weights=average_spectrum[g-npixels:g+npixels] - np.median(average_spectrum))
                                        for g in guessed_xvals]
                print("improved_xval_guesses :", improved_xval_guesses)
                result += f"improved_xval_guesses:{improved_xval_guesses}\n"


                fig, axs = plt.subplots(3, 1, 
                                        figsize=(18, 10), 
                                        sharex=False, sharey=False, 
                                        gridspec_kw={'height_ratios': [2, 1, 2]},
                                        )

                im0 = axs[0].plot(average_spectrum, label="Direct Average")
                axs[0].plot(trace_avg_spectrum, label='Trace-weighted average')
                axs[0].plot(gaussian_trace_avg_spectrum, label='Gaussian-model-Trace-weighted average', alpha=0.5, linewidth=0.5, color='r')
                axs[0].set_xlabel("pixel")
                axs[0].set_xlim((0,hdul[0].data.shape[1]))
                axs[0].xaxis.set_major_locator(ticker.AutoLocator())
                axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
                axs[0].grid(which='both', axis='both', ls='dotted')
                axs[0].set_title(f'Spectrum (raw)')

                axs[0].plot(guessed_xvals, [100]*len(guessed_xvals), 'x')
                axs[0].vlines(guessed_xvals, 0, 40000, 'b', alpha=0.5)
                axs[0].set_xlim((0,hdul[0].data.shape[1]))
                axs[0].annotate(f"improved_xval_guesses: {improved_xval_guesses}\nH-Balmer: {HBalmer}", 
                        xy=(0, -0.5), xycoords='axes fraction', #     xy=(0, -15), xycoords='axes pixels'
                        )
                axs[0].xaxis.set_major_locator(ticker.AutoLocator())
                axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
                axs[0].set_title(f'Spectrum for calibration H-Balmer (raw)')

                ##############################################################################

                linfitter = LinearLSQFitter()
                wlmodel = Linear1D()
                linfit_wlmodel = linfitter(model=wlmodel, x=improved_xval_guesses, y=HBalmer.value)
                wavelengths = linfit_wlmodel(xaxis) * u.nm
                print(f"linfit_wlmodel: {linfit_wlmodel}")
                print(f"wavelengths: {wavelengths}")

                minwave = wavelengths.min()
                maxwave = wavelengths.max()

                im1 = axs[1].plot(improved_xval_guesses, HBalmer, 'r+', markersize=8)
                axs[1].plot(xaxis, wavelengths, '-')
                axs[1].set_ylabel(r"$\lambda(x)$")
                axs[1].set_xlabel(r"x (pixels)")
                axs[1].set_xlim((0,hdul[0].data.shape[1]))
                axs[1].xaxis.set_major_locator(ticker.AutoLocator())
                axs[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
                axs[1].grid(which='both', axis='both', ls='dotted')

                axs[1].annotate(f"linfit_wlmodel: {linfit_wlmodel}", 
                        xy=(0, -0.5), xycoords='axes fraction', #     xy=(0, -15), xycoords='axes pixels'
                        )

                ##############################################################################
                im2 = axs[2].plot(wavelengths, average_spectrum, label="Direct Average")
                axs[2].plot(wavelengths, trace_avg_spectrum, label='Trace-weighted average')
                axs[2].plot(wavelengths, gaussian_trace_avg_spectrum, label='Gaussian-model-Trace-weighted average', alpha=0.5, linewidth=0.5, color='r')
                axs[2].vlines(HBalmer, 0, 40000, 'r', alpha=0.5)
                axs[2].set_title(f'Spectrum with wavelength (raw)')
                axs[2].set_xlabel("wavelength (nm)")
                axs[2].legend(loc='best')
                axs[2].xaxis.set_major_locator(ticker.AutoLocator())
                axs[2].xaxis.set_minor_locator(ticker.AutoMinorLocator())
                axs[2].grid(which='both', axis='both', ls='dotted')

                plt.suptitle(f'filename : {fpath.name}')

                plt.tight_layout(pad=2.0)
                plt.savefig(f"{TRACERESULTDIR}/{fpath.stem}_calibration_using_HBalmer.png")
                # plt.show()
                # plt.closed()

                print(type(fitted_polymodel))
                print(dir(fitted_polymodel))
                result += f"fitted_polymodel.parameters:{fitted_polymodel.parameters}\n"
                result += f"wavelengths(unit):{wavelengths.unit.to_string()}\n"
                result += f"wavelengths:{wavelengths.value.tolist()}\n"
                result += f"Direct average:{average_spectrum.tolist()}\n"
                result += f"Trace-weighted average:{trace_avg_spectrum.tolist()}\n"
                result += f"Gaussian-model-Trace-weighted average:{gaussian_trace_avg_spectrum.tolist()}\n"
                # result += f"Direct average_cal:{average_spectrum.tolist()}\n"
                # result += f"Trace-weighted average_cal:{trace_avg_spectrum.tolist()}\n"
                # result += f"Gaussian-model-Trace-weighted average_cal:{gaussian_trace_avg_spectrum.tolist()}\n"

                with open(f"{TRACERESULTDIR}/{fpath.stem}_spctrum_result_using_Hbalmer.txt", 'w') as f:
                    f.write(result) 
                                                

            except Exception as err :
                print("X"*60)
                # _Python_utilities.write_log(err_log_file, str(err), verbose=verbose)
                print(f"err messgae : {err}")
            