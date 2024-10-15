# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 23:15:44 2018
@author: user

"""
#%%
import os
from warnings import warn
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats

from astroquery.jplhorizons import Horizons

from photutils.aperture import CircularAperture as CAp
from photutils.aperture import CircularAnnulus as CAn
from photutils.detection import DAOStarFinder
from photutils.psf.groupstars import DAOGroup
from astropy.nddata import CCDData, Cutout2D

from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec

import _astro_utilities
import _Python_utilities
import _tool_visualization

from astropy.table import Table, vstack
from astroquery.vizier import Vizier
from astroquery.mast import Catalogs
from scipy.interpolate import UnivariateSpline
from xarray import Coordinate

import ysfitsutilpy as yfu
import ysphotutilpy as ypu

import warnings

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

BASEDIR = "../Post_processing/M35_Light_-_2018-10-31_-_TMB130ss_STF-8300M_-_1bin/"
BASEDIR = "../RnE_2022/KLEOPATRA_Light_-_2022-10-07_-_GSON300_STF-8300M_-_1bin/"

### make all fits file list...
fullnames = _Python_utilities.getFullnameListOfallFiles("{}/input".format(BASEDIR))
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))

c_method = 'median'
master_dir = "master_files/"
reduced_dir = "readuced_files/"
result_dir = "result_files/"

if not os.path.exists('{0}'.format("{}{}".format(BASEDIR, result_dir))):
    os.makedirs("{}{}".format(BASEDIR, result_dir))
    print("{}{}is created".format(BASEDIR, result_dir))

fullnames_light = [w for w in fullnames \
            if ("_bias_" not in w.lower()) \
                and ("_dark_" not in w.lower()) \
                    and ("_flat_" not in w.lower())]
print ("len(fullnames_light): {}".format(len(fullnames_light)))

#%%
n = 0
for fullname in fullnames_light[:]:
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames_light), 
                                            (n/len(fullnames_light))*100, os.path.basename(__file__)))
    print ("Starting...\nfullname: {}".format(fullname))

fullname = fullnames_light[2]
    
fullname_el = fullname.split("/")
hdul = fits.open(fullname)
hdr = hdul[0].header
data = hdul[0].data

img = hdul[0].data
print("img: {}".format(img))
print("img.shape: {}".format(img.shape))

# Set WCS and print for your information
w = WCS(hdr)
print("WCS: {}".format(w))

#%%
# Query object
objname = "Kleopatra"
observat = {"lon": 127.0, "lat": 37.5, "elevation": 101}
# Observatory code: see https://en.wikipedia.org/wiki/List_of_observatory_codes
# dict(lon=, lat=, elevation=)

#t_obs = Time(hdr["DATE-OBS"]) + hdr["EXPTIME"] * u.s / 2  # middle of observation time
t_obs = Time(hdr["DATE-OBS"]) + hdr["EXPOSURE"] * u.s / 2  # middle of observation time

obj = Horizons(id=objname, location=observat, epochs=t_obs.jd)
obj_q = obj.ephemerides()

print("obj_q: {}".format(obj_q))

pos_sky = SkyCoord(obj_q["RA"][0], obj_q["DEC"][0], unit='deg')
pos_pix = pos_sky.to_pixel(wcs=w)

print("pos_sky: {}".format(pos_sky))
print("pos_pix: {}".format(pos_pix))

#%%
# Position of the telescope FOV center 
# (RA/DEC of the pixel at the center)
cent_coord = yfu.center_radec(ccd_or_header=hdr, center_of_image=True)
print("cent_coord: {}".format(cent_coord))

# Get the radius of the smallest circle which encloses all the pixels
rad = yfu.fov_radius(header=hdr, unit=u.deg)
print("rad: {}".format(rad))

# Initialize PanSTARRS1 class
q = ypu.PanSTARRS1(
    ra=cent_coord.ra.value, 
    dec=cent_coord.dec.value, 
    radius=rad,
    column_filters={"gmag":"13.0..20.0", "e_gmag":"<0.10"}
)

# Query to the website (VizieR)
# This is where the most of the time is spent.
q.query()

# Only select the stars within 50-pixel bezel in the FOV.
q.select_xyinFOV(header=hdr, bezel=50)

# Remove objects not suitable for differential photometry (see description below)
q.drop_for_diff_phot(drop_by_Kron=True)

# Remove redundant columns, remove objects with too few observations:
q.select_filters(filter_names=['g', 'r', 'i'], n_mins=5)
# You can try a list of ``n_mins``:
# q.select_filters(filter_names=['g', 'r', 'i'], n_mins=[10, 3, 5])

q_stars_orig = q.queried.copy()
pos_stars_orig = np.array([q_stars_orig["x"], q_stars_orig["y"]]).T
q_stars_orig

fwhm = 4
q.drop_star_groups(crit_separation=6*fwhm)
q_stars_diropped = q.queried.copy()  # This will be overridden: see below

# %%
'''
fwhm = 4
# Query sidereal objects (PS1)
cent_coord = yfu.center_radec(ccd_or_header=hdr, center_of_image=True)

rad = yfu.fov_radius(header=hdr, unit=u.deg)

# Initialize PanSTARRS1 class
ps1 = ypu.PanSTARRS1(
    ra=cent_coord.ra.value, 
    dec=cent_coord.dec.value,
    radius=rad,
    column_filters={"gmag":"10.0..18.0", "e_gmag":"<0.10"}
)

_ = ypu.organize_ps1_and_isnear(
    ps1=ps1,
    header=hdr,
    bezel=50,
    group_crit_separation=6*fwhm,
    select_filter_kw=dict(filter_names=['g', 'r', 'i'], n_mins=5), # 'g, 'r', 'i'
    drop_by_Kron=True,
    calc_JC=True  # This uses g-r color, so both g and r must not be removed.
)

ps1.queried

q_stars = ps1.queried.copy()
'''
q_stars = q.queried.copy()  # This will be overridden: see below

pos_sky, pos_pix

pos_stars = np.array([q_stars["x"], q_stars["y"]]).T

#%%
ap_stars_orig = CAp(positions=pos_stars_orig, r=15)
ap_stars = CAp(positions=pos_stars, r=20)

fig, axs = plt.subplots(1, 1, figsize=(20, 15), sharex=False, sharey=False, gridspec_kw=None)
_tool_visualization.norm_imshow(axs, hdul[0].data, zscale=True)
ap_stars_orig.plot(color='w', lw=2)
ap_stars.plot(color='r', lw=2)

axs.set_title("PS1 Query: Dropping Nearby Stars")
plt.tight_layout()

# %%
avg, med, std = sigma_clipped_stats(data) # default is 3-sigma, 5 iters
thresh = 5 * std
finder = DAOStarFinder(
    fwhm=4, threshold=thresh,   # In reality, FWHM must be measured a priori using, e.g., ``ginga``
    sharplo=0.2, sharphi=1.0,   # default values 0.2 and 1.0
    roundlo=-1.0, roundhi=1.0,  # default values -1 and +1
    sigma_radius=1.5,           # default values 1.5
    ratio=1.0,                  # 1.0: circular gaussian
    exclude_border=True         # To exclude sources near edges
)

# The DAOStarFinder object ``finder`` gets at least one input: the image.
# Then it returns the astropy table which contains the aperture photometry results:
found = finder(data)

# Use ``found`` for aperture photometry:
coords_SF = np.array([found['xcentroid'], found['ycentroid']]).T
ap_found = CAp(coords_SF, r=25)  

# Plot all
fig, axs = plt.subplots(1, 1, figsize=(20, 15), sharex=False, sharey=False, gridspec_kw=None)
_tool_visualization.norm_imshow(axs, data, zscale=True)
ap_found.plot(color='k', lw=2, alpha=0.7)
ap_stars.plot(color='red', lw=2, alpha=0.7)

axs.set
plt.tight_layout()


# %%
# Set the maximum allowable distance
match_distance = 20 # pixel  5

# Initialize some columns
for c in ["xcentroid", "ycentroid", "peak", "pixel_shift"]:
    q_stars[c] = np.nan

for i, coo in enumerate(q_stars):
    dx = np.abs(found['xcentroid'] - coo['x'])
    dy = np.abs(found['ycentroid'] - coo['y'])
    distances = np.sqrt(dx**2 + dy**2)
    accepted_SF = found[distances < match_distance]
    if len(accepted_SF) == 0:
        continue
    elif len(accepted_SF) > 1:
        raise ValueError(f"More than 1 star match for the {i}-th row!? Reduce match_distance.")
    else:
        q_stars[i]["pixel_shift"] = distances.min()
        q_stars[i]["xcentroid"] = found[distances.argmin()]["xcentroid"]
        q_stars[i]["ycentroid"] = found[distances.argmin()]["ycentroid"]
        q_stars[i]["peak"] = found[distances.argmin()]["peak"]

# Select only those which are matched
matched = q_stars[q_stars["pixel_shift"] < match_distance]
matched

# %%
ccd = yfu.load_ccd(fullname)
err = yfu.errormap(ccd, gain_epadu=hdr["GAIN"], rdnoise_electron=hdr["RDNOISE"])
# %%
pos_sky, pos_pix

found_targ = found[(498 < found["xcentroid"]) & (found["xcentroid"] < 502)
                   & (498 < found["ycentroid"]) & (found["ycentroid"] < 502)]
found_targ

# %%
pos_star_match = np.array((matched["xcentroid"], matched["ycentroid"])).T
pos_targ = np.array((found_targ["xcentroid"], found_targ["ycentroid"])).T

fwhm = 4
r_ap = 2*fwhm
r_in = 4*fwhm
r_out = 6*fwhm

star_ap = CAp(positions=pos_star_match, r=r_ap)
star_an = CAn(positions=pos_star_match, r_in=r_in, r_out=r_out)

targ_ap = CAp(positions=pos_targ, r=r_ap)
targ_an = CAn(positions=pos_targ, r_in=r_in, r_out=r_out)

phot_star = ypu.apphot_annulus(ccd=ccd,
                               aperture=star_ap,
                               annulus=star_an,
                               error=err)

phot_targ = ypu.apphot_annulus(ccd=ccd,
                               aperture=targ_ap,
                               annulus=targ_an,
                               error=err)
# %%
fig, axs = plt.subplots(1, 1, figsize=(8, 8), 
                        sharex=False, sharey=False, gridspec_kw=None)
_tool_visualization.norm_imshow(axs, ccd, zscale=True)
star_ap.plot(axs, color='w')
star_an.plot(axs, color='r')
targ_ap.plot(axs, color='r')
targ_an.plot(axs, color='y')

# inset axes....
axins = axs.inset_axes([0.6, 0.4, 0.3, 0.3])
_tool_visualization.norm_imshow(axins, ccd, zscale=True)
# sub region of the original image
x1, x2, y1, y2 = 450, 550, 450, 550
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels('')
axins.set_yticklabels('')
targ_ap.plot(axins, color='r')
targ_an.plot(axins, color='y')
axs.indicate_inset_zoom(axins, edgecolor='r')

plt.tight_layout()
fig.align_ylabels(axs)
fig.align_xlabels(axs)
plt.show()
# %%
# Photometry result of the target:
phot_targ
# %%
# photometry of stars
phot_star

# %%
def linf(x, a, b):
    return a*x + b

fig = plt.figure(figsize=(12, 6))
gs = gridspec.GridSpec(3, 3)
# NOTE: Z = R - R_inst = zp + 1st ext
ax_l = fig.add_subplot(gs[:2, 0])                # Linearity plot
ax_r = fig.add_subplot(gs[2, 0], sharex=ax_l)    # residual plot
ax_c = fig.add_subplot(gs[:, 1:])                # C(g-r) VS Z plot
errbfmt = dict(marker='x', ls='', capsize=5, elinewidth=0.5)

m_true = matched["Vmag"]
m_inst = phot_star["mag"]
e_m_true = matched["e_Vmag"]
e_m_inst = phot_star["merr"]
color = matched["grcolor"]
e_color = matched["e_grcolor"]
Z = m_true - m_inst
e_Z = np.sqrt(e_m_true**2 + e_m_inst**2)
s_Z = np.std(Z, ddof=1)  # standard deviation
m_Z = np.mean(Z)         # Simple mean
Z_fit = np.average(Z, weights=(1/e_Z**2))
e_Z_fit = np.sqrt(1 / np.sum(1/e_Z**2))

mm = np.linspace(m_true.min(), m_true.max(), 2)  # grid
cc = np.linspace(color.min(), color.max(), 2)  # grid
p_ll, _ = curve_fit(linf, m_true, m_inst, sigma=e_m_inst, absolute_sigma=True)
p_cz, _ = curve_fit(linf, color, Z, sigma=e_Z, absolute_sigma=True)

ax_l.errorbar(m_true, m_inst, xerr=e_m_true, yerr=e_m_inst, **errbfmt)
ax_l.plot(mm, linf(mm, *p_ll), 'k-', label=f"y={p_ll[0]:.4f}x + {-p_ll[1]:.3f}")
ax_l.plot(mm, mm - Z_fit, 'r-', label=f"y=x + {Z_fit:.3f} (see right)")
ax_l.set(ylabel="V$_\mathrm{inst}$", title="Linearity Curve")
ax_l.legend(loc=2, framealpha=0, fontsize=10)
for j, (x_i, y_i) in enumerate(zip(m_true, m_inst)):
    ax_l.text(x=x_i+0.1, y=y_i-0.1, s=j+1)
    
ax_r.errorbar(m_true, Z, yerr=e_Z, **errbfmt)
ax_r.axhline(Z_fit, color='r', ls='-', lw=1)
ax_r.axhline(Z_fit + e_Z_fit, color='r', ls='--', lw=1, label=None) 
ax_r.axhline(Z_fit - e_Z_fit, color='r', ls='--', lw=1, label=None)
ax_r.axhline(m_Z + s_Z, color='b', ls=':', lw=1, label=None) 
ax_r.axhline(m_Z - s_Z, color='b', ls=':', lw=1, label=None)
ax_r.set(xlabel="V from PS1 g/r (Tonry+ 2012)",
         ylabel="V - V$_\mathrm{inst}$",
         ylim=(m_Z - 5*s_Z, m_Z + 5*s_Z)
        )
for j, (x_i, y_i) in enumerate(zip(m_true, Z)):
    ax_r.text(x=x_i, y=y_i, s=j+1)

ax_c.errorbar(color, Z, yerr=e_Z, xerr=e_color, **errbfmt)
ax_c.axhline(Z_fit, color='r', ls='-', lw=1, label=f"wieghted avg = {Z_fit:.3f}")
ax_c.axhline(Z_fit + e_Z_fit, color='r', ls='--', lw=1, label=f"err = {e_Z_fit:.3f}") 
ax_c.axhline(Z_fit - e_Z_fit, color='r', ls='--', lw=1, label=None)
ax_c.axhline(m_Z, color='b', ls='-', lw=1, label=f"simple avg = {m_Z:.3f}")
ax_c.axhline(m_Z + s_Z, color='b', ls=':', lw=1, label=f"std = {s_Z:.3f}") 
ax_c.axhline(m_Z - s_Z, color='b', ls=':', lw=1, label=None)
ax_c.plot(cc, linf(cc, *p_cz), 'k-')
ax_c.legend(loc=2, framealpha=0, ncol=2)
ax_c.set(xlabel="g - r from PS1",
         title="Z = V - V$_\mathrm{inst}$ = (k + k''X)C + (zero + k'X)")

ax_c2 = ax_c.twinx()
ax_c2.plot(np.nan, np.nan, 'k-',
           label=f"y={p_cz[0]:.3f}x + {p_cz[1]:.3f}")
ax_c2.legend(loc=4)
ax_c2.axis('off')

for j, (x_i, y_i) in enumerate(zip(color, Z)):
    ax_c.text(x=x_i, y=y_i, s=j+1)

_tool_visualization.linticker(
    [ax_l, ax_r, ax_c],
    xmajlockws=[1  , 1  , 0.2],
    xminlockws=[0.2, 0.2, 0.1],
    ymajlockws=[1  , 0.2, 0.1],
    yminlockws=[0.2, 0.1, 0.05]
)

plt.tight_layout()

# fig, zp_fit, dzp_fit, mzp, szp, linSlope, zpSlope = zpplot(matched, phot_star)

#%%
phot_targ["mstd"] = phot_targ["mag"] + Z_fit
phot_targ["mstd_err"] = np.sqrt(phot_targ["merr"]**2 + e_Z_fit**2)

tab = phot_targ.to_pandas()
tab.round(3)