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
cd ~/Downloads/ysphttutilpy && git pull && pip install -e . 
cd ~/Downloads/SNUO1Mpy && git pull && pip install -e . 

"""
#%%
from glob import glob
from pathlib import Path
import numpy as np
import os
import astropy.units as u
from ccdproc import CCDData, ccd_process
import Python_utilities

from astropy.time import Time
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from photutils.aperture import CircularAnnulus, CircularAperture
import pandas as pd
import matplotlib.pyplot as plt

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

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
base_dir = "../RnE_2022/"
#base_dir = Path("..\RnE_2022\KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin\reduced\solved\")
base_dir = Path("../RnE_2022/KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin/reduced/solved/")
base_dir = "../RnE_2022/KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin/reduced/solved/"

c_method = 'median'
master_dir = "master_files_ys"
reduced_dir = "reduced"
solved_dir = "solved"

#%%
base_dirs = sorted(Python_utilities.getFullnameListOfsubDir(base_dir))
base_dirs = [w for w in base_dirs \
        if not (w.endswith("{}/".format(master_dir)) \
            or w.endswith("{}/".format(reduced_dir)) \
            or w.endswith("{}/".format(solved_dir)) \
                or w.endswith(".fits"))]

print ("base_dirs: {}".format(base_dirs))

#%%
# for base_dir in base_dirs :
#     print ("Starting...\n{}".format(base_dir))

#     base_dir = Path(base_dir)

#     if not (base_dir/reduced_dir).exists():
#         os.makedirs(str((base_dir/reduced_dir)))

#%%
#####################################################################
# Our object (will be queried to JPL HORIZONS)
OBJID = '216' # Kleopatra

# Observed location
LOCATION = dict(lon=127.0, lat=37.3, elevation=130)

# It is used as a rough estimate, so no need to be accurate:
PIX2ARCSEC = 1.24*u.arcsec

# Used for any `astropy.SkyCoord` object:
SKYC_KW = dict(unit=u.deg, frame='icrs')

# Initial guess of FWHM in pixel
FWHM_INIT = 6

# Photometry parameters
R_AP = 1.5*FWHM_INIT # Aperture radius
R_IN = 4*FWHM_INIT   # Inner radius of annulus
R_OUT = 6*FWHM_INIT  # Outer radius of annulus
#####################################################################

#%%
base_dir = "../RnE_2022/KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin"
print ("Starting...\n{}".format(base_dir))
base_dir = Path(base_dir)

summary = yfu.make_summary(base_dir/reduced_dir/solved_dir/"*.fits")
#print(summary)
print(summary["file"][0])

ccd = yfu.load_ccd(summary["file"][0], 
                   unit="adu")

#%%
_, eph, _ = ypu.horizons_query(OBJID, epochs=Time(ccd.header["DATE-OBS"]).jd, location=LOCATION)
eph

#%%
pos_targ_init = SkyCoord(eph["RA"], eph["DEC"], **SKYC_KW).to_pixel(ccd.wcs)
ap = CircularAperture([pos_targ_init[0][0], pos_targ_init[1][0]], r=R_AP)
an = CircularAnnulus([pos_targ_init[0][0], pos_targ_init[1][0]], r_in=R_IN, r_out=R_OUT)

pos_targ_init


#%%
phot_targ = ypu.apphot_annulus(ccd, ap, an, error=yfu.errormap(ccd))
phot_targ

#%%
r_fov = yfu.fov_radius(ccd.header+ccd.wcs.to_header())
print(r_fov)
ps1 = ypu.PanSTARRS1(ccd.wcs.wcs.crval[0]*u.deg, ccd.wcs.wcs.crval[1]*u.deg, radius=r_fov,
                     column_filters={"rmag":"10.0..14.5", "e_rmag":"<0.10", "nr":">5"})
isnear = ypu.organize_ps1_and_isnear(
    ps1, 
    header=ccd.header+ccd.wcs.to_header(), 
    bezel=5*FWHM_INIT*PIX2ARCSEC.value,
    nearby_obj_minsep=5*FWHM_INIT*PIX2ARCSEC.value,
    group_crit_separation=6*FWHM_INIT
)
df_stars = ps1.queried.to_pandas()

#%%
df_stars

#%%
fig, axs = plt.subplots(1, 1, figsize=(8, 5), sharex=False, sharey=False, gridspec_kw=None)

yvu.norm_imshow(axs, ccd, zscale=True)
ap = CircularAperture([pos_targ_init[0][0], pos_targ_init[1][0]], r=R_AP)
an = CircularAnnulus([pos_targ_init[0][0], pos_targ_init[1][0]], r_in=R_IN, r_out=R_OUT)
ap.plot(axs, color="r")
an.plot(axs, color="b")

_phot_stars = []

for i, row in df_stars.iterrows():
    pos_star = SkyCoord(row["RAJ2000"], row["DEJ2000"], **SKYC_KW).to_pixel(ccd.wcs)
    ap = CircularAperture([pos_star[0], pos_star[1]], r=R_AP)
    an = CircularAnnulus([pos_star[0], pos_star[1]], r_in=R_IN, r_out=R_OUT)
    _phot_star = ypu.apphot_annulus(ccd, ap, an, error=yfu.errormap(ccd))
    _phot_star["Rmag"] = row["Rmag"]
    _phot_star["e_Rmag"] = row["e_Rmag"]
    _phot_star["grcolor"] = row["grcolor"]
    _phot_star["e_grcolor"] = row["e_grcolor"]
    _phot_star["id"] = i
    _phot_star["objID"] = int(row["objID"])
    _phot_stars.append(_phot_star)
    axs.text(pos_star[0]+10, pos_star[1]+10, f"star {i}", fontsize=8)
    ap.plot(axs, color="orange")
    an.plot(axs, color="w")


plt.tight_layout()
plt.show();

#%%
phot_stars = pd.concat(_phot_stars)
# phot_stars = phot_stars.loc[phot_stars["objID"] != 110823405221754720].copy()  # star 15
# SEE THE LAST CELL IN THIS FILE FOR DESCRIPTION
phot_stars

# %%
fig, axs = plt.subplots(1, 1, figsize=(5, 5), sharex=False, sharey=False, gridspec_kw=None)

_xx = np.linspace(13, 15)
axs.plot(phot_stars["Rmag"], phot_stars["mag"], '+')
axs.axhline(phot_targ["mag"].values, label="Kleopatra, instrumental mag")
axs.plot(_xx, _xx + np.median(phot_stars["mag"] - phot_stars["Rmag"]))

for _, row in phot_stars.iterrows():
    axs.text(row["Rmag"], row["mag"], int(row["id"]), fontsize=8)

axs.set(
    xlabel="R magnitude (PS1 to R_C filter by Tonry+2012)",
    ylabel="R_inst"
)
axs.legend()

plt.tight_layout()
plt.show();
   
# %%
fig, axs = plt.subplots(1, 2, figsize=(9, 5), sharex=False, sharey=False, gridspec_kw=None)

axs[0].plot(phot_stars["Rmag"], phot_stars["mag"] - phot_stars["Rmag"], '+')
axs[1].plot(phot_stars["grcolor"], phot_stars["mag"] - phot_stars["Rmag"], '+')
for _, row in phot_stars.iterrows():
    axs[0].text(row["Rmag"], row["mag"] - row["Rmag"], int(row["id"]), fontsize=8)
    axs[1].text(row["grcolor"], row["mag"] - row["Rmag"], int(row["id"]), fontsize=8)
    
axs[0].set(
    xlabel="R magnitude (PS1 to R_C filter by Tonry+2012)",
    ylabel="R_inst - R"
)
axs[1].set(
    xlabel="g - r (PS1)",
    ylabel="R_inst - R"
)

plt.tight_layout()
plt.show();

#%%
print(f'{int(df_stars.iloc[15]["f_objID"]):031b}')

#%%