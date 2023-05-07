# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 23:15:44 2018

@author: user
"""
#%%
from glob import glob
import numpy as np
import os
import astropy.units as u
from ccdproc import CCDData, ccd_process
import _Python_utilities

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))

BASEDIR = "../Post_processing/M35_Light_-_2018-10-31_-_TMB130ss_STF-8300M_-_1bin/"

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
###
fullname = fullnames_light[0]
fullname_el = fullname.split("/")

img = fits.getdata(fullname)
print("img.shape: {}".format(img.shape))
print("img: {}".format(img))

#hdul = fits.open(fullname)
#img1 = hdul[0].data
#print(img == img1)


#%%
from astropy.stats import mad_std, sigma_clipped_stats

sky_th = 5    # sky_th * sky_sigma will be used for detection lower limit : sky_th = 5
sky_s  = mad_std(img)
thresh = sky_s * sky_th
print('sky_s x sky_th = threshold')
print('{0:8.6f} x {1:4d}   =   {2:8.6f}\n'.format(sky_s, sky_th, thresh))

# What if we do "sigma-clip" than MAD?
sky_mean, sky_median, sky_std = sigma_clipped_stats(img,
                        sigma=3.0) # default is 3-sigma, 5 iters
print(sky_mean, sky_median, sky_std)

thresh_sc = sky_std * sky_th
print('3 sigma 5 iters clipped case:')
print('{0:8.6f} x {1:4d}   =   {2:8.6f}\n'.format(sky_std, sky_th, thresh_sc))

#%%
from photutils import detect_threshold

#thresh = detect_threshold(data=img.data, snr=3)
thresh = detect_threshold(data=img.data, nsigma=3)
thresh = thresh[0][0]

print('detect_threshold', thresh)

#%%
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

minorLocator = AutoMinorLocator()

def radial_profile(data, center):
    x, y = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 

hdul = fits.open('{}'.format(fullname))
img = hdul[0].data
img[np.isnan(img)] = 0

# Calculate the indices from the image
y, x = np.indices(img.shape)
center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])

#center = np.unravel_index(img.argmax(), img.shape)
rad_profile = radial_profile(img, center)

fig, ax = plt.subplots()
plt.plot(rad_profile[:], 'x-')

ax.xaxis.set_minor_locator(minorLocator)

plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4, color='r')
plt.grid()
ax.set_ylabel("ADU")
ax.set_xlabel("Pixels")
plt.grid(which="minor")
plt.show()



#%%
import numpy as np
#import pyfits
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def azimuthalAverage(image, center=None):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).

    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)


    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[1]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin

    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr
    print("center: {}".format(center))
    print("i_sorted".format(i_sorted))
    print("radial_prof: {}".format(radial_prof))
    return radial_prof

#read in image
#hdulist = pyfits.open('cit6ndf2fitsexample.fits')
#scidata = np.array(hdulist[0].data)[0,:,:]
center = None
radi = 5
rad = azimuthalAverage(img, center)
#rad = azimuthalAverage(img, [horizontal, vertical])

plt.xlabel('radius(pixels?)', fontsize=12)
plt.ylabel('image intensity', fontsize=12)
plt.xlim(0,10)
plt.ylim(0, 3.2)
plt.plot(rad[radi:])
plt.savefig('testfig1.png')
plt.show()

#%%

import numpy as np
from scipy.interpolate import UnivariateSpline

def profiles(image):
    amp = np.max(image)
    ypix, xpix = np.where(image==amp)   
    x = np.take(image, ypix[0], axis=0)
    y = np.take(image, xpix[0], axis=1)

    return x, y #these are the horizontal and vertical profiles through the star's centroid

def interpolate_width(axis):
    half_max = 1/2
    x = np.linspace(0, len(axis), len(axis))

    # Do the interpolation
    spline = UnivariateSpline(x, axis-half_max, s=0)
    r1, r2 = spline.roots()

    return r2-r1 #this is the FWHM along the specified axis

horizontal, vertical = profiles(img)
print("horizontal: {}".format(horizontal))
print("vertocal: {}".format(vertical))

fwhm_x = interpolate_width(horizontal)
fwhm_y = interpolate_width(vertical)

print("fwhm_x: {}".format(fwhm_x))
print("fwhm_y: {}".format(fwhm_y))
#%%
from photutils import StarFinder
from astropy.convolution import Gaussian2DKernel

gaussian_2D_kernel = Gaussian2DKernel(10)

STARfind = StarFinder(threshold = thresh, kernel = gaussian_2D_kernel)

STARfound = STARfind(img)
print('{} stars founded by STARFinder...'.format(len(STARfound)))

#%%
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder

FWHM   = 2.5

DAOfind = DAOStarFinder(threshold = thresh, 
                        fwhm = FWHM 
                        #sharplo=0.5, sharphi=2.0,  # default values: sharplo=0.2, sharphi=1.0,
                        #roundlo=0.0, roundhi=0.2,  # default values: roundlo=-1.0, roundhi=1.0,
                        #sigma_radius=1.5,          # default values: sigma_radius=1.5,
                        #ratio=0.5,                 # 1.0: circular gaussian:  ratio=1.0,
                        #sky=None, exclude_border=False)       # To exclude sources near edges : exclude_border=True
                        )
# The DAOStarFinder object ("find") gets at least one input: the image.
# Then it returns the astropy table which contains the aperture photometry results:
DAOfound = DAOfind(img)
print('{} stars founded by DAOStarFinder...'.format(len(DAOfound)))

#%%
if len(DAOfound)==0 :
    print ('No star founded using DAOStarFinder\n'*3)
else : 
    # Use the object "found" for aperture photometry:
    # save XY coordinates:
    DAOfound.pprint(max_width=1800)
    
    DAOfound.write("{}{}{}_DAOStarfinder.csv".\
                        format(BASEDIR, result_dir, fullname_el[-1][:-4]), 
                        overwrite = True,
                        format='ascii.fast_csv')
    
    #DAOcoord = (DAOfound['xcentroid'], DAOfound['ycentroid']) 
    DAOcoord = zip(DAOfound['xcentroid'], DAOfound['ycentroid'])   

    # Save apertures as circular, 4 pixel radius, at each (X, Y)
    DAOapert = CircAp(DAOcoord, r=4.)  
    print('DAOapert\n ', DAOapert)
    
    DAOimgXY = np.array(DAOcoord)
    print('DAOimgXY \n', DAOimgXY)
    
    plt.figure(figsize=(12,12))
    ax = plt.gca()
    im = plt.imshow(img, vmax=65536, origin='lower')
    DAOapert.plot(color='red', lw=2., alpha=0.7)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    plt.colorbar(im, cax=cax)
    
    plt.show()
    
#%%
import numpy as np
import matplotlib.pyplot as plt
from photutils import IRAFStarFinder
from photutils import CircularAperture as CircAp

IRAFfind = IRAFStarFinder(threshold = thresh,
                            fwhm = FWHM
                            #sigma_radius=1.5, minsep_fwhm=1.0,  # default values: sigma_radius=1.5, minsep_fwhm=2.5,
                            #sharplo=0.5, sharphi=2.0,   # default values: sharplo=0.5, sharphi=2.0,
                            #roundlo=0.0, roundhi=0.2,   # default values: roundlo=0.0, roundhi=0.2,
                            #sky=None, exclude_border=False)  # default values: sky=None, exclude_border=False)
                            )

# The DAOStarFinder object ("find") gets at least one input: the image.
# Then it returns the astropy table which contains the aperture photometry results:
IRAFfound = IRAFfind(img)
print('{} stars founded by IRAFStarFinder...'.format(len(IRAFfound)))

#%%
if len(IRAFfound)==0 :
    print ('No star founded using IRAFStarFinder')
else : 
    IRAFfound.pprint(max_width=1800)
    IRAFfound.write("{}{}{}_IRAFStarfinder.csv".\
                        format(BASEDIR, result_dir, fullname_el[-1][:-4]), 
                        overwrite =True,
                        format='ascii.fast_csv')
    
    # Use the object "found" for aperture photometry:
    # save XY coordinates:
    print (len(IRAFfound), 'stars founded')
    
    #IRAFcoord = (IRAFfound['xcentroid'], IRAFfound['ycentroid']) 
    IRAFcoord = zip(IRAFfound['xcentroid'], IRAFfound['ycentroid'])   
    
    # Save apertures as circular, 4 pixel radius, at each (X, Y)
    IRAFapert = CircAp(IRAFcoord, r=4.)  
    #print('IRAFapert\n ', IRAFapert)
    
    IRAFimgXY = np.array(IRAFcoord)
    #print('IRAFimgXY \n', IRAFimgXY)
            
    plt.figure(figsize=(12,12))
    ax = plt.gca()
    im = plt.imshow(img, vmax=65536, origin='lower')
    IRAFapert.plot(color='red', lw=2., alpha=0.7)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    plt.colorbar(im, cax=cax)
    
    plt.show()