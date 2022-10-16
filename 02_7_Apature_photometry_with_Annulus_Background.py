# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 23:15:44 2018

@author: user
"""
#%%
from cmath import log
from glob import glob
import numpy as np
import pandas as pd
import os
import Python_utilities
import matplotlib.pyplot as plt
from photutils import detect_threshold
from photutils import DAOStarFinder
from photutils import CircularAnnulus as CircAn
from photutils import CircularAperture as CircAp
from photutils import aperture_photometry as APPHOT
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from astropy.wcs import WCS

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))

base_dr = "../Post_processing/M35_Light_-_2018-10-31_-_TMB130ss_STF-8300M_-_1bin/"

### make all fits file list...
fullnames = Python_utilities.getFullnameListOfallFiles("{}/input".format(base_dr))
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))

c_method = 'median'
master_dr = "master_files/"
reduced_dr = "readuced_files/"
result_dr = "result_files/"

if not os.path.exists('{0}'.format("{}{}".format(base_dr, result_dr))):
    os.makedirs("{}{}".format(base_dr, result_dr))
    print("{}{}is created".format(base_dr, result_dr))

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

    thresh = detect_threshold(data=img, nsigma=3)
    thresh = thresh[0][0]
    print('detect_threshold', thresh)

    from photutils import DAOStarFinder
    FWHM   = 3.5
    DAOfind = DAOStarFinder(threshold=thresh, 
                        fwhm=FWHM 
                        #sharplo=0.5, sharphi=2.0,  # default values: sharplo=0.2, sharphi=1.0,
                        #roundlo=0.0, roundhi=0.2,  # default values: roundlo=-1.0, roundhi=1.0,
                        #sigma_radius=1.5,          # default values: sigma_radius=1.5,
                        #ratio=0.5,                 # 1.0: circular gaussian:  ratio=1.0,
                        #sky=None, exclude_border=False)       # To exclude sources near edges : exclude_border=True
                        )
    # The DAOStarFinder object ("DAOfind") gets at least one input: the image.
    # Then it returns the astropy table which contains the aperture photometry results:
    DAOfound = DAOfind(img)
    print('{} star(s) founded by DAOStarFinder...'.format(len(DAOfound)))
    
    if len(DAOfound)==0 :
        print ('No star was founded by DAOStarFinder...\n'*3)
    else : 
        # Use the object "found" for aperture photometry:
        N_stars = len(DAOfound)
        print('{} star(s) founded by DAOStarFinder...'.format(N_stars))
        DAOfound.pprint(max_width=1800)
        
        # save XY coordinates:
        DAOfound.write("{}{}{}_DAOStarfinder_fwhm{}.csv".\
                        format(base_dr, result_dr, fullname_el[-1][:-4], FWHM), 
                        overwrite = True,
                        format='ascii.fast_csv')
        
        print('type(DAOfound): {}'.format(type(DAOfound)))
        print('DAOfound: {}'.format(DAOfound))

        #DAOcoord = (DAOfound['xcentroid'], DAOfound['ycentroid']) 
        DAOcoord = (DAOfound['xcentroid'], DAOfound['ycentroid']) 
        print('type(DAOcoord): {}'.format(type(DAOcoord)))
        print('DAOcoord: {}'.format(DAOcoord))

        DAOcoord_zip = list(zip(DAOfound['xcentroid'], DAOfound['ycentroid']))
        print('type(DAOcoord_zip): {}'.format(type(DAOcoord_zip)))
        print('DAOcoord_zip: {}'.format(DAOcoord_zip))
        print('DAOcoord_zip[0]: {}'.format(DAOcoord_zip[0]))
        print('DAOcoord_zip[1]: {}'.format(DAOcoord_zip[1]))
        
        # Save apertures as circular, 4 pixel radius, at each (X, Y)
        DAOapert = CircAp((DAOcoord_zip), r=4.)  
        print('type(DAOapert): {}'.format(type(DAOapert)))
        print('DAOapert: {}'.format(DAOapert))
        print('dir(DAOapert): {}'.format(dir(DAOapert)))
        
        DAOimgXY = np.array(DAOcoord)
        print('type(DAOimgXY): {}'.format(type(DAOimgXY)))
        print('DAOimgXY: {}'.format(DAOimgXY))

        DAOannul = CircAn(positions = (DAOcoord_zip), r_in = 4*FWHM, r_out = 6*FWHM) 
        print('type(DAOannul): {}'.format(type(DAOannul)))
        print('DAOannul: {}'.format(DAOannul))

       
        
        plt.figure(figsize=(24,16))
        ax = plt.gca()
        im = plt.imshow(img, vmax=thresh*4, origin='lower')
        DAOannul.plot(color='red', lw=2., alpha=0.4)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        plt.colorbar(im, cax=cax)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        
        plt.colorbar(im, cax=cax)

        ###########################################################
        # input some text for explaination. 
        plt.title("Result of DAOStarfinder", fontsize = 28, 
            ha='center')
        plt.annotate('filename: {}'.format(fullname_el[-1]), fontsize=10,
            xy=(1, 0), xytext=(-500, -40), va='top', ha='left',
            xycoords='axes fraction', textcoords='offset points')
                    
        plt.annotate('FWHM: {}'.format(FWHM), fontsize=10,
            xy=(1, 0), xytext=(-1200, -30), va='top', ha='left',
            xycoords='axes fraction', textcoords='offset points')
            
        plt.annotate('Sky threshold: {:02f}'.format(thresh), fontsize=10,
            xy=(1, 0), xytext=(-1200, -40), va='top', ha='left',
            xycoords='axes fraction', textcoords='offset points')

        plt.annotate('Number of star(s): {}'.format(len(DAOfound)), fontsize=10,
            xy=(1, 0), xytext=(-1200, -50), va='top', ha='left',
            xycoords='axes fraction', textcoords='offset points')
        
        plt.savefig("{}{}{}_DAOStarfinder_fwhm{}.png".\
                        format(base_dr, result_dr, fullname_el[-1][:-4], FWHM))
        #plt.show()
        plt.close()
        
        
        ####################################
        if not 'RDNOISE' in hdul[0].header : 
            ronoise = 10.0 #STF-8300M
        else :
            ronoise = hdul[0].header['RDNOISE']
        
        if 'EGAIN' in hdul[0].header : 
            gain = hdul[0].header['EGAIN']
        elif 'GAIN' in hdul[0].header : 
            gain = hdul[0].header['GAIN']
        else :
            gain = 0.5
        
        ####################################
        N_star = len(DAOfound)
        
        mag_ann  = np.zeros(N_star)
        merr_ann = np.zeros(N_star)
        
        # aperture sum
        apert_sum = APPHOT(img, DAOapert, method='exact')['aperture_sum']
        ap_area   = DAOapert.area()
        print(apert_sum)
        
        apert_result = 'ID, Msky, sky_std, Sky count Pixel_N, Sky reject Pixel_N, mag_ann, merr_ann\n'
        
        for star_ID in range(0, N_stars):

            mask_annul = (DAOannul.to_mask(method='center'))[star_ID]
            mask_apert = (DAOapert.to_mask(method='center'))[star_ID]
            # CAUTION!! YOU MUST USE 'center', NOT 'exact'!!!
    
            cutimg = mask_annul.cutout(img)
            df_cutimg = pd.DataFrame(cutimg)
            df_cutimg.to_csv("{}{}{}_DAOstarfinder__starID_{:04}_Star_Area_pixel_value.csv".\
                        format(base_dr, result_dr, fullname_el[-1][:-4], star_ID)) 
            break
            
            apert_apply = mask_apert.multiply(img)  # change from 'sky_apply  = mask_annul.apply(img)'
            df_apert_apply = pd.DataFrame(apert_apply)
            df_apert_apply.to_csv("{}{}{}_DAOstarfinder__starID_{:04}_Apeture_area_pixel_value.csv".\
                        format(base_dr, result_dr, fullname_el[-1][:-4], star_ID)) 
            
            apert_non0   = np.nonzero(apert_apply)
            apert_pixel  = apert_apply[apert_non0]
            apert_nan = apert_apply.copy()
            apert_nan[apert_nan == 0] = np.nan
            #print(apert_nan)
            
            sky_apply = mask_annul.multiply(img)  # change from 'sky_apply  = mask_annul.apply(img)'
            df_sky_apply = pd.DataFrame(sky_apply)
            df_sky_apply.to_csv("{}{}{}_DAOstarfinder_starID_{:04}_Sky_Annulus_pixel_value.csv".\
                        format(base_dr, result_dr, fullname_el[-1][:-4], star_ID)) 
            
            sky_non0   = np.nonzero(sky_apply)
            sky_pixel  = sky_apply[sky_non0]
            sky_nan = sky_apply.copy()
            sky_nan[sky_nan == 0] = np.nan
            
            msky, sky_std, nsky, nrej = sky_fit(sky_pixel, method='mode', mode_option='sex')
            
            # sky estimation
            flux_star = apert_sum[star_ID] - msky * ap_area  # total - sky
            flux_err  = np.sqrt(apert_sum[star_ID] * gain    # Poissonian (star + sky)
                                + ap_area * ronoise**2 # Gaussian
                                + (ap_area * (gain * sky_std))**2 / nsky ) 
            mag_ann[star_ID], merr_ann[star_ID] = mag_inst(flux_star, flux_err)
            
            print('{0:7d}: {1:.5f} {2:.5f} {3:4d} {4:3d} {5:.3f} {6:.3f}'.format(i, msky, sky_std, nsky, nrej, mag_ann[i], merr_ann[i]))
            apert_result += '{0:04}, {1:.5f}, {2:.5f}, {3:4d}, {4:3d}, {5:.3f}, {6:.3f}\n'\
            .format(star_ID, msky, sky_std, nsky, nrej, mag_ann[star_ID], merr_ann[star_ID])
            
            fig = plt.figure(figsize=(12,12))
            fig.add_subplot(2,3,1)
            plt.imshow(cutimg, vmin=np.mean(cutimg/5.0), vmax=np.mean(cutimg*2.0), origin='lower')
            plt.ylabel('pixels')
            plt.grid(ls=':')
            plt.xlabel('DAO Star area image\n sum: {0} \n mean: {1:.3f}\n std: {2:.3f} \n max: {3} \n min: {4}'\
                       .format(np.sum(cutimg), np.mean(cutimg), np.std(cutimg), np.max(cutimg), np.min(cutimg)))
            
            fig.add_subplot(2,3,2)
            plt.imshow(apert_apply, vmin=np.mean(cutimg/5.0), vmax=np.mean(cutimg*2.0), origin='lower')
            plt.ylabel('pixels')
            plt.grid(ls=':')
            plt.xlabel('DAO aperture area \n (r=2.0*FWHM)\n sum: {0:.0f}  /  {1:.3f} \n max: {2:.0f}\n min: {3:.0f} \n mean: {4:.3f} \n std: {5:.3f} \n Number of Pixel: {6}  /  {7:.5f}'\
                       .format(np.sum(apert_pixel), apert_sum[star_ID], np.max(apert_pixel), np.min(apert_pixel), np.mean(apert_pixel), np.std(apert_pixel), len(apert_pixel), ap_area))

            fig.add_subplot(2,3,3)
            plt.imshow(sky_apply, vmin=np.mean(cutimg/5.0), vmax=np.mean(cutimg*2.0), origin='lower')
            plt.ylabel('pixels')
            plt.grid(ls=':')
            plt.xlabel('DAO annulus area \n (r_in=4*FWHM, r_out=6*FWHM)\n sum: {0:.0f} \n max: {1:.0f} \n min: {2:.0f} \n mean: {3:.3f}  /  {4:.3f}\n std: {5:.3f}  /  {6:.3f} \n Number of sky pixels: {7} \n Number of reject pixels: {8}'\
                       .format(np.sum(sky_pixel), np.max(sky_pixel), np.min(sky_pixel), np.mean(sky_pixel), msky, np.std(sky_pixel), sky_std, nsky, nrej))
            
            #fig.add_subplot(2,3,4)
            #plt.imshow(cutimg, vmin=np.mean(cutimg/5.0), vmax=np.mean(cutimg*2.0), origin='lower')
            #plt.ylabel('pixels')
            #plt.grid(ls=':')
            #plt.xlabel('DAO Star area image\n sum: {0} \n mean: {1:.3f}\n std: {2:.3f} \n max: {3} \n min: {4}'\
                        #.format(np.sum(cutimg), np.mean(cutimg), np.std(cutimg), np.max(cutimg), np.min(cutimg)))
            
            fig.add_subplot(2,3,5)
            plt.imshow(apert_nan, vmin=np.mean(cutimg/5.0), vmax=np.mean(cutimg*2.0), origin='lower')
            plt.ylabel('pixels')
            plt.grid(ls=':')
            plt.xlabel('DAO aperture area \n (r=2.0*FWHM)\n sum: {0:.0f}  /  {1:.3f} \n max: {2:.0f}\n min: {3:.0f} \n mean: {4:.3f} \n std: {5:.3f} \n Total pixel number: {6} \n Pixel number of aperture: {7}  /  {8:.5f}'\
                       .format(np.nansum(apert_nan), apert_sum[star_ID], np.nanmax(apert_nan), np.nanmin(apert_nan), np.nanmean(apert_nan), np.nanstd(apert_nan), np.shape(apert_nan)[0]*np.shape(apert_nan)[1], np.count_nonzero(~np.isnan(apert_nan)), ap_area))

            fig.add_subplot(2,3,6)
            plt.imshow(sky_nan, vmin=np.mean(cutimg/5.0), vmax=np.mean(cutimg*2.0), origin='lower')
            plt.ylabel('pixels')
            plt.grid(ls=':')
            plt.xlabel('DAO annulus area \n (r_in=4*FWHM, r_out=6*FWHM)\n sum: {0:.0f} \n max: {1:.0f} \n min: {2:.0f} \n mean: {3:.3f}  /  {4:.3f}\n std: {5:.3f}  /  {6:.3f} \n Pixel number of annulus: {7} \n Pixel number of accept: {8} \n Pixel number of reject: {9}'\
                       .format(np.nansum(sky_nan), np.nanmax(sky_nan), np.nanmin(sky_nan), np.nanmean(sky_nan), msky, np.nanstd(sky_nan), sky_std, np.count_nonzero(~np.isnan(sky_nan)), nsky, nrej))
            
            
            plt.gcf().suptitle('DAOstarfinder Annulus result: {0!s} \n \
                   FWHM: {1:.2f}, Sky threshold: {2:.4f}\n\
                   star ID: {3:04},  read out noise: {4:.2f}, gain: {5:.3f}\n \
                   flux_star: {6:.3f}, flux_err: {7:.3f} \n \
                   instrumental magnitude: {8:.3f}, mag_err: {9:.3f}'\
                   .format(fullname_el[-1], FWHM, thresh, \
                           star_ID, ronoise, gain, \
                           flux_star, flux_err, \
                           mag_ann[star_ID], merr_ann[star_ID]), fontsize=12)
            #plt.colorbar(size="5%", pad=0.05)                        
            
            plt.savefig("{}{}{}_DAOstarfinder__starID_{:04}_Annulus_result.png".\
                            format(base_dr, result_dr, fullname_el[-1][:-4], star_ID),
                            overwrite=True)

            print('{0!s}_DAOstarfinder_starID_{1:04}_Annulus_result.png is saved'\
                  .format(fullname_el[-1], star_ID))
            #plt.show()
        print(apert_result)
        #with open(f_name[:-4]+'_DAOstarfinder_all_AP_Annulus_result.csv', 'w') as f:
        #    f.write(apert_result)
        
        
        
        
        break