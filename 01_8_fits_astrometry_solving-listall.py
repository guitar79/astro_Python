# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

No module named 'ccdproc'
conda install -c conda-forge ccdproc
"""
import os
import subprocess
import Python_utilities
import astro_utilities

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))

base_dir = "../CCD_obs_raw/"
save_dir = "../astrometry_solved"
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
fullnames = Python_utilities.getFullnameListOfallFiles(base_dir)
print ("fullnames: {}".format(fullnames))

n = 0
for fullname in fullnames[:] :
#fullname = fullnames[5]
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames))*100, os.path.basename(__file__)))
    print ("Starting...   fullname: {}".format(fullname))

    if fullname[-4:].lower() == ".fit" \
        and fullname[-7:].lower() != "wcs.fit":
        fullname_el = fullname.split("/")
        filename_el = fullname_el[-1].split("_")

        if filename_el[1].lower() == "light" :

            with subprocess.Popen(['solve-field',
                    # '-O', #--overwrite: overwrite output files if they already exist
                    # '--scale-units', 'arcsecperpix', #pixel scale
                    # '--scale-low', '0.1', '--scale-high', '0.40', #pixel scale
                    '-g',  # --guess-scale: try to guess the image scale from the FITS headers
                    # '-p', # --no-plots: don't create any plots of the results
                    '-D', '{0}'.format(save_dir),
                    '{0}'.format(fullname)],
                    stdout=subprocess.PIPE) as proc:
                print(proc.stdout.read())

'''
SIMPLE  =                    T                                                  
BITPIX  =                   16 /8 unsigned int, 16 & 32 int, -32 & -64 real     
NAXIS   =                    2 /number of axes                                  
NAXIS1  =                 1676 /fastest changing axis                           
NAXIS2  =                 1266 /next to fastest changing axis                   
BSCALE  =   1.0000000000000000 /physical = BZERO + BSCALE*array_value           
BZERO   =   32768.000000000000 /physical = BZERO + BSCALE*array_value           
DATE-OBS= '2018-12-23T14:57:39' /YYYY-MM-DDThh:mm:ss observation start, UT      
EXPTIME =   60.000000000000000 /Exposure time in seconds                        
EXPOSURE=   60.000000000000000 /Exposure time in seconds                        
SET-TEMP=  -30.000000000000000 /CCD temperature setpoint in C                   
CCD-TEMP=  -29.724637204887038 /CCD temperature at start of exposure in C       
XPIXSZ  =   10.800000000000001 /Pixel Width in microns (after binning)          
YPIXSZ  =   10.800000000000001 /Pixel Height in microns (after binning)         
XBINNING=                    2 /Binning factor in width                         
YBINNING=                    2 /Binning factor in height                        
XORGSUBF=                    0 /Subframe X position in binned pixels            
YORGSUBF=                    0 /Subframe Y position in binned pixels            
READOUTM= 'Raw     ' /          Readout mode of image                           
FILTER  = 'R       ' /          Filter used when taking image                   
IMAGETYP= 'Light frame'        / Type of image                                  
TRAKTIME=   2.0000000000000000 /Exposure time used for autoguiding              
FOCALLEN=   910.00000000000000 /Focal length of telescope in mm                 
APTDIA  =   130.00000000000000 /Aperture diameter of telescope in mm            
APTAREA =   13273.229330778122 /Aperture area of telescope in mm^2              
EGAIN   =  0.37999999523162842 /Electronic gain in e-/ADU                       
SBSTDVER= 'SBFITSEXT Version 1.0' /Version of SBFITSEXT standard in effect      
SWCREATE= 'MaxIm DL Version 6.05 150830 1Y265' /Name of software                
SWSERIAL= '1Y265-XV72V-AKE3M-01NQJ-KJ2M7-2J' /Software serial number            
FOCUSPOS=                11186 /Focuser position in steps                       
FOCUSSSZ=  0.00000000000000000 /Focuser step size in microns                    
FOCUSTEM=  0.81250000000000000 /Focuser temperature in deg C                    
OBJECT  = '46P_Wirtanen'                                                        
OBJCTRA = '05 19 39' /          Nominal Right Ascension of center of image      
OBJCTDEC= '+45 03 34' /         Nominal Declination of center of image          
OBJCTALT= ' 81.9799' /          Nominal altitude of center of image             
OBJCTAZ = '341.9778' /          Nominal azimuth of center of image              
OBJCTHA = '  0.2336' /          Nominal hour angle of center of image           
PIERSIDE= 'WEST    ' /          Side of pier telescope is on                    
SITELAT = '37 30 00' /          Latitude of the imaging location                
SITELONG= '127 00 00' /         Longitude of the imaging location               
JD      =   2458476.1233680556 /Julian Date at start of exposure                
JD-HELIO=   2458476.1289197467 /Heliocentric Julian Date at exposure midpoint   
AIRMASS =   1.0097004620005927 /Relative optical path length through atmosphere 
TELESCOP= '        ' /          telescope used to acquire this image            
INSTRUME= 'SBIG STF-8300 CCD Camera'                                            
OBSERVER= 'Kiehyun Park'                                                        
NOTES   = 'created by guitar79@naver.com'                                       
FLIPSTAT= '        '                                                            
SWOWNER = 'Inhyun Cho' /        Licensed owner of software                      
END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
'''

'''
SIMPLE  =                    T                                                  
BITPIX  =                   16 /8 unsigned int, 16 & 32 int, -32 & -64 real     
NAXIS   =                    2 /number of axes                                  
NAXIS1  =                 1676 /fastest changing axis                           
NAXIS2  =                 1266 /next to fastest changing axis                   
BSCALE  =   1.0000000000000000 /physical = BZERO + BSCALE*array_value           
BZERO   =   32768.000000000000 /physical = BZERO + BSCALE*array_value           
DATE-OBS= '2018-12-23T14:57:39' /YYYY-MM-DDThh:mm:ss observation start, UT      
EXPTIME =   60.000000000000000 /Exposure time in seconds                        
EXPOSURE=   60.000000000000000 /Exposure time in seconds                        
SET-TEMP=  -30.000000000000000 /CCD temperature setpoint in C                   
CCD-TEMP=  -29.724637204887038 /CCD temperature at start of exposure in C       
XPIXSZ  =   10.800000000000001 /Pixel Width in microns (after binning)          
YPIXSZ  =   10.800000000000001 /Pixel Height in microns (after binning)         
XBINNING=                    2 /Binning factor in width                         
YBINNING=                    2 /Binning factor in height                        
XORGSUBF=                    0 /Subframe X position in binned pixels            
YORGSUBF=                    0 /Subframe Y position in binned pixels            
READOUTM= 'Raw     ' /          Readout mode of image                           
FILTER  = 'R       ' /          Filter used when taking image                   
IMAGETYP= 'Light frame'        / Type of image                                  
TRAKTIME=   2.0000000000000000 /Exposure time used for autoguiding              
FOCALLEN=   910.00000000000000 /Focal length of telescope in mm                 
APTDIA  =   130.00000000000000 /Aperture diameter of telescope in mm            
APTAREA =   13273.229330778122 /Aperture area of telescope in mm^2              
EGAIN   =  0.37999999523162842 /Electronic gain in e-/ADU                       
SBSTDVER= 'SBFITSEXT Version 1.0' /Version of SBFITSEXT standard in effect      
SWCREATE= 'MaxIm DL Version 6.05 150830 1Y265' /Name of software                
SWSERIAL= '1Y265-XV72V-AKE3M-01NQJ-KJ2M7-2J' /Software serial number            
FOCUSPOS=                11186 /Focuser position in steps                       
FOCUSSSZ=  0.00000000000000000 /Focuser step size in microns                    
FOCUSTEM=  0.81250000000000000 /Focuser temperature in deg C                    
OBJECT  = '46P_Wirtanen'                                                        
OBJCTRA = '05 19 39' /          Nominal Right Ascension of center of image      
OBJCTDEC= '+45 03 34' /         Nominal Declination of center of image          
OBJCTALT= ' 81.9799' /          Nominal altitude of center of image             
OBJCTAZ = '341.9778' /          Nominal azimuth of center of image              
OBJCTHA = '  0.2336' /          Nominal hour angle of center of image           
PIERSIDE= 'WEST    ' /          Side of pier telescope is on                    
SITELAT = '37 30 00' /          Latitude of the imaging location                
SITELONG= '127 00 00' /         Longitude of the imaging location               
JD      =   2458476.1233680556 /Julian Date at start of exposure                
JD-HELIO=   2458476.1289197467 /Heliocentric Julian Date at exposure midpoint   
AIRMASS =   1.0097004620005927 /Relative optical path length through atmosphere 
TELESCOP= '        ' /          telescope used to acquire this image            
INSTRUME= 'SBIG STF-8300 CCD Camera'                                            
OBSERVER= 'Kiehyun Park'                                                        
NOTES   = 'created by guitar79@naver.com'                                       
FLIPSTAT= '        '                                                            
SWOWNER = 'Inhyun Cho' /        Licensed owner of software                      
COMMENT Original key: "END"                                                     
COMMENT                                                                         
COMMENT --Start of Astrometry.net WCS solution--                                
COMMENT --Put in by the new-wcs program--                                       
COMMENT                                                                         
WCSAXES =                    2 / no comment                                     
CTYPE1  = 'RA---TAN-SIP' / TAN (gnomic) projection + SIP distortions            
CTYPE2  = 'DEC--TAN-SIP' / TAN (gnomic) projection + SIP distortions            
EQUINOX =               2000.0 / Equatorial coordinates definition (yr)         
LONPOLE =                180.0 / no comment                                     
LATPOLE =                  0.0 / no comment                                     
CRVAL1  =        79.9797535396 / RA  of reference point                         
CRVAL2  =        45.1896148715 / DEC of reference point                         
CRPIX1  =        616.484283447 / X reference pixel                              
CRPIX2  =        526.408245087 / Y reference pixel                              
CUNIT1  = 'deg     ' / X pixel scale units                                      
CUNIT2  = 'deg     ' / Y pixel scale units                                      
CD1_1   =   -0.000685708286005 / Transformation matrix                          
CD1_2   =    6.11132326012E-05 / no comment                                     
CD2_1   =   -6.12471179978E-05 / no comment                                     
CD2_2   =   -0.000685877062746 / no comment                                     
IMAGEW  =                 1676 / Image width,  in pixels.                       
IMAGEH  =                 1266 / Image height, in pixels.                       
A_ORDER =                    2 / Polynomial order, axis 1                       
A_0_0   =                    0 / no comment                                     
A_0_1   =                    0 / no comment                                     
A_0_2   =   -1.20624858871E-06 / no comment                                     
A_1_0   =                    0 / no comment                                     
A_1_1   =   -1.88088919802E-07 / no comment                                     
A_2_0   =   -9.00166605547E-08 / no comment                                     
B_ORDER =                    2 / Polynomial order, axis 2                       
B_0_0   =                    0 / no comment                                     
B_0_1   =                    0 / no comment                                     
B_0_2   =   -4.56034470461E-08 / no comment                                     
B_1_0   =                    0 / no comment                                     
B_1_1   =    1.54493659934E-07 / no comment                                     
B_2_0   =   -5.48765540292E-07 / no comment                                     
AP_ORDER=                    2 / Inv polynomial order, axis 1                   
AP_0_0  =   -4.49372854066E-05 / no comment                                     
AP_0_1  =    3.51227393882E-07 / no comment                                     
AP_0_2  =    1.20633114157E-06 / no comment                                     
AP_1_0  =    -4.1883593719E-08 / no comment                                     
AP_1_1  =    1.88661495925E-07 / no comment                                     
AP_2_0  =    9.02394682508E-08 / no comment                                     
BP_ORDER=                    2 / Inv polynomial order, axis 2                   
BP_0_0  =    -3.9232026875E-05 / no comment                                     
BP_0_1  =   -4.69457676758E-08 / no comment                                     
BP_0_2  =    4.58272972396E-08 / no comment                                     
BP_1_0  =    1.59223501687E-07 / no comment                                     
BP_1_1  =   -1.54104023126E-07 / no comment                                     
BP_2_0  =    5.48803133579E-07 / no comment                                     
HISTORY Created by the Astrometry.net suite.                                    
HISTORY For more details, see http://astrometry.net.                            
HISTORY Git URL https://github.com/dstndstn/astrometry.net                      
HISTORY Git revision 0.73                                                       
HISTORY Git date Thu_Nov_16_08:30:44_2017_-0500                                 
HISTORY This WCS header was created by the program "blind".                     
DATE    = '2020-03-12T23:42:40' / Date this file was created.                   
COMMENT -- blind solver parameters: --                                          
COMMENT Index(0): /usr/local/astrometry/data/index-5000-09.fits                 
COMMENT Index(1): /usr/local/astrometry/data/index-5000-08.fits                 
COMMENT Index(2): /usr/local/astrometry/data/index-5000-07.fits                 
COMMENT Index(3): /usr/local/astrometry/data/index-5000-06.fits                 
COMMENT Index(4): /usr/local/astrometry/data/index-5000-05.fits                 
COMMENT Index(5): /usr/local/astrometry/data/index-5000-04.fits                 
COMMENT Index(6): /usr/local/astrometry/data/index-5000-03.fits                 
COMMENT Index(7): /usr/local/astrometry/data/index-5000-02.fits                 
COMMENT Index(8): /usr/local/astrometry/data/index-5000-01.fits                 
COMMENT Index(9): /usr/local/astrometry/data/index-5000-00.fits                 
COMMENT Index(10): /usr/local/astrometry/data/index-4119.fits                   
COMMENT Index(11): /usr/local/astrometry/data/index-4118.fits                   
COMMENT Index(12): /usr/local/astrometry/data/index-4117.fits                   
COMMENT Index(13): /usr/local/astrometry/data/index-4116.fits                   
COMMENT Index(14): /usr/local/astrometry/data/index-4115.fits                   
COMMENT Index(15): /usr/local/astrometry/data/index-4114.fits                   
COMMENT Index(16): /usr/local/astrometry/data/index-4113.fits                   
COMMENT Index(17): /usr/local/astrometry/data/index-4112.fits                   
COMMENT Index(18): /usr/local/astrometry/data/index-4111.fits                   
COMMENT Index(19): /usr/local/astrometry/data/index-4110.fits                   
COMMENT Index(20): /usr/local/astrometry/data/index-4109.fits                   
COMMENT Index(21): /usr/local/astrometry/data/index-4108.fits                   
COMMENT Index(22): /usr/local/astrometry/data/index-4107.fits                   
COMMENT Field name:                                                             
COMMENT   ./46P-Wirtanen_Light_R_2018-12-23-14-57-39_60.0sec_TMB130s        
COMMENT   s_STF-8300M_-29C_2x2bin.axy                                           
COMMENT Field scale lower: 0.214797 arcsec/pixel                                
COMMENT Field scale upper: 386.635 arcsec/pixel                                 
COMMENT X col name: X                                                           
COMMENT Y col name: Y                                                           
COMMENT Start obj: 10                                                           
COMMENT End obj: 20                                                             
COMMENT Solved_in: (null)                                                       
COMMENT Solved_out:                                                             
COMMENT   ./46P-Wirtanen_Light_R_2018-12-23-14-57-39_60.0sec_TMB130s            
COMMENT   s_STF-8300M_-29C_2x2bin.solved                                        
COMMENT Solvedserver: (null)                                                    
COMMENT Parity: 2                                                               
COMMENT Codetol: 0.01                                                           
COMMENT Verify pixels: 1 pix                                                    
COMMENT Maxquads: 0                                                             
COMMENT Maxmatches: 0                                                           
COMMENT Cpu limit: 300.000000 s                                                 
COMMENT Time limit: 0 s                                                         
COMMENT Total time limit: 0 s                                                   
COMMENT Total CPU limit: 300.000000 s                                           
COMMENT Tweak: yes                                                              
COMMENT Tweak AB order: 2                                                       
COMMENT Tweak ABP order: 2                                                      
COMMENT --                                                                      
COMMENT -- properties of the matching quad: --                                  
COMMENT index id: 4109                                                          
COMMENT index healpix: -1                                                       
COMMENT index hpnside: 0                                                        
COMMENT log odds: 137.823                                                       
COMMENT odds: 7.17355e+59                                                       
COMMENT quadno: 289052                                                          
COMMENT stars: 698612,698750,698616,698744                                      
COMMENT field: 3,10,16,13                                                       
COMMENT code error: 0.000793638                                                 
COMMENT nmatch: 17                                                              
COMMENT nconflict: 0                                                            
COMMENT nfield: 2191                                                            
COMMENT nindex: 18                                                              
COMMENT scale: 2.47866 arcsec/pix                                               
COMMENT parity: 1                                                               
COMMENT quads tried: 1221                                                       
COMMENT quads matched: 209                                                      
COMMENT quads verified: 208                                                     
COMMENT objs tried: 17                                                          
COMMENT cpu time: 0.021728                                                      
COMMENT --                                                                      
COMMENT                                                                         
COMMENT --End of Astrometry.net WCS--                                           
COMMENT --(Put in by the new-wcs program)--                                     
COMMENT                                                                         
    END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
'''    
                
        