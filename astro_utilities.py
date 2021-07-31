# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

2019.09.29  modify - missing 'IMAGETYP' on APT

ModuleNotFoundError: No module named 'ccdproc'
conda install -c condaforge ccdproc
"""

from datetime import datetime
#from astropy.io import fits
 
def write_log(log_file, log_str):
    import os
    with open(log_file, 'a') as log_f:
        log_f.write("{}, {}\n".format(os.path.basename(__file__), log_str))
    return print ("{}, {}\n".format(os.path.basename(__file__), log_str))

# =============================================================================
# for checking time
# =============================================================================
cht_start_time = datetime.now()
def print_working_time(cht_start_time):
    working_time = (datetime.now() - cht_start_time) #total days for downloading
    return print('working time ::: %s' % (working_time))

master_file_dir_name = 'master_file_Python/'
processing_dir_name = 'processing_Python/'
integration_dir_name = 'integration_Python/'
alignment_dir_name = 'alignment_Python/'

# =============================================================================
#     
# =============================================================================
def get_new_filename(fullname, **kargs):
    print('Starting get_new_filename ...\n{0}'.format(fullname))
    from astropy.io import fits
    hdul = fits.open(fullname)
    
    if hdul[0].header['NAXIS1'] == 2048 \
        and hdul[0].header['NAXIS2'] == 2048 :
        for binning in ['XBINNING', 'YBINNING'] :
            if not binning in hdul[0].header :
                with fits.open('{0}'.format(fullname), mode="append") as hdul1 :
                    hdul1[0].header.append(binning, '2', 'Binning factor in ')
                    hdul1.flush()
            elif hdul[0].header[binning] is None :
                with fits.open('{0}'.format(fullname), mode="update") as hdul1 :
                    hdul1[0].header[binning] = '2'
                    hdul1.flush()
            hdul[0].header[binning] = '2'

        hdul[0].header['INSTRUME'] = 'STX-16803' 
        hdul[0].header['TELESCOPE'] = 'RiLA600' 
        hdul[0].header['OPTIC'] = 'RiLA600' 
        
    if hdul[0].header['NAXIS1'] == 4096 \
        and hdul[0].header['NAXIS2'] == 4096 :
        for binning in ['XBINNING', 'YBINNING'] :
            if not binning in hdul[0].header :
                with fits.open('{0}'.format(fullname), mode="append") as hdul1 :
                    hdul1[0].header.append(binning, '1', 'Binning factor in ')
                    hdul1.flush()
            elif hdul[0].header[binning]  is None :
                with fits.open('{0}'.format(fullname), mode="update") as hdul1 :
                    hdul1[0].header[binning] = '1'
                    hdul1.flush()
            hdul[0].header[binning] = '1'
                
        hdul[0].header['INSTRUME'] = 'STX-16803' 
        hdul[0].header['TELESCOPE'] = 'RiLA600' 
        hdul[0].header['OPTIC'] = 'RiLA600' 
    
    if not 'INSTRUME' in hdul[0].header : 
        instrument = 'UNKNOWN'
    elif  'qsi' in hdul[0].header['INSTRUME'].lower() :     
        instrument = 'QSI683ws'
    elif  'st-8300' in hdul[0].header['INSTRUME'].lower() :     
        instrument = 'ST-8300M'
    elif  'stf-8300' in hdul[0].header['INSTRUME'].lower() :     
        instrument = 'STF-8300M'
    elif  'stl-11000' in hdul[0].header['INSTRUME'].lower() :     
        instrument = 'STL-11000M'
    else :
        instrument = hdul[0].header['INSTRUME']
    instrument = instrument.replace(" ","+")
    
    if 'CCD-TEMP' in hdul[0].header :     
        ccd_temp_el = str(hdul[0].header['CCD-TEMP']).split('.')
    else : 
        ccd_temp_el = 'NAN'
    
    if 'TIME-OBS' in hdul[0].header : 
        obs_date  = hdul[0].header['DATE-OBS'][:10]+'-'+hdul[0].header['TIME-OBS']
    elif 'DATE-OBS' in hdul[0].header :
        obs_date = hdul[0].header['DATE-OBS'][:19]
    else :
        obs_date = "No-obsdate"
    obs_date = obs_date.replace("T", "-")
    obs_date = obs_date.replace(":", '-')

    if 'EXPOSURE' in hdul[0].header : 
        esposure = "{:03d}".format(int(hdul[0].header['EXPOSURE']))
    elif 'EXPTIME' in hdul[0].header : 
        esposure = "{:03d}".format(int(hdul[0].header['EXPTIME']))
    else : 
        esposure = 'No_exptime' 
   
    if not 'OBJECT' in hdul[0].header : 
        object_name = '-'
    
    elif 'dark ' in hdul[0].header['OBJECT'] : 
        image_type = 'Dark'
        filter_name = '-'
        object_name = '-'
        optic = '-'
    elif 'bias ' in hdul[0].header['OBJECT'] : 
        image_type = 'Bias'
        filter_name = '-'
        object_name = '-'        
        optic = '-'
    elif 'Bias ' in hdul[0].header['OBJECT'] : 
        image_type = 'Bias'
        filter_name = '-'
        object_name = '-'        
        optic = '-'
    elif 'flat ' in hdul[0].header['OBJECT'] : 
        image_type = 'Flat'
        object_name = '-'
    elif hdul[0].header['OBJECT'] =='' : 
        object_name = '-'
    else : 
        object_name = hdul[0].header['OBJECT']
    
    if not 'IMAGETYP' in hdul[0].header : 
        image_type = '-'
        if 'FILTER' in hdul[0].header:
            filter_name = hdul[0].header['FILTER']
        else :
            filter_name = "-"

        if 'OBJECT' in hdul[0].header:
            object_name = hdul[0].header['OBJECT']
        else:
            object_name = "-"
    elif hdul[0].header['IMAGETYP'][:1].lower() == 'b' \
        or hdul[0].header['IMAGETYP'][0:1].lower() == 'z':
        image_type = 'Bias'
        filter_name = '-'
        object_name = '-'
        optic = '-'
    elif hdul[0].header['IMAGETYP'][0:1] == 'D' or hdul[0].header['IMAGETYP'][0:1] == 'd':
        image_type = 'Dark'
        filter_name = '-'
        object_name = '-'
        optic = '-'
    elif hdul[0].header['IMAGETYP'][0:1] == 'F' or hdul[0].header['IMAGETYP'][0:1] == 'f':
        image_type = 'Flat'
        object_name = '-'
    elif hdul[0].header['IMAGETYP'][0:1] == 'L' :
        image_type = 'Light'
        #filter_name = hdul[0].header['FILTER'] 
        object_name = hdul[0].header['OBJECT']   
    elif hdul[0].header['IMAGETYP'][0:1] == 'o' :
        image_type = 'Light'
        filter_name = hdul[0].header['FILTER'] 
        object_name = hdul[0].header['OBJECT']
        
    if not 'FILTER' in hdul[0].header : 
        filter_name = '-'
    elif hdul[0].header['FILTER'] == 'Ha' :
        filter_name = 'H'
    elif hdul[0].header['FILTER'] == 'S2' :
        filter_name = 'S'
    elif hdul[0].header['FILTER'] == 'O3' :
        filter_name = 'O'
    elif hdul[0].header['FILTER'] == 'Luminance' :
        filter_name = 'L'
    elif hdul[0].header['FILTER'] == 'Blue' :
        filter_name = 'B'
    elif hdul[0].header['FILTER'] == 'Green' :
        filter_name = 'L'
    elif hdul[0].header['FILTER'] == 'Red' :
        filter_name = 'R'
    else : 
        filter_name = hdul[0].header['FILTER'] 

    if not 'XBINNING' in hdul[0].header \
        or not 'YBINNING' in hdul[0].header : 
        xbin = '-'
        ybin = '-'
    else : 
        xbin = hdul[0].header['XBINNING']
        ybin = hdul[0].header['YBINNING']

    if isinstance(xbin, float) : 
        xbin = int(xbin)
        xbin = str(xbin)
    if isinstance(ybin, float) : 
        ybin = int(ybin)
        ybin = str(ybin)
    
    if not 'CTYPE1' in hdul[0].header \
        or not 'CTYPE2' in hdul[0].header : 
        wcs = "-"
    elif hdul[0].header['CTYPE1'] == 'RA---TAN' \
        and hdul[0].header['CTYPE2'] == 'DEC--TAN':
        wcs = "wcs"
    
    object_name = object_name.replace('_', '-')
    object_name = object_name.replace(':', '-')
    object_name = object_name.replace('.', '-')
    object_name = object_name.replace(' ', '')
    object_name = object_name.replace('NGC', 'N')
    object_name = object_name.replace('ngc', 'N')
    object_name = object_name.replace('bias', '-')
    object_name = object_name.replace('Bias', '-')
    object_name = object_name.replace('dark', '-')
    object_name = object_name.replace('Dark', '-')
    object_name = object_name.replace('flat', '-')
    object_name = object_name.replace('Flat', '-')
    
    if not 'OPTIC' in hdul[0].header : 
        optic = 'OPTIC'
    else :
        optic = hdul[0].header['OPTIC']
        
    new_filename = '{0}_{1}_{2}_{3}_{4}sec_{5}_{6}_{7}C_{8}bin_{9}.fit'\
        .format(object_name.upper(),
        image_type,
        filter_name,
        obs_date,
        esposure,
        optic,
        instrument,
        ccd_temp_el[0],
        xbin,
        wcs)
    hdul.close()
    return new_filename


def get_new_foldername_from_filename(filename):
    log_file = 'get_new_foldername.log'
    print('Starting get_new_foldername ...\n{0}'.format(filename))
    
    filename_el = filename[:-4].split("_")
    from datetime import datetime, timedelta
    timez = 9
    if int(filename_el[3][17:19])>=60 :
        obs_UT = datetime.strptime("{}59".format(filename_el[3][:17]), '%Y-%m-%d-%H-%M-%S')
    else:
         obs_UT = datetime.strptime(filename_el[3], '%Y-%m-%d-%H-%M-%S')
    obs_LST = obs_UT + timedelta(hours = timez)
    if obs_LST.hour < 12 :
        obs_LST = obs_LST - timedelta(days = 1)
    filename_el[3] = obs_LST.strftime('%Y-%m-%d-%H-%M-%S')
    if filename_el[1] == 'Bias':
        new_foldername = '{6}_{8}/Cal/-_{3}_-_{1}_-_{4}_-_{6}_-_{8}/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    elif filename_el[1] == 'Dark' :
        new_foldername = '{6}_{8}/Cal/-_{3}_-_{1}_-_{4}_-_{6}_-_{8}/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    elif filename_el[1] == 'Flat' :
        new_foldername = '{6}_{8}/Cal_{5}/-_{3}_-_{1}_-_{5}_{6}_-_{8}/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    else : 
        new_foldername = '{6}_{8}/Light_{5}/{0}_{1}_-_{3}_-_{5}_{6}_-_{8}/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    write_log(log_file, 
                '{1} ::: \nNew foldername is {0} ...'\
                .format(new_foldername, datetime.now()))    
    return new_foldername


def get_new_foldername(filename):
    log_file = 'get_new_foldername.log'
    print('Starting get_new_foldername ...\n{0}'.format(filename))
    
    filename_el1 = filename.split("bin")
    filename_el = filename_el1[0].split("_")
    
    if filename_el[1] == 'Bias':
        new_foldername = '{6}_{8}bin/Cal/-_{3}_-_{1}_-_{4}_-_{6}_-_{8}bin/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    elif filename_el[1] == 'Dark' :
        new_foldername = '{6}_{8}bin/Cal/-_{3}_-_{1}_-_{4}_-_{6}_-_{8}bin/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    elif filename_el[1] == 'Flat' :
        new_foldername = '{6}_{8}bin/Cal_{5}/-_{3}_-_{1}_-_{5}_{6}_-_{8}bin/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    else : 
        new_foldername = '{6}_{8}bin/Light_{5}/{0}_{1}_-_{3}_-_{5}_{6}_-_{8}bin/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    write_log(log_file, 
                '{1} ::: \nNew foldername is {0} ...'\
                .format(new_foldername, datetime.now()))    
    return new_foldername

def getFullnameListOfallFiles(dirName):
    ##############################################3
    import os
    # create a list of file and sub directories 
    # names in the given directory 
    listOfFile = sorted(os.listdir(dirName))
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getFullnameListOfallFiles(fullPath)
        else:
            allFiles.append(fullPath)
                
    return allFiles


def getFullnameListOfallsubDirs1(dirName):
    ##############################################3
    import os
    allFiles = list()
    for file in sorted(os.listdir(dirName)):
        d = os.path.join(dirName, file)
        allFiles.append(d)
        if os.path.isdir(d):
            allFiles.extend(getFullnameListOfallsubDirs1(d))

    return allFiles

def getFullnameListOfallsubDirs(dirName):
    ##############################################3
    import os
    allFiles = list()
    for it in os.scandir(dirName):
        if it.is_dir():
            allFiles.append(it.path)
            allFiles.extend(getFullnameListOfallsubDirs(it))

    return allFiles


                                
def connectMariaDB():
    #import pymysql
    import pymysql.cursors
    #conda install pymysql
    
    #mariaDB info
    db_host = 'parksparks.iptime.org'
    db_user = 'root'
    db_pass = 'rlgusl01'
    db_name = 'CCD_obs'
    db_port = 3307
        
    conn = pymysql.connect(host = db_host,
                          port = db_port,
                          user = db_user, password = db_pass,
                          db = db_name, charset = 'utf8mb4',
                          cursorclass = pymysql.cursors.DictCursor)
    
    return conn

import numpy as np
from astropy.io import fits
from ccdproc import combine
from datetime import datetime
#ModuleNotFoundError: No module named 'ccdproc'
#conda install -c conda-forge ccdproc


def print_subworking_time(sub_start_time):

    from datetime import datetime
    working_time = (datetime.now() - cht_start_time) #total days for downloading
    return print('working time ::: %s' % (working_time))



def sub_p_solve_field(fullname, save_dir_name, sub_start_time): 
    import subprocess
    print('-'*60)
    print(fullname)
    with subprocess.Popen(['solve-field', 
                           '-O', #--overwrite: overwrite output files if they already exist
                           #'--scale-units', 'arcsecperpix', #pixel scale
                           #'--scale-low', '0.1', '--scale-high', '0.40', #pixel scale
                           '-g', #--guess-scale: try to guess the image scale from the FITS headers
                           #'-p', # --no-plots: don't create any plots of the results
                           '-D', '{0}'.format(save_dir_name), 
                           '{0}'.format(fullname)], 
                          stdout=subprocess.PIPE) as proc :
        print(proc.stdout.read())
        print(print_subworking_time(sub_start_time))
        '''
        solve-field -O fullname
       '''
    return 0

def align_image(im1, im2):
    import cv2
    #conda install -c conda-forge opencv

    #code from https://www.learnopencv.com/image-alignment-ecc-in-opencv-c-python/
    # Convert images to grayscale
    #im1_gray = cv2.cvtColor(im1,cv2.COLOR_BGR2GRAY)
    #im2_gray = cv2.cvtColor(im2,cv2.COLOR_BGR2GRAY)

    im1_gray = im1
    im2_gray = im2
    
    im1_32f_gray = np.array(im1_gray/65536.0, dtype=np.float32)
    im2_32f_gray = np.array(im2_gray/65536.0, dtype=np.float32)
    
    # Find size of image1
    sz = im1.shape
    # Define the motion model
    #warp_mode = cv2.MOTION_TRANSLATION
    warp_mode = cv2.MOTION_EUCLIDEAN
    # Define 2x3 or 3x3 matrices and initialize the matrix to identity
    if warp_mode == cv2.MOTION_HOMOGRAPHY :
        warp_matrix = np.eye(3, 3, dtype=np.float32)
    else :
        warp_matrix = np.eye(2, 3, dtype=np.float32)
    # Specify the number of iterations.
    number_of_iterations = 1000  #5000
    # Specify the threshold of the increment
    # in the correlation coefficient between two iterations
    termination_eps = 1e-7   #1e-10
    # Define termination criteria
    criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, number_of_iterations, termination_eps)
    
   # Run the ECC algorithm. The results are stored in warp_matrix.
    (cc, warp_matrix) = cv2.findTransformECC (im1_32f_gray, im2_32f_gray, warp_matrix, warp_mode, criteria)
    if warp_mode == cv2.MOTION_HOMOGRAPHY :
        # Use warpPerspective for Homography 
        im2_aligned = cv2.warpPerspective (im2, warp_matrix, (sz[1],sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP)
    else :
        # Use warpAffine for Translation, Euclidean and Affine
        im2_aligned = cv2.warpAffine(im2, warp_matrix, (sz[1],sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP);
    # Show final results
    return im2_aligned


def combine_BiasDark(file_list, c_method, 
        base_dir_name, master_file_dir_name, current_dir_name) :
        
    try :
        
        combine_result = combine(file_list,       # ccdproc does not accept numpy.ndarray, but only python list.
               method = c_method,  # default is average so I specified median.
               unit='adu')  
        
        combine_result.data = np.array(combine_result.data, dtype=np.float32)

        combine_result.write('{0}/{1}{2}_master_{3}_float32.fit'\
                  .format(base_dir_name, master_file_dir_name, 
                  current_dir_name, c_method), overwrite =True, format='fits')
        
        ##### fits header update
        with fits.open('{0}/{1}{2}_master_{3}_float32.fit'\
              .format(base_dir_name, master_file_dir_name, current_dir_name, c_method),
              mode='update') as hdul:
            hdul[0].header.append(('COMMENT', ', '.join(file_list), 'combine file list'))
            hdul[0].header.append('COMMENT', 
                                  '{0}'.format(len(file_list)), 
                                  'combine file number')
    
    except Exception as err :
        print('{5} ::: {4} with {0}/{1}{2}_master_{3}_float32.fit ...'\
            .format(base_dir_name, master_file_dir_name, 
            current_dir_name, c_method, err, datetime.now()))

    return 0

def combine_Flat(file_list, c_method, 
        base_dir_name, master_file_dir_name, current_dir_name, chl) :
        
    try :
        
        combine_result = combine(file_list,       # ccdproc does not accept numpy.ndarray, but only python list.
               method = c_method,  # default is average so I specified median.
               unit='adu')  
        
        combine_result.data = np.array(combine_result.data, dtype=np.float32)

        combine_result.write('{0}/{1}{2}_master_{3}_{4}_float32.fit'\
                  .format(base_dir_name, master_file_dir_name, 
                  current_dir_name, c_method, chl), overwrite =True, format='fits')
        
        ##### fits header update
        with fits.open('{0}/{1}{2}_master_{3}_{4}_float32.fit'\
              .format(base_dir_name, master_file_dir_name, current_dir_name, c_method, chl),
              mode='update') as hdul:
            hdul[0].header.append(('COMMENT', ', '.join(file_list), 'combine file list'))
            hdul[0].header.append('COMMENT', 
                                  '{0}'.format(len(file_list)), 
                                  'combine file number')
    
    except Exception as err :
        print('{6} ::: {5} with {0}/{1}{2}_master_{3}_{4}_float32.fit ...'\
            .format(base_dir_name, master_file_dir_name, 
            current_dir_name, c_method, chl, err, datetime.now()))

    return 0


def combine_master_file(file_list, c_method, 
        base_dir_name, master_file_dir_name, current_dir_name) :
        
    try :
        combine_result = combine(file_list,       # ccdproc does not accept numpy.ndarray, but only python list.
                       method = c_method,  # default is average so I specified median.
                       sigma_clip = True, sigma_clip_low_thresh=3, sigma_clip_high_thresh=3,
                       unit = 'adu')              # unit is required: it's ADU in our case.
        
        combine_result.data = np.array(combine_result.data, dtype=np.float32)

        combine_result.write('{0}/{1}{2}_master_{3}_float32.fit'\
                  .format(base_dir_name, master_file_dir_name, 
                  current_dir_name, c_method), overwrite =True, format='fits')
        
        ##### fits header update
        with fits.open('{0}/{1}{2}_master_{3}_float32.fit'\
              .format(base_dir_name, master_file_dir_name, current_dir_name, c_method),
              mode='update') as hdul:
            hdul[0].header.append(('COMMENT', ', '.join(file_list), 'combine file list'))
            hdul[0].header.append('COMMENT', 
                                  '{0}'.format(len(file_list)), 
                                  'combine file number')
    
    except Exception as err :
        print('{5} ::: {4} with {0}/{1}{2}_master_{3}_float32.fit ...'\
            .format(base_dir_name, master_file_dir_name, 
            current_dir_name, c_method, err, datetime.now()))

    return 0

def combine_master_flat_file(file_list, c_method, 
         base_dir_name, master_file_dir_name, current_dir_name, chl) :
        
    try :
        combine_result = combine(file_list,       # ccdproc does not accept numpy.ndarray, but only python list.
                       method = c_method,  # default is average so I specified median.
                       sigma_clip = True, sigma_clip_low_thresh=3, sigma_clip_high_thresh=3,
                       unit = 'adu')              # unit is required: it's ADU in our case.
        
        combine_result.data = np.array(combine_result.data, dtype=np.float32)
        
        combine_result.write('{0}{1}{2}_master_{3}_{4}_float32.fit'\
                  .format(base_dir_name, master_file_dir_name, 
                  current_dir_name, c_method, chl), overwrite =True, format='fits')
        ##### fits header update
        with fits.open('{0}{1}{2}_master_{3}_{4}_float32.fit'\
                  .format(base_dir_name, master_file_dir_name, 
                  current_dir_name, c_method, chl),
                  mode='update') as hdul:
            hdul[0].header.append(('COMMENT', ', '.join(file_list), 'combine file list'))
            hdul[0].header.append('COMMENT', 
                                  '{0}'.format(len(file_list)), 
                                  'combine file number')


    except Exception as err :
        print('{6} ::: {5} with {0}{1}{2}_master_{3}_{4}_float32.fit ...'\
                  .format(base_dir_name, master_file_dir_name, 
                  current_dir_name, c_method, chl, err, datetime.now()))
                
    return 0