# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 01:00:19 2018
@author: user

2019.09.29  modify - missing 'IMAGETYP' on APT

ModuleNotFoundError: No module named 'ccdproc'
conda install -c condaforge ccdproc
"""
#%%
from pathlib import Path
from astropy.io import fits
import subprocess
from datetime import datetime, timedelta
import os
from pathlib import Path
import numpy as np
from ccdproc import combine
import shutil

import Python_utilities

#%%
#########################################
#directory variables
#########################################
c_method = "median"
#CCD_obs_dir = "../CCD_obs_raw/"
#base_dir = "R:\CCD_obs\CCD_obs_raw"
#CCD_obs_raw_dir = "/mnt/Rdata/CCD_obs/CCD_obs_raw/"
CCD_obs_raw_dir = "CCD_obs_raw"

#base_dir = "../RnE_2022/RiLA600_STX-16803_2bin/"
#base_dir = "/mnt/Rdata/CCD_obs/RnE_2022/RiLA600_STX-16803_2bin/"
#base_dir = "R:\CCD_obs\RnE_2022\RiLA600_STX-16803_2bin"
#base_dir = "/mnt/Rdata/CCD_obs/RiLA600_2022"
base_dir = "R:\CCD_obs\RiLA600_2022"
CCD_NEW_dir = "CCD_new_files"
CCD_duplicate_dir = "CCD_duplicate_files"


master_dir = "master_files_ys"
reduced_dir = "reduced"
reduced_dir2 = "reduced2"
solved_dir = "solved"
solved_dir2 = "solved2"
DAOfinder_result_dir = "DAOfinder_result"
IRAFfinder_result_dir = "IRAFfinder_result"
APh_result_dir = "APh_result"
Asteroid_result_dir = "Asteroid_result"

#
master_file_dir = 'master_file_Python/'
processing_dir = 'processing_Python/'
integration_dir = 'integration_Python/'
alignment_dir = 'alignment_Python/'


#######################################################
# OBS instruments information 
#######################################################
gain = 0
rdnoise = 0
binning = 1

GAINDIC = {"STF-8300M": 0.37, 
        "STX-16803": 1.27, 
        "STL-11000M": 0.8, 
        "QSI683ws": 0.13 } 

RDNOISEDIC = {"STF-8300M": 9.3, 
            "STX-16803": 9.0, 
            "STL-11000M": 9.6, 
            "QSI683ws": 8.0 } 

PIXSIZEEDIC = {"STF-8300M": 5.4, 
            "STX-16803": 9.0, 
            "STL-11000M": 9.0, 
            "QSI683ws": 5.4 } 


PIXSCALEDIC = {"FSQ106ED_STF-8300M": 2.1, 
            "FSQ106ED-x72_STF-8300M": 2.1/0.72, 
            "FSQ106ED-x73_STF-8300M": 2.1/0.73, 
            "FS60CB_STF-8300M": 3.14, 
            "RiLA600_STX-16803": 0.62, 
            "FSQ106ED_STL-11000": 3.5, 
            "FSQ106ED-x72_STL-11000": 3.5/0.72, 
            "FSQ106ED-x73_STL-11000": 3.5/0.73, 
            "FSQ106ED_QSI683ws": 2.1,
            "FSQ106ED-x72_QSI683ws": 2.1/0.72,
            "FSQ106ED-x73_QSI683ws": 2.1/0.73
             } 

FOCALLENDIC = {"TMB130ss": 910, 
            "TMB130ss-x75": 910*0.75, 
            "RiLA600": 3000, 
            "GSON300": 1200, 
            "FS60CB": 355, 
            "SVX080T": 480, 
            "SVX080T-x80": 480*0.8, 
            "FSQ106ED": 530, 
            "FSQ106ED-x72": 530*0.72,
            "FSQ106ED-x73": 530*0.73} 

#CCDNAME, PIXSIZE, GAIN, RENOISE    
CCDDIC = {"STF-8300M": [5.4, 0.37, 9.3], 
        "QSI683ws": [5.4, 0.13, 8.0],
        "STL-11000M": [9.0, 0.8, 9.6],
        "STX-16803": [9.0, 1.27, 9.0]}
        
# OPTICDIC = {"TMB130ss": [910], 
#             "RiLA600": [3000], 
#             "GSON300": [1200], 
#             "FS-60CB": [355], 
#             "SVX080T": [480], 
#             "FSQ106ED": [530]} 

#OPTICNAME, Focal_length, Aperature
OPTICDIC = {"TMB130ss": [910, 4], 
            "RiLA600": [3000, 5], 
            "GSON300": [1200, 4], 
            "FS60CB": [355, 355/60], 
            "SVX080T": [480, 6], 
            "FSQ106ED": [530, 5]}

#Optic Acc name
OptAccDIC = {"x80": 0.8, 
            "x72": 0.72, 
            "x73": 0.73, 
            "x75": 0.75} 
#######################################################

def calPixScale (F_length, Opt_acc, Pix_size) :
    # Pixel Size   /   Telescope Focal Length   )   X 206.265  
    # Pixel Size : um, 
    # Telescope Focal Length: mm
    # PixScale: arcsec / pixel
    PIXScale = Pix_size / (F_length * Opt_acc) *  206.265
    return PIXScale


#%%
KevinFitsHeader():
    def __init__(self, fpath):
        self.fpath = Path(fpath)
        self.checkKEYs = ["OBJECT", "TELESCOP", "OPTIC", "CCDNAME", 'FILTER',
            "GAIN", "EGAIN", "RDNOISE", "FOCALLEN", "PIXSCALE",
            "XBINNING", "YBINNING", "FLIPSTAT"]
        """
        Parameters
        ----------
        fpath : string
            The fullname of input file...

        """
    def append_header(self):
        with fits.open(str(self.fpath), mode="append") as self.hdul :
            for self.checkKEY in self.checkKEYs: 
                if not self.checkKEY in self.hdul[0].header :
                    self.hdul[0].header.append(self.checkKEY, 
                                    '', 
                                    f"The keyword '{self.checkKEY}' is added.") 
                print(f"{self.checkKEY}: ", self.hdul[0].header[self.checkKEY])

            self.hdul.flush()  # changes are written back to original.fits
        return self.hdul

    def update_header(self): 
        self.foldername_el = self.fpath.parts[-2].split('_')
        self.fname_el = self.fpath.parts[-1].split('_')
        print("foldername_el", self.foldername_el)
        print("fname_el", self.fname_el)
        self.object_name = self.foldername_el[0]
        self.filter_name = self.fname_el[2]
        self.optic_name = self.foldername_el[5]
        self.ccd_name = self.foldername_el[6]
        print("object_name", self.object_name)
        print("filter_name", self.filter_name)
        print("optic_name", self.optic_name)
        print("ccd_name", self.ccd_name)
   
        # Change something in hdul.
        with fits.open(str(self.fpath), mode="update") as self.hdul :
            
            #if object_name != hdul[0].header["OBJECT"] : 
            self.hdul[0].header["OBJECT"] = self.object_name.upper()
            print(f"The 'OBJECT' is set {self.object_name.upper()}")
        
            if len(self.hdul[0].header['DATE-OBS']) == 10 \
                and 'TIME-OBS' in self.hdul[0].header : 
                self.hdul[0].header['DATE-OBS'] += 'T' + self.hdul[0].header['TIME-OBS']
                print(f"The 'DATE-OBS' is set {self.hdul[0].header['DATE-OBS']}")

            if "ze" in self.hdul[0].header["IMAGETYP"].lower() \
                    or "bi" in self.hdul[0].header["IMAGETYP"].lower() :
                self.hdul[0].header["IMAGETYP"] = "BIAS"
                print(f"The 'IMAGETYP' is set {self.hdul[0].header['IMAGETYP']}")
            elif "da" in self.hdul[0].header["IMAGETYP"].lower() :
                self.hdul[0].header["IMAGETYP"] = "DARK"
                print(f"The 'IMAGETYP' is set {self.hdul[0].header['IMAGETYP']}")
            elif "fl" in self.hdul[0].header["IMAGETYP"].lower() :
                self.hdul[0].header["IMAGETYP"] = "FLAT"
                print(f"The 'IMAGETYP' is set {self.hdul[0].header['IMAGETYP']}")
            elif "da" in self.hdul[0].header["IMAGETYP"].lower() \
                    or "lig" in self.hdul[0].header["IMAGETYP"].lower() :
                self.hdul[0].header["IMAGETYP"] = "LIGHT"
                print(f"The 'IMAGETYP' is set {self.hdul[0].header['IMAGETYP']}")
            
            if "BIAS" in self.hdul[0].header["IMAGETYP"] \
                or "DARK" in self.hdul[0].header["IMAGETYP"] :
                self.hdul[0].header["FILTER"] = "-"
                print(f"The 'FILTER' is set {self.hdul[0].header['FILTER']}")
                self.hdul[0].header['OPTIC'] = "-"
                print(f"The 'OPTIC' is set {self.hdul[0].header['OPTIC']}")

            if "FLAT" in self.hdul[0].header["IMAGETYP"] \
                or "LIGHT" in self.hdul[0].header["IMAGETYP"] :
                self.hdul[0].header["FILTER"] = self.filter_name.upper()
                print(f"The 'FILTER' is set {self.hdul[0].header['FILTER']}")
                self.hdul[0].header["OPTIC"] = self.optic_name
                print(f"The 'OPTIC' is set {self.hdul[0].header['OPTIC']}")

            try : 
                if 'qsi' in self.hdul[0].header['INSTRUME'] : 
                    self.CCDNAME = 'QSI683ws'
                elif 'st-8300' in self.hdul[0].header['INSTRUME'] : 
                    self.CCDNAME = 'ST-8300M'
                elif 'stf-8300' in self.hdul[0].header['INSTRUME'] : 
                    self.CCDNAME = 'STF-8300M'
                elif '11000' in hdul[0].header['INSTRUME'] : 
                    self.CCDNAME = 'STL-11000M'
                elif '16803' in self.hdul[0].header['INSTRUME'] : 
                    self.CCDNAME = 'STX-16803'
                elif "SBIG" in self.hdul[0].header['INSTRUME'] :
                    if self.hdul[0].header['XPIXSZ'] == 5.4 \
                            or self.hdul[0].header['XPIXSZ'] == 10.8 :
                        self.CCDNAME = 'STF-8300M'
                    elif self.hdul[0].header['XPIXSZ'] == 9.0 \
                            or self.hdul[0].header['XPIXSZ'] == 18.0 :
                        if self.hdul[0].header['NAXIS1'] == 2048 \
                                or self.hdul[0].header['NAXIS1'] == 4096 :
                            self.CCDNAME = 'STX-16803'
                        elif self.hdul[0].header['NAXIS1'] == 4008 \
                                or  self.hdul[0].header['NAXIS1'] == 2672 \
                                or  self.hdul[0].header['NAXIS1'] == 2004 \
                                or  self.hdul[0].header['NAXIS1'] == 1336 :
                            self.CCDNAME = 'STL-11000M'
                else :
                    self.CCDNAME = self.ccd_name    
            except :
                self.CCDNAME = self.ccd_name
            print("CCDNAME", self.CCDNAME)

            self.hdul[0].header["CCDNAME"] = self.CCDNAME
            print(f"The 'CCDNAME' is set {self.CCDNAME}...")

            if not "CCD-TEMP" in self.hdul[0].header :
                self.hdul[0].header['CCD-TEMP'] = 'N'
                print(f"The 'CCD-TEMP' is set {self.hdul[0].header['CCD-TEMP']}...")

            if not "EXPOSURE" in self.hdul[0].header :
                self.hdul[0].header["EXPOSURE"] = self.hdul[0].header["EXPTIME"]
                print(f"The 'EXPOSURE' is set {self.hdul[0].header['EXPTIME']}...")

            self.hdul[0].header['GAIN'] = GAINDIC[self.CCDNAME]
            print(f"The 'GAIN' is set {self.hdul[0].header['GAIN']}...")
            self.hdul[0].header['EGAIN'] = GAINDIC[self.CCDNAME]
            print(f"The 'EGAIN' is set {self.hdul[0].header['EGAIN']}...")
            self.hdul[0].header['RDNOISE'] = RDNOISEDIC[self.CCDNAME]
            print(f"The 'RDNOISE' is set {self.hdul[0].header['RDNOISE']}...")
            self.hdul.flush()  # changes are written back to original.fits
            print('*'*30)
            print(f"The header of {self.fpath.name} is updated..")
        return self.hdul

#%%
#########################################
#single  KevinPSolver
#########################################
class KevinSolver():
    def __init__(self, fullname, solved_dir):
        self.fullname = fullname
        self.solved_dir = solved_dir
        """
        Parameters
        ----------
        fullname : string
            The fullname of input file...

        solved dir: string
            The directory where the output file              
        """
    def astap(self):
        print("Starting... \n{}".format(self.fullname))
        self.fullname_el = self.fullname.split("/")
        self.filename_el = self.fullname_el[-1].split("_")

        #Path(os.path.dirname(str(f_path))).parents[0]
        if os.path.exists('{}/{}/{}'.format((os.path.dirname(self.fullname)), 
                                            self.solved_dir, self.fullname_el[-1])):
            print("{} is already solved ...".format(self.fullname_el[-1]))

        else : 
            try:    
                with subprocess.Popen(['astap', 
                            '-f', '{0}'.format(self.fullname), 
                            '-o', 
                            '{}/{}/{}.tmp'.format((os.path.dirname(self.fullname)), 
                                self.solved_dir, self.fullname_el[-1][:-5]), 
                            '-wcs',
                            '-analyse2',
                            '-update',],
                            stdout=subprocess.PIPE) as proc :
                    print(proc.stdout.read())
                                
            except Exception as err :
                print('{1} ::: {2} with {0} ...'\
                            .format(self.fullname, datetime.now(), err))

    def astap(self):
        print("Starting... \n{}".format(self.fullname))
        self.fullname_el = self.fullname.split("/")
        self.filename_el = self.fullname_el[-1].split("_")

        #Path(os.path.dirname(str(f_path))).parents[0]
        if os.path.exists('{}/{}/{}'.format((os.path.dirname(self.fullname)), 
                                            self.solved_dir, self.fullname_el[-1])):
            print("{} is already solved ...".format(self.fullname_el[-1]))

        else : 
            try:
                with subprocess.Popen(['solve-field', 
                                        '-O', #--overwrite: overwrite output files if they already exist
                                        #'--scale-low', '0.1', '--scale-high', '0.40', #pixel scale
                                        '-g', #--guess-scale: try to guess the image scale from the FITS headers
                                        '--cpulimit', '120',  #will make it give up after 30 seconds.
                                        '--nsigma', '15',
                                        #'--downsample', '4',
                                        '-u', 'app', #'--scale-units', 'arcsecperpix', #pixel scale
                                        '-L', '1.2', '-U', '1.3',
                                        #'-N',  '{}'.format(self.fullname[-1]), #--new-fits <filename>: output filename of the new FITS file containingthe WCS header; "none" to not create this file
                                        '-p', 
                                        '--no-plots',#: don't create any plots of the results
                                        '-D', '{}/{}/'.format((os.path.dirname(self.fullname)), self.solved_dir),
                                        '{0}'.format(self.fullname)], 
                                        stdout=subprocess.PIPE) as proc :
                    print(proc.stdout.read())
            except Exception as err :
                print('{1} ::: {2} with {0} ...'\
                            .format(self.fullname, datetime.now(), err))
#########################################

#%%
#########################################
# Astrometry Solver
#########################################
def AstrometrySolver(fullname, 
                    solved_dir
                    ): 
    """
    Parameters
    ----------
    fullname : string
        The fullname of input file...

    solved dir: string
        The directory where the output file


    """

    fpath = Path(fullname)
    SOLVEDDIR = Path(solved_dir)

    print("Starting... \n{}".format(str(fullname)))
    #fullname_el = fullname.split("/")
    #filename_el = fullname_el[-1].split("_")

    print("str(SOLVEDDIR):", str(SOLVEDDIR))
    print('str(SOLVEDDIR/fpath.name)', SOLVEDDIR/fpath.name)

    if (SOLVEDDIR/fpath.name).exists():
        print(f"{str((SOLVEDDIR/fpath.name))} is already exist...")
    
    else: 
        try : 
            # solve command.
            # solve-field fullname.fit -O --cpulimit 120 --nsigma 15 -u app -L 1.2 -U 1.3 -N new_filename.fits -p --no-plots -D output_directory {0}
            with subprocess.Popen(['solve-field', 
                                    '-O', #--overwrite: overwrite output files if they already exist
                                    #'--scale-low', '0.1', '--scale-high', '0.40', #pixel scale
                                    #'-g', #--guess-scale: try to guess the image scale from the FITS headers
                                    '--cpulimit', '10',  #will make it give up after 30 seconds.
                                    '--nsigma', '15',
                                    '--downsample', '4',
                                    '-u', 'app', #'--scale-units', 'arcsecperpix', #pixel scale
                                    #'-L', '1.2', '-U', '1.3', #16803_2bin
                                    '-L', '0.60', '-U', '0.63',   #16803_1bin
                                    #'-N',  '{}'.format(fullname[-1]), #--new-fits <filename>: output filename of the new FITS file containingthe WCS header; "none" to not create this file
                                    '-p', 
                                    #'--no-plots',#: don't create any plots of the results
                                    '-D', str(SOLVEDDIR),
                                    str(fpath)
                                    ], 
                                    stdout=subprocess.PIPE) as proc :
                print(proc.stdout.read())
            
            if (SOLVEDDIR/f'{fpath.stem}.new').exists():
                print(f"{fpath.name} is solved successfully ...")
                shutil.move(str((SOLVEDDIR/f'{fpath.stem}.new')), \
                            str((SOLVEDDIR/f'{fpath.stem}.fit')))
                print(f"{str(SOLVEDDIR/f'{fpath.stem}.new')} is renamed to fits ...")
            
            else : 
                print(f"solving {fpath.name} is failed...")
            
        except Exception as err :
            print('{1} ::: {2} with {0} ...'\
                .format(fullname, datetime.now(), err))


# #%%
# #########################################
# #single calss Astrometry Solver
# #########################################
# class AstrometrySolver():
#     def __init__(self, fullname, solved_dir):
        
#         """
#         Parameters
#         ----------
#         fullname : string
#             The fullname of input file...

#         solved dir: string
#             The directory where the output file              
#         """

#         self.fullname = fullname
#         self.solved_dir = solved_dir

#         print("Starting... \n{}".format(self.fullname))
#         self.fullname_el = self.fullname.split("/")
#         self.filename_el = self.fullname_el[-1].split("_")

#         print("self.solved_dir:", self.solved_dir)
#         print('{}/{}'.format(self.solved_dir, self.fullname_el[-1]))

#         if os.path.exists('{}/{}'.format(self.solved_dir, self.fullname_el[-1])):
#             print("{} is already solved ...".format(self.fullname_el[-1]))
        
#         else: 

#             try : 
#                 # solve command.
#                 # solve-field fullname.fit -O --cpulimit 120 --nsigma 15 -u app -L 1.2 -U 1.3 -N new_filename.fits -p --no-plots -D output_directory {0}
#                 with subprocess.Popen(['solve-field', 
#                                         '-O', #--overwrite: overwrite output files if they already exist
#                                         #'--scale-low', '0.1', '--scale-high', '0.40', #pixel scale
#                                         '-g', #--guess-scale: try to guess the image scale from the FITS headers
#                                         '--cpulimit', '120',  #will make it give up after 30 seconds.
#                                         '--nsigma', '15',
#                                         #'--downsample', '4',
#                                         '-u', 'app', #'--scale-units', 'arcsecperpix', #pixel scale
#                                         '-L', '1.2', '-U', '1.3',
#                                         #'-N',  '{}'.format(self.fullname[-1]), #--new-fits <filename>: output filename of the new FITS file containingthe WCS header; "none" to not create this file
#                                         '-p', 
#                                         '--no-plots',#: don't create any plots of the results
#                                         '-D', '{}/'.format(solved_dir),
#                                         '{0}'.format(self.fullname)], 
#                                         stdout=subprocess.PIPE) as proc :
#                     print(proc.stdout.read())
                
#                 if os.path.exists('{}/{}.new'.format(self.solved_dir, self.fullname_el[-1][:-5])):
#                     print("{} is solved successful ...".format(self.fullname_el[-1]))
                    
#                     shutil.move('{}/{}.new'.format(self.solved_dir, self.fullname_el[-1][:-5]), \
#                                 '{}/{}'.format(self.solved_dir, self.fullname_el[-1]))
#                     print("{} is renamed to fits ...".format(self.fullname_el[-1]))
                
#                 else : 
#                     print("{} solving fail ...".format(self.fullname_el[-1]))
                
#             except Exception as err :
#                     print('{1} ::: {2} with {0} ...'\
#                             .format(self.fullname, datetime.now(), err))
#########################################


# #%%
# #########################################
# #single  Astrometry Solver1
# #########################################
# class AstrometrySolver1():
#     def __init__(self, fullname, solved_dir):
        
#         """
#         Parameters
#         ----------
#         fullname : Path-like
#             The fullname of input file...

#         solved dir: string
#             The directory where the output file              
#         """

#         self.fullname = fullname
#         self.solved_dir = solved_dir

#         print("Starting... \n{}".format(self.fullname))
#         self.fullname_el = self.fullname.split("/")
#         self.filename_el = self.fullname_el[-1].split("_")

#         print("self.solved_dir:", self.solved_dir)
#         print('{}/{}'.format(self.solved_dir, self.fullname_el[-1]))

#         if os.path.exists('{}/{}'.format(self.solved_dir, self.fullname_el[-1])):
#             print("{} is already solved ...".format(self.fullname_el[-1]))
        
#         else: 

#             try : 
#                 # solve command.
#                 # solve-field fullname.fit -O --cpulimit 120 --nsigma 15 -u app -L 1.2 -U 1.3 -N new_filename.fits -p --no-plots -D output_directory {0}
#                 with subprocess.Popen(['solve-field', 
#                                         '-O', #--overwrite: overwrite output files if they already exist
#                                         #'--scale-low', '0.1', '--scale-high', '0.40', #pixel scale
#                                         '-g', #--guess-scale: try to guess the image scale from the FITS headers
#                                         '--cpulimit', '120',  #will make it give up after 30 seconds.
#                                         '--nsigma', '15',
#                                         #'--downsample', '4',
#                                         '-u', 'app', #'--scale-units', 'arcsecperpix', #pixel scale
#                                         '-L', '1.2', '-U', '1.3',
#                                         #'-N',  '{}'.format(self.fullname[-1]), #--new-fits <filename>: output filename of the new FITS file containingthe WCS header; "none" to not create this file
#                                         '-p', 
#                                         '--no-plots',#: don't create any plots of the results
#                                         '-D', '{}/'.format(solved_dir),
#                                         '{0}'.format(self.fullname)], 
#                                         stdout=subprocess.PIPE) as proc :
#                     print(proc.stdout.read())
                
#                 if os.path.exists('{}/{}.new'.format(self.solved_dir, self.fullname_el[-1][:-5])):
#                     print("{} is solved successful ...".format(self.fullname_el[-1]))
                    
#                     shutil.move('{}/{}.new'.format(self.solved_dir, self.fullname_el[-1][:-5]), \
#                                 '{}/{}'.format(self.solved_dir, self.fullname_el[-1]))
#                     print("{} is renamed to fits ...".format(self.fullname_el[-1]))
                
#                 else : 
#                     print("{} solving fail ...".format(self.fullname_el[-1]))
                    
#             except Exception as err :
#                     print('{1} ::: {2} with {0} ...'\
#                             .format(self.fullname, datetime.now(), err))
# #########################################
#%%
# #########################################
# #single ASTAPSolver
# #########################################
# class ASTAPSolver():
#     def __init__(self, fullname, solved_dir):
#         self.fullname = fullname
#         self.solved_dir = solved_dir
#         """
#         Parameters
#         ----------
#         fullname : string
#             The fullname of input file...

#         solved dir: string
#             The directory where the output file              
#         """

#         print("Starting... \n{}".format(self.fullname))
#         self.fullname_el = self.fullname.split("/")
#         self.filename_el = self.fullname_el[-1].split("_")

#         print("self.solved_dir:", self.solved_dir)
#         print('{}/{}'.format(self.solved_dir, self.fullname_el[-1]))

#         if os.path.exists('{}/{}'.format(self.solved_dir, self.fullname_el[-1])):
#             print("{} is already solved ...".format(self.fullname_el[-1]))
        
#         else : 

#             try:
#                 # solve command.
#                 # astap -f fullname.fit -o output_file.fits -wcs -analyse2 -update
#                 #astap -f ../RnE_2022/KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin/KLEOPATRA_Light_v_2022-11-04-11-48-17_160sec_RiLA600_STX-16803_-20C_2bin.fit -o ../RnE_2022/KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin/solved1/11.fit -update
#                 with subprocess.Popen(['astap', 
#                             '-f', '{0}'.format(self.fullname), 
#                             '-o', 
#                             '{}/{}'.format(self.solved_dir, 
#                                     self.fullname_el[-1]), 
#                             '-wcs',
#                             '-analyse2',
#                             '-update',],
#                                 stdout=subprocess.PIPE) as proc :
#                     print(proc.stdout.read())
                
#             except Exception as err :
#                     print('{1} ::: {2} with {0} ...'\
#                             .format(self.fullname, datetime.now(), err))
               

#########################################

#########################################
#single ASTAPSolver
#########################################
def ASTAPSolver(fullname, solved_dir):

    """
    Parameters
    ----------
    fullname : string
        The fullname of input file...

    solved dir: string
        The directory where the output file              
    """

    print("Starting... \n{}".format(fullname))
    fullname_el = fullname.split("/")
    filename_el = fullname_el[-1].split("_")

    print("solved_dir:", solved_dir)
    print('{}/{}'.format(solved_dir, fullname_el[-1]))

    if os.path.exists('{}/{}'.format(solved_dir, fullname_el[-1])):
        print("{} is already solved ...".format(fullname_el[-1]))
    
    else : 

        try:
            # solve command.
            # astap -f fullname.fit -o output_file.fits -wcs -analyse2 -update
            #astap -f ../RnE_2022/KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin/KLEOPATRA_Light_v_2022-11-04-11-48-17_160sec_RiLA600_STX-16803_-20C_2bin.fit -o ../RnE_2022/KLEOPATRA_Light_-_2022-11-04_-_RiLA600_STX-16803_-_2bin/solved1/11.fit -update
            with subprocess.Popen(['astap', 
                        '-f', '{0}'.format(fullname), 
                        '-o', 
                        '{}/{}'.format(solved_dir, 
                                fullname_el[-1]), 
                        '-wcs',
                        '-analyse2',
                        '-update',],
                            stdout=subprocess.PIPE) as proc :
                print(proc.stdout.read())
            
        except Exception as err :
            print('{1} ::: {2} with {0} ...'\
                        .format(fullname, datetime.now(), err))

#%%        
# =============================================================================
# for checking time
# =============================================================================
cht_start_time = datetime.now()
def print_working_time(cht_start_time):
    working_time = (datetime.now() - cht_start_time) #total days for downloading
    return print('working time ::: %s' % (working_time))


# =============================================================================
#     
# =============================================================================
#%%
def get_new_filename(fullname, **kargs):
    print('Starting get_new_filename ...\n{0}'.format(fullname))
    from astropy.io import fits
    hdul = fits.open(fullname)
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
    

    if hdul[0].header['NAXIS1'] == 1024 \
        and hdul[0].header['NAXIS2'] == 1024 :
        for binning in ['XBINNING', 'YBINNING'] :
            if not binning in hdul[0].header :
                with fits.open('{0}'.format(fullname), mode="append") as hdul1 :
                    hdul1[0].header.append(binning, '3', 'Binning factor in ')
                    hdul1.flush()
            elif hdul[0].header[binning] is None :
                with fits.open('{0}'.format(fullname), mode="update") as hdul1 :
                    hdul1[0].header[binning] = '3'
                    hdul1.flush()
            hdul[0].header[binning] = '3'
        hdul[0].header['INSTRUME'] = 'STX-16803' 
        hdul[0].header['TELESCOPE'] = 'RiLA600' 
        hdul[0].header['OPTIC'] = 'RiLA600' 


    if not 'INSTRUME' in hdul[0].header : 
        if hdul[0].header['CCDNAME'].lower() :     
            instrument = hdul[0].header['CCDNAME']
        else:
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
    
    if 'DATE-OBS' in hdul[0].header and 'TIME-OBS' in hdul[0].header : 
        if len(hdul[0].header['DATE-OBS']) == 10 :
            with fits.open('{0}'.format(fullname), mode="append") as hdul1 :        
                hdul[0].header['DATE-OBS'] += 'T{}'.format(hdul[0].header['TIME-OBS'])
                hdul1.flush()
                
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
        or hdul[0].header['IMAGETYP'][:1].lower() == 'z':
        image_type = 'Bias'
        filter_name = '-'
        object_name = '-'
        optic = '-'
    elif hdul[0].header['IMAGETYP'][:1].lower() == 'd':
        image_type = 'Dark'
        filter_name = '-'
        object_name = '-'
        optic = '-'
    elif hdul[0].header['IMAGETYP'][:1].lower() == 'f':
        image_type = 'Flat'
        object_name = '-'
    elif hdul[0].header['IMAGETYP'][:1].lower() == 'l' :
        image_type = 'Light'
        #filter_name = hdul[0].header['FILTER'] 
        object_name = hdul[0].header['OBJECT']   
    elif hdul[0].header['IMAGETYP'][0:1] == 'o' :
        image_type = 'Light'
        filter_name = hdul[0].header['FILTER'] 
        object_name = hdul[0].header['OBJECT']
        
    
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
        
    new_filename = '{0}_{1}_{2}_{3}_{4}sec_{5}_{6}_{7}C_{8}bin.fit'\
        .format(object_name.upper(),
        image_type,
        filter_name,
        obs_date,
        esposure,
        optic,
        instrument,
        ccd_temp_el[0],
        xbin)
    hdul.close()
    return new_filename


def get_new_foldername_from_filename(filename):
    #print('Starting get_new_foldername ...\n{0}'.format(filename))   
    filename_el = filename[:-4].split("_")
    #print("filename_el: ", filename_el)
    timez = 9
    obs_UT = datetime.strptime(filename_el[3], '%Y-%m-%d-%H-%M-%S')
    obs_LST = obs_UT + timedelta(hours = timez)
    if obs_LST.hour < 12 :
        obs_LST = obs_LST - timedelta(days = 1)
    filename_el[3] = obs_LST.strftime('%Y-%m-%d-%H-%M-%S')
    if filename_el[1].lower() == 'bias':
        new_foldername = '{6}_{8}/Cal/-_{1}_-_{3}_-_-_{6}_-_{8}/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    elif filename_el[1].lower() == 'dark' :
        new_foldername = '{6}_{8}/Cal/-_{1}_-_{3}_{4}_-_{6}_-_{8}/'\
        .format(filename_el[0],
        filename_el[1],
        filename_el[2],
        filename_el[3][:10],
        filename_el[4],
        filename_el[5],
        filename_el[6],
        filename_el[7],
        filename_el[8])
    elif filename_el[1].lower() == 'flat' :
        new_foldername = '{6}_{8}/Cal_{5}/-_{1}_-_{3}_-_{5}_{6}_-_{8}/'\
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
    #write_log(log_file, 
    #            '{1} ::: \nNew foldername is {0} ...'\
    #            .format(new_foldername, datetime.now()))    
    return new_foldername


def get_new_foldername(filename):
    #log_file = 'get_new_foldername.log'
    print('Starting get_new_foldername ...\n{0}'.format(filename))
    
    filename_el1 = filename.split("bin")
    filename_el = filename_el1[0].split("_")
    
    if filename_el[1].lower() == 'bias':
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
    elif filename_el[1].lower() == 'dark' :
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
    elif filename_el[1].lower() == 'flat' :
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
    #write_log(log_file, 
    #            '{1} ::: \nNew foldername is {0} ...'\
    #            .format(new_foldername, datetime.now()))    
    return new_foldername

#%%
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

#%%
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


#%%                                
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


#%%
def print_subworking_time(sub_start_time):
    from datetime import datetime
    working_time = (datetime.now() - cht_start_time) #total days for downloading
    return print('working time ::: %s' % (working_time))

#%%
def subp_solve_field(fullname, save_dir_name, sub_start_time): 
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

#%%
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
    criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT,