"""
# -*- coding: utf-8 -*-

Created on Thu Nov 22 01:00:19 2018
@author: guitar79@naver.com

"""
#%%
from pathlib import Path
import shutil
import ysfitsutilpy as yfu

import _Python_utilities
import _astro_utilities
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
# Set directory variables.
#######################################################

# BASEDIR = Path(r"r:\OBS_data")   # for windows
BASEDIR = Path("/mnt/Rdata/ASTRO_data")  # for ubuntu
# BASEDIR = Path("/Volumes/OBS_data")  # for mac OS
 
DOINGDIR = BASEDIR/ _astro_utilities.CCD_NEW_dir

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfallsubDirs(DOINGDIR))
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))
#######################################################            

#%%
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    # print("DOINGDIR", DOINGDIR)
    # print(f"Starting: {str(DOINGDIR.parts[-1])}")
    try: 
        summary = yfu.make_summary(DOINGDIR/"*.fit*",)
        if summary is not None : 
            # print("summary: ", summary)
            print("len(summary)", len(summary))

            for _, row in summary.iterrows():
                # 파일명 출력
                # print (row["file"])
                fpath = Path(row["file"])
                try:
                    hdul = _astro_utilities.KevinFitsUpdater(fpath,
                                                    # checkKEYs = ["OBJECT", "TELESCOP", "OPTIC", "CCDNAME", 'FILTER',
                                                    #             #"GAIN", "EGAIN", "RDNOISE", 
                                                    #             "PIXSCALE", "FOCALLEN", "APATURE", "CCD-TEMP",
                                                    #             'XPIXSZ', 'YPIXSZ',
                                                    #             "XBINNING", "YBINNING", "FLIPSTAT", "EXPTIME", "EXPOSURE"],
                                                    # imgtype_update=True,  
                                                    # fil_update=False,
                                                    )
                    print("hdul: ", hdul)

                except Exception as err :
                    print("X"*60)
                    print(err)
                    pass
                 
        shutil.move(str(DOINGDIR), str(BASEDIR/_astro_utilities.CCD_NEWUP_dir / DOINGDIR.stem))
        print(str(DOINGDIR), str(BASEDIR/_astro_utilities.CCD_NEWUP_dir / DOINGDIR.stem))
    except Exception as err :
        print("X"*60)
        print(err)
        pass