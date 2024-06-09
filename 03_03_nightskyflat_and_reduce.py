#%%
import os
from glob import glob
from astropy.nddata import CCDData
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
# read all files in base directory for processing
BASEDIR = Path("/mnt/Rdata/OBS_data") 
DOINGDIR = Path(BASEDIR/ "2023-Asteroid" / "RiLA600_STX-16803_-_1bin")
# DOINGDIR = Path(BASEDIR/ "2023-Asteroid" / "RiLA600_STX-16803_-_2bin")
# DOINGDIR = Path(BASEDIR/ "2022-Asteroid" / "GSON300_STF-8300M_-_1bin")

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(DOINGDIR))
DOINGDIRs = sorted([x for x in DOINGDIR.iterdir() if x.is_dir()])
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

mas1 = [x for x in DOINGDIRs if "CAL-BDF" in str(x)]
mas1 = mas1[0]/ _astro_utilities.master_dir
print ("mas1: ", format(mas1))

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_".lower() in str(x).lower()])
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

# filter_str = '2023-12-'
# DOINGDIRs = [x for x in DOINGDIRs if filter_str in x]
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
for DOINGDIR in DOINGDIRs[:] :
    DOINGDIR = Path(DOINGDIR)
    print (DOINGDIR.parts[-1])
    
    fits_in_dir = sorted(list(DOINGDIR.glob('*.fit*')))
    #print("fits_in_dir", fits_in_dir)
    print("len(fits_in_dir)", len(fits_in_dir))

    if len(fits_in_dir) == 0 :
        print(f"There is no fits fils in {DOINGDIR.parts[-1]}")
        pass
    else :     
        MASTERDIR = DOINGDIR / _astro_utilities.master_dir
        REDUCEDDIR = DOINGDIR / _astro_utilities.reduced_dir
        REDUC_nightsky = DOINGDIR / _astro_utilities.REDUC_nightsky_dir

        # if not MASTERDIR.exists():
        shutil.copytree(mas1, MASTERDIR, dirs_exist_ok=True)

        if not REDUC_nightsky.exists():
            os.makedirs("{}".format(str(REDUC_nightsky)))
            print("{} is created...".format(str(REDUC_nightsky)))

        summary = yfu.make_summary(DOINGDIR/"*.fit*",
                                    # keywords = ["FILTER", 'SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                                    #     'EXTEND', 'BZERO', 'IMAGETYP', 'EXPOSURE', 'EXPTIME', 'DATE-LOC',
                                    #     'DATE-OBS', 'XBINNING', 'YBINNING', 'EGAIN', 'XPIXSZ', 'YPIXSZ',
                                    #     'INSTRUME', 'SET-TEMP', 'CCD-TEMP', 'TELESCOP', 'FOCALLEN', 'FOCRATIO',
                                    #     'OBJECT', 'OBJCTRA', 'OBJCTDEC', 'OBJCTROT', 'ROWORDER', 'EQUINOX',
                                    #     'SWCREATE', 'NOTES']
                                    )

        print("len(summary):", len(summary))
        print("summary:", summary)
        #print(summary["file"][0])   

        summary_light = summary.loc[summary["IMAGETYP"] == "LIGHT"].copy()
        summary_light = summary_light.reset_index(drop=True) 

        for filt in ["b", "v", "r", "L", "R", "G", "B"]:
        #for filt in ["V"]:
            summary_light_filt = summary_light.loc[summary_light["FILTER"] == filt].copy()
            
            if summary_light_filt.empty:
                print("The dataframe(summary_light_filt) is empty")
                pass
            else:
                try:
                    print("len(summary_light_filt):", len(summary_light_filt))
                    print("summary_light_filt:", summary_light_filt)

                    ccd = yfu.imcombine(
                        summary_light_filt["file"].tolist(), 
                        combine="med",
                        scale="avg", 
                        scale_to_0th=False, 
                        reject="sc", 
                        sigma=2.5,
                        verbose=True,
                        memlimit = 2.e+11,
                        )
                    ccd.write(MASTERDIR / f"nightskyflat-{filt}.fits", overwrite=True)
                except Exception as err: 
                    print ('Error messgae .......')
                    #_Python_utilities.write_log(err_log_file, err)

        for _, row in summary_light.iterrows():
            try: 
                fpath = Path(row["file"])
                filt = row["FILTER"]
                ccd = yfu.ccdred(
                    fpath, 
                    mflatpath=str(MASTERDIR / f"nightskyflat-{filt}.fits"),
                    output=REDUC_nightsky/fpath.name
                )
            except FileNotFoundError: 
                _Python_utilities.write_log(err_log_file, "FileNotFoundError")