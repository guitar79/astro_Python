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
BASEDIR = Path("/mnt/Rdata/OBS_data") 
PROJECDIR = Path("/mnt/Rdata/OBS_data/2024-EXO")
TODODIR = PROJECDIR / "_-_-_2024-05_-_GSON300_STF-8300M_-_1bin"
TODODIR = PROJECDIR / "_-_-_2024-06_-_GSON300_STF-8300M_-_1bin"

PROJECDIR = Path("/mnt/Rdata/OBS_data/2022-Asteroid")
TODODIR = PROJECDIR / "GSON300_STF-8300M_-_1bin"
TODODIR = PROJECDIR / "RiLA600_STX-16803_-_1bin"
TODODIR = PROJECDIR / "RiLA600_STX-16803_-_2bin"

# PROJECDIR = Path("/mnt/Rdata/OBS_data/2023-Asteroid")
# TODODIR = PROJECDIR / "GSON300_STF-8300M_-_1bin"
# TODODIR = PROJECDIR / "RiLA600_STX-16803_-_1bin"
# TODODIR = PROJECDIR / "RiLA600_STX-16803_-_2bin"

PROJECDIR = Path("/mnt/Rdata/OBS_data/2016-Variable")
TODODIR = PROJECDIR / "-_-_-_2016-_-_RiLA600_STX-16803_-_2bin"

# PROJECDIR = Path("/mnt/Rdata/OBS_data/2017-Variable")
# TODODIR = PROJECDIR / "-_-_-_2017-_-_RiLA600_STX-16803_-_2bin"

DOINGDIRs = sorted(_Python_utilities.getFullnameListOfsubDirs(TODODIR))
print ("DOINGDIRs: ", format(DOINGDIRs))
print ("len(DOINGDIRs): ", format(len(DOINGDIRs)))

CALDIR = [x for x in DOINGDIRs if "CAL-BDF" in str(x)]
MASTERDIR = Path(CALDIR[0]) / _astro_utilities.master_dir
if not MASTERDIR.exists():
    os.makedirs("{}".format(str(MASTERDIR)))
    print("{} is created...".format(str(MASTERDIR)))

print ("MASTERDIR: ", format(MASTERDIR))

DOINGDIRs = sorted([x for x in DOINGDIRs if "_LIGHT_" in str(x)])
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
    sMASTERDIR = DOINGDIR / _astro_utilities.master_dir
    REDUCEDDIR = DOINGDIR / _astro_utilities.reduced_dir
    REDUC_nightsky = DOINGDIR / _astro_utilities.reduced_nightsky_dir

    # if not MASTERDIR.exists():
    # shutil.copytree(MASTERDIR, MASTERDIR, dirs_exist_ok=True)

    if not sMASTERDIR.exists():
        os.makedirs(str(sMASTERDIR))
        print("{} is created...".format(str(sMASTERDIR)))

    if not REDUC_nightsky.exists():
        os.makedirs("{}".format(str(REDUC_nightsky)))
        print("{} is created...".format(str(REDUC_nightsky)))
    
    summary = yfu.make_summary(DOINGDIR/"*.fit*")
    if summary is not None :
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
                if (not (sMASTERDIR / f"nightskyflat-{filt}.fits").exists()) :
                    print("len(summary_light_filt):", len(summary_light_filt))
                    print("summary_light_filt:", summary_light_filt)
                    
                    File_Num = 80
                    if len(summary_light_filt["file"]) > File_Num :
                        combine_lst = summary_light_filt["file"].tolist()[:File_Num]
                    else : 
                        combine_lst = summary_light_filt["file"].tolist()
                    try : 
                        ccd = yfu.imcombine(
                                            combine_lst, 
                                            combine="med",
                                            scale="avg", 
                                            scale_to_0th=False, 
                                            reject="sc", 
                                            sigma=2.5,
                                            verbose=True,
                                            memlimit = 2.e+11,
                                            )
                        ccd.write(sMASTERDIR / f"nightskyflat-{filt}.fits", overwrite=True)
                        print (f"Create Create nightskyflat-{filt}.fits +++...")
                    except :
                        ccd = yfu.imcombine(
                                            combine_lst, 
                                            combine="med",
                                            scale="avg", 
                                            scale_to_0th=False, 
                                            reject="sc", 
                                            # sigma=2.5,
                                            verbose=True,
                                            memlimit = 2.e+11,
                                            )
                        ccd.write(sMASTERDIR / f"nightskyflat-{filt}.fits", overwrite=True)
                        print (f"Create Create nightskyflat-{filt}.fits +++...")

        for _, row in summary_light.iterrows():
            try: 
                fpath = Path(row["file"])
                filt = row["FILTER"]
                ccd = yfu.ccdred(
                                fpath, 
                                mflatpath=str(sMASTERDIR / f"nightskyflat-{filt}.fits"),
                                output=REDUC_nightsky/fpath.name
                               )
            except FileNotFoundError: 
                _Python_utilities.write_log(err_log_file, "FileNotFoundError")