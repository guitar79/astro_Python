#%%
# 
import os
from glob import glob
from astropy.nddata import CCDData
from pathlib import Path

import ysfitsutilpy as yfu

import Python_utilities
import astro_utilities


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
BASEDIR = "../RnE_2022/"
BASEDIR = "../RnE_2022/RiLA600_STX-16803_2bin/"
#BASEDIR = "../RnE_2022/GSON300_STF-8300M/"

#%%
BASEDIRs = sorted(Python_utilities.getFullnameListOfsubDir(BASEDIR))
print ("BASEDIRs: {}".format(BASEDIRs))
print ("len(BASEDIRs): {}".format(len(BASEDIRs)))

#%%
for BASEDIR in BASEDIRs[7:] :
    print ("Starting...\n{}".format(BASEDIR))

    BASEDIR = Path(BASEDIR)
    
    MASTERDIR = BASEDIR / astro_utilities.master_dir
    REDUCEDDIR = BASEDIR / astro_utilities.reduced_dir
    SOLVEDDIR = BASEDIR / astro_utilities.solved_dir
    RESULTDIR = BASEDIR / astro_utilities.DAOfinder_result_dir
    OBSRAWDIR = BASEDIR / astro_utilities.CCD_obs_dir
    REDUCEDDIR2 = BASEDIR / "reduced2"
    
    #REDUCEDDIR = BASEDIR / "reduced2"

    if not REDUCEDDIR.exists():
        os.makedirs("{}".format(str(REDUCEDDIR)))
        print("{} is created...".format(str(REDUCEDDIR)))
        
    #%%
    summary = yfu.make_summary(BASEDIR/"*.fit*")
    if summary.empty:
        pass
    else:
        print(summary)
        print("len(summary):", len(summary))
        print(summary["file"][0])

        # %%
        fpaths = list((REDUCEDDIR).glob("*.fits"))
        fpaths.sort()
        summary = yfu.make_summary(fpaths)
        #%%
        for filt in ["b", "v", "r"]:
            df = summary.loc[summary["FILTER"] == filt].copy()
            ccd = yfu.imcombine(
                df["file"].tolist(), 
                combine="med",
                scale="avg", 
                scale_to_0th=False, 
                #reject="sc", 
                #sigma=2.5,
                verbose=2
                )
            ccd.write(MASTERDIR / f"nightskyflat-{filt}.fits", overwrite=True)


        # %%
        for _, row in summary.iterrows():
            fpath = Path(row["file"])
            filt = row["FILTER"]
            ccd = yfu.ccdred(
                fpath, 
                mflatpath=str(MASTERDIR / f"nightskyflat-{filt}.fits"),
                output=REDUCEDDIR2/fpath.name
            )
