#%%
# 
import os
from glob import glob
from astropy.nddata import CCDData
from pathlib import Path

import ysfitsutilpy as yfu
import ysphotutilpy as ypu
import ysvisutilpy as yvu

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
BASEDIR = astro_utilities.base_dir
OBSRAWDIR = astro_utilities.CCD_obs_dir
BASEDIRs = sorted(Python_utilities.getFullnameListOfsubDir(BASEDIR))
print ("BASEDIRs: {}".format(BASEDIRs))
print ("len(BASEDIRs): {}".format(len(BASEDIRs)))

#%%
for BASEDIR in BASEDIRs[:] :
    print ("Starting...\n{}".format(BASEDIR))

    BASEDIR = Path(BASEDIR)

    MASTERDIR = BASEDIR / astro_utilities.master_dir
    REDUCEDDIR = BASEDIR / astro_utilities.reduced_dir
    REDUCEDDIR2 = BASEDIR / astro_utilities.reduced_dir2

    if not REDUCEDDIR2.exists():
        os.makedirs("{}".format(str(REDUCEDDIR2)))
        print("{} is created...".format(str(REDUCEDDIR2)))
        
    #%%
    summary = yfu.make_summary(REDUCEDDIR/"*.fit*")
    if summary.empty:
        print("The dataframe(summary) is empty")
        pass
    else:
        print(summary)
        print("len(summary):", len(summary))

        # %%
        for filt in ["b", "v", "r"]:
            summary_filt = summary.loc[summary["FILTER"] == filt].copy()
            
            if summary_filt.empty:
                print("The dataframe(summary_filt) is empty")
                pass
            else:
                print("len(summary_filt):", len(summary_filt))
                print("summary_filt:", summary_filt)

                ccd = yfu.imcombine(
                    summary_filt["file"].tolist(), 
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
