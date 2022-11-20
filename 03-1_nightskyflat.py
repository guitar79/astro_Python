#%%
# 
import ysfitsutilpy as yfu
from astropy.nddata import CCDData
from pathlib import Path

#######################################################
# read all files in base directory for processing
BASEDIR = "../RnE_2022/RiLA600_STX-16803_2bin/KLEOPATRA_Light_-_2022-11-02_-_RiLA600_STX-16803_-_2bin"

c_method = 'median'
master_dir = "master_files_ys"
reduced_dir = "reduced"

# %%
fpaths = list((Path(BASEDIR)/reduced_dir).glob("*.fits"))
fpaths.sort()
summary = yfu.make_summary(fpaths)
#%%
for filt in ["b"]:
    df = summary.loc[summary["FILTER"] == filt].copy()
    ccd = yfu.imcombine(
        df["file"].tolist(), 
        combine="med",
        scale="avg", 
        scale_to_0th=False, 
        # reject="sc", 
        # sigma=2.5,
        verbose=2
    )
    ccd.write(Path(BASEDIR)/f"nightskyflat-{filt}.fits", overwrite=True)


# %%
for _, row in summary.iterrows():
    fpath = Path(row["file"])
    filt = row["FILTER"]
    ccd = yfu.ccdred(
        fpath, 
        mflatpath=str(Path(BASEDIR)/f"nightskyflat-{filt}.fits"),
        output=Path(BASEDIR)/"reduced2"/fpath.name
    )


# %%
