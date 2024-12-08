#!/bin/sh
# PROJECT="C3-EXO";
# YRMN="2024-08" ;
# OPTIC="GSON300";
# CCDNAME="STF-8300M";
# BIN="1bin" ; 

# # for BDF in "BIAS" "DARK" ;
# # do 
# #     find /mnt/Rdata/ASTRO_data/A3_CCD_obs_raw/${CCDNAME}_${BIN}/Cal -type f -name "-_${BDF}*${YRMN}*.fit*" -exec cp -f "{}" "/mnt/Rdata/ASTRO_data/${PROJECT}/-_-_-_${YRMN}_-_${OPTIC}_${CCDNAME}_-_${BIN}/-_CAL-BDF_-_${YRMN}_-_${OPTIC}_${CCDNAME}_-_${BIN}/" \;
# # done ;
# # find /mnt/Rdata/ASTRO_data/A3_CCD_obs_raw/${CCDNAME}_${BIN}/Cal_${OPTIC} -type f -name "-_FLAT*${YRMN}*.fit*" -exec cp -f "{}" "/mnt/Rdata/ASTRO_data/${PROJECT}/-_-_-_${YRMN}_-_${OPTIC}_${CCDNAME}_-_${BIN}/-_CAL-BDF_-_${YRMN}_-_${OPTIC}_${CCDNAME}_-_${BIN}/" \;

# find /mnt/Rdata/ASTRO_data/A3_CCD_obs_raw/${CCDNAME}_${BIN}/Cal_${OPTIC} -type f -name "-_FLAT*${YRMN}*.fit*" -exec cp -f "{}" "/mnt/Rdata/ASTRO_data/${PROJECT}/-_-_-_2024-09_-_${OPTIC}_${CCDNAME}_-_${BIN}/-_CAL-BDF_-_2024-09_-_${OPTIC}_${CCDNAME}_-_${BIN}/" \;

# rsync -avuz --progress --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/ASTRO_data/C1-Variable' '/mnt/Rdata/ASTRO_data/' 
# rsync -avuz --progress --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/ASTRO_data/C2-Asteroid' '/mnt/Rdata/ASTRO_data/' 
# rsync -avuz --progress --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/ASTRO_data/C3-EXO' '/mnt/Rdata/ASTRO_data/'
# rsync -avuz --progress --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/ASTRO_data/C4-Spectra' '/mnt/Rdata/ASTRO_data/'
# rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/ASTRO_data/2024-OA' '/mnt/Rdata/ASTRO_data/'


# rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/ASTRO_data/C1-Variable' '/mnt/Rdata/ASTRO_data/' 
# rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/ASTRO_data/C2-Asteroid' '/mnt/Rdata/ASTRO_data/' 
# rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/ASTRO_data/C3-EXO' '/mnt/Rdata/ASTRO_data/'
# rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/ASTRO_data/C4-Spectra' '/mnt/Rdata/ASTRO_data/'
# rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/ASTRO_data/2024-OA' '/mnt/Rdata/ASTRO_data/'

# BASEDIR = "/mnt/Rdata/ASTRO_data"
# PROJECTDIR = "C1-Variable

# for dir_list in "-_-_-_2016-_-_RiLA600_STX-16803_-_2bin"
# do
#         rsync -avuz --progress --rsh='ssh -p2022' "guitar79@parksparks.iptime.org:/volume1/Rdata/ASTRO_data/C1-Variable/$dir_list" "/mnt/Rdata/ASTRO_data/C1-Variable/"
# done

# find /mnt/Rdata/ASTRO_data/A3_CCD_obs_raw/ -type f -name '*.wcs' -delete
# find /mnt/Rdata/ASTRO_data/A3_CCD_obs_raw/ -type f -name '*.rdls' -delete
# find /mnt/Rdata/ASTRO_data/A3_CCD_obs_raw/ -type f -name '*.asy' -delete
# find /mnt/Rdata/ASTRO_data/A3_CCD_obs_raw/ -type f -name '*.corr' -delete
# find /mnt/Rdata/ASTRO_data/A3_CCD_obs_raw/ -type f -name '*.match' -delete
# find /mnt/Rdata/ASTRO_data/A3_CCD_obs_raw/ -type f -name '*.solved' -delete

# find /mnt/Rdata/ASTRO_data/A3_CCD_obs_raw/ -type f \( -name '*.wcs' -o -name '*.rdls' -o -name '*.asy' -o -name '*.corr' -o -name '*.match' -o -name '*.soled' \) -delete

# do 
#     # rsync -avuz --progress /mnt/Rdata/ASTRO_data/2024-OA/_측광과제-교사참고/1반/${NNAME}/ /mnt/Rdata/ASTRO_data/2024-OA/_측광과제-학생수행/1반/${NNAME}/ ;
#     rsync -avuz --progress /mnt/Rdata/ASTRO_data/2024-OA/_측광과제-교사참고/2반/${NNAME}/ /mnt/Rdata/ASTRO_data/2024-OA/_측광과제-학생수행/2반/${NNAME}/ ;
# done ;

# rsync -avuz --progress /mnt/Rdata/ASTRO_data/2024-OA/_측광과제-교사참고/2반/ /mnt/Rdata/ASTRO_data/2024-OA/_측광과제-학생수행/2반/ ;


