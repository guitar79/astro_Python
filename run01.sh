#!/bin/sh
# PROJECT="03-EXO";
# YRMN="2024-08" ;
# OPTIC="GSON300";
# CCDNAME="STF-8300M";
# BIN="1bin" ; 

# # for BDF in "BIAS" "DARK" ;
# # do 
# #     find /mnt/Rdata/OBS_data/CCD_obs_raw/${CCDNAME}_${BIN}/Cal -type f -name "-_${BDF}*${YRMN}*.fit*" -exec cp -f "{}" "/mnt/Rdata/OBS_data/${PROJECT}/-_-_-_${YRMN}_-_${OPTIC}_${CCDNAME}_-_${BIN}/-_CAL-BDF_-_${YRMN}_-_${OPTIC}_${CCDNAME}_-_${BIN}/" \;
# # done ;
# # find /mnt/Rdata/OBS_data/CCD_obs_raw/${CCDNAME}_${BIN}/Cal_${OPTIC} -type f -name "-_FLAT*${YRMN}*.fit*" -exec cp -f "{}" "/mnt/Rdata/OBS_data/${PROJECT}/-_-_-_${YRMN}_-_${OPTIC}_${CCDNAME}_-_${BIN}/-_CAL-BDF_-_${YRMN}_-_${OPTIC}_${CCDNAME}_-_${BIN}/" \;

# find /mnt/Rdata/OBS_data/CCD_obs_raw/${CCDNAME}_${BIN}/Cal_${OPTIC} -type f -name "-_FLAT*${YRMN}*.fit*" -exec cp -f "{}" "/mnt/Rdata/OBS_data/${PROJECT}/-_-_-_2024-09_-_${OPTIC}_${CCDNAME}_-_${BIN}/-_CAL-BDF_-_2024-09_-_${OPTIC}_${CCDNAME}_-_${BIN}/" \;

# rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/01-Variable' '/mnt/Rdata/OBS_data/' 
# rsync -avuz --progress --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/02-Asteroid' '/mnt/Rdata/OBS_data/' 
# rsync -avuz --progress --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/03-EXO' '/mnt/Rdata/OBS_data/'
# rsync -avuz --progress --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/04-Spectra' '/mnt/Rdata/OBS_data/'
# rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/2024-OA' '/mnt/Rdata/OBS_data/'

# rsync -avuz --progress --rsh='ssh -p2022' '/mnt/Rdata/OBS_data/01-Variable' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/' 
# rsync -avuz --progress --rsh='ssh -p2022' '/mnt/Rdata/OBS_data/02-Asteroid' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/' 
# rsync -avuz --progress --delete --rsh='ssh -p2022' '/mnt/Rdata/OBS_data/03-EXO' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/' 
# rsync -avuz --progress --rsh='ssh -p2022' '/mnt/Rdata/OBS_data/04-Spectra' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/' 
# rsync -avuz --progress --delete --rsh='ssh -p2022' '/mnt/Rdata/OBS_data/2024-OA' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/' 

# rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/01-Variable' '/mnt/Rdata/OBS_data/' 
# rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/02-Asteroid' '/mnt/Rdata/OBS_data/' 
# rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/03-EXO' '/mnt/Rdata/OBS_data/'
# rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/04-Spectra' '/mnt/Rdata/OBS_data/'
rsync -avuz --progress --delete --rsh='ssh -p2022' 'guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/2024-OA' '/mnt/Rdata/OBS_data/'


# BASEDIR = "/mnt/Rdata/OBS_data"
# PROJECTDIR = "01-Variable

# for dir_list in "-_-_-_2016-_-_RiLA600_STX-16803_-_2bin"
# do
#         rsync -avuz --progress --rsh='ssh -p2022' "guitar79@parksparks.iptime.org:/volume1/Rdata/OBS_data/01-Variable/$dir_list" "/mnt/Rdata/OBS_data/01-Variable/"
# done

# find /mnt/Rdata/OBS_data/CCD_obs_raw/ -type f -name '*.wcs' -delete
# find /mnt/Rdata/OBS_data/CCD_obs_raw/ -type f -name '*.rdls' -delete
# find /mnt/Rdata/OBS_data/CCD_obs_raw/ -type f -name '*.asy' -delete
# find /mnt/Rdata/OBS_data/CCD_obs_raw/ -type f -name '*.corr' -delete
# find /mnt/Rdata/OBS_data/CCD_obs_raw/ -type f -name '*.match' -delete
# find /mnt/Rdata/OBS_data/CCD_obs_raw/ -type f -name '*.solved' -delete

# find /mnt/Rdata/OBS_data/CCD_obs_raw/ -type f \( -name '*.wcs' -o -name '*.rdls' -o -name '*.asy' -o -name '*.corr' -o -name '*.match' -o -name '*.soled' \) -delete

# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/23116최현준
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/23067신재헌
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/23097임윤준
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/23054박하람
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/22029박주영
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/22115최준서
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/22032박지환
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/22028박성환
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/22106조형준
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/22027박성민
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/22030박지원
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/22098정이찬
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/22024문준원
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/22072이재우
# mkdir /mnt/Rdata/OBS_data/2024-OA/1반/22080이혁준

# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/23074오서준
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22018김한준
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22004권민우
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22095정은재
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22022노현우
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22073이재욱
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/21100정영우
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22045양현서
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22050오태원
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22125홍은찬
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22093정우현
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22110최석원
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22035박홍준
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22053용승주
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22107지민기
# mkdir /mnt/Rdata/OBS_data/2024-OA/2반/22118최현진

# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/23075오승민
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/23027김재우
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/23108조형석
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/23069안선우
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/22005권순민
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/22008김도현
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/22039손희원
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/22088장태훈
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/22012김수아
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/22034박현수
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/22082임비건
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/22048오은총
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/22069이은우
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/22103조연우
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/22121함석규
# mkdir /mnt/Rdata/OBS_data/2024-OA/3반/22108차무겸