---
marp: Flase
paginate: true
---

# astro_Python
anaconda environment

# Windows
## 가상환경 리스트
conda env list

## 가상환경 만들기
conda create -n astro_Python_env

## activate 가상환경 시작
conda activate astro_Python_env

## deactivate 가상환경 종료
deactivate

## install module
conda install astropy pymysql seaborn
conda install -c conda-forge ccdproc xarray
pip install opencv-python photutils astroquery
pip install ysfitsutilpy ysphotutilpy chardet


# 가상환경 내보내기 (export)
conda env export > astro_Python_env.yaml

# .yaml 파일로 새로운 가상환경 만들기
conda env create -f astro_Python_env.yaml

## .yaml 파일로 가상환경 업데이트(activate 되어있을 때)
conda env update --file astro_Python_env.yaml

## .yaml 파일로 가상환경 업데이트(deactivate 되어있을 때)
conda env update --astro_Python_env envname --file astro_Python_env.yaml

## 가상환경 제거하기
conda env remove -n astro_Python_env
