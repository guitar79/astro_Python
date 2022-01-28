# astro_Python
anaconda environment

# Windows
## 가상환경 리스트
conda env list

## 가상환경 만들기
conda create -n astro_Python_win_env

## activate 가상환경 시작
activate astro_Python_win_env

## deactivate 가상환경 종료
deactivate

## install module
conda install astropy pymysql spyder
conda install -c conda-forge ccdproc
pip install opencv-python


## 가상환경 내보내기 (export)
conda env export > astro_Python_win_env.yaml

## .yaml 파일로 새로운 가상환경 만들기
conda env create -f astro_Python_win_env.yaml

## 가상환경 제거하기
conda env remove -n astro_Python_win_env


# ubuntu

## 가상환경 리스트 출력
conda env list

## 가상환경 만들기 
conda create -n astro_Python_ubuntu_env

## activate 가상환경 시작
conda activate astro_Python_ubuntu_env

# deactivate 가상환경 종료
conda deactivate

# install module
conda install astropy pymysql spyder
conda install -c conda-forge ccdproc
pip install opencv-python

# 가상환경 내보내기 (export)
conda env export > astro_Python_ubuntu_env.yaml

# .yaml 파일로 새로운 가상환경 만들기
conda env create -f astro_Python_ubuntu_env.yaml

# 가상환경 제거하기
conda env remove -n astro_Python_ubuntu_env