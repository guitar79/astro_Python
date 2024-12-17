
# astro_Python
anaconda environment

## 가상환경 리스트
conda env list

## 가상환경 만들기
conda create -n astro_Python_env

## activate 가상환경 시작
conda activate astro_Python_env

## deactivate 가상환경 종료
conda deactivate

## install module
conda install -c conda-forge jupyter-book
conda install astropy pymysql
pip install opencv-python astroquery sep
pip install ysfitsutilpy ysphotutilpy chardet ccdproc xarray
pip install photutils==1.12

## 가상환경 내보내기 (export)
conda env export > astro_Python_env.yaml

## .yaml 파일로 새로운 가상환경 만들기
conda env create -f astro_Python_env.yaml

## .yaml 파일로 가상환경 업데이트(activate 되어있을 때)
conda env update --file astro_Python_env_241208.yaml

## .yaml 파일로 가상환경 업데이트(deactivate 되어있을 때)
conda env update --astro_Python_env --file astro_Python_env.yaml

## 가상환경 제거하기
conda remove -n astro_Python_env --all
