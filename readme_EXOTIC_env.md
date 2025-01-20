## 가상환경 리스트
conda env list

## 가상환경 만들기
conda create -n EXOTIC_env python=3.10

## activate 가상환경 시작
conda activate EXOTIC_env

## deactivate 가상환경 종료
deactivate

## install module
pip install exotic

## 가상환경 내보내기 (export)
conda env export > EXOTIC_env.yaml

## .yaml 파일로 새로운 가상환경 만들기
conda env create -f EXOTIC_env.yaml

## .yaml 파일로 가상환경 업데이트(activate 되어있을 때)
conda env update --file EXOTIC_env.yaml

## .yaml 파일로 가상환경 업데이트(deactivate 되어있을 때)
conda env update --EXOTIC_env envname --file EXOTIC_env.yaml

## 가상환경 제거하기
conda env remove -n EXOTIC_env --all


https://github.com/rzellem/EXOTIC



