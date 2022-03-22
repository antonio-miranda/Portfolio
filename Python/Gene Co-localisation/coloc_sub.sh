#!/bin/bash

#PBS -N Coloc
#PBS -l select=1:ncpus=40:ompthreads=40:mem=480gb
#PBS -lwalltime=24:00:00
#PBS -o /rdsgpfs/general/user/aalmeid2/home/coloc.out

export PATH=$Home/rds/general/user/aalmeid2/home/.conda/envs/personal/bin/python:$PATH
cd $Home

module load anaconda3/personal
unset PYTHONPATH
source activate personal

python dcm_coloc_all2.py