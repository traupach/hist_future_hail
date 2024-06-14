#!/bin/bash
# Calculate convective parameters from basic_params* files.

#PBS -q normal
#PBS -P li18
#PBS -l storage=gdata/up6+gdata/hh5+gdata/w42
#PBS -l ncpus=28
#PBS -l walltime=06:00:00
#PBS -l mem=192GB
#PBS -j oe
#PBS -W umask=0022
#PBS -l wd
#PBS -l jobfs=1GB
#PBS -N process_WRF_conv_job

module use /g/data3/hh5/public/modules
module load conda/analysis3-22.04

time python3 ~/git/hist_future_hail/scripts/process_WRF_conv.py
