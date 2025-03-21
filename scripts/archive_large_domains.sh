#!/bin/bash

#PBS -q normal
#PBS -P w42
#PBS -l storage=gdata/up6+gdata/hh5+scratch/up6+massdata/up6
#PBS -l walltime=10:00:00
#PBS -l mem=10GB       
#PBS -l ncpus=40
#PBS -j oe
#PBS -l wd
#PBS -W umask=0022
#PBS -N archive_job
#PBS -l jobfs=2GB

module load parallel

# Run from e.g. /g/data/up6/tr2908/hist_future_hail/WRF_v4.4/simulations/remote/
# And remember to increase ncpus for parallel archiving.

parallel tar cvf /scratch/up6/tr2908/d03_{= s/\\//_/g =}.tar {}/wrfout_d0[1,2,4]* --remove-files ::: */*/WRF
