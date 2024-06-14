#!/bin/bash

#PBS -q copyq
#PBS -P up6
#PBS -l storage=gdata/up6+gdata/hh5+scratch/up6+massdata/up6
#PBS -l walltime=10:00:00
#PBS -l mem=2GB       
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -W umask=0022
#PBS -N archive_job_ssp
#PBS -l jobfs=2GB

module load parallel

# Run from e.g. /g/data/up6/tr2908/hist_future_hail/WRF_v4.4/simulations
# And remember to increase ncpus for parallel archiving.
# parallel tar cvf /scratch/up6/tr2908/{= s/\\//_/g =}.tar {}/wrfout_d0[1,2,4]* ::: */*/WRF

mdss -P up6 get 'tr2908/hist_future_hail/WRF_v4.4/simulations/ssp*.tar'
