#!/bin/bash

#PBS -q normal
#PBS -P li18
#PBS -l storage=gdata/up6+gdata/hh5
#PBS -l walltime=06:00:00
#PBS -l mem=1900GB       
#PBS -l ncpus=480
#PBS -j oe
#PBS -l wd
#PBS -W umask=0022
#PBS -N real_job
#PBS -l jobfs=100GB

module load openmpi
ulimit -s unlimited
limit stacksize unlimited

echo 'Running in directory:' `pwd`
env > run_environment_real.txt

echo 'Running real.exe using $PBS_NCPUS mpi nodes...'
time mpirun -np $PBS_NCPUS -report-bindings ./real.exe 
