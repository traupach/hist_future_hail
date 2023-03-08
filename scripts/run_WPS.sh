#!/bin/bash
# Run WPS to prepare a WRF run.

#PBS -q normal
#PBS -P li18
#PBS -l storage=gdata/up6+gdata/hh5+gdata/rt52+gdata/zz93+gdata/sx70
#PBS -l ncpus=48
#PBS -l walltime=12:00:00
#PBS -l mem=190GB
#PBS -j oe
#PBS -W umask=0022
#PBS -l wd
#PBS -l jobfs=100GB
#PBS -N WPS_job

exec_dir=/g/data/up6/tr2908/hist_future_hail/xu_code/nc2wrf/

module use /g/data3/hh5/public/modules
module load conda/analysis3

# Process downloaded netcdf files from Xu et al., 2021 
# (https://www.nature.com/articles/s41597-021-01079-3) into FILE* files.
bash mk.inputlist.sh                     # Generate file list to process.
$exec_dir/GCM_gfortran.exe input.GCM2WRF # Process files.

# Run geogrid.
mpirun -np $PBS_NCPUS -report-bindings geogrid/geogrid.exe

# Run metgrid to interpolate input data.
mpirun -np $PBS_NCPUS -report-bindings metgrid/metgrid.exe

# Remove intermediate files.
rm FILE*
rm SST*
