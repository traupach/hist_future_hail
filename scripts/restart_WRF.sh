#!/bin/bash

if [ ! -e runs ]; then
    mkdir runs
fi

if [ ! -e rsl.out.0000 ]; then
    echo 'no rsl.out.0000 file.'
    exit
fi

run_dir=runs/`date -ur rsl.out.0000 +%Y%m%d_%H%MZ`
if [ -e $run_dir ]; then
    echo Run directory ${run_dir} already exists.
    
else
    mkdir $run_dir
    cp namelist.input $run_dir
    mv namelist.output $run_dir
    mv rsl* $run_dir
    mv run_envir* $run_dir
    mv WRF_job* $run_dir
    mv real_job* $run_dir
    
    cd $run_dir
    tar cfz rsl_logs.tar.gz rsl* --remove-file
    cd ../../

    # Remove all but the most recent restart file (for each of 7 domains).
    for i in 1 2 3 4 5 6 7; do
	if ls wrfrst_d0${i}* 1> /dev/null 2> /dev/null; then
            if [ `ls wrfrst_d0${i}* | wc -l` > 1 ]; then
		rm `ls -t wrfrst_d0${i}* | awk 'NR>1'`
            fi
	fi
    done
    
    year=`ls wrfrst_d01* | cut -d_ -f3 | cut -d- -f1`
    month=`ls wrfrst_d01* | cut -d_ -f3 | cut -d- -f2`
    day=`ls wrfrst_d01* | cut -d_ -f3 | cut -d- -f3`
    
    sed -i s/start_year.*$/"start_year = ${year}, ${year}, ${year}, ${year}, ${year}, ${year}, ${year},"/g namelist.input
    sed -i s/start_month.*$/"start_month = ${month}, ${month}, ${month}, ${month}, ${month}, ${month}, ${month},"/g namelist.input
    sed -i s/start_day.*$/"start_day = ${day}, ${day}, ${day}, ${day}, ${day}, ${day}, ${day},"/g namelist.input
    sed -i s/^\ restart\ .*$/" restart = .true.,"/g namelist.input
    
    # Add the reset_simulation_start option if it is not there already.
    if [ `grep "write_hist_at_0h_rst" namelist.input | wc -l` == 0 ]; then
        sed -i '/^\ restart\ .*/a \ write_hist_at_0h_rst = .true.,' namelist.input
    fi
    
    # Remove AFWA options that don't work with OpenMP compiled WRF.
    sed -i '/^ &afwa.*$/Q' namelist.input
fi
