import pandas as pd
import numpy as np
import re
import os

def gen_download_script(years, out_file, fut_ssp='ssp245', link_list='data/xu_file_list.csv'):
    """
    Write a script to wget all the boundary condition files required for a run over given years.
    
    Arguments:
        years: Years to run over. For year x data will be downloaded for Sep x to Feb x+1.
        out_files: The script file to write.
        fut_ssp: Which SSP to target.
        link_list: A file containing links to all archive files possible to download.
    """
    
    # Open the list of all download links, to subset from.
    files = pd.read_csv(link_list, header=None, names=['file'])
    files['name'] = files.file.str[97:None]
    
    # Loop through seasons, selecting the convective months.
    for season in years:
        season_from = season
        season_to = season+1
        season_string = f'{season_from}-{season_to}'
        files.loc[files.file.str.contains(f'atm_hist_{season_from}_(09|10|11|12)'), 'season'] = season_string
        files.loc[files.file.str.contains(f'atm_hist_{season_to}_(01|02)'), 'season'] = season_string
        files.loc[files.file.str.contains(f'atm_{fut_ssp}_{season_to}_(09|10|11|12)'), 'season'] = season_string   
        files.loc[files.file.str.contains(f'atm_{fut_ssp}_{season_to}_(01|02)'), 'season'] = season_string
    files = files.dropna()

    # Future files are listed twice with different download links, but link to the same files. Remove duplicates.
    files = files.loc[~files.duplicated('name'),:]
    files = files.sort_values(by='name')
    
    # Write a list of wget commands.
    cmds = [f'wget -c "{x}" -O "{y}"' for x, y in zip(files['file'], files['name'])]
    fp = open(out_file, 'w')
    fp.write('#!/bin/bash\n')
    fp.write('#PBS -l ncpus=1\n')
    fp.write('#PBS -l mem=2GB\n')
    fp.write('#PBS -l jobfs=2GB\n')
    fp.write('#PBS -q copyq\n')
    fp.write('#PBS -P or60\n')
    fp.write('#PBS -l walltime=02:00:00\n')
    fp.write('#PBS -l storage=gdata/up6\n')
    fp.write('#PBS -l wd\n')
    for cmd in cmds:
        fp.write(cmd + '\n')
    fp.close()
    
def set_up_WRF(year, template_dir, sims_dir, exp):
    """
    Set up directories ready for WPS and WRF runs for a given season, including 
    updating namelist files.
    
    Arguments:
        year: Event year (includes Sep year to Feb year+1).
        template_dir: Directory with template WPS and WRF setups.
        sims_dir: Directory in which simulations will be run.
        exp: Experiment to run (hist or ssp245).
    """

    sim_dir = f'{sims_dir}/{year}-{year+1}/'
    if not os.path.exists(sim_dir):
        os.mkdir(sim_dir)
        
    start_time = f'{year}-09-30_00:00:00'
    end_time = f'{year+1}-02-28_23:00:00'
        
    # WPS setup. Link executables + data files, copy and update namelist.
    if not os.path.exists(f'{sim_dir}/WPS'):
        os.system(f'cp -ar {template_dir}/WPS {sim_dir}/')
        os.system(f'sed -i s/year=.*$/"year={year}"/g {sim_dir}/WPS/mk.inputlist.sh')
        os.system(f'sed -i s/exp=.*$/"exp={exp}"/g {sim_dir}/WPS/mk.inputlist.sh')
        os.system(f'sed -i s/start_date.*$/"start_date = \'{start_time}\',"/g {sim_dir}/WPS/namelist.wps')
        os.system(f'sed -i s/end_date.*$/"end_date = \'{end_time}\',"/g {sim_dir}/WPS/namelist.wps')
    else:
        print('Skipping existing WPS...')

    # WRF setup. Link executables + data files, copy and update namelist.
    if not os.path.exists(f'{sim_dir}/WRF/'):
        os.mkdir(f'{sim_dir}/WRF')
        os.system(f'cp -ar {template_dir}/WRF {sim_dir}/')

        os.system(f'sed -i s/start_year.*$/"start_year = {start_time[0:4]}, {start_time[0:4]}, {start_time[0:4]},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/start_month.*$/"start_month = {start_time[5:7]}, {start_time[5:7]}, {start_time[5:7]},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/start_day.*$/"start_day = {start_time[8:10]}, {start_time[8:10]}, {start_time[8:10]},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/start_hour.*$/"start_hour = {start_time[11:13]}, {start_time[11:13]}, {start_time[11:13]},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/end_year.*$/"end_year = {end_time[0:4]}, {end_time[0:4]}, {end_time[0:4]},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/end_month.*$/"end_month = {end_time[5:7]}, {end_time[5:7]}, {end_time[5:7]},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/end_day.*$/"end_day = {end_time[8:10]}, {end_time[8:10]}, {end_time[8:10]},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/end_hour.*$/"end_hour = {end_time[11:13]}, {end_time[11:13]}, {end_time[11:13]},/g" {sim_dir}/WRF/namelist.input')
    else:
        print(f'Skipping existing WRF...')