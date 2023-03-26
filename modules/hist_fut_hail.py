import re
import os
import glob
import xarray
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from shapely.geometry.polygon import LinearRing
from matplotlib.patches import Rectangle

# Coordinates for cities of interest.
cities = {'Perth': (115.8606, -31.9559),
          'Sydney': (151.21, -33.867778),
          'Melbourne': (144.963056, -37.814167),
          'Brisbane': (153.028056, -27.467778),
          'Canberra': (149.126944, -35.293056)}

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
        files.loc[files.file.str.contains(f'atm_{fut_ssp}_{season_from}_(09|10|11|12)'), 'season'] = season_string   
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
    comma_q = "', '"
    comma = ', '
        
    # WPS setup. Link executables + data files, copy and update namelist.
    if not os.path.exists(f'{sim_dir}/WPS'):
        os.mkdir(f'{sim_dir}/WPS')
        os.system(f'ln -s {template_dir}/WPS/* {sim_dir}/WPS/')        
        os.system(f'rm {sim_dir}/WPS/namelist.wps {sim_dir}/WPS/mk.inputlist.sh')
        os.system(f'cp {template_dir}/WPS/namelist.wps {sim_dir}/WPS/')
        os.system(f'cp {template_dir}/WPS/mk.inputlist.sh {sim_dir}/WPS/')
        os.system(f'sed -i s/year=.*$/"year={year}"/g {sim_dir}/WPS/mk.inputlist.sh')
        os.system(f'sed -i s/exp=.*$/"exp={exp}"/g {sim_dir}/WPS/mk.inputlist.sh')
        os.system(f'sed -i s/start_date.*$/"start_date = \'{comma_q.join(np.repeat(start_time, 7))}\',"/g {sim_dir}/WPS/namelist.wps')
        os.system(f'sed -i s/end_date.*$/"end_date = \'{comma_q.join(np.repeat(end_time, 7))}\',"/g {sim_dir}/WPS/namelist.wps')
    else:
        print('Skipping existing WPS...')

    # WRF setup. Link executables + data files, copy and update namelist.
    if not os.path.exists(f'{sim_dir}/WRF/'):
        os.mkdir(f'{sim_dir}/WRF')
        os.system(f'ln -s {template_dir}/WRF/* {sim_dir}/WRF/')
        os.system(f'rm {sim_dir}/WRF/namelist.input')
        os.system(f'cp -f {template_dir}/WRF/namelist.input {sim_dir}/WRF/')
        os.system(f'sed -i s/start_year.*$/"start_year = {comma.join(np.repeat(start_time[0:4], 7))},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/start_month.*$/"start_month = {comma.join(np.repeat(start_time[5:7], 7))},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/start_day.*$/"start_day = {comma.join(np.repeat(start_time[8:10], 7))},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/start_hour.*$/"start_hour = {comma.join(np.repeat(start_time[11:13], 7))},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/end_year.*$/"end_year = {comma.join(np.repeat(end_time[0:4], 7))},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/end_month.*$/"end_month = {comma.join(np.repeat(end_time[5:7], 7))},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/end_day.*$/"end_day = {comma.join(np.repeat(end_time[8:10], 7))},/g" {sim_dir}/WRF/namelist.input')
        os.system(f'sed -i s/end_hour.*$/"end_hour = {comma.join(np.repeat(end_time[11:13], 7))},/g" {sim_dir}/WRF/namelist.input')
    else:
        print(f'Skipping existing WRF...')
        
def plot_wrf_domains(wps_dir, pts=None):
    """
    Plot WRF domains for a WPS setup. Requires geogrid to have been run first.
    
    Arguments:
        wps_dir: The WPS directory to find domain information in.
        pts: Points to add as a dictionary with keys a name and values an (x,y) pair.
    """
    
    dom_files = sorted(glob.glob(wps_dir+'/geo*.nc'))
    doms = [xarray.open_dataset(x) for x in dom_files] 
    
    max_map_fac = np.array([np.max(np.abs(x.MAPFAC_M.values-1) / 1 * 100) for x in doms])
    if np.any(max_map_fac > 5):
        print('WARNING: Map factor is too large or small; reconfigure domains.')
    
    # Custom functions for WRF domain plotting.
    def add_box(x, y, ax, colour='white', linestyle='solid'):
        geom = LinearRing(list(zip(list(x.isel(south_north=-1).values) + list(x.isel(south_north=0).values)[::-1],
                               list(y.isel(south_north=-1).values) + list(y.isel(south_north=0).values)[::-1])))
        ax.add_geometries([geom], crs=ccrs.PlateCarree(), linewidth=1.5, edgecolor=colour, facecolor='none', linestyle=linestyle)

    def add_doms(dom_list, ax, colour='white', linestyle='solid'):
        for dom in dom_list:
            x = dom.XLONG_M.isel(Time=0)
            y = dom.XLAT_M.isel(Time=0)
            add_box(x, y, ax, colour=colour, linestyle=linestyle)
    
    # Do the plotting.
    plt.rcParams['font.size'] = 18
    fig, ax = plt.subplots(figsize=(10,8), subplot_kw={'projection': ccrs.PlateCarree()})
    z = doms[0].isel(Time=0).HGT_M.values
    x = doms[0].isel(Time=0).XLONG_M.values
    y = doms[0].isel(Time=0).XLAT_M.values
    ax.pcolormesh(x,y,z,cmap='terrain', transform=ccrs.PlateCarree())
    ax.coastlines(linewidth=1.5)
    add_doms(doms[1:(len(doms))], ax=ax, colour='red')
    
    # Add optional points.
    for i in pts.keys():
        ax.scatter(pts[i][0], pts[i][1], color='red')
    
    plt.show()
    
def plot_wrf_domain_def(parent_id, i_parent_start, j_parent_start, ref_lon, ref_lat, e_we, e_sn, dx, dy, num_doms, scale_factor, pts=None):
    """
    Plot a WRF domain definition, assuming the projection is 'lat-lon' (ie rotated cylindrical equidistant) in WRF.
    
    Arguments:
        parent_id, i_parent_start, j_parent_start, ref_lon, 
        ref_lat, e_we, e_sn, dx, dy, num_doms, scale_factor: WRF parameters from namelist.wps.
        pts: Points to plot as dots on the map.
    """
    
    fig, ax = plt.subplots(figsize=(10,8), subplot_kw={'projection': ccrs.PlateCarree()})

    # Define the WRF projection.
    ref_lon_orig = ref_lon
    ref_lat_orig = ref_lat
    proj = ccrs.RotatedPole(pole_longitude=-(180 - ref_lon_orig),  
                            pole_latitude=-(90 + ref_lat_orig),
                            central_rotated_longitude=180.0,
                            globe=ccrs.Globe(ellipse=None,
                                             semimajor_axis=6370000,
                                             semiminor_axis=6370000))
    
    # Reproject the reference lat/long to the WRF projection (should be 0,0 really).
    (ref_lon, ref_lat) = proj.transform_point(x=ref_lon_orig, y=ref_lat_orig, src_crs=ccrs.Geodetic())
    
    # Add optional points.
    if not pts is None:
        for i in pts.keys():
            ax.scatter(pts[i][0], pts[i][1], color='red', transform=ccrs.PlateCarree())

    domains = {1: [ref_lon + e_we[0]*dx/2, ref_lat + e_sn[0]*dy/2, dx, dy]}
    for d in np.arange(num_doms):
        parent = parent_id[d]

        s = scale_factor
        if d == 0:
            s = 1

        p_dx = domains[parent][2]
        p_dy = domains[parent][3]
        d_dx = p_dx/s
        d_dy = p_dy/s

        x = domains[parent][0] - p_dx * i_parent_start[d]
        y = domains[parent][1] - p_dy * j_parent_start[d]
        w = -d_dx * e_we[d]
        h = -d_dy * e_sn[d]

        dom = Rectangle(xy=(x, y), width=w, height=h, fill=False, edgecolor='red', transform=proj)
        ax.add_patch(dom)
        if not d+1 in domains:
            domains[d+1] = [x, y, d_dx, d_dy]

    ax.coastlines()
    ax.scatter(ref_lon_orig, ref_lat_orig, marker='o', s=100)
    ax.set_xlim(ref_lon_orig - e_we[0]*dx/2 - 13, ref_lon_orig + e_we[0]*dx/2 + 13)
    ax.set_ylim(ref_lat_orig - e_sn[0]*dy/2 - 2, ref_lat_orig + e_sn[0]*dy/2 + 6)
    plt.show()
    
    print('parent_id = ' + ', '.join([str(x) for x in parent_id]) + ',')
    print('i_parent_start = ' + ', '.join([str(x) for x in i_parent_start]) + ',')
    print('j_parent_start = ' + ', '.join([str(x) for x in j_parent_start]) + ',')
    print('e_we = ' + ', '.join([str(x) for x in e_we]) + ',')
    print('e_sn = ' + ', '.join([str(x) for x in e_sn]) + ',')
    print('dx = ' + str(dx))
    print('dy = ' + str(dy))
    print('ref_lat = ' + str(ref_lat_orig))
    print('ref_lon = ' + str(ref_lon_orig))
    print('pole_lat = ' + str(90 + ref_lat_orig))
    print('pole_lon = 0')
    print('stand_lon = ' + str(round(180 - ref_lon_orig, 1)))