"""Module providing support functions for historic vs future hail comparisons."""

import os
import glob
import xarray
import geopandas
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
from skimage import morphology
import matplotlib.pyplot as plt
from cartopy.io import shapereader
import matplotlib.ticker as mticker
from matplotlib import colors
from matplotlib.patches import Rectangle
from shapely.geometry.polygon import LinearRing
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Coordinates for cities of interest.
cities = {
    'Perth': (115.8606, -31.9559),
    'Sydney': (151.21, -33.867778),
    'Melbourne': (144.963056, -37.814167),
    'Brisbane': (153.028056, -27.467778),
    'Canberra': (149.126944, -35.293056),
}

# Coordinates for remote regions of interest.
remote = {
    'Adelaide': (138.5999, -34.9287),
    'Burketown': (139.5423, -17.7383),
    'Kalgoorlie': (121.4656, -30.7582),
}


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
        season_to = season + 1
        season_string = f'{season_from}-{season_to}'
        files.loc[files.file.str.contains(f'atm_hist_{season_from}_(09|10|11|12)'), 'season'] = (
            season_string
        )
        files.loc[files.file.str.contains(f'atm_hist_{season_to}_(01|02)'), 'season'] = (
            season_string
        )
        files.loc[
            files.file.str.contains(f'atm_{fut_ssp}_{season_from}_(09|10|11|12)'), 'season'
        ] = season_string
        files.loc[files.file.str.contains(f'atm_{fut_ssp}_{season_to}_(01|02)'), 'season'] = (
            season_string
        )
    files = files.dropna()

    # Future files are listed twice with different download links, but link to the same files. Remove duplicates.
    files = files.loc[~files.duplicated('name'), :]
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


def set_up_WRF(year, template_dir, namelist_dir, sims_dir, exp):
    """
    Set up directories ready for WPS and WRF runs for a given season, including
    updating namelist files.

    Arguments:
        year: Event year (includes Sep year to Feb year+1).
        template_dir: Directory with template WPS and WRF setups.
        namelist_dir: Directory from which to source WPS and WRF namelists.
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
        os.system(f'rm {sim_dir}/WPS/mk.inputlist.sh')
        os.system(f'cp {namelist_dir}/namelist.wps {sim_dir}/WPS/')
        os.system(f'cp {template_dir}/WPS/mk.inputlist.sh {sim_dir}/WPS/')
        os.system(f'sed -i s/year=.*$/"year={year}"/g {sim_dir}/WPS/mk.inputlist.sh')
        os.system(f'sed -i s/exp=.*$/"exp={exp}"/g {sim_dir}/WPS/mk.inputlist.sh')
        os.system(
            f'sed -i s/start_date.*$/"start_date = \'{comma_q.join(np.repeat(start_time, 7))}\',"/g {sim_dir}/WPS/namelist.wps'
        )
        os.system(
            f'sed -i s/end_date.*$/"end_date = \'{comma_q.join(np.repeat(end_time, 7))}\',"/g {sim_dir}/WPS/namelist.wps'
        )
    else:
        print('Skipping existing WPS...')

    # WRF setup. Link executables + data files, copy and update namelist.
    if not os.path.exists(f'{sim_dir}/WRF/'):
        os.mkdir(f'{sim_dir}/WRF')
        os.system(f'ln -s {template_dir}/WRF/* {sim_dir}/WRF/')
        os.system(f'cp -f {namelist_dir}/namelist.input {sim_dir}/WRF/')
        os.system(
            f'sed -i s/start_year.*$/"start_year = {comma.join(np.repeat(start_time[0:4], 7))},"/g {sim_dir}/WRF/namelist.input'
        )
        os.system(
            f'sed -i s/start_month.*$/"start_month = {comma.join(np.repeat(start_time[5:7], 7))},"/g {sim_dir}/WRF/namelist.input'
        )
        os.system(
            f'sed -i s/start_day.*$/"start_day = {comma.join(np.repeat(start_time[8:10], 7))},"/g {sim_dir}/WRF/namelist.input'
        )
        os.system(
            f'sed -i s/start_hour.*$/"start_hour = {comma.join(np.repeat(start_time[11:13], 7))},"/g {sim_dir}/WRF/namelist.input'
        )
        os.system(
            f'sed -i s/end_year.*$/"end_year = {comma.join(np.repeat(end_time[0:4], 7))},"/g {sim_dir}/WRF/namelist.input'
        )
        os.system(
            f'sed -i s/end_month.*$/"end_month = {comma.join(np.repeat(end_time[5:7], 7))},"/g {sim_dir}/WRF/namelist.input'
        )
        os.system(
            f'sed -i s/end_day.*$/"end_day = {comma.join(np.repeat(end_time[8:10], 7))},"/g {sim_dir}/WRF/namelist.input'
        )
        os.system(
            f'sed -i s/end_hour.*$/"end_hour = {comma.join(np.repeat(end_time[11:13], 7))},"/g {sim_dir}/WRF/namelist.input'
        )
    else:
        print('Skipping existing WRF...')


def plot_wrf_domains(
    wps_dir, pts=None, figsize=(10, 8), proj=ccrs.PlateCarree(), fontsize=14, labels=None, file=None
):
    """
    Plot WRF domains for a WPS setup. Requires geogrid to have been run first.

    Arguments:
        wps_dir: The WPS directory to find domain information in.
        pts: Points to add as a dictionary with keys a name and values an (x,y) pair.
        figsize: The figure size (width x height).
        proj: Projection to display map in.
        labels: Point name: (label, offset) for each point to display, with offset in map coordinates.
        file: Output file for figure.
    """

    dom_files = sorted(glob.glob(wps_dir + '/geo*.nc'))
    doms = [xarray.open_dataset(x) for x in dom_files]

    max_map_fac = np.array([np.max(np.abs(x.MAPFAC_M.values - 1) / 1 * 100) for x in doms])
    if np.any(max_map_fac > 5):
        print('WARNING: Map factor is large or small (below 0.95 or over 1.05).')
        facts = [np.max(x.MAPFAC_M.values) for x in doms]
        print(f'Max factor is {np.max(facts):.3} in domain {np.argmax(facts)+1}.')
        print(f'Min factor is {np.min(facts):.3} in domain {np.argmin(facts)+1}.')

    # Custom functions for WRF domain plotting.
    def add_box(x, y, ax, colour='white', linestyle='solid'):
        geom = LinearRing(
            list(
                zip(
                    list(x.isel(south_north=-1).values) + list(x.isel(south_north=0).values)[::-1],
                    list(y.isel(south_north=-1).values) + list(y.isel(south_north=0).values)[::-1],
                )
            )
        )
        ax.add_geometries(
            [geom],
            crs=ccrs.PlateCarree(),
            linewidth=1.5,
            edgecolor=colour,
            facecolor='none',
            linestyle=linestyle,
        )

    def add_doms(dom_list, ax, colour='white', linestyle='solid'):
        for dom in dom_list:
            x = dom.XLONG_M.isel(Time=0)
            y = dom.XLAT_M.isel(Time=0)
            add_box(x, y, ax, colour=colour, linestyle=linestyle)

    # Do the plotting.
    plt.rcParams['font.size'] = fontsize
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': proj})
    # z = doms[0].isel(Time=0).HGT_M.values
    x = doms[0].isel(Time=0).XLONG_M.values
    y = doms[0].isel(Time=0).XLAT_M.values
    # pc = ax.pcolormesh(x,y,z,cmap='terrain', transform=ccrs.PlateCarree())
    ax.coastlines(linewidth=1.5)
    # fig.colorbar(pc)
    add_doms([doms[0]], ax=ax, colour='orange')
    add_doms([doms[x] for x in [1, 3]], ax=ax, colour='darkgreen')
    add_doms([doms[x] for x in [2, 4, 5, 6]], ax=ax, colour='blue')

    ax.set_xlim(np.min(x) - 2, np.max(x) + 2)
    ax.set_ylim(np.min(y) - 2, np.max(y) + 2)

    gl = ax.gridlines(crs=proj, draw_labels=True, alpha=0.3)
    gl.top_labels = gl.right_labels = False

    # Add optional points.
    for i in pts.keys():
        ax.scatter(pts[i][0], pts[i][1], color='red', transform=ccrs.PlateCarree())

        if labels is not None and i in labels:
            ax.annotate(
                labels[i][0],
                (pts[i][0], pts[i][1]),
                textcoords='offset points',
                xytext=labels[i][1],
                ha='center',
                color='black',
            )

    if file is not None:
        plt.savefig(fname=file, dpi=300, bbox_inches='tight')

    plt.show()


def plot_wrf_domain_def(
    parent_id,
    i_parent_start,
    j_parent_start,
    ref_lon,
    ref_lat,
    e_we,
    e_sn,
    dx,
    dy,
    num_doms,
    scale_factor,
    pts=None,
):
    """
    Plot a WRF domain definition, assuming the projection is 'lat-lon' (ie rotated cylindrical equidistant) in WRF.

    Arguments:
        parent_id, i_parent_start, j_parent_start, ref_lon,
        ref_lat, e_we, e_sn, dx, dy, num_doms, scale_factor: WRF parameters from namelist.wps.
        pts: Points to plot as dots on the map.
    """

    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

    # Define the WRF projection.
    ref_lon_orig = ref_lon
    ref_lat_orig = ref_lat
    proj = ccrs.RotatedPole(
        pole_longitude=-(180 - ref_lon_orig),
        pole_latitude=-(90 + ref_lat_orig),
        central_rotated_longitude=180.0,
        globe=ccrs.Globe(ellipse=None, semimajor_axis=6370000, semiminor_axis=6370000),
    )

    # Reproject the reference lat/long to the WRF projection (should be 0,0 really).
    (ref_lon, ref_lat) = proj.transform_point(
        x=ref_lon_orig, y=ref_lat_orig, src_crs=ccrs.Geodetic()
    )

    # Add optional points.
    if pts is not None:
        for i in pts.keys():
            ax.scatter(pts[i][0], pts[i][1], color='red', transform=ccrs.PlateCarree())

    domains = {1: [ref_lon + e_we[0] * dx / 2, ref_lat + e_sn[0] * dy / 2, dx, dy]}
    for d in np.arange(num_doms):
        parent = parent_id[d]

        s = scale_factor
        if d == 0:
            s = 1

        p_dx = domains[parent][2]
        p_dy = domains[parent][3]
        d_dx = p_dx / s
        d_dy = p_dy / s

        x = domains[parent][0] - p_dx * i_parent_start[d]
        y = domains[parent][1] - p_dy * j_parent_start[d]
        w = -d_dx * e_we[d]
        h = -d_dy * e_sn[d]

        dom = Rectangle(xy=(x, y), width=w, height=h, fill=False, edgecolor='red', transform=proj)
        ax.add_patch(dom)
        if d + 1 not in domains:
            domains[d + 1] = [x, y, d_dx, d_dy]

    ax.coastlines()
    ax.scatter(ref_lon_orig, ref_lat_orig, marker='o', s=100)
    ax.set_xlim(ref_lon_orig - e_we[0] * dx / 2 - 13, ref_lon_orig + e_we[0] * dx / 2 + 13)
    ax.set_ylim(ref_lat_orig - e_sn[0] * dy / 2 - 2, ref_lat_orig + e_sn[0] * dy / 2 + 6)
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


def mean_temp_per_year(
    infiles={
        'SSP2-4.5': '/scratch/w42/tr2908/xu_data/lnd.ssp245.*.nc',
        'SSP5-8.5': '/scratch/w42/tr2908/xu_data/lnd.ssp585.*.nc',
    },
    outfile='/g/data/up6/tr2908/hist_future_hail/xu_data/ssp245/mean_temperature_by_year.nc',
):
    """
    Calculate mean near-surface temperatue from Xu et al lnd* files, and write a timeseries per year to a netcdf file.

    Arguments:
        infiles: Dictionary with name: path_string for the lnd* files to use.
        outfile: The output file to write.
    """

    if not os.path.exists(outfile):
        mean_temps = []
        for sim, path in infiles.items():
            d = xarray.open_mfdataset(path, parallel=True)
            d['year'] = d.time.dt.year
            d = d.drop('height')

            # Weight by cosine of latitude.
            d['area_weight'], _ = xarray.broadcast(np.cos(np.deg2rad(d.lat)), d.tas)
            d['weighted_temp'] = d.tas * d.area_weight

            # Calculate means and write to disk.
            mean_temp = d.weighted_temp.groupby(d.year).sum(
                ['lat', 'lon', 'time']
            ) / d.area_weight.groupby(d.year).sum(['lat', 'lon', 'time'])
            mean_temp.attrs['description'] = (
                f"Average near-surface temperature calculated using d.tas.groupby(d.year).mean(['lat', 'lon', 'time']) where d is all Xu et al. {sim} data."
            )
            mean_temp.name = 'tas'
            mean_temp = mean_temp.expand_dims({'simulation': [sim]})
            mean_temps.append(mean_temp)

        mean_temps = xarray.merge(mean_temps)
        mean_temps.to_netcdf(outfile)

    return xarray.open_dataset(outfile)


def conv_stats(
    basic_vars=['hailcast_diam_max', 'longitude', 'latitude'],
    exclude_vars=['shear_u', 'shear_v', 'positive_shear'],
    cache_dir='/g/data/up6/tr2908/hist_future_hail/WRF_v4.4/results_cache',
    cities={'Perth': 3, 'Melbourne': 5, 'Brisbane': 6, 'Sydney + Canberra': 7},
    sims={
        'Historical': '/g/data/up6/tr2908/hist_future_hail/WRF_v4.4/simulations/hist/',
        'SSP2-4.5': '/g/data/up6/tr2908/hist_future_hail/WRF_v4.4/simulations/ssp245/',
    },
    city_cells=None,
):
    """
    Process means and quantiles per timestep in basic* and conv* files.

    Arguments:
        basic_vars: basic*.nc variables to read.
        exclude_vars: Exclude these variables.
        cache_dir: Output cache directory.
        cities: WRF domain numbers for each city.
        sims: name: directory dictionary for simulations to process.
        city_cells: isel definition to select city pixels if required.

    Returns: Means and quantiles for each variable.
    """

    res = []
    comp = dict(zlib=True, shuffle=True, complevel=5)

    for city in cities:
        for sim in sims:
            print(f'{sim} {city}')

            outfile = f'{cache_dir}/{city}_{sim}_hourly_stats.nc'
            if not os.path.exists(outfile):
                basic_files = sorted(glob.glob(f'{sims[sim]}/*/WRF/basic_params_d0{cities[city]}*'))
                conv_files = sorted(glob.glob(f'{sims[sim]}/*/WRF/conv_params_d0{cities[city]}*'))

                print('Opening datasets...')
                b = xarray.open_mfdataset(
                    basic_files, parallel=True, combine='nested', concat_dim='time'
                )[basic_vars]
                c = xarray.open_mfdataset(
                    conv_files, parallel=True, combine='nested', concat_dim='time'
                )
                dat = xarray.merge([b, c])

                dat = dat.drop(exclude_vars)
                dat['latitude'] = dat.latitude.isel(time=0)
                dat['longitude'] = dat.longitude.isel(time=0)
                dat = dat.chunk({'time': 500})

                # We are only interested in non-zero hail sizes.
                dat['hailcast_diam_max'] = dat.hailcast_diam_max.where(dat.hailcast_diam_max > 0)

                if city in city_cells:
                    dat = dat.isel(city_cells[city])

                means = dat.mean(['south_north', 'west_east'])
                quants = dat.quantile(dim=['south_north', 'west_east'], q=[0, 0.05, 0.5, 0.95, 1])
                means = means.rename({n: f'{n}_mean' for n in list(means.data_vars)})
                quants = quants.rename({n: f'{n}_quantile' for n in list(quants.data_vars)})
                stats = xarray.merge([means, quants])

                stats = stats.drop(
                    ['longitude_mean', 'latitude_mean', 'longitude_quantile', 'latitude_quantile']
                )
                stats = stats.expand_dims({'city': [city], 'sim': [sim]})

                stats['latitude'] = dat.latitude
                stats['longitude'] = dat.longitude

                print('Writing to netcdf...')
                encoding = {var: comp for var in stats.data_vars}
                stats.to_netcdf(outfile, encoding=encoding)

            d = xarray.open_dataset(outfile).load()
            d['latitude'] = d.latitude.expand_dims({'city': d.city})
            d['longitude'] = d.longitude.expand_dims({'city': d.city})
            res.append(d)

    res = xarray.merge(res)
    res = res.assign_coords({'sim': ['Historical', '+2.8 C']})
    return res


def plot_map_to_ax(
    dat,
    ax,
    coastlines=True,
    grid=True,
    dat_proj=ccrs.PlateCarree(),
    disp_proj=ccrs.PlateCarree(),
    title=None,
    colour_scale=None,
    cmap=None,
    norm=None,
    cbar_ticks=None,
    tick_labels=None,
    contour=False,
    stippling=None,
    stipple_size=3,
    colourbar=True,
    ticks_left=True,
    ticks_bottom=True,
    cbar_aspect=25,
    cbar_fraction=0.07,
    cbar_shrink=0.4,
    cbar_pad=0.015,
    cbar_label=None,
    cbar_orientation='vertical',
    coastlines_colour='black',
    xlims=None,
    ylims=None,
    num_ticks=None,
    divergent=False,
    cbar_inset=False,
    title_inset=False,
    discrete=False,
    log_scale=False,
    nan_colour='#eeeeee',
    axis_off=False,
    country=None,
    annotations=None,
    polygons=None,
    polygon_colour='black',
    hatch='.',
    pts=None,
    left_title=None,
):
    """
    Plot data on a map to a specified plot axis object.

    Arguments:

        - dat: DataSet to plot or list of datasets to plot.
        - ax: GeoAxes object to plot to.
        - dat_proj, dist_proj: Data and display projections.
        - figsize: Figure size width x height.
        - coastlines: Show coastlines?
        - grid: Show grid?
        - ncol/nrows: Number of columns/rows to plot.
        - title: Title for the plot.
        - colour_scale: None for default, or a tuple of min/max values for the scale.
        - cmap: The matplotlib colour map to use.
        - norm: A norm object for colours (e.g. colors.BoundaryNorm).
        - cbar_ticks: Colour bar ticks.
        - tick_labels: Colour bar tick labels.
        - contour: Plot using xarray's contourf function?
        - stippling: True where trippling should appear.
        - stipple_size: Size for stippling points.
        - colourbar: Include a colourbar?
        - ticks_left: Include ticks on left of plot?
        - ticks_bottom: Include ticks on bottom of plot?
        - cbar_aspect: colorbar aspect ratio?
        - cbar_fraction: fraction argument for colorbar().
        - cbar_shrink: shrink argument for colorbar().
        - cbar_pad: pad argument for colorbar().
        - cbar_label: Overwrite label?
        - cbar_orientation: orientation argument for colorbar().
        - coastlines_colour: Colour for coastlines.
        - xlims, ylims: x and y plot limits.
        - num_ticks: Number of ticks for x and y axes (None for auto).
        - divergent: Is the colour scale divergent? If so make zero central.
        - cbar_inset: Inset the colorbar in lower left?
        - title_inset: Inset the title in the upper left?
        - discrete: Make the colour bar discrete?
        - log_scale: Make the colour scale log-scaled?
        - nan_colour: Colour for missing values.
        - axis_off: Turn off all axes.
        - country: Plot coastlines only for a specific country.
        - annotations: Add annotations to the map - dictionary of {'Text': [x, y, xadj, yadj, ha]} where
                       x, y are position of text, xadj, yadj give offsets to the label, ha is 'left' or
                       'right' for horizontal anchor.
        - polygons: If specified, draw each polygon onto the plot.
        - polygon_colour: Colour for polygons.
        - hatch: Hatching to use for stippling ('.' for scatterplot, '//' etc for contour hatching).

    """

    col_min = None
    col_max = None

    # Rasterize elements with zorder below 0.
    ax.set_rasterization_zorder(0)

    if colour_scale is not None:
        col_min = colour_scale[0]
        col_max = colour_scale[1]

    if divergent:
        if col_min is None or col_max is None:
            col_min = dat.min()
            col_max = dat.max()

        col_min = -1 * np.max(np.abs([col_min, col_max]))
        col_max = np.max(np.abs([col_min, col_max]))

    cbar_spacing = 'proportional'

    if discrete:
        assert cbar_ticks is not None, 'Discrete colorbar requires cbar_ticks'
        assert cbar_ticks == sorted(cbar_ticks), 'cbar_ticks must be sorted for discrete plot.'
        cbar_ticks = np.array(cbar_ticks + [np.max(cbar_ticks) + 1])
        norm = colors.BoundaryNorm(cbar_ticks, ncolors=len(cbar_ticks) - 1)
        cbar_ticks = (cbar_ticks[0:-1] + cbar_ticks[1:]) / 2
        cbar_spacing = 'uniform'

    if log_scale:
        norm = colors.LogNorm(vmin=col_min, vmax=col_max)

    fmt = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    fmt.set_powerlimits((-4, 6))
    cbar_args = {
        'spacing': cbar_spacing,
        'fraction': cbar_fraction,
        'ticks': cbar_ticks,
        'aspect': cbar_aspect,
        'shrink': cbar_shrink,
        'pad': cbar_pad,
        'orientation': cbar_orientation,
        'format': '%g',
    }

    if cbar_inset:
        cax = inset_axes(
            ax,
            width='50%',
            height='3%',
            loc='lower left',
            bbox_to_anchor=(0.05, 0.15, 1, 1),
            bbox_transform=ax.transAxes,
        )
        cbar_args['cax'] = cax

    if cbar_label is not None:
        cbar_args['label'] = cbar_label

    if colourbar is False:
        cbar_args = None

    cmap = plt.get_cmap(cmap).copy()
    cmap.set_bad(nan_colour)

    if not contour:
        res = dat.plot(
            ax=ax,
            transform=dat_proj,
            vmin=col_min,
            vmax=col_max,
            cmap=cmap,
            norm=norm,
            cbar_kwargs=cbar_args,
            add_colorbar=colourbar,
            zorder=-10,
        )
    else:
        res = dat.plot.contourf(
            ax=ax,
            transform=dat_proj,
            vmin=col_min,
            vmax=col_max,
            cmap=cmap,
            norm=norm,
            cbar_kwargs=cbar_args,
            add_colorbar=colourbar,
        )

    if stippling is not None:
        ax.autoscale(False)
        if hatch == '.':
            pts = stippling.where(stippling).to_dataframe().dropna().reset_index()
            ax.scatter(
                x=pts[stippling.dims[1]],
                y=pts[stippling.dims[0]],
                marker='.',
                color='black',
                transform=dat_proj,
                s=stipple_size,
            )
        else:
            stippling.plot.contourf(
                hatches=['', hatch], levels=[0, 0.5, 1], colors='none', ax=ax, add_colorbar=False
            )

    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)
    if tick_labels is not None:
        assert len(tick_labels) == len(cbar_ticks), 'Labels and ticks must have same length'
        res.colorbar.ax.set_yticklabels(tick_labels)
    if left_title is not None:
        if title_inset:
            title = f'{left_title} {title}'
        else:
            ax.set_title(left_title, fontsize=plt.rcParams['font.size'], loc='left')
    if title is not None:
        if title_inset:
            ax.annotate(
                text=title,
                xy=(0.05, 0.9),
                xycoords='axes fraction',
                fontweight='bold',
                fontsize=plt.rcParams['font.size'],
            )
        else:
            ax.set_title(title, fontsize=plt.rcParams['font.size'])
    if polygons is not None:
        poly = geopandas.GeoSeries(polygons).unary_union
        ax.add_geometries(
            poly, crs=ccrs.PlateCarree(), facecolor='none', edgecolor=polygon_colour, linewidth=1.75
        )
    if coastlines:
        if country is not None:
            shpfilename = shapereader.natural_earth(
                resolution='10m', category='cultural', name='admin_0_countries'
            )
            df = geopandas.read_file(shpfilename)
            poly = df.loc[df['ADMIN'] == country]['geometry'].values[0]
            ax.add_geometries(
                poly,
                crs=ccrs.PlateCarree(),
                facecolor='none',
                edgecolor=coastlines_colour,
                linewidth=0.75,
            )
        else:
            ax.coastlines(color=coastlines_colour)
    if grid:
        locator = None
        if num_ticks is not None:
            locator = mticker.MaxNLocator(nbins=num_ticks + 1)
        gl = ax.gridlines(crs=disp_proj, draw_labels=True, alpha=0.5, xlocs=locator, ylocs=locator)
        gl.top_labels = gl.right_labels = False
        gl.left_labels = ticks_left
        gl.bottom_labels = ticks_bottom
    if axis_off:
        ax.axis('off')

    if pts is not None:
        for i in pts.keys():
            ax.scatter(pts[i][0], pts[i][1], color='black')

    if annotations is not None:
        for text, [x, y, xadj, yadj, ha] in annotations.items():
            if np.abs(xadj) >= 1 or np.abs(yadj) >= 1:
                if ha == 'right' or ha == 'left':
                    ax.plot(
                        [x, x + xadj - (0.2 * np.sign(xadj))], [y, y + yadj + 0.2], color='black'
                    )
                elif ha == 'center':
                    ax.plot(
                        [x, x + xadj - (0.2 * np.sign(xadj))], [y, y + yadj - 0.2], color='black'
                    )
                    if yadj < 0:
                        print('Warning: ha=center and negative y adjustment are not supported.')
                else:
                    assert 1 == 0, 'Invalid value of ha.'
            ax.annotate(xy=(x + xadj, y + yadj), text=text, ha=ha)

    return res


def plot_map(
    dat,
    dat_proj=ccrs.PlateCarree(),
    disp_proj=ccrs.PlateCarree(),
    figsize=(12, 8),
    grid=True,
    ncols=1,
    nrows=1,
    title=None,
    share_scale=False,
    colour_scale=None,
    cbar_ticks=None,
    tick_labels=None,
    file=None,
    scale_label='',
    share_axes=False,
    ticks_left=True,
    ticks_bottom=True,
    wspace=0.05,
    hspace=0.05,
    stippling=None,
    cbar_adjust=0.862,
    cbar_pad=0.015,
    col_labels=None,
    row_labels=None,
    xlims=None,
    ylims=None,
    show=True,
    shared_scale_quantiles=(0, 1),
    polygons=None,
    letter_labels=False,
    **kwargs,
):
    """
    Plot data on a map.

    Arguments:

        - dat: DataSet to plot or list of datasets to plot.
        - dat_proj, dist_proj: Data and display projections.
        - figsize: Figure size width x height.
        - grid: Show grid?
        - ncols/nrows: Number of columns/rows to plot.
        - title: Title(s) for the plot(s).
        - share_scale: Make the range of values in each plot the same?
        - colour_scale: Tuple with min/max values to use on scale. Overwritten by share_scale.
        - cbar_ticks: Ticks for the colourbar.
        - tick_labels: Colour bar tick labels.
        - file: If specified save to 'file' instead of showing onscreen.
        - scale_label: The label for a shared scale.
        - share_axes: Share left/bottom axes?
        - ticks_left, ticks_bottom: Display ticks on the left/bottom of plots?
        - wspace, hspace: gridspec wspace and hspace arguments.
        - stippling: Stippling per axis.
        - cbar_adjust: Amount to shrink plots by to add cbar for shared scale.
        - cbar_pad: Padding between figure and colour bar.
        - col_labels/row_labels: Labels for each column/row; overwrites individial plot titles.
        - xlims, ylims: x and y limits.
        - show: Show the map?
        - shared_scale_quantiles: Quantiles for a shared scale.
        - polygon: If specified, polygons to put on each map.
        - letter_labels: Use a letter to label each subplot?
        - kwargs: Extra arguments to plot_map_to_ax.

    Return:
        - The axis plotted to.

    """

    fig, ax = plt.subplots(
        figsize=figsize,
        ncols=ncols,
        nrows=nrows,
        subplot_kw={'projection': disp_proj},
        gridspec_kw={'wspace': wspace, 'hspace': hspace},
    )
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    letters = [
        'a)',
        'b)',
        'c)',
        'd)',
        'e)',
        'f)',
        'g)',
        'h)',
        'i)',
        'j)',
        'k)',
        'l)',
        'm)',
        'n)',
        'o)',
        'p)',
        'q)',
        'r)',
        's)',
        't)',
        'u)',
        'v)',
        'w)',
        'x)',
        'y)',
        'z)',
    ]

    if not isinstance(dat, list):
        im = plot_map_to_ax(
            dat=dat,
            ax=ax,
            grid=grid,
            dat_proj=dat_proj,
            disp_proj=disp_proj,
            title=title,
            stippling=stippling,
            colour_scale=colour_scale,
            cbar_ticks=cbar_ticks,
            cbar_pad=cbar_pad,
            tick_labels=tick_labels,
            ticks_left=ticks_left,
            ticks_bottom=ticks_bottom,
            xlims=xlims,
            ylims=ylims,
            polygons=polygons,
            **kwargs,
        )
    else:
        assert ncols * nrows >= len(dat), 'Not enough cols/rows to fit all plots.'

        if share_scale:
            all_vals = np.array([])

            for d in dat:
                all_vals = np.concatenate([all_vals, np.array(d.values.flat)])

            colour_scale = (
                np.nanquantile(all_vals, shared_scale_quantiles[0]),
                np.nanquantile(all_vals, shared_scale_quantiles[1]),
            )
            assert not (
                np.isnan(colour_scale[0]) or np.isnan(colour_scale[1])
            ), 'share_scale cannot be used with subplots missing data.'

        for i, d in enumerate(dat):
            ax_title = None
            if title is not None:
                ax_title = title[i]

            ax_poly = None
            if polygons is not None:
                ax_poly = polygons[i]

            tb = ticks_bottom
            tl = ticks_left
            if share_axes:
                if i < (ncols * nrows) - ncols:
                    tb = False
                if i % ncols != 0:
                    tl = False

            stipple = None if stippling is None else stippling[i]
            proj = dat_proj if not isinstance(dat_proj, list) else dat_proj[i]
            xlim = xlims if not isinstance(xlims, list) else xlims[i]
            ylim = ylims if not isinstance(ylims, list) else ylims[i]

            left_title = None
            if letter_labels:
                left_title = letters.pop(0)

            im = plot_map_to_ax(
                dat=d,
                ax=ax.flat[i],
                grid=grid,
                dat_proj=proj,
                disp_proj=disp_proj,
                title=ax_title,
                colour_scale=colour_scale,
                cbar_pad=cbar_pad,
                cbar_ticks=cbar_ticks,
                tick_labels=tick_labels,
                colourbar=(not share_scale),
                stippling=stipple,
                xlims=xlim,
                ylims=ylim,
                ticks_left=tl,
                ticks_bottom=tb,
                polygons=ax_poly,
                left_title=left_title,
                **kwargs,
            )

        while i + 1 < len(ax.flat):
            fig.delaxes(ax.flat[i + 1])
            i = i + 1

        if share_scale:
            fig.subplots_adjust(right=cbar_adjust)
            cbar_ax = fig.add_axes([cbar_adjust + cbar_pad, 0.23, 0.02, 0.55])
            fmt = mticker.ScalarFormatter(useOffset=False, useMathText=True)
            fmt.set_powerlimits((-4, 6))
            _ = fig.colorbar(
                im, ax=ax, cax=cbar_ax, ticks=cbar_ticks, label=scale_label, format=fmt
            )

        if col_labels is not None or row_labels is not None:
            for a in ax.flat:
                a.set_title('')

            if col_labels is not None:
                axes = ax if ax.ndim == 1 else ax[0, :]
                for a, lab in zip(axes, col_labels):
                    a.set_title(lab, fontsize=plt.rcParams['font.size'])

            if row_labels is not None:
                fig.subplots_adjust(left=0.02)
                lab_ax = fig.add_axes([0, 0.11, 0.02, 0.78], autoscale_on=True)
                lab_ax.axis('off')

                for i, lab in enumerate(row_labels):
                    p = 0.03 + (i / len(row_labels)) * 1.33
                    lab_ax.annotate(
                        lab, xy=(0.5, 1 - p), rotation=90, xycoords='axes fraction', ha='center'
                    )

    if file is not None:
        plt.savefig(fname=file, dpi=300, bbox_inches='tight')

        if show:
            plt.show()
        plt.close()
    else:
        if show:
            plt.show()

    return ax


def plot_maxes(
    maxes,
    lm,
    domain,
    figsize=(12, 2),
    variable='hailcast_diam_max',
    lab='Max hail diam [mm]',
    file_dir=None,
):
    """
    Plot maximum hail diameters by location for historical and future
    epochs both nonmasked and masked to land points only.

    Arguments:
        maxes: Max hail diameters per location per epoch.
        lm: The landmask.
        domain: The domain to use in titles.
        figsize: Figure size.
        variable: Variable to plot.
        lab: Label for the variable.
        file_dir: Optional directory to save plot to.
    """

    x = maxes.longitude
    y = maxes.latitude
    z_hist = maxes[variable].sel(epoch='historical')
    z_futu = maxes[variable].sel(epoch='ssp245')

    plt.rcParams['font.size'] = 12
    fig, axs = plt.subplots(
        ncols=4,
        figsize=figsize,
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace': 0.3},
    )
    pc_hist = axs[0].pcolormesh(x, y, z_hist, cmap='Spectral_r', transform=ccrs.PlateCarree())
    pc_futu = axs[1].pcolormesh(x, y, z_futu, cmap='Spectral_r', transform=ccrs.PlateCarree())
    pc_hist_land = axs[2].pcolormesh(
        x, y, z_hist.where(lm == 1), cmap='Spectral_r', transform=ccrs.PlateCarree()
    )
    pc_futu_land = axs[3].pcolormesh(
        x, y, z_futu.where(lm == 1), cmap='Spectral_r', transform=ccrs.PlateCarree()
    )
    plt.colorbar(pc_hist, label=lab)
    plt.colorbar(pc_futu, label=lab)
    plt.colorbar(pc_hist_land, label=lab)
    plt.colorbar(pc_futu_land, label=lab)

    for ax in axs:
        ax.coastlines(linewidth=1.5)

    domain_str = domain.replace('_', '/')
    axs[0].set_title(f'{domain_str}\nhist')
    axs[1].set_title(f'{domain_str}\nssp245')
    axs[2].set_title(f'{domain_str}\nhist')
    axs[3].set_title(f'{domain_str}\nssp245')
    
    if file_dir is not None:
        plt.savefig(fname=f'{file_dir}/maxes_{domain}_{variable}.pdf', dpi=300, bbox_inches='tight')

def process_maxima(
    sim_dir,
    plot=True,
    domains={'Perth': 3, 'Melbourne': 5, 'Brisbane': 6, 'Sydney_Canberra': 7},
    time_adjust={'Perth': 8, 'Melbourne': 11, 'Brisbane': 10, 'Sydney_Canberra': 11},
    results_dir='/g/data/up6/tr2908/hist_future_hail/results/',
    variables=['hailcast_diam_max', 'wind_10m'],
    file_dir='paper/figures/'
):
    """
    Find block maxima for selected variables. On the way save the landmask
    and plot max hail sizes by location.

    Argumnets:
        sim_dir: Simulation base directory.
        plot: Plot max hail sizes?
        figsize: Figure size.
        domains: The name: number of domains to consider.
        time_adjust: Time adjustments to make from UTC time (in hours).
        results_dir: Output directory for results.
        variables: Variables to process.
        file_dir: Output directory for hailcost maxima plots.
    """

    for domain, d in domains.items():
        print(f'Processing d0{d}...')

        # Define cache files.
        dmn = domain.replace(' + ', '_')
        lm_file = f'{results_dir}/{dmn}_landmask.nc'
        maxes_file = f'{results_dir}/{dmn}_domain_maximums.nc'
        h_maxima_file = f'{results_dir}/{dmn}_hist_block_maxima.feather'
        f_maxima_file = f'{results_dir}/{dmn}_ssp245_block_maxima.feather'

        # Get land mask.
        if not os.path.exists(lm_file):
            print('Computing dilated landmask...')
            lm = xarray.open_dataset(glob.glob(f'{sim_dir}/hist/*/WRF/wrfout_d0{d}*')[0])
            lm = lm.LANDMASK.isel(Time=0).astype(bool)
            lm.values = morphology.remove_small_holes(lm.values, area_threshold=9, connectivity=2)
            lm.values = morphology.binary_erosion(lm.values, footprint=np.ones((2, 2)))
            lm.to_netcdf(lm_file)
        lm = xarray.open_dataset(lm_file).LANDMASK

        # Open hist and future simulations.
        if (
            not os.path.exists(maxes_file)
            or not os.path.exists(h_maxima_file)
            or not os.path.exists(f_maxima_file)
        ):
            hist_files = sorted(glob.glob(f'{sim_dir}/hist/*/WRF/basic*d0{d}*.nc'))
            fut_files = sorted(glob.glob(f'{sim_dir}/ssp245/*/WRF/basic*d0{d}*.nc'))

            assert len(hist_files) == 3040, 'Missing files for historical period.'
            assert len(fut_files) == 3040, 'Missing files for future period.'

            hist = xarray.open_mfdataset(hist_files, parallel=True)
            futu = xarray.open_mfdataset(fut_files, parallel=True)

            hist = hist.chunk({'time': -1, 'south_north': 50, 'west_east': 50})
            futu = futu.chunk({'time': -1, 'south_north': 50, 'west_east': 50})

            # Adjust to local time.
            hist = hist.assign_coords(
                {'time': hist.time + np.timedelta64(time_adjust[domain], 'h')}
            )
            futu = futu.assign_coords(
                {'time': futu.time + np.timedelta64(time_adjust[domain], 'h')}
            )

            # Remove spin up time and trim so timeseries are the same length.
            hist = hist.where(hist.time.dt.month != 9, drop=True)
            hist = hist.where(hist.time.dt.month != 3, drop=True)
            futu = futu.where(futu.time.dt.month != 9, drop=True)
            futu = futu.where(futu.time.dt.month != 3, drop=True)
            assert np.all(hist.time == futu.time)

            # Subset to only required variables.
            lats = hist.isel(time=0).latitude
            lons = hist.isel(time=0).longitude
            hist = hist[variables]
            futu = futu[variables]

            # Subset to only those points where surface hail was simulated.
            hist = hist.where(hist.hailcast_diam_max > 0)
            futu = futu.where(futu.hailcast_diam_max > 0)

            # If getting wind speed, don't worry about wind direction.
            if 'wind_10m' in variables:
                hist['wind_10m'] = hist.wind_10m.sel(wspd_wdir='wspd')
                futu['wind_10m'] = futu.wind_10m.sel(wspd_wdir='wspd')
                hist = hist.drop_vars('wspd_wdir').reset_coords()
                futu = futu.drop_vars('wspd_wdir').reset_coords()

            # Subset to land points only.
            hist_land = hist.where(lm == 1).reset_coords(drop=True)
            futu_land = futu.where(lm == 1).reset_coords(drop=True)

        # Calculate maximums over two epochs for full domains.
        if not os.path.exists(maxes_file):
            h_max = hist.max('time').expand_dims({'epoch': ['historical']})
            f_max = futu.max('time').expand_dims({'epoch': ['ssp245']})
            maxes = xarray.merge([h_max, f_max])
            maxes['latitude'] = lats
            maxes['longitude'] = lons
            maxes.to_netcdf(maxes_file)
        maxes = xarray.open_dataset(maxes_file)

        if plot:
            if 'hailcast_diam_max' in variables:
                plot_maxes(
                    maxes=maxes,
                    variable='hailcast_diam_max',
                    lab='Max hail diam [mm]',
                    lm=lm,
                    domain=domain,
                    file_dir=file_dir
                )
            if 'wind_10m' in variables:
                plot_maxes(
                    maxes=maxes, variable='wind_10m', 
                    lab='Max wind [m/s]', lm=lm, domain=domain
                )

        # Output block maxima for land points only.
        if not os.path.exists(h_maxima_file):
            h_maxima = hist_land.max(['south_north', 'west_east']).resample(time='1D').max()
            h_maxima = h_maxima.expand_dims(
                {'domain': [domain.replace('_', ' + ')], 'epoch': ['historic']}
            )
            h_maxima.to_dataframe().dropna(how='all').to_feather(h_maxima_file)

        if not os.path.exists(f_maxima_file):
            f_maxima = futu_land.max(['south_north', 'west_east']).resample(time='1D').max()
            f_maxima = f_maxima.expand_dims(
                {'domain': [domain.replace('_', ' + ')], 'epoch': ['ssp245']}
            )
            f_maxima.to_dataframe().dropna(how='all').to_feather(f_maxima_file)
