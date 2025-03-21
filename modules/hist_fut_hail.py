"""Module providing support functions for historic vs future hail comparisons."""

import glob
import itertools
import os

import cartopy.crs as ccrs
import geopandas
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import xarray
from cartopy.io import shapereader
from matplotlib import colors
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from shapely.geometry import Polygon
from shapely.geometry.polygon import LinearRing
from skimage import morphology

# Coordinates for cities of interest.
cities = {
    'Perth': (115.8606, -31.9559),
    'Sydney': (151.21, -33.867778),
    'Melbourne': (144.963056, -37.814167),
    'Brisbane': (153.028056, -27.467778),
    'Canberra': (149.126944, -35.293056),
    'Adelaide': (138.5999, -34.9287),
    'Kalgoorlie': (121.4656, -30.7582),
}

# Letters for plot labels.
letters = [
    'a',
    'b',
    'c',
    'd',
    'e',
    'f',
    'g',
    'h',
    'i',
    'j',
    'k',
    'l',
    'm',
    'n',
    'o',
    'p',
    'q',
    'r',
    's',
    't',
    'u',
    'v',
    'w',
    'x',
    'y',
    'z',
]


def gen_download_script(years, out_file, fut_ssp='ssp245', link_list='data/xu_file_list.csv'):
    """Write a script to wget all the boundary condition files required for a run over given years.

    Arguments:
        years: Years to run over. For year x data will be downloaded for Sep x to Feb x+1.
        out_file: The script file to write.
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
        files.loc[files.file.str.contains(f'atm_hist_{season_from}_(09|10|11|12)'), 'season'] = season_string
        files.loc[files.file.str.contains(f'atm_hist_{season_to}_(01|02)'), 'season'] = season_string
        files.loc[
            files.file.str.contains(f'atm_{fut_ssp}_{season_from}_(09|10|11|12)'),
            'season',
        ] = season_string
        files.loc[files.file.str.contains(f'atm_{fut_ssp}_{season_to}_(01|02)'), 'season'] = season_string
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
    """Set up directories ready for WPS and WRF runs for a given season.

    Including updating namelist files.

    Arguments:
        year: Event year (includes Sep year to Feb year+1).
        template_dir: Directory with template WPS and WRF setups.
        namelist_dir: Directory from which to source WPS and WRF namelists.
        sims_dir: Directory in which simulations will be run.
        exp: Experiment to run (hist or ssp245).

    """
    sim_dir = f'{sims_dir}/{year}-{year + 1}/'
    if not os.path.exists(sim_dir):
        os.mkdir(sim_dir)

    start_time = f'{year}-09-30_00:00:00'
    end_time = f'{year + 1}-02-28_23:00:00'
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
            f'sed -i s/start_date.*$/"start_date = \'{comma_q.join(np.repeat(start_time, 7))}\',"/g {sim_dir}/WPS/namelist.wps',
        )
        os.system(
            f'sed -i s/end_date.*$/"end_date = \'{comma_q.join(np.repeat(end_time, 7))}\',"/g {sim_dir}/WPS/namelist.wps',
        )
    else:
        print('Skipping existing WPS...')

    # WRF setup. Link executables + data files, copy and update namelist.
    if not os.path.exists(f'{sim_dir}/WRF/'):
        os.mkdir(f'{sim_dir}/WRF')
        os.system(f'ln -s {template_dir}/WRF/* {sim_dir}/WRF/')
        os.system(f'cp -f {namelist_dir}/namelist.input {sim_dir}/WRF/')
        os.system(
            f'sed -i s/start_year.*$/"start_year = {comma.join(np.repeat(start_time[0:4], 7))},"/g {sim_dir}/WRF/namelist.input',
        )
        os.system(
            f'sed -i s/start_month.*$/"start_month = {comma.join(np.repeat(start_time[5:7], 7))},"/g {sim_dir}/WRF/namelist.input',
        )
        os.system(
            f'sed -i s/start_day.*$/"start_day = {comma.join(np.repeat(start_time[8:10], 7))},"/g {sim_dir}/WRF/namelist.input',
        )
        os.system(
            f'sed -i s/start_hour.*$/"start_hour = {comma.join(np.repeat(start_time[11:13], 7))},"/g {sim_dir}/WRF/namelist.input',
        )
        os.system(
            f'sed -i s/end_year.*$/"end_year = {comma.join(np.repeat(end_time[0:4], 7))},"/g {sim_dir}/WRF/namelist.input',
        )
        os.system(
            f'sed -i s/end_month.*$/"end_month = {comma.join(np.repeat(end_time[5:7], 7))},"/g {sim_dir}/WRF/namelist.input',
        )
        os.system(
            f'sed -i s/end_day.*$/"end_day = {comma.join(np.repeat(end_time[8:10], 7))},"/g {sim_dir}/WRF/namelist.input',
        )
        os.system(
            f'sed -i s/end_hour.*$/"end_hour = {comma.join(np.repeat(end_time[11:13], 7))},"/g {sim_dir}/WRF/namelist.input',
        )
    else:
        print('Skipping existing WRF...')


def plot_wrf_domains(
    wps_files,
    pts=None,
    figsize=(10, 8),
    proj=ccrs.PlateCarree(),
    labels=None,
    file=None,
):
    """Plot WRF domains for a WPS setup. Requires geogrid to have been run first.

    Arguments:
        wps_files: A dictionary of WPS files and the (linestyle, colour) they should be plotted as.
        pts: Points to add as a dictionary with keys a name and values an (x,y) pair.
        figsize: The figure size (width x height).
        proj: Projection to display map in.
        labels: Point name: (label, offset) for each point to display, with offset in map coordinates.
        file: Output file for figure.

    """
    dom_files = list(wps_files.keys())
    doms = [xarray.open_dataset(x) for x in dom_files]

    max_map_fac = np.array([np.max(np.abs(x.MAPFAC_M.values - 1) / 1 * 100) for x in doms])
    if np.any(max_map_fac > 5):  # noqa: PLR2004
        print('WARNING: Map factor is large or small (below 0.95 or over 1.05).')
        facts = [np.max(x.MAPFAC_M.values) for x in doms]
        print(f'Max factor is {np.max(facts):.3} in domain {np.argmax(facts) + 1}.')
        print(f'Min factor is {np.min(facts):.3} in domain {np.argmin(facts) + 1}.')

    # Custom functions for WRF domain plotting.
    def add_box(x, y, ax, colour='white', linestyle='solid'):
        geom = LinearRing(
            list(
                zip(
                    list(x.isel(south_north_stag=-1).values) + list(x.isel(south_north_stag=0).values)[::-1],
                    list(y.isel(south_north_stag=-1).values) + list(y.isel(south_north_stag=0).values)[::-1],
                ),
            ),
        )
        ax.add_geometries(
            [geom],
            crs=ccrs.PlateCarree(),
            linewidth=2,
            edgecolor=colour,
            facecolor='none',
            linestyle=linestyle,
        )

    def add_doms(dom_list, ax, colour='white', linestyle='solid'):
        for dom in dom_list:
            x = dom.XLONG_C.isel(Time=0)
            y = dom.XLAT_C.isel(Time=0)
            add_box(x, y, ax, colour=colour, linestyle=linestyle)

    # Do the plotting.
    _, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': proj})
    ax.coastlines(linewidth=1.5)

    min_x = min_y = np.inf
    max_x = max_y = -np.inf

    for i, (ls, col) in enumerate([wps_files[x] for x in wps_files]):
        add_doms([doms[i]], ax=ax, colour=col, linestyle=ls)
        x = doms[i].isel(Time=0).XLONG_C.values
        y = doms[i].isel(Time=0).XLAT_C.values
        min_x = min(min_x, np.min(x))
        max_x = max(max_x, np.max(x))
        min_y = min(min_y, np.min(y))
        max_y = max(max_y, np.max(y))

    ax.set_xlim(min_x - 2, max_x + 2)
    ax.set_ylim(min_y - 2, max_y + 2)

    gl = ax.gridlines(crs=proj, draw_labels=True, alpha=0.3)
    gl.top_labels = gl.right_labels = False

    # Add optional points.
    for i in pts:
        ax.scatter(pts[i][0], pts[i][1], color='red', transform=ccrs.PlateCarree())

        if labels is not None and i in labels:
            ax.annotate(
                labels[i][0],
                (pts[i][0], pts[i][1]),
                textcoords='offset points',
                xytext=labels[i][1],
                ha='center',
                color='black',
                fontweight='bold',
                fontsize=plt.rcParams['font.size'] + 1,
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
    """Plot a WRF domain definition.

    Assuming the projection is 'lat-lon' (ie rotated cylindrical equidistant) in WRF.

    Arguments:
        parent_id, i_parent_start, j_parent_start, ref_lon,
        ref_lat, e_we, e_sn, dx, dy, num_doms, scale_factor: WRF parameters from namelist.wps.
        pts: Points to plot as dots on the map.
        dx: X resolution.
        dy: Y resolution.
        e_sn: Pixels in S-N direction.
        e_we: Pixels in W-E direction.
        i_parent_start: Starting pixel inside parent domain (i).
        j_parent_start: Starting pixel inside parent domain (j).
        num_doms: Number of domains.
        parent_id: ID of each parent domain.
        ref_lat: Reference (centre) lat.
        ref_lon: Reference (centre) lon.
        scale_factor: Scale factor(s) between parent and child domains.

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
        x=ref_lon_orig,
        y=ref_lat_orig,
        src_crs=ccrs.Geodetic(),
    )

    # Add optional points.
    if pts is not None:
        for i in pts:
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
    infiles=None,
    outfile='/g/data/up6/tr2908/hist_future_hail/xu_data/ssp245/mean_temperature_by_year.nc',
):
    """Calculate mean near-surface temperatue from Xu et al lnd* files.

    Write a timeseries per year to a netcdf file.

    Arguments:
        infiles: Dictionary with name: path_string for the lnd* files to use.
        outfile: The output file to write.

    """
    if infiles is None:
        infiles = {
            'SSP2-4.5': '/scratch/w42/tr2908/xu_data/lnd.ssp245.*.nc',
            'SSP5-8.5': '/scratch/w42/tr2908/xu_data/lnd.ssp585.*.nc',
        }
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
            mean_temp = d.weighted_temp.groupby(d.year).sum(['lat', 'lon', 'time']) / d.area_weight.groupby(d.year).sum(['lat', 'lon', 'time'])
            mean_temp.attrs['description'] = (
                ("Average near-surface temperature calculated using d.tas.groupby(d.year).mean(['lat', 'lon', 'time'])",
                 f'where d is all Xu et al. {sim} data.'),
            )
            mean_temp.name = 'tas'
            mean_temp = mean_temp.expand_dims({'simulation': [sim]})
            mean_temps.append(mean_temp)

        mean_temps = xarray.merge(mean_temps)
        mean_temps.to_netcdf(outfile)

    return xarray.open_dataset(outfile)


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
    """Plot data on a map to a specified plot axis object.

    Arguments:
        dat: DataSet to plot or list of datasets to plot.
        ax: GeoAxes object to plot to.
        dat_proj, dist_proj: Data and display projections.
        figsize: Figure size width x height.
        coastlines: Show coastlines?
        grid: Show grid?
        ncol/nrows: Number of columns/rows to plot.
        title: Title for the plot.
        colour_scale: None for default, or a tuple of min/max values for the scale.
        cmap: The matplotlib colour map to use.
        norm: A norm object for colours (e.g. colors.BoundaryNorm).
        cbar_ticks: Colour bar ticks.
        tick_labels: Colour bar tick labels.
        contour: Plot using xarray's contourf function?
        stippling: True where trippling should appear.
        stipple_size: Size for stippling points.
        colourbar: Include a colourbar?
        ticks_left: Include ticks on left of plot?
        ticks_bottom: Include ticks on bottom of plot?
        cbar_aspect: colorbar aspect ratio?
        cbar_fraction: fraction argument for colorbar().
        cbar_shrink: shrink argument for colorbar().
        cbar_pad: pad argument for colorbar().
        cbar_label: Overwrite label?
        cbar_orientation: orientation argument for colorbar().
        coastlines_colour: Colour for coastlines.
        xlims: x plot limit.
        ylims: y plot limits.
        num_ticks: Number of ticks for x and y axes (None for auto).
        divergent: Is the colour scale divergent? If so make zero central.
        cbar_inset: Inset the colorbar in lower left?
        title_inset: Inset the title in the upper left?
        discrete: Make the colour bar discrete?
        log_scale: Make the colour scale log-scaled?
        nan_colour: Colour for missing values.
        axis_off: Turn off all axes.
        country: Plot coastlines only for a specific country.
        annotations: Add annotations to the map - dictionary of {'Text': [x, y, xadj, yadj, ha]} where
                     x, y are position of text, xadj, yadj give offsets to the label, ha is 'left' or
                     'right' for horizontal anchor.
        polygons: If specified, draw each polygon onto the plot.
        polygon_colour: Colour for polygons.
        hatch: Hatching to use for stippling ('.' for scatterplot, '//' etc for contour hatching).
        dat_proj: The data projection.
        disp_proj: Projection to plot in.
        left_title: Put titles to left instead of centre?
        pts: Extra points to plot.

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
        cbar_ticks = np.array([*cbar_ticks, np.max(cbar_ticks) + 1])
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
        ax.autoscale(False)  # noqa: FBT003
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
            stippling.plot.contourf(hatches=['', hatch], levels=[0, 0.5, 1], colors='none', ax=ax, add_colorbar=False)

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
        ax.add_geometries(poly, crs=ccrs.PlateCarree(), facecolor='none', edgecolor=polygon_colour, linewidth=1.75)
    if coastlines:
        if country is not None:
            shpfilename = shapereader.natural_earth(resolution='10m', category='cultural', name='admin_0_countries')
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
        for i in pts:
            ax.scatter(pts[i][0], pts[i][1], color='black')

    if annotations is not None:
        for text, [x, y, xadj, yadj, ha] in annotations.items():
            if np.abs(xadj) >= 1 or np.abs(yadj) >= 1:
                if ha in ('right', 'left'):
                    ax.plot([x, x + xadj - (0.2 * np.sign(xadj))], [y, y + yadj + 0.2], color='black')
                elif ha == 'center':
                    ax.plot([x, x + xadj - (0.2 * np.sign(xadj))], [y, y + yadj - 0.2], color='black')
                    if yadj < 0:
                        print('Warning: ha=center and negative y adjustment are not supported.')
                else:
                    assert 1 == 0, 'Invalid value of ha.'  # noqa: PLR0133
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
    """Plot data on a map.

    Arguments:
        dat: DataSet to plot or list of datasets to plot.
        dat_proj, dist_proj: Data and display projections.
        figsize: Figure size width x height.
        grid: Show grid?
        ncols/nrows: Number of columns/rows to plot.
        title: Title(s) for the plot(s).
        share_scale: Make the range of values in each plot the same?
        colour_scale: Tuple with min/max values to use on scale. Overwritten by share_scale.
        cbar_ticks: Ticks for the colourbar.
        tick_labels: Colour bar tick labels.
        file: If specified save to 'file' instead of showing onscreen.
        scale_label: The label for a shared scale.
        share_axes: Share left/bottom axes?
        ticks_left, ticks_bottom: Display ticks on the left/bottom of plots?
        wspace, hspace: gridspec wspace and hspace arguments.
        stippling: Stippling per axis.
        cbar_adjust: Amount to shrink plots by to add cbar for shared scale.
        cbar_pad: Padding between figure and colour bar.
        col_labels/row_labels: Labels for each column/row; overwrites individial plot titles.
        xlims: x limits.
        ylims: y limits.
        show: Show the map?
        shared_scale_quantiles: Quantiles for a shared scale.
        polygon: If specified, polygons to put on each map.
        letter_labels: Use a letter to label each subplot?
        kwargs: Extra arguments to plot_map_to_ax.
        col_labels: Labels for each column.
        dat_proj: Data projection.
        disp_proj: Display projection.
        hspace: Height spacing.
        wspace: Width spacing.
        ncols: Number of columns to use.
        nrows: Number of rows.
        polygons: Polygons to overplot.
        row_labels: Labels for each row.
        ticks_bottom: Include bottom ticks?
        ticks_left: Include left ticks?

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
            assert not (np.isnan(colour_scale[0]) or np.isnan(colour_scale[1])), 'share_scale cannot be used with subplots missing data.'

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
            _ = fig.colorbar(im, ax=ax, cax=cbar_ax, ticks=cbar_ticks, label=scale_label, format=fmt)

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
                    lab_ax.annotate(lab, xy=(0.5, 1 - p), rotation=90, xycoords='axes fraction', ha='center')

    if file is not None:
        plt.savefig(fname=file, dpi=300, bbox_inches='tight')

        if show:
            plt.show()
        plt.close()
    elif show:
        plt.show()

    return ax


def plot_maxes(
    all_maxes,
    subset_maxes,
    domain,
    figsize=(12, 2),
    variable='hailcast_diam_max',
    lab='Max hail diam [mm]',
    file_dir=None,
):
    """Plot maximum values by location for historical and future epochs both nonmasked and masked.

    Arguments:
        all_maxes: Maximums per location per epoch.
        subset_maxes: As for maxes but subset to study regions.
        domain: The domain to use in titles.
        figsize: Figure size.
        variable: Variable to plot.
        lab: Label for the variable.
        file_dir: Optional directory to save plot to.

    """
    x = all_maxes.longitude.values
    y = all_maxes.latitude.values

    plt.rcParams['font.size'] = 12
    _, axs = plt.subplots(
        ncols=4,
        figsize=figsize,
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace': 0.3},
    )
    pc_hist = axs[0].pcolormesh(
        x,
        y,
        all_maxes[variable].sel(epoch='historical'),
        cmap='Spectral_r',
        transform=ccrs.PlateCarree(),
        rasterized=True,
    )
    pc_futu = axs[1].pcolormesh(
        x,
        y,
        all_maxes[variable].sel(epoch='ssp245'),
        cmap='Spectral_r',
        transform=ccrs.PlateCarree(),
        rasterized=True,
    )
    pc_hist_land = axs[2].pcolormesh(
        x,
        y,
        subset_maxes[variable].sel(epoch='historical'),
        cmap='Spectral_r',
        transform=ccrs.PlateCarree(),
        rasterized=True,
    )
    pc_futu_land = axs[3].pcolormesh(
        x,
        y,
        subset_maxes[variable].sel(epoch='ssp245'),
        cmap='Spectral_r',
        transform=ccrs.PlateCarree(),
        rasterized=True,
    )
    plt.colorbar(pc_hist, label=lab)
    plt.colorbar(pc_futu, label=lab)
    plt.colorbar(pc_hist_land, label=lab)
    plt.colorbar(pc_futu_land, label=lab)

    for ax in axs:
        ax.coastlines(linewidth=1.5)

    domain_str = domain.replace('_', '/')
    axs[0].set_title(f'{domain_str}\nHistorical')
    axs[1].set_title(f'{domain_str}\nFuture')
    axs[2].set_title(f'{domain_str}\nHistorical')
    axs[3].set_title(f'{domain_str}\nFuture')

    if file_dir is not None:
        plt.savefig(fname=f'{file_dir}/maxes_{domain}_{variable}.pdf', dpi=300, bbox_inches='tight')


def process_maxima_set(
    dat,
    lm,
    time_adjust,
    epoch,
    domain,
    time_adjust_mins=0,
    drop_months=None,
    max_hailsize=180,
):
    """Process block (daily) maxima for a dataset.

    Args:
        dat: The raw data.
        lm: Landmask for this dataset.
        time_adjust: Adjustment in hours to make to UTC times.
        epoch: Epoch descriptor to add to outputs.
        domain: The domain descriptor.
        time_adjust_mins: Additional minutes to add to the time.
        drop_months: Months to ignore. Defaults to [9, 3].
        max_hailsize: Maximum hail size to allow [mm].

    """
    if drop_months is None:
        drop_months = [9, 3]

    dat = dat.chunk({'time': 3000, 'south_north': -1, 'west_east': -1})

    # Adjust to local time.
    dat = dat.assign_coords({'time': dat.time + np.timedelta64(time_adjust, 'h') + np.timedelta64(time_adjust_mins, 'm')})

    # Remove spin up time and trim so all timeseries are the same length.
    # Note adjustment to local time can create leap year days (Feb 29s) which
    # are not removed here, so must be removed later before analysis is done.
    # In R code these are removed.
    for m in drop_months:
        dat = dat.where(dat.time.dt.month != m, drop=True)

    # Subset to only required variables.
    lats = dat.isel(time=0).latitude
    lons = dat.isel(time=0).longitude

    # Subset to only those points where surface hail was simulated.
    dat = dat.where(dat.hailcast_diam_max > 0)

    # Don't worry about wind direction.
    dat['wind_10m'] = dat.wind_10m.sel(wspd_wdir='wspd')
    dat = dat.drop_vars('wspd_wdir').reset_coords()

    # Subset to land points only.
    dat_land = dat.where(lm == 1).reset_coords(drop=True)

    # Remove too-large hail.
    num_large_hail = (dat_land.hailcast_diam_max > max_hailsize).sum().values
    num_any_hail = (dat_land.hailcast_diam_max > 0).sum().values
    perc_large_hail = num_large_hail / num_any_hail * 100
    dat_land = dat_land.where(dat_land.hailcast_diam_max <= max_hailsize)

    maxes_all = dat.max('time').expand_dims({'epoch': [epoch]})
    maxes_land = dat_land.max('time')

    # Save percentage of hail that is too large and removed.
    maxes_land['perc_large_hail'] = perc_large_hail
    maxes_land.perc_large_hail.attrs['description'] = f'Percentage of hail removed because it is over {max_hailsize} mm.'

    maxes_land = maxes_land.expand_dims({'epoch': [epoch]})
    for i in [maxes_all, maxes_land]:
        i['latitude'] = (('south_north', 'west_east'), lats.values)
        i['longitude'] = (('south_north', 'west_east'), lons.values)

    # Persist for speed.
    dat_land = dat_land.chunk({'time': -1, 'south_north': 10, 'west_east': 10}).persist()

    # Calculate daily maxima for land points only.
    maxima = dat_land.resample(time='1D').max(dim=['time', 'south_north', 'west_east'])
    maxima = maxima.expand_dims({'domain': [domain.replace('_', ' + ')], 'epoch': [epoch]})
    maxima = maxima.to_dataframe().dropna(how='all')

    # Calculate daily mean for land points only.
    means = dat_land.resample(time='1D').mean(dim=['time', 'south_north', 'west_east'])
    means = means.expand_dims({'domain': [domain.replace('_', ' + ')], 'epoch': [epoch]})
    means = means.to_dataframe().dropna(how='all')

    return maxima, means, maxes_all, maxes_land


def process_maxima(
    sim_dir,
    domains=None,
    time_adjust=None,
    time_adjust_mins=None,
    results_dir='/g/data/up6/tr2908/hist_future_hail/results/',
    variables=None,
    file_dir=None,
    drop_vars_basic=None,
    drop_vars_conv=None,
    **kwargs,
):
    """Process maxima of hail and wind values.

    Args:
        sim_dir: Directory where simulations are stored. Defaults to '/g/data/up6/tr2908/hist_future_hail/WRF_v4.4/simulations/cities/'.
        domains: Domain IDs in sims. Defaults to {'Perth': 3, 'Melbourne': 5, 'Brisbane': 6, 'Sydney_Canberra': 7}.
        time_adjust: Time adjustments to get to local time. Defaults to {'Perth': 8, 'Melbourne': 11, 'Brisbane': 10, 'Sydney_Canberra': 11}.
        results_dir: Where to write results. Defaults to '/g/data/up6/tr2908/hist_future_hail/results/'.
        variables: Variables to process. Defaults to ['hailcast_diam_max', 'wind_10m'].
        file_dir: Figure directory. Defaults to 'paper/figures/'.
        drop_vars_basic: Basic variables to not include.
        drop_vars_conv: Conv variables to not include.
        time_adjust_mins: Adjust time if required, per domain.
        kwargs: Extra arguments to process_maxima_set().

    """
    if domains is None:
        domains = {'Perth': 3, 'Melbourne': 5, 'Brisbane': 6, 'Sydney_Canberra': 7}
    if time_adjust is None:
        time_adjust = {'Perth': 8, 'Melbourne': 11, 'Brisbane': 10, 'Sydney_Canberra': 11}
    if time_adjust_mins is None:
        time_adjust_mins = {'Perth': 0, 'Melbourne': 0, 'Brisbane': 0, 'Sydney_Canberra': 0}
    if variables is None:
        variables = ['hailcast_diam_max', 'wind_10m']
    if drop_vars_basic is None:
        drop_vars_basic = ['pressure', 'temperature', 'u', 'v', 'z', 'z_agl', 'mixing_ratio', 'specific_humidity', 'bottom_top']
    if drop_vars_conv is None:
        drop_vars_conv = ['shear_u', 'shear_v', 'positive_shear']

    maxima = []

    for domain, d in domains.items():
        print(f'Processing d0{d}...')

        # Define cache files.
        dmn = domain.replace(' + ', '_')
        lm_file = f'{results_dir}/{dmn}_landmask.nc'
        all_maxes_file = f'{results_dir}/{dmn}_domain_maximums_all.nc'
        land_maxes_file = f'{results_dir}/{dmn}_domain_maximums_subset.nc'
        h_maxima_file = f'{results_dir}/{dmn}_hist_block_maxima.feather'
        f_maxima_file = f'{results_dir}/{dmn}_ssp245_block_maxima.feather'
        h_means_file = f'{results_dir}/{dmn}_hist_block_means.feather'
        f_means_file = f'{results_dir}/{dmn}_ssp245_block_means.feather'

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
            not os.path.exists(all_maxes_file)
            or not os.path.exists(land_maxes_file)
            or not os.path.exists(h_maxima_file)
            or not os.path.exists(f_maxima_file)
            or not os.path.exists(h_means_file)
            or not os.path.exists(f_means_file)
        ):
            hist_basic = open_set(pattern=f'{sim_dir}/hist/*/WRF/basic*d0{d}*.nc', drop_vars=drop_vars_basic)
            futu_basic = open_set(pattern=f'{sim_dir}/ssp245/*/WRF/basic*d0{d}*.nc', drop_vars=drop_vars_basic)
            hist_conv = open_set(pattern=f'{sim_dir}/hist/*/WRF/conv*d0{d}*.nc', drop_vars=drop_vars_conv)
            futu_conv = open_set(pattern=f'{sim_dir}/ssp245/*/WRF/conv*d0{d}*.nc', drop_vars=drop_vars_conv)

            hist = xarray.merge([hist_basic, hist_conv])
            futu = xarray.merge([futu_basic, futu_conv])

            hist_maxima, hist_means, hist_maxes_all, hist_maxes_land = process_maxima_set(
                dat=hist,
                lm=lm,
                time_adjust=time_adjust[domain],
                epoch='historical',
                domain=domain,
                time_adjust_mins=time_adjust_mins[domain],
                **kwargs,
            )

            futu_maxima, futu_means, futu_maxes_all, futu_maxes_land = process_maxima_set(
                dat=futu,
                lm=lm,
                time_adjust=time_adjust[domain],
                epoch='ssp245',
                domain=domain,
                time_adjust_mins=time_adjust_mins[domain],
                **kwargs,
            )

            hist_maxima.to_feather(h_maxima_file)
            futu_maxima.to_feather(f_maxima_file)
            hist_means.to_feather(h_means_file)
            futu_means.to_feather(f_means_file)
            del hist_maxima, futu_maxima

            all_maxes = xarray.merge([hist_maxes_all, futu_maxes_all])
            all_maxes.to_netcdf(all_maxes_file)
            del all_maxes

            land_maxes = xarray.merge([hist_maxes_land, futu_maxes_land])
            land_maxes.to_netcdf(land_maxes_file)
            del land_maxes

        all_maxes = xarray.open_dataset(all_maxes_file)
        land_maxes = xarray.open_dataset(land_maxes_file)

        if 'hailcast_diam_max' in variables:
            plot_maxes(
                all_maxes=all_maxes,
                subset_maxes=land_maxes,
                variable='hailcast_diam_max',
                lab='Max hail diam [mm]',
                domain=domain,
                file_dir=file_dir,
            )
        if 'wind_10m' in variables:
            plot_maxes(
                all_maxes=all_maxes,
                subset_maxes=land_maxes,
                variable='wind_10m',
                lab='Max wind [m/s]',
                domain=domain,
            )

        maxima.append(land_maxes.expand_dims({'domain': [domain.replace('_', '/')]}))

    return xarray.merge(maxima)


def plot_maxima(
    maxima,
    variable,
    scale_label,
    factor=1,
    cbar_max=None,
    cbar_min=None,
    file=None,
    figsize=(12, 8.2),
    cbar_adjust=0.862,
    cbar_pad=0.015,
    nrows=3,
    ncols=4,
    locator_base=None,
    title_xs=None,
    title_ys=None,
    city_polygons_file='data/SUA_2021_AUST_GDA2020_SHP.zip',
    polygons_colour='fuchsia',
    start_letter=0,
):
    """Plot maxima for each epoch.

    Args:
        maxima: The maxima to plot.
        variable: Variable to plot.
        scale_label: Label for the shared scale bar.
        factor: Factor by which to multiply values (for unit conversion).
        cbar_max: Maximum value to plot, None if data max.
        cbar_min: Minimum value to plot, None if data min.
        file: Output file for plot.
        figsize: Figure suze. 
        cbar_adjust: cbar adjustment factor. Defaults to 0.862.
        cbar_pad: cbar padding factor. Defaults to 0.015.
        nrows: Number of rows.
        ncols: Number of columns.
        locator_base: 'base' argument for locator, per domain, to override default of 2.
        title_xs: x position for inset title, to override default per domain.
        title_ys: y position for inset title, to override default per domain.
        city_polygons_file: File containin urban regions polygons to include.
        polygons_colour: The colour to use to plot polygons.
        start_letter: Letter to start from (numeric, 1=a).

    """
    if locator_base is None:
        locator_base = {'Sydney/Canberra': 1.5, 'Brisbane': 1.5}
    if title_xs is None:
        title_xs = {'Kalgoorlie': 0.02, 'Adelaide': 0.05}
    if title_ys is None:
        title_ys = {'Kalgoorlie': 0.9, 'Melbourne': 0.9}

    maxima = maxima.copy(deep=True)
    maxima[variable] = maxima[variable] * factor

    if city_polygons_file is not None:
        cities = geopandas.read_file(city_polygons_file)
        cities = cities[['Not in any' not in x for x in cities.SUA_NAME21]]

    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols + 1,
        figsize=figsize,
        subplot_kw={'projection': ccrs.PlateCarree()},
        gridspec_kw={'wspace': 0.1, 'hspace': 0.3, 'width_ratios': [1, 1, 0.22, 1, 1]},
    )

    zmin = maxima[variable].min()
    zmax = maxima[variable].max()

    plot_letters = letters.copy()[start_letter:][::-1]
    pnum = 0
    for i, (d, e) in enumerate(itertools.product(maxima.domain, maxima.epoch)):
        x = maxima.sel(epoch=e, domain=d).longitude
        y = maxima.sel(epoch=e, domain=d).latitude
        z = maxima.sel(epoch=e, domain=d)[variable]

        x = x.dropna('west_east', how='all').dropna('south_north', how='all')
        y = y.dropna('west_east', how='all').dropna('south_north', how='all')
        z = z.sel(south_north=x.south_north, west_east=x.west_east)

        m = xarray.DataArray(
            z.values,
            dims=['y', 'x'],
            coords={'lat': (('y', 'x'), y.values), 'lon': (('y', 'x'), x.values)},
        )

        im = m.plot(
            ax=axs.flat[pnum],
            transform=ccrs.PlateCarree(),
            add_colorbar=False,
            x='lon',
            y='lat',
            cmap='Spectral_r',
            vmin=cbar_min if cbar_min is not None else zmin,
            vmax=cbar_max if cbar_max is not None else zmax,
            rasterized=True,
        )

        axs.flat[pnum].coastlines()
        dom = d.values
        if dom == "Sydney/Canberra":
            dom = "Syd./Canb."
        axs.flat[pnum].set_title(f'{dom} ({"H" if e.values == "historical" else "F"})')

        if city_polygons_file is not None:
            bbox = Polygon([(np.min(x), np.min(y)), (np.min(x), np.max(y)), (np.max(x), np.max(y)), (np.max(x), np.min(y))])
            city_polygons = cities.clip(bbox).geometry
            city_polygons = geopandas.GeoSeries(city_polygons).unary_union

            axs.flat[pnum].add_geometries(city_polygons, crs=ccrs.PlateCarree(), facecolor='none', edgecolor=polygons_colour, linewidth=0.75)

        title_x = title_xs.get(str(d.values), 0.9)
        title_y = title_ys.get(str(d.values), 0.05)
        axs.flat[pnum].annotate(
            text=f'{plot_letters.pop()}',
            xy=(title_x, title_y),
            xycoords='axes fraction',
            fontweight='bold',
            fontsize=plt.rcParams['font.size'],
        )
        base = locator_base.get(str(d.values), 2)
        locator = mticker.MultipleLocator(base=base)
        gl = axs.flat[pnum].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, xlocs=locator, ylocs=locator)
        gl.top_labels = gl.right_labels = False
        if (i + 1) % int(ncols / 2) == 0:
            gl.left_labels = False
            if (i + 1) % int(ncols) != 0:
                pnum = pnum + 1
                axs.flat[pnum].set_visible(False)

        pnum = pnum + 1

    fig.subplots_adjust(right=cbar_adjust)
    cbar_ax = fig.add_axes([cbar_adjust + cbar_pad, 0.23, 0.02, 0.55])
    fmt = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    fmt.set_powerlimits((-4, 6))

    extend = 'neither'
    if cbar_max is not None and cbar_max < zmax:
        extend = 'max'
    if cbar_min is not None and cbar_min > zmin:
        extend = 'min'

    _ = fig.colorbar(im, ax=axs.flat[0], cax=cbar_ax, ticks=None, label=scale_label, format=fmt, extend=extend)

    if file is not None:
        plt.savefig(fname=file, dpi=300, bbox_inches='tight')


def open_set(pattern, expected_n=3040, drop_vars=None):
    """Open a dataset.

    Args:
        pattern: The pattern to match files for.
        expected_n: Check there are n files being opened. Defaults to 3040.
        drop_vars: Optionally drop variables.

    """
    files = sorted(glob.glob(pattern))
    assert len(files) == expected_n, f'Missing files for {pattern}'
    dat = xarray.open_mfdataset(files, parallel=True, combine='nested', concat_dim='time')
    if drop_vars is not None:
        dat = dat.drop_vars(drop_vars)
    return dat
