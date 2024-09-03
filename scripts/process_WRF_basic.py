## Use wrf-python to extract basic properties from wrfout_d03 files in the current directory
## and write output files basic_properties_*.nc

import sys
import os

import wrf
import dask
import metpy
from metpy import calc
import xarray
import numpy as np
import pandas as pd
from glob import glob
from netCDF4 import Dataset

# Match d03 (Perth), d05 (Melbourne), d06 (Brisbane), d07 (Syd + Canberra)
files = sorted(glob('wrfout_d0[3,5,6,7]*'))
for filename in files:
    print(filename)
    outfile = filename.replace('wrfout', 'basic_params') + '.nc'

    if os.path.exists(outfile):
        print(f'skipping {outfile}...')
        continue

    nc = Dataset(filename)
    dat = xarray.Dataset({'pressure': wrf.getvar(nc, 'pressure', timeidx=wrf.ALL_TIMES, squeeze=False),
                          'temperature': wrf.getvar(nc, 'tk', timeidx=wrf.ALL_TIMES, squeeze=False),
                          #'rh': wrf.getvar(nc, 'rh', timeidx=wrf.ALL_TIMES, squeeze=False),
                          'u': wrf.getvar(nc, 'ua', timeidx=wrf.ALL_TIMES, squeeze=False), 
                          'v': wrf.getvar(nc, 'va', timeidx=wrf.ALL_TIMES, squeeze=False),
                          #'w': wrf.getvar(nc, 'wa', timeidx=wrf.ALL_TIMES, squeeze=False),
                          'z': wrf.getvar(nc, 'height', timeidx=wrf.ALL_TIMES, squeeze=False),
                          'z_agl': wrf.getvar(nc, 'height_agl', timeidx=wrf.ALL_TIMES, squeeze=False),
                          'mixing_ratio': wrf.getvar(nc, 'QVAPOR', timeidx=wrf.ALL_TIMES, squeeze=False),
                          'wind_10m': wrf.getvar(nc, 'uvmet10_wspd_wdir', timeidx=wrf.ALL_TIMES, squeeze=False),
                          #'updraft_helicity': wrf.getvar(nc, 'updraft_helicity', timeidx=wrf.ALL_TIMES, squeeze=False),
                          #'slp': wrf.getvar(nc, 'slp', timeidx=wrf.ALL_TIMES, squeeze=False),
                          #'ter': wrf.getvar(nc, 'ter', timeidx=wrf.ALL_TIMES, squeeze=False),
                          #'td': wrf.getvar(nc, 'td', timeidx=wrf.ALL_TIMES, squeeze=False),
                          'hailcast_diam_max': wrf.getvar(nc, 'HAILCAST_DIAM_MAX', timeidx=wrf.ALL_TIMES, squeeze=False)})
                          
    assert not np.any(dat.pressure == dat.pressure.attrs['_FillValue'])
    assert not np.any(dat.pressure == dat.pressure.attrs['missing_value'])
    del dat.pressure.attrs['_FillValue']
    del dat.pressure.attrs['missing_value']
    
    dat['specific_humidity'] = metpy.calc.specific_humidity_from_mixing_ratio(dat.mixing_ratio).metpy.dequantify()
    dat = dat.rename({'Time': 'time', 'XLONG': 'longitude', 'XLAT': 'latitude'})
    dat = dat.reset_coords().drop('XTIME')
    dat = dat.assign_coords({'south_north': dat.south_north,
                             'west_east': dat.west_east,
                             'bottom_top': dat.bottom_top})

    #assert np.all(dat.ter == dat.ter.max('time')), 'Terrain is not constant.'
    #dat['ter'] = dat.ter.max('time', keep_attrs=True)

    dat.attrs['projection'] = str(dat.u.attrs['projection'])

    for k in dat.keys():
        if 'projection' in dat[k].attrs:
            del dat[k].attrs['projection']
        
    comp = dict(zlib=True, shuffle=True, complevel=5)
    encoding = {var: comp for var in dat.data_vars}
    dat.to_netcdf(outfile, encoding=encoding)
    del dat
    nc.close()
    del nc
