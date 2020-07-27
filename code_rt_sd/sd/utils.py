#!/usr/bin/env python

"""utils.py: module is dedicated to utility methods to simulate algorithms."""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import matplotlib
matplotlib.use("Agg")
import os
import datetime as dt
import pandas as pd
from netCDF4 import Dataset, num2date
import numpy as np
import calendar
from pysolar.solar import get_altitude

def create_folder_structures(dn, stn):
    """ Create full folder structures for the simulations """
    folder = "data/sim/{dn}/{stn}".format(dn=dn.strftime("%Y.%m.%d.%H.%M"), stn=stn)
    if not os.path.exists(folder): os.system("mkdir -p " + folder) 
    return

def extrap1d(x,y,kind="linear"):
    """ This method is used to extrapolate 1D paramteres """
    interpolator = interp1d(x,y,kind=kind)
    xs = interpolator.x
    ys = interpolator.y
    def pointwise(x):
        if x < xs[0]: return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]: return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else: return interpolator(x)
    def ufunclike(xs):
        return array(list(map(pointwise, array(xs))))
    return ufunclike

def get_sd_radar(rad):
    """ Get lat, lon, and bearing angles of the radar """
    sd = pd.read_csv("config/radars.csv")
    r = sd[sd.rad==rad]
    lat, lon, bearing = r["lat"].tolist()[0], r["lon"].tolist()[0], r["ray_bearing"].tolist()[0]
    return lat, lon, bearing

def get_geolocate_range_cells(rad, beam=None):
    """
    Fetch geolocate range cells
    rad: Radar code
    """
    fname = "config/geo/{rad}.geolocate.data.nc".format(rad=rad)
    gzfname = "config/geo/{rad}.geolocate.data.nc.gz".format(rad=rad)
    os.system("gzip -d " + gzfname)
    rootgrp = netCDF4.Dataset(fname)
    os.system("gzip " + fname)
    lat, lon = rootgrp.variables["lat"][:], rootgrp.variables["lon"][:]
    if beam is not None: lat, lon = lat[beam, :], lon[beam, :]
    return lat, lon

def calculate_sza(d, lat, lon, alt=300):
    """
    This method is used to estimate the solar zenith angle for a specific date and
    sepcific location in space. Note that this method uses skyfield api to estimate
    solar zenith angle. This has been validated against NOAA website values.
    """
    d = d.replace(tzinfo=dt.timezone.utc)
    sza = 90. - get_altitude(lat, lon, d)
    return sza

def download_goes_data(dn, sat=15, v=True):
    """ Download GOES data """
    def _get_month_bounds_(start_time):
        """ This method is used to get the first and last date of the month """
        month_start = start_time.replace(day = 1).strftime("%Y%m%d")
        _, month_end = calendar.monthrange(start_time.year, start_time.month)
        month_end = (start_time.replace(day = 1) + dt.timedelta(days=month_end-1)).strftime("%Y%m%d")
        return month_start, month_end
    fname = "data/sim/{dnx}/goes.csv".format(dnx=dn.strftime("%Y.%m.%d.%H.%M"))
    if not os.path.exists(fname+".gz"):
        month_start, month_end = _get_month_bounds_(dn)
        url = "https://satdat.ngdc.noaa.gov/sem/goes/data/avg/{year}/{month}/goes{sat}/netcdf/"\
                "g{sat}_xrs_1m_{mstart}_{mend}.nc".format(year=dn.year, month="%02d"%dn.month, sat=sat,
                        mstart=month_start, mend=month_end)
        if v: print("\n Download file -from- " + url)
        tag_vars = ["A_AVG","B_AVG"]
        fn = fname.replace(".csv",".nc")
        os.system("wget -O {fn} {url}".format(fn=fn, url=url))
        if os.path.exists(fn):
            nc = Dataset(fn)
            tt = nc.variables["time_tag"]
            jd = np.array(num2date(tt[:],tt.units))
            data = {}
            for var in tag_vars:  data[var] = nc.variables[var][:]
            data["date"] = jd
            data_dict = pd.DataFrame(data)
            data_dict = data_dict[(data_dict.date >= dn-dt.timedelta(minutes=30)) 
                    & (data_dict.date < dn+dt.timedelta(minutes=2))]
            data_dict.to_csv(fname, index=False, header=True)
            os.system("gzip {fname}".format(fname=fname))
            if v: print("\n File saved  -to- " + fname)
            os.remove(fn)
        else: print("\n Unable to download file.")
    return

def read_goes(dn):
    """ This method is used to fetch GOES x-ray data for a given day """
    download_goes_data(dn)
    gzfname = "data/sim/{dnx}/goes.csv.gz".format(dnx=dn.strftime("%Y.%m.%d.%H.%M"))
    fname = "data/sim/{dnx}/goes.csv".format(dnx=dn.strftime("%Y.%m.%d.%H.%M"))
    os.system("gzip -d " + gzfname)
    _o = pd.read_csv(fname,parse_dates=["date"])
    os.system("gzip {fname}".format(fname=fname))
    return _o

def get_rtime(dn, sig="B_AVG", th=3.e-6):
    """ This method is used to calculate rise time in minutes """
    rtime = np.nan
    _o = read_goes(dn)
    sig = np.array(_o[sig])
    ix = np.min([i for i,v in enumerate(sig) if v>th and v<np.max(sig)])
    dx = _o.date.tolist()[ix]
    print("\n start time - ", dx)
    rtime = ((dn-dx).total_seconds())/60. 
    return rtime
