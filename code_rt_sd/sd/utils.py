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
