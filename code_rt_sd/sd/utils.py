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
import pytz
from timezonefinder import TimezoneFinder


def smooth(x, window_len=51, window="hanning"):
    if x.ndim != 1: raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len: raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<3: return x
    if not window in ["flat", "hanning", "hamming", "bartlett", "blackman"]: raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    s = np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    if window == "flat": w = numpy.ones(window_len,"d")
    else: w = eval("np."+window+"(window_len)")
    y = np.convolve(w/w.sum(),s,mode="valid")
    d = window_len - 1
    y = y[int(d/2):-int(d/2)]
    return y

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

def get_sd_data(fname, bmnum):
    """ Get SD data for a bmnum """
    os.system("gzip -d " + fname)
    x = pd.read_csv(fname.replace(".gz", ""), parse_dates=["time"])
    x = x[x.bmnum==bmnum]
    os.system("gzip " + fname.replace(".gz", ""))
    return x

def calculate_LT(d, lat, lon):
    """ Calculate local time from latitude and longitude """
    tf = TimezoneFinder()
    tzf = tf.timezone_at(lng=lon, lat=lat)
    dn = d.replace(tzinfo=pytz.utc)
    lt = dn.astimezone(pytz.timezone(tzf))
    return lt

def get_gridded_parameters(q, xparam="time", yparam="slist", zparam="v"):
    """
    Method converts scans to "beam" and "slist" or gate
    """
    plotParamDF = q[ [xparam, yparam, zparam] ]
    plotParamDF[xparam] = plotParamDF[xparam].tolist()
    plotParamDF[yparam] = plotParamDF[yparam].tolist()
    plotParamDF = plotParamDF.groupby( [xparam, yparam] ).mean().reset_index()
    plotParamDF = plotParamDF[ [xparam, yparam, zparam] ].pivot( xparam, yparam )
    x = plotParamDF.index.values
    y = plotParamDF.columns.levels[1].values
    X, Y  = np.meshgrid( x, y )
    # Mask the nan values! pcolormesh can't handle them well!
    Z = np.ma.masked_where(
            np.isnan(plotParamDF[zparam].values),
            plotParamDF[zparam].values)
    print(Z.shape, X.shape, Y.shape)
    return X,Y,Z

def medfilt2D_weight(X, kernel=np.array([[1,3,1],[3,5,3],[1,3,1]]), tau=0.7):
    """
    Weighted median filter
    """
    Y = np.zeros_like(X) * np.nan
    for i in range(1, X.shape[0]-2):
        for j in range(1, X.shape[1]-2):
            x = X[i-1:i+2, j-1:j+2]
            w = np.repeat(x.filled(np.nan).ravel(), kernel.ravel())
            if np.sum(kernel[np.logical_not(np.ma.getmask(x))])/np.sum(kernel) >= tau: Y[i,j] = np.nanmedian(w)
    Y = np.ma.masked_invalid(Y)
    return Y

class InterpolateData(object):
    """ Interpolate data from lev to Zh """

    def __init__(self):
        return

    def extrap1d(self,x,y,kind="linear"):
        """ This method is used to extrapolate 1D paramteres """
        from scipy.interpolate import interp1d
        from scipy import array
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

    def _get_ij_(self, lats, lons, lat, lon):
        """ Get (lat, lon) index """
        _ij_ = (np.argmin(np.abs(lats-lat)), np.argmin(np.abs(lons-lon)))
        return _ij_

    def _intp_heights_(self, h, hx, param, scale="linear", kind="cubic"):
        from scipy.interpolate import interp1d
        if scale == "linear": pnew = interp1d(h, param, kind=kind)(hx)
        if scale == "log": pnew = 10**interp1d(h, np.log10(param), kind=kind)(hx)
        return pnew

    def _intp_latlon_(self, lats, lons, latx, lonx, param, scale="linear", kind="cubic"):
        from scipy.interpolate import interp2d
        if scale == "linear": pnew = interp2d(lats, lons, param.T, kind=kind)(latx, lonx).T
        if scale == "log": pnew = 10**interp2d(lats, lons, np.log10(param.T), kind=kind)(latx, lonx).T
        return pnew

    def _intp_(self, h, lats, lons, param, hd=[50,350,1], dlat=0.5, dlon=1, scale="log", kind="cubic", v=True):
        lonx = np.arange(np.min(lons), np.max(lons)+dlon, dlon)
        latx = np.arange(np.min(lats), np.max(lats)+dlat, dlat)
        param_intm = np.zeros((len(h), len(latx), len(lonx)))
        h_intm = np.zeros((len(h), len(latx), len(lonx)))
        for k,_ in enumerate(h):
            param_intm[k, :, :] = self._intp_latlon_(lats, lons, latx, lonx, param[k, :, :], scale, kind)
            h_intm[k, :, :] = self._intp_latlon_(lats, lons, latx, lonx, h[k, :, :], scale, kind)
        if v: print("\tLatlon convt.")
        hx = np.arange(hd[0],hd[1],hd[2])
        pnew = np.zeros((len(hx), len(latx), len(lonx)))
        for i,_ in enumerate(latx):
            for j,_ in enumerate(lonx):
                pnew[:,i,j] = 10**(self.extrap1d(h_intm[:, i, j], np.log10(param_intm[:, i, j]))(hx))
        if v: print("\tHeight convt.")
        return pnew, hx, latx, lonx

    def _intrp_(self, h, lats, lons, param, hd=300, v=True):
        param_x = np.zeros((len(lats), len(lons))) * np.nan
        for i in range(len(lats)):
            for j in range(len(lons)):
                param_x[i,j] = param[np.argmin(np.abs(h[:, i, j]-hd)), i, j]
        return param_x

    def _intp_by_beam_(self, h, lats, lons, param, beam, hd=200, v=True):
        param_x = np.zeros((len(lats), len(lons))) * np.nan
        for i in range(len(lats)):
            for j in range(len(lons)):
                param_x[i,j] = param[np.argmin(np.abs(h[:, i, j]-hd)), i, j]
        from scipy.io import loadmat
        x = loadmat("data/op/2015.05.05.22.11/waccmx/bks/bm.%02d/bearing.mat"%beam)
        xlats, xlons = x["lat"][0,:], x["lon"][0,:]
        dat = []
        for xlat, xlon in zip(xlats, xlons):
            dat.append(param_x[np.argmin(np.abs(lats-xlat)), np.argmin(np.abs(lons-xlon))])
        return np.max(dat)
