#!/usr/bin/env python

"""wxproc.py: module is dedicated to download and process the WACCMX data."""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import os
os.system("source ncar_tgcm_waccm_proc/ncar/bin/activate")
from netCDF4 import Dataset, num2date
import datetime as dt
import time
import numpy as np
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-ev", "--event")
    parser.add_argument("-dn", "--date")
    args = parser.parse_args()
    event, date = dt.datetime.strptime(args.event, "%Y.%m.%d.%H.%M"), dt.datetime.strptime(args.date, "%Y.%m.%d.%H.%M")
    ldat = Dataset("tmp/waccmx.nc", "w", format="NETCDF4")
    dates = [date, date-dt.timedelta(minutes=1)]
    typs = ["flr", "flr"]
    cases = ["flr","dly"]
    fname = "/glade/scratch/shibaji/archive/R.{ev}.{case}-FXHIST-f19_f19/atm/hist/R.{ev}.{case}-FXHIST-f19_f19.cam.h1.{dn}-{tm}.nc"
    t = []
    zg = {"dly": np.zeros((1, 126, 96, 144)), "flr": np.zeros((1, 126, 96, 144))}
    ne = {"dly": np.zeros((1, 126, 96, 144)), "flr": np.zeros((1, 126, 96, 144))}
    sp = {"dly": np.zeros((1, 126, 96, 144)), "flr": np.zeros((1, 126, 96, 144))}
    sh = {"dly": np.zeros((1, 126, 96, 144)), "flr": np.zeros((1, 126, 96, 144))}
    ui = {"dly": np.zeros((1, 126, 96, 144)), "flr": np.zeros((1, 126, 96, 144))}
    vi = {"dly": np.zeros((1, 126, 96, 144)), "flr": np.zeros((1, 126, 96, 144))}
    wi = {"dly": np.zeros((1, 126, 96, 144)), "flr": np.zeros((1, 126, 96, 144))}
    for case, date, ty in zip(cases, dates, typs):
        tm = int((date - date.replace(hour=0, minute=0, second=0)).total_seconds()/(24*60))*(24*60)
        f = fname.format(ev=event.strftime("%Y.%m.%d.%H.%M"), case=ty, dn=event.strftime("%Y-%m-%d"), tm="%05d"%tm)    
        rdat = Dataset(f)
        print(" Fname - ",f)
        lat, lon, lev = rdat.variables["lat"], rdat.variables["lon"], rdat.variables["lev"]
        tx = num2date(rdat.variables["time"][:], units=rdat.variables["time"].units)
        tx = np.array([x._to_real_datetime() for x in tx]).astype("datetime64[ns]")
        tx = [dt.datetime.utcfromtimestamp(x.astype(int) * 1e-9) for x in tx]
        tm = tx.index(date)
        zg[case][0,:,:,:] = rdat.variables["Z3"][tm,:,:,:]
        ne[case][0,:,:,:] = rdat.variables["EDens"][tm,:,:,:]
        sp[case][0,:,:,:] = rdat.variables["SIGMAPED"][tm,:,:,:]
        sh[case][0,:,:,:] = rdat.variables["SIGMAHAL"][tm,:,:,:]
        ui[case][0,:,:,:] = rdat.variables["UI"][tm,:,:,:]
        vi[case][0,:,:,:] = rdat.variables["VI"][tm,:,:,:]
        wi[case][0,:,:,:] = rdat.variables["WI"][tm,:,:,:]
    ldat.createDimension("lat", len(lat[:]))
    ldat.createDimension("lon", len(lon[:]))
    ldat.createDimension("lev", len(lev[:]))
    ldat.createDimension("time", 1)
    rlat, rlon, rlev = ldat.createVariable("lat","f4",("lat",)), ldat.createVariable("lon","f4",("lon",)),\
            ldat.createVariable("lev","f4",("lev",))
    rlat[:], rlon[:], rlev[:] = lat[:], np.mod( (lon[:] + 180), 360 ) - 180, lev[:]
    rZGd, rZGf, rNEd, rNEf = ldat.createVariable("ZGd","f4",("time", "lev", "lat","lon")),\
            ldat.createVariable("ZGf","f4",("time", "lev", "lat","lon")),\
            ldat.createVariable("NEd","f4",("time", "lev", "lat","lon")),\
            ldat.createVariable("NEf","f4",("time", "lev", "lat","lon"))
    rSPd, rSPf, rSHd, rSHf = ldat.createVariable("SPd","f4",("time", "lev", "lat","lon")),\
            ldat.createVariable("SPf","f4",("time", "lev", "lat","lon")),\
            ldat.createVariable("SHd","f4",("time", "lev", "lat","lon")),\
            ldat.createVariable("SHf","f4",("time", "lev", "lat","lon"))
    rUId, rUIf, rVId, rVIf = ldat.createVariable("UId","f4",("time", "lev", "lat","lon")),\
            ldat.createVariable("UIf","f4",("time", "lev", "lat","lon")),\
            ldat.createVariable("VId","f4",("time", "lev", "lat","lon")),\
            ldat.createVariable("VIf","f4",("time", "lev", "lat","lon"))
    rWId, rWIf = ldat.createVariable("WId","f4",("time", "lev", "lat","lon")),\
            ldat.createVariable("WIf","f4",("time", "lev", "lat","lon"))
    rZGd[:], rZGf[:], rNEd[:], rNEf[:] = zg["dly"][:], zg["flr"][:], ne["dly"][:], ne["flr"][:]
    rSPd[:], rSPf[:], rSHd[:], rSHf[:] = sp["dly"][:], sp["flr"][:], sh["dly"][:], sh["flr"][:]
    rUId[:], rUIf[:], rVId[:], rVIf[:] = ui["dly"][:], ui["flr"][:], vi["dly"][:], vi["flr"][:]
    rWId[:], rWIf[:] = wi["dly"][:], wi["flr"][:]
    ldat.close()
    rdat.close()
    os.system("gzip tmp/waccmx.nc")
