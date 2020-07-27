#!/usr/bin/env python

"""proc.py: module is dedicated to download and process the TIME-GCM data."""

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
    parser.add_argument("-rn", "--run")
    args = parser.parse_args()
    event, date = dt.datetime.strptime(args.event, "%Y.%m.%d.%H.%M"), dt.datetime.strptime(args.date, "%Y.%m.%d.%H.%M")
    ldat = Dataset("tmp/tgcm.nc", "w", format="NETCDF4")
    cases = ["1","2"]
    fname = "/glade/work/shibaji/timegcm/cases/{dn}-{r}.{c}/timegcm_trunk/timegcm-ch/timegcm.s_{doy}.nc"
    t = []
    zg = {"dly": np.zeros((1, 49, 36, 72)), "flr": np.zeros((1, 49, 36, 72))}
    ne = {"dly": np.zeros((1, 49, 36, 72)), "flr": np.zeros((1, 49, 36, 72))}
    cds = {"1":"dly","2":"flr"}
    for case in cases:
        tm = int((date - date.replace(hour=0, minute=0, second=0)).total_seconds()/(24*60))*(24*60)
        f = fname.format(dn=event.strftime("%Y.%m.%d.%H.%M"), c=case, r=args.run, doy=date.strftime("%j"))    
        rdat = Dataset(f)
        print(" Fname - ",f)
        lat, lon, lev = rdat.variables["lat"], rdat.variables["lon"], rdat.variables["lev"]
        tx = num2date(rdat.variables["time"][:], units=rdat.variables["time"].units)
        tx = np.array([x._to_real_datetime() for x in tx]).astype("datetime64[ns]")
        tx = [dt.datetime.utcfromtimestamp(x.astype(int) * 1e-9) - dt.timedelta(days=1) for x in tx]
        tm = tx.index(date)
        zg[cds[case]][0,:,:,:] = rdat.variables["ZG"][tm,:,:,:]
        ne[cds[case]][0,:,:,:] = rdat.variables["NE"][tm,:,:,:]
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
    rZGd[:], rZGf[:], rNEd[:], rNEf[:] = zg["dly"][:], zg["flr"][:], ne["dly"][:], ne["flr"][:]
    ldat.close()
    rdat.close()
    os.system("gzip tmp/tgcm.nc")
