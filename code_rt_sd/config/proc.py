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
from netCDF4 import Dataset
import time
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fname", help="Source location of file")
    parser.add_argument("-s", "--start", type=int, help="Source location of file")
    parser.add_argument("-e", "--end", type=int, help="Source location of file")
    args = parser.parse_args()
    fname, start, end = args.fname, args.start, args.end
    ldat = Dataset("tmp/timegcm.nc", "w", format="NETCDF4")
    rdat = Dataset(fname)
    ZG = rdat.variables["ZG"]
    NE = rdat.variables["NE"]
    lat, lon, time, lev = rdat.variables["lat"], rdat.variables["lon"], rdat.variables["time"], rdat.variables["lev"]
    ldat.createDimension("lat", len(lat[:]))
    ldat.createDimension("lon", len(lon[:]))
    ldat.createDimension("time", len(time[start:end]))
    ldat.createDimension("lev", len(lev[:]))
    rlat, rlon, rtime, rlev = ldat.createVariable("lat","f8",("lat",)), ldat.createVariable("lon","f8",("lon",)),\
            ldat.createVariable("time","f8",("time",)), ldat.createVariable("lev","f8",("lev",))
    rlat[:], rlon[:], rtime[:], rlev[:] = lat[:], lon[:], time[start:end], lev[:]
    rtime.units = time.units
    rZG, rNE = ldat.createVariable("ZG","f8",("time", "lev", "lat","lon")), ldat.createVariable("NE","f8",("time", "lev", "lat","lon"))
    rZG[:], rNE[:] = ZG[start:end,:,:,:], NE[start:end,:,:,:]
    ldat.close()
    rdat.close()
