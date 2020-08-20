#!/usr/bin/env python

"""optim.py: optimize experiments python program for Doppler shift"""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
import sys
sys.path.append("sd/")
import datetime as dt
import pandas as pd
import argparse
from dateutil import parser as dparser
import numpy as np
from netCDF4 import Dataset
import time

from superdarn import Senstitivity
from plotlib import SensitivityAnalysis as SAS

from SALib.sample import saltelli
from SALib.analyze import sobol

if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("-r", "--rad", default="bks", help="Radar code (default bks)")
        parser.add_argument("-ev", "--event", default=dt.datetime(2015,5,5,22,10), help="Start date (default 2015-05-05T22:10)",
                type=dparser.isoparse)
        parser.add_argument("-b", "--bmnum", type=int, default=0, help="Beam num starting from 0 in center and +, - on right-left")
        parser.add_argument("-dra", "--d_ratio", type=float, default=1., help="D region e-denisity peak times")
        parser.add_argument("-era", "--e_ratio", type=float, default=1., help="E region e-denisity peak order")
        parser.add_argument("-fra", "--f_ratio", type=float, default=1., help="F region e-denisity peak order")
        parser.add_argument("-v", "--verbose", action="store_false", help="Increase output verbosity (default True)")
        parser.add_argument("-fr", "--frequency", type=float, default=12., help="Frequency of oprrations in MHz (default 12 MHz)")
        parser.add_argument("-mr", "--mrange", type=float, default=1990., help="Max ground range in km (default 1000 km)")
        parser.add_argument("-nmr", "--nmrange", type=int, default=200, help="Number of steps in ground range (default 101)")
        parser.add_argument("-sh", "--sheight", type=float, default=50., help="Start height in km (default 50 km)")
        parser.add_argument("-nmh", "--eheight", type=float, default=350., help="End height in km (default 350 km)")
        parser.add_argument("-hinc", "--hinc", type=float, default=1., help="Step in height in km (default 1 km)")
        parser.add_argument("-es", "--selev", type=float, default=16., help="Start elevation angle deg")
        parser.add_argument("-ee", "--eelev", type=float, default=35., help="End elevation angle deg")
        parser.add_argument("-ei", "--ielev", type=float, default=1., help="Inc of elevation angle deg (default 20*)")
        parser.add_argument("-nhops", "--nhops", type=float, default=1, help="Number of hops (default 1)")
        parser.add_argument("-pl", "--plot", action="store_false", help="Analyze sensitivity (default True)")
        args = parser.parse_args()
        if args.verbose:
            print("\n Parameter list for simulation ")
            for k in vars(args).keys():
                print("     " , k , "->" , str(vars(args)[k]))
        N, n = 500, 3
        problem = {
                "num_vars": n,
                "names": ["D-Ratio", "E-Ratio", "F-Ratio"],
                "bounds": [[1, 200],
                    [1., 1.2],
                    [1, 1.1]]
                }
        if args.plot:
            dat = Dataset("config/sensitivity.nc")
            y = dat.variables["vd_mean"][:]
            sas = SAS(problem, dat)
            sas.analyze()
        else:
            samples = saltelli.sample(problem, N, calc_second_order=True)
            N = N * (2*n+2)
            os.system("rm data/sim/sensitivity.nc")
            rootgrp = Dataset("data/sim/sensitivity.nc", "w", format="NETCDF4")
            rootgrp.description = """Optimization based on 4 parameters."""
            rootgrp.history = "Created " + time.ctime(time.time())
            rootgrp.source = "SuperDARN Sudden Phase Anomaly: Doppler Flash"
            rootgrp.createDimension("samples", N)
            rootgrp.createDimension("params", 3)
            x = rootgrp.createVariable("parameters","f4",("samples","params"))
            x.description = "All the parameetrs"
            x[:] = samples
            vd, vf = np.zeros(N), np.zeros(N)
            vd[:], vf[:] = np.nan, np.nan
            recs = []
            for i in range(N):
                args.d_ratio, args.e_ratio, args.f_ratio = samples[i][0], samples[i][1], samples[i][2]
                rec = Senstitivity(args)._exe_()
                recs.append(rec)
            recs = np.array(recs)
            attrs = ["mean", "max", "min"]
            for i, a in enumerate(attrs):
                x = rootgrp.createVariable("vd_"+a,"f8",("samples"))
                x.description = "Velocity estimated from D-region density change. (in m/s)"
                x[:] = recs[:,0+i*3]
                x = rootgrp.createVariable("vf_"+a,"f8",("samples"))
                x.description = "Velocity estimated from F-region density change. (in m/s)"
                x[:] = recs[:,1+i*3]
                x = rootgrp.createVariable("vt_"+a,"f8",("samples"))
                x.description = "Velocity estimated from total density change. (in m/s)"
                x[:] = recs[:,1+i*3]
            print("")
            rootgrp.close()
            os.system("cp data/sim/sensitivity.nc config/")
        if os.path.exists("sd/__pycache__/"):
            os.system("rm -rf sd/__pycache__/")
            os.system("rm -rf py*.log")
