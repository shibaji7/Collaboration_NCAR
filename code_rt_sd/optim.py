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

from superdarn import SuperDARN


_I_, _J_, _K_, _L_ = 7, 5, 5, 5

if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("-r", "--rad", default="bks", help="Radar code (default bks)")
        parser.add_argument("-ev", "--event", default=dt.datetime(2015,5,5,22,10), help="Start date (default 2015-05-05T22:10)",
                type=dparser.isoparse)
        parser.add_argument("-b", "--bmnum", type=int, default=0, help="Beam num starting from 0 in center and +, - on right-left")
        parser.add_argument("-dra", "--d_ratio", type=float, default=100., help="D region e-denisity peak times (100)")
        parser.add_argument("-ds", "--d_start", type=float, default=10., help="D region start (60.)")
        parser.add_argument("-de", "--d_end", type=float, default=35., help="D region end (80.)")
        parser.add_argument("-drt", "--d_rtime", type=float, default=2., help="D region e-denisity rise time (2.)")
        parser.add_argument("-era", "--e_ratio", type=float, default=1., help="E region e-denisity peak order (1.)")
        parser.add_argument("-fra", "--f_ratio", type=float, default=1.01, help="F region e-denisity peak order")
        parser.add_argument("-fs", "--f_start", type=float, default=130., help="F region start (180.)")
        parser.add_argument("-fe", "--f_end", type=float, default=190, help="F region end (240.)")
        parser.add_argument("-frt", "--f_rtime", type=float, default=0.5, help="F region e-denisity rise time (0.5)")
        parser.add_argument("-rd", "--save_radar", action="store_false", help="Save riometer data (default True)")
        parser.add_argument("-ps", "--plot_summary", action="store_true", help="Plot summary report (default False)")
        parser.add_argument("-plt", "--plot", action="store_false", help="Plot summary report (default True)")
        parser.add_argument("-sr", "--save_result", action="store_false", help="Save results (default True)")
        parser.add_argument("-c", "--clear", action="store_true", help="Clear pervious stored files (default False)")
        parser.add_argument("-v", "--verbose", action="store_false", help="Increase output verbosity (default True)")
        parser.add_argument("-fr", "--frequency", type=float, default=12., help="Frequency of oprrations in MHz (default 12 MHz)")
        parser.add_argument("-mr", "--mrange", type=float, default=1990., help="Max ground range in km (default 1000 km)")
        parser.add_argument("-nmr", "--nmrange", type=int, default=200, help="Number of steps in ground range (default 101)")
        parser.add_argument("-sh", "--sheight", type=float, default=50., help="Start height in km (default 50 km)")
        parser.add_argument("-nmh", "--eheight", type=float, default=350., help="End height in km (default 350 km)")
        parser.add_argument("-hinc", "--hinc", type=float, default=1., help="Step in height in km (default 1 km)")
        parser.add_argument("-es", "--selev", type=float, default=30., help="Start elevation angle deg")
        parser.add_argument("-ee", "--eelev", type=float, default=35., help="End elevation angle deg")
        parser.add_argument("-ei", "--ielev", type=float, default=1., help="Inc of elevation angle deg (default 20*)")
        parser.add_argument("-nhops", "--nhops", type=float, default=1, help="Number of hops (default 1)")
        args = parser.parse_args()
        if args.verbose:
            print("\n Parameter list for simulation ")
            for k in vars(args).keys():
                print("     " , k , "->" , str(vars(args)[k]))
        if args.plot:
            dat = Dataset("config/optim.nc")
            d_ratios = dat.variables["d_ratio"][:]
            d_rtimes = dat.variables["d_rtime"][:]
            f_ratios = dat.variables["f_ratio"][:]
            f_rtimes = dat.variables["f_rtime"][:]
            v_d, v_f = dat.variables["vd"][:], dat.variables["vf"][:]
            u = []
            for i, dr in enumerate(d_ratios):
                for j, dtr in enumerate(d_rtimes):
                    for k, fr in enumerate(f_ratios):
                        for l, ftr in enumerate(f_rtimes):
                            x = i+j+k+l
                            u.append({"dratio": d_ratios[i], "drtime": d_rtimes[j], "fratio": f_ratios[k], "frtime": f_rtimes[l],
                                "vd": v_d[i,j,k,l], "vf": v_f[i,j,k,l]})
            df = pd.DataFrame.from_records(u)
            import statsmodels.api as sm
            from statsmodels.formula.api import ols
            models = {}
            models["m1"] = ols("vd~dratio+fratio+drtime+frtime", data=df).fit()
            table = sm.stats.anova_lm(models["m1"], typ=2)
            print(table)
        else:
            d_ratios = np.linspace(2., 400., _I_)
            d_rtimes = np.linspace(0.5, 3.5, _J_)
            f_ratios = np.linspace(1., 1.2, _K_)
            f_rtimes = np.linspace(0.5, 3., _L_)
            os.system("rm data/sim/optim.nc")
            rootgrp = Dataset("data/sim/optim.nc", "w", format="NETCDF4")
            rootgrp.description = """Optimization based on 4 parameters."""
            rootgrp.history = "Created " + time.ctime(time.time())
            rootgrp.source = "SuperDARN Sudden Phase Anomaly: Doppler Flash"
            rootgrp.createDimension("ndratio", len(d_ratios))
            rootgrp.createDimension("ndrtime", len(d_rtimes))
            rootgrp.createDimension("nfratio", len(f_ratios))
            rootgrp.createDimension("nfrtime", len(f_rtimes))
            x = rootgrp.createVariable("d_ratio","f4",("ndratio",))
            x.description = "Fractor by which D-region density goes up."
            x[:] = d_ratios
            x = rootgrp.createVariable("f_ratio","f4",("nfratio",))
            x.description = "Fractor by which F-region density goes up."
            x[:] = f_ratios
            x = rootgrp.createVariable("d_rtime","f4",("ndrtime",))
            x.description = "Time taken by the D region e-denisty rise. (in mins)"
            x[:] = d_rtimes
            x = rootgrp.createVariable("f_rtime","f4",("nfrtime",))
            x.description = "Time taken by the F region e-denisty rise. (in mins)"
            x[:] = f_rtimes
            vd, vf = np.zeros((len(d_ratios),len(d_rtimes),len(f_ratios),len(f_rtimes))),\
                    np.zeros((len(d_ratios),len(d_rtimes),len(f_ratios),len(f_rtimes)))
            vd[:], vf[:] = np.nan, np.nan
            for i, dr in enumerate(d_ratios):
                for j, dtr in enumerate(d_rtimes):
                    for k, fr in enumerate(f_ratios):
                        for l, ftr in enumerate(f_rtimes):
                            args.d_ratio, args.d_rtime, args.f_ratio, args.f_rtime = dr, dtr, fr, ftr
                            vd[i,j,k,l], vf[i,j,k,l], _ = SuperDARN(args)._exe_()
            x = rootgrp.createVariable("vd","f8",("ndratio","ndrtime","nfratio","nfrtime"))
            x.description = "Velocity estimated from D-region density change. (in m/s)"
            x[:] = vd
            x = rootgrp.createVariable("vf","f8",("ndratio","ndrtime","nfratio","nfrtime"))
            x.description = "Velocity estimated from F-region density change. (in m/s)"
            x[:] = vf
            print("")
            rootgrp.close()
            os.system("cp data/sim/optim.nc config/")
        if os.path.exists("sd/__pycache__/"):
            os.system("rm -rf sd/__pycache__/")
            os.system("rm -rf py*.log")
