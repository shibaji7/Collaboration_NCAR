#!/usr/bin/env python

"""model_sensitivity.py: Sensitivity test experiments python program for Doppler shift"""

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

from superdarn import ModelSensitivity
from plotlib import ModelSensitivity as SAS


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_false", help="Increase output verbosity (default True)")
    parser.add_argument("-mr", "--mrange", type=float, default=1990., help="Max ground range in km (default 1000 km)")
    parser.add_argument("-nmr", "--nmrange", type=int, default=200, help="Number of steps in ground range (default 101)")
    parser.add_argument("-sh", "--sheight", type=float, default=50., help="Start height in km (default 50 km)")
    parser.add_argument("-nmh", "--eheight", type=float, default=350., help="End height in km (default 350 km)")
    parser.add_argument("-hinc", "--hinc", type=float, default=1., help="Step in height in km (default 1 km)")
    parser.add_argument("-es", "--selev", type=float, default=30., help="Start elevation angle deg")
    parser.add_argument("-ee", "--eelev", type=float, default=35., help="End elevation angle deg")
    parser.add_argument("-ei", "--ielev", type=float, default=1., help="Inc of elevation angle deg (default 20*)")
    parser.add_argument("-nhops", "--nhops", type=float, default=1, help="Number of hops (default 1)")
    parser.add_argument("-pl", "--plot", action="store_false", help="Analyze sensitivity (default True)")
    parser.add_argument("-reg", "--regs", action="store_true", help="Analyze sensitivity (default False)")
    args = parser.parse_args()
    if args.verbose:
        print("\n Parameter list for simulation ")
        for k in vars(args).keys():
            print("     " , k , "->" , str(vars(args)[k]))
    model = "waccmx"
    if not args.plot:
        os.system("rm data/sim/model_sensitivity.nc")
        d_ratio, e_ratio, f_ratio, d_rate, e_rate, f_rate, frequency, sza, LT = [], [], [], [], [], [], [], [], []
        vn, vh, vn_max, vn_min, vh_max, vh_min = [], [], [], [], [], []
        rads = ["bks", "wal", "fhe", "fhw", "cve", "cvw", "gbr", "kap", "sas", "pgr"]
        events = pd.read_csv("config/events.csv", parse_dates=["dn"])
        events = events.drop(4)
        for dn in events.dn.tolist():
            for rad in rads:
                for bm in range(0,24):
                    for j in range(15, 19):
                        fname = "data/op/{dn}/{model}/{rad}/bm.{bm}/ne.ti({jx}).d.mat".format(dn=dn.strftime("%Y.%m.%d.%H.%M"),
                                bm="%02d"%bm, model=model, jx=j, rad=rad)
                        if os.path.exists(fname):
                            print(" Path -> ", fname)
                            v1, v2, v1_max, v1_min, v2_max, v2_min, dr, er, fr, drr,\
                                    err, frr, frq, chi, lt = ModelSensitivity(dn, rad, args, model=model, bm=bm, jx=j).exe()
                            vn.append(v1)
                            vh.append(v2)
                            vn_max.append(v1_max)
                            vn_min.append(v1_min)
                            vh_max.append(v2_max)
                            vh_min.append(v2_min)
                            d_ratio.append(dr)
                            e_ratio.append(er)
                            f_ratio.append(fr)
                            d_rate.append(drr)
                            e_rate.append(err)
                            f_rate.append(frr)
                            frequency.append(frq)
                            sza.append(chi)
                            LT.append(lt)
                break
        N = len(vn)
        print(" Total samples - ", N)
        rootgrp = Dataset("data/sim/model_sensitivity.nc", "w", format="NETCDF4")
        rootgrp.description = """Optimization based on driving parameters."""
        rootgrp.history = "Created " + time.ctime(time.time())
        rootgrp.source = "SuperDARN Sudden Phase Anomaly: Doppler Flash"
        rootgrp.createDimension("samples", N)
        params = {"vn": vn, "vh": vh, "d_ratio": d_ratio, "e_ratio": e_ratio, "f_ratio": f_ratio, 
                "frequency": frequency, "sza": sza, "d_rate":d_rate, "e_rate":e_rate, "f_rate":f_rate, 
                "vt": np.array(vn)+np.array(vh), "vn_max":vn_max, "vn_min":vn_min, "vh_max":vh_max, "vh_min":vh_min,
                "vt_min":np.array(vn_min)+np.array(vh_min), "vt_max":np.array(vn_max)+np.array(vh_max)}
        for k in params.keys():
            x = rootgrp.createVariable(k,"f8",("samples"))
            x[:] = params[k]
        rootgrp.close()
        os.system("cp data/sim/model_sensitivity.nc config/")
    else: 
        dat = Dataset("config/model_sensitivity.nc")
        sas = SAS(dat)
        sas.analyze(args.regs)
    if os.path.exists("sd/__pycache__/"):
        os.system("rm -rf sd/__pycache__/")
        os.system("rm -rf py*.log")
