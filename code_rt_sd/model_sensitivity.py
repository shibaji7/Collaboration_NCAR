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
import glob
from scipy.integrate import trapz
from scipy import signal

import plotlib

INT_F = 300
def get_freq(dn, rad):
    f = 12.
    fname = "data/op/{dn}/waccmx/sd_{rad}_data.csv.gz".format(dn=dn.strftime("%Y.%m.%d.%H.%M"),rad=rad)
    if os.path.exists(fname):
        os.system("gzip -d " + fname)
        du = pd.read_csv(fname.replace(".gz", ""))
        os.system("gzip " + fname.replace(".gz", ""))
        if len(du) > 0: f = np.median(du.tfreq)/1e3
    return f

def get_vdeta(d, rlim, freq):
    d = d[(d.height>=rlim[0]) & (d.height<rlim[1])]
    f = trapz(signal.resample(d.dop,INT_F))
    v = (0.5 * f * 3e8 / (freq * 1e6))
    return v

def _estimate_dop_delh_(x, y, freq, phi=0):
    dh = (np.max(x.height) - np.max(y.height)) * 1000.
    xf = (-2.*freq*1e6/3e8) * (dh/60.) * np.cos(np.deg2rad(phi))
    xd = 0.5 * xf * 3e8 / (freq * 1e6)
    return xd

sim_fname = "finescale_simulate.csv"
t, T = True, True
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_false", help="Increase output verbosity (default True)")
    parser.add_argument("-pl", "--plot", action="store_false", help="Analyze sensitivity (default True)")
    d_reg, e_reg, f_reg = [60,90], [100,130], [140,300]
    args = parser.parse_args()
    model = "waccmx"
    vD, vE, vF, vFh, vT = [], [], [], [], []
    if not args.plot: 
        events = pd.read_csv("config/events.csv", parse_dates=["dn"])
        rads = pd.read_csv("config/radars.csv")
        for d in events.dn.tolist():
            for r in rads.rad.tolist():
                r = "bks"
                if r in ["bks"]:
                    freq = get_freq(d, r)
                    for bm in range(24):
                        dic = "data/op/{dn}/waccmx/{r}/bm.{bm}/".format(r=r,bm="%02d"%bm,dn=d.strftime("%Y.%m.%d.%H.%M"))
                        for i in range(18,21):
                            i_start, i_end = 16, 30
                            for elv in np.linspace(i_start, i_end, (i_end-i_start)*2 + 1):
                                if elv.is_integer(): fname = dic + "ti(%02d)_elv(%d)_f.csv"%(i,elv)
                                else: fname = dic + "ti(%02d)_elv(%.1f)_f.csv"%(i,elv)
                                files = glob.glob(fname)
                                for f in files:
                                    T = True
                                    print(f, f.replace("_f.","_d."))
                                    ff = pd.read_csv(f)
                                    vd = np.abs(get_vdeta(ff, d_reg, freq))
                                    ve = np.abs(get_vdeta(ff, e_reg, freq))
                                    vf = np.abs(get_vdeta(ff, f_reg, freq))
                                    bf = pd.read_csv(f.replace("_f.","_d."))
                                    vfh = np.abs(_estimate_dop_delh_(ff,bf,freq))
                                    vD.append(vd)
                                    vE.append(ve)
                                    vF.append(vf)
                                    vFh.append(vfh)
                                    vT.append(vd+ve+vf+vfh)
                if t and T: break
            if t and T: break
        x = pd.DataFrame()
        x["vD"], x["vE"], x["vF"], x["vFh"], x["vT"] = vD, vE, vF, vFh, vT
        x.to_csv("config/"+sim_fname, index=False)
    else: 
        x = pd.read_csv("config/"+sim_fname)
        x = x[(np.abs(x.vT)>30) & (np.abs(x.vT)<300)]
        vd, ve, vf = np.array(x.vD/x.vT), np.array(x.vE/x.vT), np.array((x.vF + x.vFh)/x.vT)
        vd = vd[vd>.1]
        plotlib.plot_region_distribution(vd, ve, vf)
        vdn, vdh = np.array((x.vD + x.vE + x.vF)/x.vT), np.array(x.vFh/x.vT)
        plotlib.plot_distribution(vdn, vdh)
        plotlib.plot_htstogram(vd, ve, vf, vdn, vdh)
    if os.path.exists("sd/__pycache__/"):
        os.system("rm -rf sd/__pycache__/")
        os.system("rm -rf py*.log")
