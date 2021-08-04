#!/usr/bin/env python

"""capture_stats.py: Sensitivity test experiments python program for Doppler shift"""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2021, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
import sys
sys.path.append("../sd/")
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
from scipy.io import loadmat

from pysolar.solar import get_altitude
import plotlib

def calculate_sza(d, lats, lons, alt=300):
    d = d.replace(tzinfo=dt.timezone.utc)
    szas = []
    for la, lo in zip(lats, lons):
        szas.append(90. - get_altitude(la, lo, d))
    return szas

INT_F = 300
def get_freq(dn, rad):
    f = 12.
    fname = "../data/op/{dn}/waccmx/sd_{rad}_data.csv.gz".format(dn=dn.strftime("%Y.%m.%d.%H.%M"),rad=rad)
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

sim_fname = "finescale_simulate_total_{kind}.csv"
t, T = True, True
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_false", help="Increase output verbosity (default True)")
    parser.add_argument("-pl", "--plot", action="store_true", help="Analyze sensitivity (default False)")
    parser.add_argument("-t", "--type", default="A", help="Flare type (A/X/M)")
    parser.add_argument("-n", "--number", default=-1, type=int, help="Number of radar event")
    d_reg, e_reg, f_reg = [60,90], [100,130], [140,300]
    args = parser.parse_args()
    model = "waccmx"
    vD, vE, vF, vFh, vT, SZA, E, F = [], [], [], [], [], [], [], []
    events = pd.read_csv("op/radar_event_list.csv", parse_dates=["date"])
    T, kind = 1, "A"
    if args.type == "M" or args.type == "X":
        events = events[events.type.str.contains(args.type)]
        T, kind = 0, args.type
    if args.number >= 0:
        events = events.iloc[[args.number]]
        T, kind = 0, events.rad.tolist()[0] + "_" + events.date.tolist()[0].strftime("%Y-%m-%d-%H-%M")
    if not args.plot: 
        ix = 0
        for d, r in zip(events.date.tolist(), events.rad.tolist()):
            print("Events - ", r, d)
            if ix >= T: 
                freq = get_freq(d, r)
                for bm in range(24):
                    bearing_file = "../data/op/%s/waccmx/%s/bm.%02d/bearing.mat"%(d.strftime("%Y.%m.%d.%H.%M"), r, bm)
                    print(bearing_file)
                    obj = loadmat(bearing_file)
                    lat, lon = obj["lat"], obj["lon"]
                    sza = np.mean(calculate_sza(d, lat, lon))
                    dic = "../data/op/{dn}/waccmx/{r}/bm.{bm}/".format(r=r,bm="%02d"%bm,dn=d.strftime("%Y.%m.%d.%H.%M"))
                    for i in range(18,19):
                        i_start, i_end = 16, 30
                        for elv in np.linspace(i_start, i_end, (i_end-i_start)*2 + 1):
                            if elv.is_integer(): fname = dic + "ti(%02d)_elv(%d)_f.csv"%(i,elv)
                            else: fname = dic + "ti(%02d)_elv(%.1f)_f.csv"%(i,elv)
                            files = glob.glob(fname)
                            for f in files:
                                try:
                                    T = False
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
                                    SZA.append(sza)
                                    E.append(elv)
                                    F.append(freq)
                                except: pass
            #if t and T: break
            ix += 1
        x = pd.DataFrame()
        x["vD"], x["vE"], x["vF"], x["vFh"], x["vT"], x["sza"], x["elv"], x["freq"] = vD, vE, vF, vFh, vT, SZA, E, F
        x = x.round(3)
        x.to_csv("op/" + sim_fname.format(kind=kind), index=False)
    else:
        print(" Compile plots...")
        print("op/"+sim_fname.format(kind=kind))
        x0 = pd.read_csv("op/"+sim_fname.format(kind=kind))
        x0 = x0[(np.abs(x0.vT)>30) & (np.abs(x0.vT)<300)]
        vd0, ve0, vf0 = np.array(x0.vD/x0.vT), np.array(x0.vE/x0.vT), np.array((x0.vF + x0.vFh)/x0.vT)
        vdn0, vdh0 = np.array((x0.vD + x0.vE + x0.vF)/x0.vT), np.array(x0.vFh/x0.vT)
        x = pd.read_csv("op/finescale_simulate.csv")
        x = x[(np.abs(x.vT)>30) & (np.abs(x.vT)<300)]
        vd, ve, vf = np.array(x.vD/x.vT), np.array(x.vE/x.vT), np.array((x.vF + x.vFh)/x.vT)
        vdn, vdh = np.array((x.vD + x.vE + x.vF)/x.vT), np.array(x.vFh/x.vT)
        plotlib.plot_mastogram(vd, ve, vf, vdn, vdh, vd0, ve0, vf0, vdn0, vdh0)
    if os.path.exists("__pycache__/"):
        os.system("rm -rf __pycache__/")
        os.system("rm -rf py*.log")
