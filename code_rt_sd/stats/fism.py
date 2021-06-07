#!/usr/bin/env python

"""fism.py: module is dedicated to plottting FISM related EUV, XUV and X-rays."""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2021, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import os
import requests
import datetime as dt
import pandas as pd
import argparse
from scipy import signal
from dateutil import parser as dparser

def compare_TS_by_bin(dns, bins):
    """ 
    Compare FISM TS on several days for a number of bins
    parameters:
    -----------
    dns <list[(datetime, datetime)]>: List of touple of datetime needed to be compared (start, end)
    bins [list[float]]: Comparisons on wavebands
    """
    fig, axes = plt.subplots(figsize=(6,7), dpi=120, nrows=2, ncols=1, sharex=True)
    colors = ["darkred", "darkblue"]
    lws = [1.2, .8]
    secs_to_min = 60
    steps_sec = 3
    frac = int(secs_to_min/steps_sec)
    for y, b in enumerate(bins):
        ax = axes[y]
        ax.set_xlabel(r"Time to Flare Peak")
        ax.set_ylabel(r"Intensity ($\Phi_0$), $Wm^{-2}nm^{-1}$")
        for z, dn in enumerate(dns):
            start, end = dn - dt.timedelta(minutes=20), dn + dt.timedelta(minutes=10)
            fname = "op/tmpfism.csv"
            uri = "https://lasp.colorado.edu/lisird/latis/dap/fism_flare_hr.csv?&"+\
                    "time>={:d}-{:02d}-{:02d}T{:02d}:{:02d}:00.000Z&time<={:d}-{:02d}-{:02d}T{:02d}:{:02d}:00.000Z&".format(start.year,
                            start.month, start.day, start.hour, start.minute, end.year,
                            end.month, end.day, end.hour, end.minute)+\
                    "wavelength~{:f}".format(b)
            resp = requests.get(uri)
            with open(fname, "w") as f: f.write(resp.text)
            data = pd.read_csv(fname)["irradiance (W/m^2/nm)"]
            wv = pd.read_csv(fname)["wavelength (nm)"].tolist()[0]
            ax.semilogy(np.arange(len(data))-20, data, color=colors[z], alpha=0.8, lw=lws[z],
                    label=dn.strftime("%Y-%m-%d %H:%M UT"))
            ax.axvline(np.argmax(np.diff(data))-20, color=colors[z], alpha=0.8, ls="--")
            ax.text(1.05, 0.5, "$\lambda$=%.2f nm"%wv, ha="center", va="center",
                    rotation=90, transform=ax.transAxes)
            os.remove(fname)
        ax.legend(loc=2)
    fig.savefig("op/compare_TS_by_bin.png", bbox_inches="tight")
    return

def compare_FISM_diff_spect_by_dn(dns, _f107):
    """
    Compare FISM spectrum for several dates
    paraeters:
    ----------
    dns <List[datetime]>: List of datetime of the events
    """
    fig = plt.figure(figsize=(6,4), dpi=120)
    ax = fig.add_subplot(111)
    ax.set_xlabel(r"Waveband ($\lambda$), nm")
    ax.set_ylabel(r"Intensity ($\Phi_0$), $Wm^{-2}nm^{-1}$")
    colors = ["darkred", "darkblue"]
    lws = [1.2, .8]
    bins = np.arange(0.05,190.05,0.1)
    for z,dn in enumerate(dns):
        fname = "op/tmpfism.csv"
        spec = {}
        for ix, dly in enumerate([20,0]):
            spec[ix] = []
            start, end = dn - dt.timedelta(minutes=1+dly), dn + dt.timedelta(minutes=1+dly)
            for b in bins:
                uri = "https://lasp.colorado.edu/lisird/latis/dap/fism_flare_hr.csv?&"+\
                        "time>={:d}-{:02d}-{:02d}T{:02d}:{:02d}:00.000Z&time<={:d}-{:02d}-{:02d}T{:02d}:{:02d}:00.000Z&".format(start.year,
                                start.month, start.day, start.hour, start.minute, end.year,
                                end.month, end.day, end.hour, end.minute)+\
                                        "wavelength~{:.02f}".format(b)
                resp = requests.get(uri)
                print(uri)
                with open(fname, "w") as f: f.write(resp.text)
                data = pd.read_csv(fname)
                data["time"] = [dt.datetime(1970,1,1)+dt.timedelta(seconds=x) for x in data["time (seconds since 1970-01-01)"]]
                spec[ix].append(data["irradiance (W/m^2/nm)"].tolist()[1])
                os.remove(fname)
        spec = np.array(spec[1]) - np.array(spec[0])
        ax.loglog(bins, spec, color=colors[z], alpha=0.8, lw=lws[z], 
                label=dn.strftime("%Y-%m-%d %H:%M UT, F107=") + "%.1f"%_f107[z])
        ax.set_xlim(bins[0], bins[-1])
        if z==0: 
            ax.axhline(2.2e-4,color="k",ls="--")
            ax.axvline(0.1,color="k",ls="--")
            ax.axvline(0.8,color="k",ls="--")
    ax.legend(loc=3)
    fig.savefig("op/compare_FISM_diff_spect.png", bbox_inches="tight")
    return

def compare_FISM_spect_by_dn(dns, _f107):
    """
    Compare FISM spectrum for several dates
    paraeters:
    ----------
    dns <List[datetime]>: List of datetime of the events
    """
    fig = plt.figure(figsize=(6,4), dpi=120)
    ax = fig.add_subplot(111)
    ax.set_xlabel(r"Waveband ($\lambda$), nm")
    ax.set_ylabel(r"Intensity ($\Phi_0$), $Wm^{-2}nm^{-1}$")
    colors = ["darkred", "darkblue"]
    lws = [1.2, .8]
    bins = np.arange(0.05,190.05,0.1)
    for z,dn in enumerate(dns):
        fname = "op/tmpfism.csv"
        start, end = dn - dt.timedelta(minutes=1), dn + dt.timedelta(minutes=1)
        spec = []
        for b in bins:
            uri = "https://lasp.colorado.edu/lisird/latis/dap/fism_flare_hr.csv?&"+\
                    "time>={:d}-{:02d}-{:02d}T{:02d}:{:02d}:00.000Z&time<={:d}-{:02d}-{:02d}T{:02d}:{:02d}:00.000Z&".format(start.year,
                            start.month, start.day, start.hour, start.minute, end.year,
                            end.month, end.day, end.hour, end.minute)+\
                                    "wavelength~{:.02f}".format(b)
            resp = requests.get(uri)
            print(uri)
            with open(fname, "w") as f: f.write(resp.text)
            data = pd.read_csv(fname)
            data["time"] = [dt.datetime(1970,1,1)+dt.timedelta(seconds=x) for x in data["time (seconds since 1970-01-01)"]]
            spec.append(data["irradiance (W/m^2/nm)"].tolist()[1])
            os.remove(fname)
        ax.loglog(bins, spec, color=colors[z], alpha=0.8, lw=lws[z], 
                label=dn.strftime("%Y-%m-%d %H:%M UT, F107=") + "%.1f"%_f107[z])
        ax.set_xlim(bins[0], bins[-1])
        if z==0: 
            ax.axhline(2.2e-4,color="k",ls="--")
            ax.axvline(0.1,color="k",ls="--")
            ax.axvline(0.8,color="k",ls="--")
    ax.legend(loc=2)
    fig.savefig("op/compare_FISM_spect.png", bbox_inches="tight")
    return

def disk_limb_histograms(dns):
    import glob
    bins = np.arange(0,1,.05)
    fig, axes = plt.subplots(figsize=(6,6), dpi=150, sharex="all", sharey="all", nrows=2, ncols=2)
    hs = ["step", "step"]
    for z, dn in enumerate(dns):
        df = pd.DataFrame()
        files = glob.glob("op/*%s.csv"%dn.strftime("%Y-%m-%d-%H-%M"))
        for f in files:
            df = pd.concat([df, pd.read_csv(f)])
        df = df[(df.vT>30) & (df.vT<300)]
        vd, ve, vf = np.array(df.vD/df.vT), np.array(df.vE/df.vT), np.array((df.vF + df.vFh)/df.vT)
        vn, vh = np.array((df.vD + df.vE + df.vF)/df.vT), np.array(df.vFh/df.vT)
        
        ax = axes[z,0]
        vnmean, vhmean = np.mean(vn), np.mean(vh)
        if vnmean > 0.9 : vnmean, vhmean = 0.9, 0.1
        ax.axvline(vnmean, ls="--", color="r", lw=0.6)
        ax.hist(vn, bins=bins, color="r", alpha=0.5, density=True, label=r"$\frac{V_{d\eta}}{V_T}[\mu=%.2f]$"%vnmean, histtype=hs[z])
        ax.axvline(vhmean, ls="--", color="b", lw=0.6)
        ax.hist(vh, bins=bins, color="b", alpha=0.5, density=True, label=r"$\frac{V_{dh}}{V_T}[\mu=%.2f]$"%vhmean, histtype=hs[z])
        ax.text(0.1,0.9, "(a.%d)"%z, horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
        ax.set_xlim(0,1)
        ax.legend(loc=1, prop={"size": 8})
        ax.set_ylabel(r"Density $\left(\frac{V_x}{V_T}\right)$")
        if z==1: ax.set_xlabel(r"$\frac{V_x}{V_T}$")
        ax = axes[z,1]
        ax.axvline(np.mean(vd), ls="--", color="r", lw=0.6)
        ax.hist(vd, bins=bins, color="r", alpha=0.5, density=True, label=r"$\frac{V_D}{V_T}[\mu=%.2f]$"%np.mean(vd), histtype=hs[z])
        ax.axvline(np.mean(ve), ls="--", color="g", lw=0.6)
        ax.hist(ve, bins=bins, color="g", alpha=0.5, density=True, label=r"$\frac{V_E}{V_T}[\mu=%.2f]$"%np.mean(ve), histtype=hs[z])
        ax.axvline(np.mean(vf), ls="--", color="b", lw=0.6)
        ax.hist(vf, bins=bins, color="b", alpha=0.5, density=True, label=r"$\frac{V_F}{V_T}[\mu=%.2f]$"%np.mean(vf), histtype=hs[z])
        ax.text(0.1,0.9, "(b.%d)"%z, horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
        ax.text(1.05,0.5, dn.strftime("%Y-%m-%d %H:%M UT"), ha="center", va="center", transform=ax.transAxes, rotation=90)
        ax.set_xlim(0,1)
        ax.legend(loc=1, prop={"size": 8})
        if z==1: ax.set_xlabel(r"$\frac{V_x}{V_T}$")
    fig.savefig("op/disk_limb_histograms.png", bbox_inches="tight")
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-ev", "--event", default=dt.datetime(2015,5,5,22,11), nargs="+", help="Start date", type=dparser.isoparse)
    parser.add_argument("-b", "--bins", type=int, default=10, help="Wavebins in nm")
    parser.add_argument("-c", "--code", default="")
    args = parser.parse_args()
    args.event = [dt.datetime(2014,6,10,11,42), dt.datetime(2015,3,11,16,20)]
    args.bins = [10,30]
    print("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     " , k , "->" , str(vars(args)[k]))
    f107 = pd.read_csv("data.csv")
    #f107.time = [dt.datetime.strptime(x, "%Y%m%d") for x in f107.time]
    _f107 = []
    for d in args.event:
        d = int(d.strftime("%Y%m%d"))
        _f107.append(f107[f107.time==d]["f107"].tolist()[0])
    if args.code == "hist":
        disk_limb_histograms(args.event)
    else:
        #compare_FISM_spect_by_dn(args.event, _f107)
        compare_FISM_diff_spect_by_dn(args.event, _f107)
        #compare_TS_by_bin(args.event, args.bins)
    pass
