import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.append("../../sd_cartopy/")
import os
import numpy as np
import datetime as dt
import pandas as pd
from scipy.stats import median_absolute_deviation as MAD

import pydarn
import rad_fov
from pysolar.solar import get_altitude

def get_SZA(lat, lon, d):
    d = d.replace(tzinfo=dt.timezone.utc)
    sza = 90. - get_altitude(lat, lon, d)
    return sza

def plot_sza_versus_height_contour_plots(rads=["bks"], dates=[dt.datetime(2015,5,5,22,11)], tis=[17,18,19], dint=3):
    elvs = np.linspace(16,30,15)
    V, H, E, S = [], [], [], []
    for rad in rads:
        print(rad)
        hdw = pydarn.read_hdw_file(rad)
        rf = rad_fov.CalcFov(hdw=hdw, ngates=60)
        lons, lats = rf.lonFull, rf.latFull
        for d in dates:
            for bm in range(hdw.beams):
                for ti in tis:
                    for e in elvs:
                        fname = "../../data/op/%s/waccmx/%s/bm.%02d/ti(%02d)_elv(%d)_f.csv"%(d.strftime("%Y.%m.%d.%H.%M"), 
                                rad, bm, ti, e)
                        if os.path.exists(fname):
                            out = pd.read_csv(fname)
                            v, h = np.max(np.abs(out.dop)), out.height.tolist()[np.argmax(np.abs(out.dop))]
                            V.append(v) 
                            H.append(h) 
                            E.append(int(e))
                            S.append(get_SZA(lats[bm,:].mean(), lons[bm,:].mean(), d))
    df = pd.DataFrame()
    df["h"] = H
    df["e"] = E
    df["s"] = S
    df["s"] = np.rint(df.s/dint)*dint
    md = df.groupby(by="s").agg([np.median]).reset_index()
    sd = df.groupby(by="s").agg([MAD]).reset_index()
    print(md.head())
    fig = plt.figure(dpi=150, figsize=(4,3))
    ax = fig.add_subplot(111)
    ax.set_xlabel(r"SZA, $\chi$ ($^o$)")
    ax.set_ylabel(r"Height of max($\delta\eta_o$), km")
    ax.set_xlim(30,100)
    ax.set_ylim(70,210)
    ax.errorbar(md.s, md.h["median"], yerr=0.15*sd.h["median_absolute_deviation"], xerr=dint/2, fmt="o", ecolor="r", color="b", 
            capthick=0.5, lw=0.5, ms=2., capsize=1)
    fig.savefig("figures/sza_dist.png", bbox_inches="tight")
    return

def plot():
    def qu(m, a=0.9):
        return np.quantile(m, a)
    dint=3
    x = pd.read_csv("../op/finescale_simulate_total_A.csv")
    x = x[(np.abs(x.vT)>30) & (np.abs(x.vT)<300)]
    x.vF = x.vF + x.vFh
    x = x[["vD", "vE", "vF", "sza"]]
    x["sza"] = np.rint(x.sza/dint)*dint
    x = x.groupby(by="sza").agg([qu, MAD]).reset_index()
    print(x)
    fig, axes = plt.subplots(dpi=150, figsize=(4,3), nrows=1, ncols=1, sharex=True)
    ax = axes
    ax.set_xlabel(r"SZA, $\chi$ ($^o$)")
    ax.set_ylabel(r"$V_x$, $ms^{-1}$")
    ax.set_ylim(0,130)
    t = 0.4
    ct, lw, cs = 0.3, 0.3, 0.7
    ax.errorbar(x.sza, x.vD["qu"], yerr=t*x.vD["median_absolute_deviation"], xerr=dint/2, fmt="o", ecolor="k", color="r",
            capthick=ct, lw=lw, ms=2., capsize=cs, label=r"$V_D$")
    ax.errorbar(x.sza, x.vE["qu"], yerr=t*x.vE["median_absolute_deviation"], xerr=dint/2, fmt="o", ecolor="k", color="g",
            capthick=ct, lw=lw, ms=2., capsize=cs, label=r"$V_E$")
    ax.errorbar(x.sza, x.vF["qu"]-40, yerr=t*x.vF["median_absolute_deviation"], xerr=dint/2, fmt="o", ecolor="k", color="b",
            capthick=ct, lw=lw, ms=2., capsize=cs, label=r"$V_F$")
    ax.legend(loc=1)
    ax.set_xlim(30, 100)
    fig.savefig("figures/sza_dist.png", bbox_inches="tight")
    return

if __name__ == "__main__":
    rads = ["bks", "cvw", "cve", "fhw", "fhe", "wal", "pgr", "kap", "sas"]
    dates = [dt.datetime(2015,5,5,22,11), dt.datetime(2015,3,11,16,22), dt.datetime(2015,9,28,14,58),
            dt.datetime(2016,4,18,0,29), dt.datetime(2014,2,25,0,49), dt.datetime(2014,6,10,11,42)]
    #plot_sza_versus_height_contour_plots(rads, dates)
    plot()
    os.system("rm -rf __pycache__")
    pass
