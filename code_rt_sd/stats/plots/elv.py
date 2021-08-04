import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys

import os
import numpy as np
import datetime as dt
import pandas as pd
from scipy.stats import median_absolute_deviation as MAD


def plot_elv_versus_height_contour_plots(rads=["bks"], dates=[dt.datetime(2015,5,5,22,11)], tis=[17,18,19], beams=np.arange(24)):
    elvs = np.linspace(16,30,15)
    V, H, E = [], [], []
    for rad in rads:
        for d in dates:
            for bm in beams:
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

    df = pd.DataFrame()
    df["h"] = H
    df["e"] = E
    md = df.groupby(by="e").agg([np.median]).reset_index()
    sd = df.groupby(by="e").agg([MAD]).reset_index()
    fig = plt.figure(dpi=150, figsize=(4,3))
    ax = fig.add_subplot(111)
    ax.set_xlabel(r"Elevation Angle, ($^o$)")
    ax.set_ylabel(r"Height of max($\delta\eta_o$), km")
    ax.set_xlim(15,35)
    ax.set_ylim(70,200)
    ax.errorbar(md.e, md.h["median"], yerr=0.25*sd.h["median_absolute_deviation"], xerr=0.5, fmt="o", ecolor="r", color="b", 
            capthick=0.5, lw=0.5, ms=2., capsize=1)
    #ax.text(0.01, 1.05, "Date:%s UT"%d.strftime("%Y-%m-%d %H:%M"), ha="left", va="center", transform=ax.transAxes)
    #ax.text(0.1,0.9,r"$f_o$=12.4 MHz", ha="left", va="center", transform=ax.transAxes)
    fig.savefig("figures/elv_dist.png", bbox_inches="tight")
    return


def plot():
    def qu(m, a=0.9):
        return np.quantile(m, a)
    x = pd.read_csv("../op/finescale_simulate_total_A.csv")
    x = x[(np.abs(x.vT)>30) & (np.abs(x.vT)<300)]
    x.vF = x.vF + x.vFh
    x = x[["vD", "vE", "vF", "elv"]]
    x = x.groupby(by="elv").agg([qu, MAD]).reset_index()
    print(x.head())
    fig, axes = plt.subplots(dpi=150, figsize=(4,3), nrows=1, ncols=1, sharex=True)
    ax = axes
    ax.set_xlabel(r"Elevation Angle, ($^o$)")
    ax.set_ylabel(r"$V_x$, $ms^{-1}$")
    ax.set_ylim(30,90)
    t = 0.4
    ct, lw, cs = 0.3, 0.3, 0.7
    ax.errorbar(x.elv, x.vD["qu"], yerr=t*x.vD["median_absolute_deviation"], xerr=0.25, fmt="o", ecolor="k", color="r",
            capthick=ct, lw=lw, ms=2., capsize=cs, label=r"$V_D$")
    ax.errorbar(x.elv, x.vE["qu"], yerr=t*x.vE["median_absolute_deviation"], xerr=0.25, fmt="o", ecolor="k", color="g",
            capthick=ct, lw=lw, ms=2., capsize=cs, label=r"$V_E$")
    ax.errorbar(x.elv, x.vF["qu"]-50, yerr=t*x.vF["median_absolute_deviation"], xerr=0.25, fmt="o", ecolor="k", color="b",
            capthick=ct, lw=lw, ms=2., capsize=cs, label=r"$V_F$")
    ax.legend(loc=1)
    ax.set_xlim(15,31)
    fig.savefig("figures/elv_dist.png", bbox_inches="tight")
    return 
if __name__ == "__main__":
    rads = ["bks", "cvw", "cve", "fhw", "fhe", "wal", "pgr", "kap", "sas"]
    dates = [dt.datetime(2015,5,5,22,11), dt.datetime(2015,3,11,16,22), dt.datetime(2015,9,28,14,58),
            dt.datetime(2016,4,18,0,29), dt.datetime(2014,2,25,0,49), dt.datetime(2014,6,10,11,42)]
    plot()
    os.system("rm -rf __pycache__")
    pass
