import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.append("../../sd")
sys.path.append("../../sd_cartopy")

import os
import glob
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy import array
import datetime as dt

from pysolar.solar import get_altitude
import sza
import pandas as pd
from scipy.stats import median_absolute_deviation as MAD

def plot_freq_versus_height_contour_plots(fig, ax, freqs=np.linspace(5,14,101)*1e6, loc=[10,-100], date=dt.datetime(2015,5,5,22,11)):
    mdata = sza.from_pickle("data_DL.pkl")
    cent, cent0 = mdata["cent"], mdata["cent0"]
    #fig = plt.figure(dpi=150, figsize=(4,3))
    #ax = fig.add_subplot(111)
    i,j = sza._get_ij_(cent["latx"], cent["lonx"], loc[0], loc[1])
    alts = cent["alts"]
    nds = np.zeros((len(freqs), len(alts)))
    for _i, f in enumerate(freqs):
        c_indx = np.argmin(np.abs(sza.calculate_fp(cent["nex"][:,i,j])-f))
        n_diff = -1*(sza.calculate_eta(cent["nex"][:,i,j], cent["alts"], loc[0], loc[1], date, f) -\
                sza.calculate_eta(cent0["nex"][:,i,j], cent0["alts"], loc[0], loc[1], date, f))
        n_diff[c_indx:] = np.nan
        nds[_i,:] = 2*f*n_diff/3e8
        if _i==-1:
            print(n_diff)
            print("Freq: - ", f/1e6)
            break
    print(alts)
    xx, yy = np.meshgrid(freqs, alts) 
    c = ax.pcolormesh(xx/1e6, yy, nds.T, cmap="Reds")
    ax.set_xlabel(r"$f_0$, MHz")
    ax.set_ylabel("Height, km")
    ax.set_ylim(100,300)
    ax.set_xlim(5,14)
    cb = fig.colorbar(c, ax=ax, shrink=0.8)
    cb.set_label(r"$V_x (ms^{-1})$")
    ax.text(0.01, 1.05, "Date:%s UT"%date.strftime("%Y-%m-%d %H:%M"), ha="left", va="center", transform=ax.transAxes)
    ax.text(0.1,0.9,"$\chi=%.1f^o$"%sza.get_SZA(loc[0], loc[1], date), ha="left", va="center", transform=ax.transAxes)
    return

def plot():
    def qu(m, a=0.9):
        return np.quantile(m, a)
    dint=0.1
    x = pd.read_csv("../op/finescale_simulate_total_A.csv")
    x = x[(np.abs(x.vT)>30) & (np.abs(x.vT)<300)]
    x.vF = x.vF + x.vFh
    x = x[["vD", "vE", "vF", "freq"]]
    x["freq"] = np.rint(x.freq/dint)*dint
    x = x.groupby(by="freq").agg([qu, MAD]).reset_index()
    fig, axes = plt.subplots(dpi=150, figsize=(8,3), nrows=1, ncols=2)
    ax = axes[0]
    ax.set_xlabel(r"$f_0$, (MHz)")
    ax.set_ylabel(r"$V_x$, $ms^{-1}$")
    ax.set_ylim(0,130)
    t = 0.7
    ct, lw, cs = 0.3, 0.3, 0.7
    ax.errorbar(x.freq, x.vD["qu"], yerr=t*x.vD["median_absolute_deviation"], xerr=dint/2, fmt="o", ecolor="k", color="r",
            capthick=ct, lw=lw, ms=2., capsize=cs, label=r"$V_D$")
    ax.errorbar(x.freq, x.vE["qu"], yerr=t*x.vE["median_absolute_deviation"], xerr=dint/2, fmt="o", ecolor="k", color="g",
            capthick=ct, lw=lw, ms=2., capsize=cs, label=r"$V_E$")
    ax.errorbar(x.freq, x.vF["qu"]-40, yerr=t*x.vF["median_absolute_deviation"], xerr=dint/2, fmt="o", ecolor="k", color="b",
            capthick=ct, lw=lw, ms=2., capsize=cs, label=r"$V_F$")
    ax.legend(loc=1)
    ax.set_xlim(10,18)
    plot_freq_versus_height_contour_plots(fig, axes[1])
    fig.subplots_adjust(wspace=0.3)
    fig.savefig("figures/freq_dist.png", bbox_inches="tight")
    return

if __name__ == "__main__":
    plot()
    os.system("rm -rf __pycache__")
    pass
