import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import os
import sys
sys.path.append("../../sd")
sys.path.append("../../sd_cartopy")
import get_sd_data as gsd
import utils

import numpy as np
import scipy.stats as sps
import datetime as dt
import pandas as pd
import pytz
from netCDF4 import Dataset

from timezonefinder import TimezoneFinder
from pysolar.solar import get_altitude

tf = TimezoneFinder(in_memory=True)

def _get_ij_(lats, lons, lat, lon):
    """ Get (lat, lon) index """
    _ij_ = (np.argmin(np.abs(lats-lat)), np.argmin(np.abs(lons-lon)))
    return _ij_

def plot_con_ui_altitude_profile(f0="../../data/op/2015.05.05.22.11/waccmx/2015.05.05.21.51.nc.gz", 
        f1="../../data/op/2015.05.05.22.11/waccmx/2015.05.05.22.09.nc.gz"):
    os.system("gzip -d " + f0)
    f0  = f0.replace(".gz","")
    ds0 = Dataset(f0)
    os.system("gzip " + f0)
    os.system("gzip -d " + f1)
    f1  = f1.replace(".gz","")
    ds1 = Dataset(f1)
    os.system("gzip " + f1)
    sdlat, sdlon = 37.93, -75.47
    isrlat, isrlon = -11.95, -76.87
    altscale = 1e-3
    date = dt.datetime(2015,5,5,22,11)
    
    lats, lons = ds1.variables["lat"][:], ds1.variables["lon"][:]
    i, j = _get_ij_(lats, lons, sdlat, sdlon)
    _i, _j = _get_ij_(lats, lons, isrlat, isrlon)
    fig = plt.figure(dpi=150, figsize=(6,3))
    ax = fig.add_subplot(121)
    vi = -1*(ds1.variables["WIf"][:][0, :, i, j] - ds0.variables["WIf"][:][0, :, i, j])
    ax.plot(vi, ds1.variables["ZGf"][0, :, i, j]*altscale, color="r", label="$\Delta\omega_I^{SD}$", ls="-", lw=1.2)
    vi = -1*(ds1.variables["WIf"][:][0, :, _i, _j] - ds0.variables["WIf"][:][0, :, _i, _j])
    ax.plot(vi, ds1.variables["ZGf"][0, :, _i, _j]*altscale, color="b", label="$\Delta\omega_I^{ISR}$", ls="-", lw=0.8)
    ax.set_ylim(90, 200)
    ax.set_xlim(0,3)
    ax.legend(loc=1)
    ax.set_xlabel(r"$\Delta\omega_I=\omega_I^f-\omega_I^b$, $ms^{-1}$")
    ax.set_ylabel(r"Height, km")
    ax.text(0.01, 1.05, "Date: %s UT"%date.strftime("%Y-%m-%d %H:%M"), ha="left", va="center", transform=ax.transAxes)

    ax = fig.add_subplot(122)
    r1 = ds1.variables["SHf"][:][0, :, i, j]/ds1.variables["SPf"][:][0, :, i, j]
    r2 = ds0.variables["SHd"][:][0, :, i, j]/ds0.variables["SPd"][:][0, :, i, j]
    ax.plot(utils.smooth(r1-r2, 21), ds1.variables["ZGf"][0, :, i, j]*altscale, color="r", label=r"$\rho^{SD}$", ls="-", lw=1.2)
    r1 = ds1.variables["SHf"][:][0, :, _i, _j]/ds1.variables["SPf"][:][0, :, _i, _j]
    r2 = ds0.variables["SHd"][:][0, :, _i, _j]/ds0.variables["SPd"][:][0, :, _i, _j]
    ax.plot(utils.smooth(r1-r2, 21), ds1.variables["ZGf"][0, :, _i, _j]*altscale, color="b", label=r"$\rho^{ISR}$", ls="-", lw=1.2)
    ax.set_xlabel(r"$\rho=\left(\frac{\sigma_H}{\sigma_P}\right)^{f}-\left(\frac{\sigma_H}{\sigma_P}\right)^{b}$")
    ax.set_ylim(90, 200)
    ax.set_xlim(-.2, .2)
    ax.axvline(0, color="k", ls="--", lw=0.8)
    ax.legend(loc=1)
    ax.text(0.99, 1.05, "Flare: X2.7", ha="right", va="center", transform=ax.transAxes)
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    fig.savefig("figures/alt_prof.png", bbox_inches="tight")
    return

def get_SZA(lat, lon, d):
    d = d.replace(tzinfo=dt.timezone.utc)
    sza = 90. - get_altitude(lat, lon, d)
    return sza

def get_SZA_array(latx, lonx, d):
    szaz = np.zeros((len(latx), len(lonx)))
    for i in range(len(latx)):
        for j in range(len(lonx)):
            szaz[i,j] = get_SZA(latx[i], lonx[j], d)
    return szaz

def global_distr(f0="../../data/op/2015.05.05.22.11/waccmx/2015.05.05.21.51.nc.gz",
        f1="../../data/op/2015.05.05.22.11/waccmx/2015.05.05.22.09.nc.gz"):
    import fov
    import matplotlib.colors as colors
    os.system("gzip -d " + f0)
    f0  = f0.replace(".gz","")
    ds0 = Dataset(f0)
    os.system("gzip " + f0)
    os.system("gzip -d " + f1)
    f1  = f1.replace(".gz","")
    ds1 = Dataset(f1)
    os.system("gzip " + f1)
    sdlat, sdlon = 37.93+2, -75.47-2
    isrlat, isrlon = -11.95, -76.87
    altscale = 1e-3
    date = dt.datetime(2015,5,5,22,11)

    lats, lons = ds1.variables["lat"][:], ds1.variables["lon"][:]
    sza = get_SZA_array(lats, lons, date)
    fig = plt.figure(dpi=150, figsize=(15,7))
    ax, fig, _to, _from = fov.get_globe(fig, 121, date)
    xx, yy = np.meshgrid(lons, lats)
    i, j = _get_ij_(lats, lons, sdlat, sdlon)
    h = np.argmin(np.abs(ds1.variables["ZGf"][0, :, i, j]*altscale-120))
    Z = -1*(ds1.variables["WIf"][:][0, h, :, :] - ds0.variables["WIf"][:][0, h, :, :])
    c = ax.pcolormesh(xx, yy, Z, cmap="Spectral")
    cb = fig.colorbar(c, ax=ax, shrink=0.4)
    cb.set_label(r"$\Delta\omega_I$, $ms^{-1}$")
    ax.scatter([sdlon],[sdlat],color="g",s=15,marker="D")
    ax.scatter([isrlon],[isrlat],color="r",s=15,marker="D")
    cs = ax.contour(lons, lats, sza, linestyles="dotted")
    ax.clabel(cs, inline=True, fontsize=6, fmt=r"$%d^{o}$")
    ax.text(-0.02, 0.5, "Geographic Coordinates", horizontalalignment="right",
            verticalalignment="center", transform=ax.transAxes, rotation=90)
    #ax.set_xlim(np.min(lons), np.max(lats))
    ax.text(0.01, 1.05, "Date: %s UT"%date.strftime("%Y-%m-%d %H:%M"), ha="left", va="center", transform=ax.transAxes)

    ax, fig, _to, _from = fov.get_globe(fig, 122, date)
    r1 = ds1.variables["SHf"][:][0, h, :, :]/ds1.variables["SPf"][:][0, h, :, :]
    r2 = ds0.variables["SHd"][:][0, h, :, :]/ds0.variables["SPd"][:][0, h, :, :]
    Z = r1-r2
    c = ax.pcolormesh(xx, yy, Z, cmap="RdBu")
    cb = fig.colorbar(c, ax=ax, shrink=0.4)
    cb.set_label(r"$\rho$")
    ax.scatter([sdlon],[sdlat],color="g",s=15,marker="D")
    ax.scatter([isrlon],[isrlat],color="r",s=15,marker="D")
    ax.text(0.99, 1.05, "Flare: X2.7", ha="right", va="center", transform=ax.transAxes)
    fig.savefig("figures/global_dist.png", bbox_inches="tight")
    return


def _add_colorbar(fig, ax, bounds, colormap, label=""):
    import matplotlib as mpl
    pos = ax.get_position()
    cpos = [pos.x1 + 0.1, pos.y0 + 0.0125,
            0.015, pos.height * 0.8]                # this list defines (left, bottom, width, height
    cax = fig.add_axes(cpos)
    norm = mpl.colors.BoundaryNorm(bounds, colormap.N)
    cb2 = mpl.colorbar.ColorbarBase(cax, cmap=colormap,
            norm=norm,
            ticks=bounds,
            spacing="uniform",
            orientation="vertical")
    cb2.set_label(label)
    return

def plot_jro_sd_goes_TS():
    st, et = dt.datetime(2005,9,7,14), dt.datetime(2005,9,7,22)
    u = pd.read_csv("../isr-jro-data/g12_xrs_1m_20050901_20050930.csv", skiprows=117, parse_dates=["time_tag"])
    u = u[(u.time_tag>=st) & (u.time_tag<et)]
    fig = plt.figure(figsize=(7,8), dpi=120)
    ax = fig.add_subplot(311)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.semilogy(u.time_tag, u["xs"], "b.", label=r"$\lambda_s(0.05-0.3 nm)$", ms=1, lw=0.8)
    ax.semilogy(u.time_tag, u["xl"], "r.", label=r"$\lambda_s(0.1-0.8 nm)$", ms=1, lw=0.8)
    ax.set_xlim(u.time_tag.tolist()[0], u.time_tag.tolist()[-1])
    ax.set_ylim(1e-8,1e-2)
    ax.set_xlim(st, et)
    ax.axvline(dt.datetime(2005,9,7,17,20), color="b", ls="--", lw=0.8)
    ax.axvline(dt.datetime(2005,9,7,17,37), color="r", ls="--", lw=0.8)
    ax.text(0.01, 1.05, "Date: 7 Sep., 2005", ha="left", va="center", transform=ax.transAxes)
    ax.text(0.9, 0.1, "(a)", ha="center", va="center", transform=ax.transAxes)
    ax.set_ylabel(r"Solar Irradiance, $Wm^{-2}$")
    ax.legend(loc=2)
    ax = fig.add_subplot(313)
    ax.text(0.9, 0.1, "(c)", ha="center", va="center", transform=ax.transAxes)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_ylabel(r"$\omega_I$, $ms^{-1}$")
    labels = ["7 Sep., 2005", "8 Sep., 2005"]
    for f, c, lb in zip(["../isr-jro-data/jul20050907_avg_150km.001.txt", "../isr-jro-data/jul20050908_avg_150km.001.txt"], ["b", "g"], 
            labels):
        dates, lat, lon, vipn2, dvipn2, vipe1, dvipe1 = [], [], [], [], [], [], []
        with open(f, "r") as f:
            lines = f.readlines()
            for l in lines[1:]:
                l = list(filter(None, l.replace("\n", "").replace("missing","NaN").split(" ")))
                d = dt.datetime(int(l[0]), int(l[1]), int(l[2]), int(l[3]), int(l[4]), int(l[5]))
                dates.append(d)
                lat.append(float(l[6]))
                lon.append(float(l[7]))
                vipn2.append(float(l[11]))
                dvipn2.append(float(l[12]))
                vipe1.append(float(l[13]))
                dvipe1.append(float(l[14]))
        u = pd.DataFrame()
        u["date"], u["lat"], u["lon"], u["vipn2"], u["dvipn2"], u["vipe1"], u["dvipe1"] = dates, lat, lon, vipn2, dvipn2, vipe1, dvipe1
        local_time_zone = tf.timezone_at(lng=u.lon.mean(), lat=u.lat.mean())
        timezone = pytz.timezone(local_time_zone)
        u.date = [timezone.localize(d).astimezone(pytz.utc) + dt.timedelta(hours=-5) for d in u.date]
        u.date = [d.replace(day=7) for d in u.date]
        ax.errorbar(u.date, u.vipn2, yerr=u.dvipn2, fmt="o", ecolor="r", color=c, capthick=0.5, lw=0.5, ms=1., capsize=1, label=lb)
    ax.set_xlim(st, et)
    ax.set_ylim(-10,30)
    ax.legend(loc=3)
    ax.axvline(dt.datetime(2005,9,7,17,20), color="b", ls="--", lw=0.8)
    ax.axvline(dt.datetime(2005,9,7,17,37), color="r", ls="--", lw=0.8)
    ax.text(0.99, 0.9, "Radar: JRO, JULIA", ha="right", va="center", transform=ax.transAxes)
    ax.set_xlabel("Time, UT")
    ax = fig.add_subplot(312)
    ax.text(0.9, 0.1, "(b)", ha="center", va="center", transform=ax.transAxes)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    fd = gsd.FetchData("wal", [st, et])
    beams, _ = fd.fetch_data(by="beams")
    u = fd.convert_to_pandas(beams)
    u = u[(u.slist>=10) & (np.abs(u.v_e)<100.)]
    print(u.head())
    X, Y, Z = utils.get_gridded_parameters(u, xparam="time", yparam="slist", zparam="v")
    bounds = list(range(-100, 100+1, 25))
    cmap = plt.cm.jet
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    Z = utils.medfilt2D_weight(Z, tau=0.99)
    #ax.pcolormesh(X, Y, Z.T, lw=0.01, edgecolors="None", cmap=cmap, norm=norm)
    ax.text(0.99, 0.9, "Radar: WAL, SD", ha="right", va="center", transform=ax.transAxes)
    #ax.set_ylabel("Range Gate")
    #ax.set_xlim(st, et)
    #_add_colorbar(fig, ax, bounds, cmap, label="")
    #ax = ax.twinx()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    M = np.nanmedian(Z.T, axis=0)
    u = pd.DataFrame(); u["date"], u["m"]= X[0,:], M
    md, sd = u.set_index("date").resample("60s").mean().reset_index(), u.set_index("date").resample("60s").std().reset_index()
    ax.errorbar(md.date, md.m, yerr=0.3*sd.m, fmt="o-", ecolor="r", color="b", 
            capthick=0.5, lw=0.8, ms=1., capsize=1, label=lb, alpha=0.2)
    ax.set_xlabel("Time, UT")
    ax.set_xlim(st, et)
    ax.set_ylim(-30,60)
    ax.axvline(dt.datetime(2005,9,7,17,20), color="b", ls="--", lw=0.8)
    ax.axvline(dt.datetime(2005,9,7,17,37), color="r", ls="--", lw=0.8)
    ax.set_ylabel(r"$V^D_{LoS}$, $ms^{-1}$")
    fig.autofmt_xdate()
    fig.savefig("figures/sd_jro.png", bbox_inches="tight")
    return


if __name__ == "__main__":
    plot_jro_sd_goes_TS()
    #plot_con_ui_altitude_profile()
    #global_distr()
    os.system("rm -rf *.log")
    pass
