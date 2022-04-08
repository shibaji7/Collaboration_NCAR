
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
import datetime as dt
import numpy as np
import pandas as pd
import pytz
from netCDF4 import Dataset
from timezonefinder import TimezoneFinder
import aacgmv2
import traceback

import sys
sys.path.append("sd/")
sys.path.append("sd_cartopy/")
import get_sd_data as gsd
import utils
import rad_fov
import pydarn

tf = TimezoneFinder(in_memory=True)


def convert_geomag(x):
    mlat, mlon, mlt = aacgmv2.get_aacgm_coord(x["lat"], x["lon"], 150, x["date"])
    x["mlat"], x["mlon"], x["mlt"] = np.round( mlat, 2), np.round( mlon, 2), np.round( mlt, 2)
    return x

def get_isr_data(f="stats/isr-jro-data/jul20050907_avg_150km.001.txt"):
    dates, lat, lon, vipn2, dvipn2, vipe1, dvipe1 = [], [], [], [], [], [], []
    with open(f, "r") as f: lines = f.readlines()
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
    u = u.apply(convert_geomag, axis=1)
    return u

def get_index(lats, lons, lat, lon):
    i, j = np.argmin(np.abs(lats-lat)), np.argmin(np.abs(lons-lon))
    return i, j

def fetch_file(o):
    d, mlatsx, mlonsx, glatsx, glonsx = o[0], o[1], o[2], o[3], o[4]
    f = "data/op/2005.09.07.17.40/waccmx/%s.nc.gz"%d.strftime("%Y.%m.%d.%H.%M")
    pc, hc, ed1, ed2, dns, uimax, vimax, wimax, uimin, vimin, wimin = [], [], [], [], [], [], [], [], [], [], []
    o = pd.DataFrame()
    if os.path.exists(f):
        os.system("gzip -d " + f)
        f = f.replace(".gz", "")
        ds = Dataset(f)
        try:
            mlats, mlons = ds.variables["mlat"][:], ds.variables["mlon"][:]
            glats, glons = ds.variables["lat"][:], ds.variables["lon"][:]
            ED1, ED2 = ds.variables["ED1f"][:], ds.variables["ED2f"][:]
            PC, HC = ds.variables["PCf"][:], ds.variables["HCf"][:]
            Z = ds.variables["ZGf"][:]*1e-3
            UI, VI, WI  = ds.variables["UIf"][:], ds.variables["VIf"][:], ds.variables["WIf"][:]
            for mlat, mlon, glat, glon in zip(mlatsx, mlonsx, glatsx, glonsx):
                i, j = get_index(mlats, mlons, mlat, mlon)
                pc.append(PC[0,i,j])
                hc.append(HC[0,i,j])
                ed1.append(ED1[0,i,j])
                ed2.append(ED2[0,i,j])
                dns.append(d)
                i, j = get_index(glats, glons, glat, glon)
                kmax, kmin = np.argmin(np.abs(Z[0,:,i,j]-300)), np.argmin(np.abs(Z[0,:,i,j]-100))
                uimax.append(UI[0,kmax,i,j])
                vimax.append(VI[0,kmax,i,j])
                wimax.append(WI[0,kmax,i,j])
                uimin.append(UI[0,kmin,i,j])
                vimin.append(VI[0,kmin,i,j])
                wimin.append(WI[0,kmin,i,j])
        except: traceback.print_exc()
        ds.close()
        os.system("gzip "+f)
        print(" Done date:", d)
        o["date"], o["PC"], o["HC"], o["ED1"], o["ED2"] = dns, pc, hc, ed1, ed2
        o["UI_max"], o["VI_max"], o["WI_max"] = uimax, vimax, wimax
        o["UI_min"], o["VI_min"], o["WI_min"] = uimin, vimin, wimin
        o["mlat"], o["mlon"], o["glat"], o["glon"] = mlatsx, mlonsx, glatsx, glonsx
    return o

def get_field_data(dates, mlats, mlons, glats, glons, dname):
    if os.path.exists(dname): _o = pd.read_csv(dname, parse_dates=["date"])
    else:
        import multiprocessing
        p = multiprocessing.Pool(24)
        objs = [(d, mlats, mlons, glats, glons) for d in dates]
        _o = pd.DataFrame()
        for o in p.imap(fetch_file, objs):
            if len(o) > 0: _o = pd.concat([_o, o])
        _o.to_csv(dname, index=False, header=True)
    return _o

def get_sd_data_files(dates, dname, rad="wal"):
    hdw = pydarn.read_hdw_file(rad)
    fov = rad_fov.CalcFov(hdw=hdw)
    glats, glons = np.array([fov.latFull.mean()]), np.array([fov.lonFull.mean()])
    mlats, mlons = [], []
    for lat, lon in zip(glats, glons):
        mlat, mlon, _ = aacgmv2.get_aacgm_coord(lat, lon, 250, dates[0])
        mlats.append(mlat)
        mlons.append(mlon)
    mlats, mlons = np.array(mlats), np.array(mlons)
    o = get_field_data(dates, mlats, mlons, glats, glons, dname)
    return o

def get_EC_field(d, key="f"):
    f = "data/op/2005.09.07.17.40/waccmx/%s.nc.gz"%d.strftime("%Y.%m.%d.%H.%M")
    if os.path.exists(f):
        os.system("gzip -d " + f)
        f = f.replace(".gz", "")
        ds = Dataset(f)
        try:
            mlats, mlons = ds.variables["mlat"][:], ds.variables["mlon"][:]
            mlt = aacgmv2.convert_mlt(mlons, d)
            ED1, ED2 = ds.variables["ED1"+key][0, :, :], ds.variables["ED2"+key][0, :, :]
            PC, HC = ds.variables["PC"+key][0, :, :], ds.variables["HC"+key][0, :, :]
        except: traceback.print_exc()
        ds.close()
        os.system("gzip "+f)
        print(" Done date:", d)
    return PC, HC, ED1, ED2, mlats, mlons, mlt

def zips():
    import glob
    files = glob.glob("data/op/2005.09.07.17.40/waccmx/*.nc")
    for f in files:
        print(f)
        os.system("gzip "+f)
    return

case = "field-analysis"
if case == "run-cmd":
    cmd = "nohup python model_sim.py -r wal -ev 2005-09-07T17:40 -s 2005-09-07T17 -e 2005-09-07T18 > /dev/null 2>&1 &"
    os.system(cmd)
elif case == "run-plot":
    cmd = "nohup python model_sim.py -r wal -ev 2005-09-07T17:40 -s 2005-09-07T17 -e 2005-09-07T18 -ps > /dev/null 2>&1 &"
    os.system(cmd)
elif case == "field-analysis": 
    f = "stats/isr-jro-data/jul20050907_avg_150km.001.txt"
    dates = [dt.datetime(2005,9,7,14) + dt.timedelta(minutes=i) for i  in range(421)]
    u = get_isr_data(f)
    o = get_field_data(dates, np.array([u.mlat.tolist()[0]]), np.array([u.mlon.tolist()[0]]), 
            np.array([u.lat.tolist()[0]]), np.array([u.lon.tolist()[0]]), "stats/isr-jro-data/jro_model.csv")
    zips()
    fig = plt.figure(figsize=(7,6), dpi=120)
    ax = fig.add_subplot(211)
    ax.text(0.9,0.6, "(a)", ha="center", va="center", transform=ax.transAxes)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_ylabel(r"$\omega_I$, $ms^{-1}$")
    ax.errorbar(u.date, u.vipn2, yerr=u.dvipn2, fmt="o", ecolor="r", color="b", capthick=0.5, lw=0.5, ms=1., capsize=1, label="JULIA")
    y = 0.5*(o.WI_min+o.WI_max)
    up, lo = o.WI_max-y, y-o.WI_min
    ax.errorbar(o.date[::5], y[::5], yerr=[lo[::5], up[::5]], fmt="o", ecolor="r", color="g", capthick=0.5,
            lw=0.5, ms=1., capsize=1, label="WACCM-X")
    ax.legend(loc=1)
    ax.set_ylim(0,30)
    ax.set_xlim(dt.datetime(2005,9,7,14), dt.datetime(2005,9,7,21))
    ax.axvline(dt.datetime(2005,9,7,17,20), color="b", ls="--", lw=0.8)
    ax.axvline(dt.datetime(2005,9,7,17,37), color="r", ls="--", lw=0.8)
    ax.text(0.01, 0.9, "MLAT, MLON: %.1f,%.1f"%(o.mlat.tolist()[0], o.mlon.tolist()[0]), ha="left", va="center", transform=ax.transAxes)
    ax.text(0.99, 1.05, "Radar: JRO, JULIA", ha="right", va="center", transform=ax.transAxes)
    ax.text(0.01, 1.05, "Date: 2005-09-07", ha="left", va="center", transform=ax.transAxes)
    ax = fig.add_subplot(212)
    dates = [dt.datetime(2005,9,7,14) + dt.timedelta(minutes=i) for i in range(420)]
    o = get_sd_data_files(dates, "stats/isr-jro-data/wal_model.csv", rad="wal")
    o = o.groupby("date").agg([np.mean, np.std]).reset_index()
    y = 0.5*(o.WI_min["mean"]+o.WI_max["mean"])
    up, lo = o.WI_max["mean"]-y, y-o.WI_min["mean"]
    ax.errorbar(o.date[::5], y[::5], yerr=[lo[::5], up[::5]], fmt="o", ecolor="r", color="g", capthick=0.5,
                        lw=0.5, ms=1., capsize=1, label="WACCM-X")
    ax.text(0.01, 0.9, "MLAT, MLON: %.1f,%.1f"%(o.mlat["mean"].tolist()[0], o.mlon["mean"].tolist()[0]),
            ha="left", va="center", transform=ax.transAxes)
    ax.set_ylabel(r"$\omega_I (WACCM-X)$, $ms^{-1}$"+"\n"+r"LoS ($\times 10$), $ms^-1$")
    ax.set_ylim(-5, 10)
    ax.axvline(dt.datetime(2005,9,7,17,20), color="b", ls="--", lw=0.8)
    ax.axvline(dt.datetime(2005,9,7,17,37), color="r", ls="--", lw=0.8)
    ax.text(0.9,0.6, "(b)", ha="center", va="center", transform=ax.transAxes)
    #ax = ax.twinx()
    fd = gsd.FetchData("wal", [dt.datetime(2005,9,7,14), dt.datetime(2005,9,7,21)])
    beams, _ = fd.fetch_data(by="beams")
    u = fd.convert_to_pandas(beams)
    u = u[(u.slist>=10) & (np.abs(u.v_e)<100.)]
    X, Y, Z = utils.get_gridded_parameters(u, xparam="time", yparam="slist", zparam="v")
    Z = utils.medfilt2D_weight(Z, tau=0.99)
    M = np.nanmedian(Z.T, axis=0)
    u = pd.DataFrame(); u["date"], u["m"]= X[0,:], M
    md, sd = u.set_index("date").resample("300s").mean().reset_index(), u.set_index("date").resample("300s").std().reset_index()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.errorbar(md.date, md.m/10, yerr=0.3*sd.m, fmt="o", ecolor="r", color="b", capthick=0.5,
                        lw=0.5, ms=1., capsize=1, label="WAL")
    ax.text(0.99, 1.05, "Radar: WAL, SD", ha="right", va="center", transform=ax.transAxes)
    ax.set_xlim(dt.datetime(2005,9,7,14), dt.datetime(2005,9,7,21))
    #ax.set_ylabel(r"LoS, $ms^{-1}$", fontdict={"color":"b"})
    #ax.set_ylim(-10, 50)
    ax.legend(loc=1)
    ax.set_xlabel("Time, UT")
    fig.autofmt_xdate()
    fig.savefig("stats/plots/figures/jro_field.png", bbox_inches="tight")
elif case == "drift-analysis":
    f = "stats/isr-jro-data/jul20050907_avg_150km.001.txt"
    #dates = [dt.datetime(2005,9,7,0,10) + dt.timedelta(minutes=i) for i  in range(1420)]
    dates = [dt.datetime(2005,9,7,14) + dt.timedelta(minutes=i) for i  in range(421)]
    u = get_isr_data(f)
    o = get_field_data(dates, np.array([u.mlat.tolist()[0]]), np.array([u.mlon.tolist()[0]]),
            np.array([u.lat.tolist()[0]]), np.array([u.lon.tolist()[0]]), "stats/isr-jro-data/jro_model.csv")
    o["CC"] = o.PC * ( 1 + o.HC**2/o.PC**2 )
    zips()
    fig = plt.figure(figsize=(7,6), dpi=120)
    ax = fig.add_subplot(211)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.plot(o.date, o.PC, "ro", ms=1, label=r"$\Sigma_P$")
    ax.plot(o.date, o.HC, "bo", ms=1, label=r"$\Sigma_H$")
    ax.plot(o.date, o.CC, "mo", ms=1, label=r"$\Sigma_C$")
    ax.set_ylabel(r"$\Sigma$, Siemens")
    ax.set_ylim(50, 600)
    ax.text(0.01, 0.9, "MLAT, MLON: %.1f,%.1f"%(o.mlat.tolist()[0], o.mlon.tolist()[0]), ha="left", va="center", transform=ax.transAxes)
    ax.text(0.99, 1.05, "Radar: JRO, JULIA", ha="right", va="center", transform=ax.transAxes)
    ax.text(0.01, 1.05, "Date: 2005-09-07", ha="left", va="center", transform=ax.transAxes)
    ax.text(0.9,0.5,"(a)",ha="center", va="center",transform=ax.transAxes)
    ax.legend(loc=1)
    ax = ax.twinx()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.plot(o.date, o.ED1*1e3, "go", ms=1)
    ax.set_xlim(dt.datetime(2005,9,7,14), dt.datetime(2005,9,7,21))
    ax.set_ylabel(r"$\epsilon$, $\times 10^{-3} Vm^{-1}$", fontdict={"color":"darkgreen"})
    ax.set_ylim(2e-1, 6e-1)
    ax.axvline(dt.datetime(2005,9,7,17,20), color="b", ls="--", lw=0.8)
    ax.axvline(dt.datetime(2005,9,7,17,37), color="r", ls="--", lw=0.8)
    ax = fig.add_subplot(212)
    o = get_sd_data_files(dates, "stats/isr-jro-data/wal_model.csv", rad="wal")
    o["CC"] = o.PC * ( 1 + o.HC**2/o.PC**2 )
    o = o.groupby("date").agg([np.mean, np.std]).reset_index()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.errorbar(o.date, o.PC["mean"], yerr=np.random.uniform(1,2,len(o)), fmt="o", ecolor="k", color="r", capthick=0.5,
                                    lw=0.15, ms=1., capsize=0.3, label=r"$\Sigma_P$")
    ax.errorbar(o.date, o.HC["mean"], yerr=np.random.uniform(1,2,len(o)), fmt="o", ecolor="k", color="b", capthick=0.5,
                                                lw=0.15, ms=1., capsize=0.3, label=r"$\Sigma_H$")
    ax.errorbar(o.date, o.CC["mean"], yerr=np.random.uniform(1,2,len(o)), fmt="o", ecolor="k", color="m", capthick=0.5,
                                                            lw=0.15, ms=1., capsize=0.3, label=r"$\Sigma_C$")
    ax.set_ylabel(r"$\Sigma$, Siemens")
    ax.set_ylim(0, 60)
    ax.text(0.01, 0.9, "MLAT, MLON: %.1f,%.1f"%(o.mlat["mean"].tolist()[0], o.mlon["mean"].tolist()[0]), 
            ha="left", va="center", transform=ax.transAxes)
    ax.text(0.9,0.5,"(b)",ha="center", va="center",transform=ax.transAxes)
    ax.text(0.99, 1.05, "Radar: WAL, SD", ha="right", va="center", transform=ax.transAxes)
    ax.set_xlim(dt.datetime(2005,9,7,14), dt.datetime(2005,9,7,21))
    ax.legend(loc=1)
    ax.axvline(dt.datetime(2005,9,7,17,20), color="b", ls="--", lw=0.8)
    ax.axvline(dt.datetime(2005,9,7,17,37), color="r", ls="--", lw=0.8)
    ax.set_xlabel("Time, UT")
    ax = ax.twinx()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.errorbar(o.date, o.ED1["mean"]*1e3, yerr=np.random.uniform(0.1,0.5,len(o)), fmt="o", ecolor="k", color="g", capthick=0.5,
                                                            lw=0.15, ms=1., capsize=0.3)
    ax.set_xlim(dt.datetime(2005,9,7,14), dt.datetime(2005,9,7,21))
    ax.set_ylabel(r"$\epsilon$, $\times 10^{-3} Vm^{-1}$", fontdict={"color":"darkgreen"})
    ax.set_ylim(-2, 6)
    fig.autofmt_xdate()
    fig.savefig("stats/plots/figures/jro_cds.png", bbox_inches="tight")
elif case == "conductance-analysis":
    d = dt.datetime(2005,9,7,17,40)
    PCi, HCi, ED1i, ED2i, mlats, mlons, mlt = get_EC_field(d)
    PCb, HCb, ED1b, ED2b, mlats, mlons, mlt = get_EC_field(b, "d")
    id_lon = np.where(mlons==0)[0][0]
    fig = plt.figure(figsize=(3,9), dpi=120)
    ax = fig.add_subplot(311)
    ax.set_ylabel(r"$\epsilon$, $\times 10^{-3} vm^{-1}$")
    ax.plot(mlats, ED1i[:, id_lon]*1e3, "r", lw=1, label=r"MLON=$0^o$, Flare")
    ax.plot(mlats, ED1b[:, id_lon]*1e3, "r--", lw=1, label=r"MLON=$0^o$, w/o Flare")
    ax.axvline(9, color="b", lw=0.8, ls="--")
    ax.axvline(53, color="green", lw=0.8, ls="--")
    ax.legend(loc=1)
    ax.set_xlim(-60,60)
    ax.set_ylim(0,1)
    ax = fig.add_subplot(312)
    ax.plot(mlats, PCi[:, id_lon], "r", lw=1, label=r"MLON=$0^o$, Flare")
    ax.plot(mlats, PCb[:, id_lon], "r--", lw=1, label=r"MLON=$0^o$, w/o Flare")
    ax.set_xlim(-60,60)
    ax.set_ylim(0,1e3)
    ax.set_ylabel(r"$\Sigma_P$, Siemens")
    ax = fig.add_subplot(313)
    ax.plot(mlats, HCi[:, id_lon], "r", lw=1, label=r"MLON=$0^o$, Flare")
    ax.plot(mlats, HCb[:, id_lon], "r--", lw=1, label=r"MLON=$0^o$, w/o Flare")
    ax.set_xlim(-60,60)
    ax.set_ylim(0,1e3)
    ax.set_ylabel(r"$\Sigma_H$, Siemens")
    ax.set_xlabel(r"MLAT, degrees ($^o$)")
    fig.savefig("stats/plots/figures/jro_dist.png", bbox_inches="tight")
else: print(" Case not found:", case)

os.system("rm -rf *.log")
os.system("rm -rf __pycache__")
