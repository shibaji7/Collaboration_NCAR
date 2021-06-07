import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os
import requests
import pandas as pd
import numpy as np
import datetime as dt
import plotlib


def get_f107(d):
    d = int(d.strftime("%Y%m%d"))
    f107 = pd.read_csv("data.csv")
    return f107[f107.time==d]["f107"].tolist()[0]

def proc_bins(df, bins, keyx, keyy, rnd=True):
    x, y, err = [], [], []
    for i in range(len(bins)-1):
        u = df[(df[keyx]>=bins[i]) & (df[keyx]<bins[i+1])]
        if len(u) > 0:
            x.append((bins[i]+bins[i+1])/2)
            y.append(u[keyy].median().round(1))
            err.append(np.abs(u[keyy]-u[keyy].median()).median().round(1))
    x, y, err = np.array(x), np.array(y), np.array(err)
    return x, y, err

case = 3
if case == 0:
    for i in range(28):
        print("python capture_stats.py -n %d"%i)
        os.system("python capture_stats.py -n %d"%i)
elif case == 1:
    ix = 0
    events = pd.read_csv("op/radar_event_list.csv", parse_dates=["date"])
    sim_fname = "finescale_simulate_total_{kind}.csv"
    df = pd.DataFrame()
    for d, r in zip(events.date.tolist(), events.rad.tolist()):
        if ix > 0:
            kind = r + "_" + d.strftime("%Y-%m-%d-%H-%M")
            u = pd.read_csv("op/"+sim_fname.format(kind=kind))
            u["f107"] = get_f107(d)
            df = pd.concat([df, u])
        ix+=1
    df = df[(np.abs(df.vT)>30) & (np.abs(df.vT)<300)]
    df["vN"] = df.vT - df.vFh
    szas = np.arange(35,95,5)
    f107s = [90, 125, 160, 200]
    fig, axes = plt.subplots(dpi=150, nrows=2, ncols=3, figsize=(9,6))
    ax = axes[0,0]
    ax.set_xlabel("SZA")
    ax.set_ylabel(r"Doppler Velocity, $v_D$")
    x, y, err = proc_bins(df, szas, "sza", "vD")
    ax.errorbar(x, y, yerr=err, color="r", marker="o", ms=1, ls="None", lw=1, label="vD")
    ax = axes[0,1]
    ax.set_xlabel("SZA")
    ax.set_ylabel(r"Doppler Velocity, $v_E$")
    x, y, err = proc_bins(df, szas, "sza", "vE")
    ax.errorbar(x, y,  yerr=err, color="g", marker="o", ms=1, ls="None", lw=0.8, label="vE")
    ax = axes[0,2]
    ax.set_xlabel("SZA")
    ax.set_ylabel(r"Doppler Velocity, $v_F$")
    x, y, err = proc_bins(df, szas, "sza", "vF")
    ax.errorbar(x, y,  yerr=err, color="b", marker="o", ms=1, ls="None", lw=0.6, label="vF")
    ax = axes[1,0]
    ax.set_xlabel("F107")
    ax.set_ylabel(r"Doppler Velocity, $v_D$")
    x, y, err = proc_bins(df, f107s, "f107", "vD")
    ax.errorbar(x, y, yerr=err, color="r", marker="o", ms=1, ls="None", lw=1, label="vD")
    ax = axes[1,1]
    ax.set_xlabel("F107")
    ax.set_ylabel(r"Doppler Velocity, $v_E$")
    x, y, err = proc_bins(df, f107s, "f107", "vE")
    ax.errorbar(x, y,  yerr=err, color="g", marker="o", ms=1, ls="None", lw=0.8, label="vE")
    ax = axes[1,2]
    ax.set_xlabel("F107")
    ax.set_ylabel(r"Doppler Velocity, $v_F$")
    x, y, err = proc_bins(df, f107s, "f107", "vF")
    ax.errorbar(x, y,  yerr=err, color="b", marker="o", ms=1, ls="None", lw=0.6, label="vF")
    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("op/stats_variations.png", bbox_inches="tight")
elif case == 2:
    dates = [dt.datetime(2014,6,10,11,42), dt.datetime(2015,3,11,16,22)]
    events = pd.read_csv("op/radar_event_list.csv", parse_dates=["date"])
    sim_fname = "finescale_simulate_total_{kind}.csv"
    for d, r in zip(events.date.tolist(), events.rad.tolist()):
        if d in dates:
            kind = r + "_" + d.strftime("%Y-%m-%d-%H-%M")
            x = pd.read_csv("op/"+sim_fname.format(kind=kind))
            x = x[(np.abs(x.vT)>30) & (np.abs(x.vT)<300)]
            print(kind, x.sza.mean())
            vd, ve, vf = np.array(x.vD/x.vT), np.array(x.vE/x.vT), np.array((x.vF + x.vFh)/x.vT)
            vdn, vdh = np.array((x.vD + x.vE + x.vF)/x.vT), np.array(x.vFh/x.vT)
            plotlib.plot_htstogram(vd, ve, vf, vdn, vdh, kind)
    pass
elif case == 3:
    wavebins = [1, 10, 30]
    events = pd.read_csv("op/radar_event_list.csv", parse_dates=["date"])
    sim_fname = "finescale_simulate_total_{kind}.csv"
    df = pd.DataFrame()
    for d, r in zip(events.date.tolist(), events.rad.tolist()):
        kind = r + "_" + d.strftime("%Y-%m-%d-%H-%M")
        x = pd.read_csv("op/"+sim_fname.format(kind=kind))
        x = x[(np.abs(x.vT)>30) & (np.abs(x.vT)<300)]
        
        fname = "op/tmpfism.csv"
        start, end = d - dt.timedelta(minutes=1), d + dt.timedelta(minutes=1)
        for b in wavebins:
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
            x["bin_"+str(b)] = data["irradiance (W/m^2/nm)"].tolist()[1]
            os.remove(fname)
        df = pd.concat([df, x])
        #break
    #print(df["bin_1"].max(), df["bin_1"].min())
    #print(df["bin_10"].max(), df["bin_10"].min())
    #print(df["bin_30"].max(), df["bin_30"].min())
    nu1 = [5e-5,1e-4,2e-4,4e-4,7e-4]
    nu10 = [2e-5,3e-5,4e-5,5e-5,6e-5]
    nu30 = [7e-6,8e-6,9e-6,11e-6]
    fig, axes = plt.subplots(dpi=150, nrows=3, ncols=3, figsize=(9,9))
    ax = axes[0,0]
    ax.set_xlabel(r"$\nu=1, 10^5$")
    ax.set_ylabel(r"Doppler Velocity, $v_D$")
    x, y, err = proc_bins(df, nu1, "bin_1", "vD")
    ax.errorbar(x*1e5, y, yerr=err, color="r", marker="o", ms=1, ls="None", lw=1, label="vD")
    ax = axes[0,1]
    ax.set_xlabel(r"$\nu=10, 10^5$")
    ax.set_ylabel(r"Doppler Velocity, $v_D$")
    x, y, err = proc_bins(df, nu10, "bin_10", "vD")
    ax.errorbar(x*1e5, y, yerr=err, color="r", marker="o", ms=1, ls="None", lw=1, label="vD")
    ax = axes[0,2]
    ax.set_xlabel(r"$\nu=30, 10^5$")
    ax.set_ylabel(r"Doppler Velocity, $v_D$")
    x, y, err = proc_bins(df, nu30, "bin_30", "vD")
    ax.errorbar(x*1e6, y, yerr=err, color="r", marker="o", ms=1, ls="None", lw=1, label="vD")

    ax = axes[1,0]
    ax.set_xlabel(r"$\nu=1, 10^5$")
    ax.set_ylabel(r"Doppler Velocity, $v_E$")
    x, y, err = proc_bins(df, nu1, "bin_1", "vE")
    ax.errorbar(x*1e5, y, yerr=err, color="g", marker="o", ms=1, ls="None", lw=1, label="vE")
    ax = axes[1,1]
    ax.set_xlabel(r"$\nu=10, 10^5$")
    ax.set_ylabel(r"Doppler Velocity, $v_E$")
    x, y, err = proc_bins(df, nu10, "bin_10", "vE")
    ax.errorbar(x*1e5, y, yerr=err, color="g", marker="o", ms=1, ls="None", lw=1, label="vE")
    ax = axes[1,2]
    ax.set_xlabel(r"$\nu=30, 10^5$")
    ax.set_ylabel(r"Doppler Velocity, $v_E$")
    x, y, err = proc_bins(df, nu30, "bin_30", "vE")
    ax.errorbar(x*1e6, y, yerr=err, color="g", marker="o", ms=1, ls="None", lw=1, label="vE")

    ax = axes[2,0]
    ax.set_xlabel(r"$\nu=1, 10^5$")
    ax.set_ylabel(r"Doppler Velocity, $v_F$")
    x, y, err = proc_bins(df, nu1, "bin_1", "vF")
    ax.errorbar(x*1e5, y, yerr=err, color="b", marker="o", ms=1, ls="None", lw=1, label="vF")
    ax = axes[2,1]
    ax.set_xlabel(r"$\nu=10, 10^5$")
    ax.set_ylabel(r"Doppler Velocity, $v_F$")
    x, y, err = proc_bins(df, nu10, "bin_10", "vF")
    ax.errorbar(x*1e5, y, yerr=err, color="b", marker="o", ms=1, ls="None", lw=1, label="vF")
    ax = axes[2,2]
    ax.set_xlabel(r"$\nu=30, 10^5$")
    ax.set_ylabel(r"Doppler Velocity, $v_F$")
    x, y, err = proc_bins(df, nu30, "bin_30", "vF")
    ax.errorbar(x*1e6, y, yerr=err, color="b", marker="o", ms=1, ls="None", lw=1, label="vF")
    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("op/stats_binpeak_variations.png", bbox_inches="tight")
    pass

if os.path.exists("__pycache__/"):
    os.system("rm -rf __pycache__/")
    os.system("rm -rf py*.log")
