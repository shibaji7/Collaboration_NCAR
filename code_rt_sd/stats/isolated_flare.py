import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.append("../sd/")
import datetime as dt
import pandas as pd
from get_sd_data import FetchData
import numpy as np

goes = pd.read_csv("isr-jro-data/g12_xrs_1m_20050901_20050930.csv", skiprows=117, parse_dates=["time_tag"])
goes = goes[(goes.time_tag >= dt.datetime(2005,9,7)) & (goes.time_tag < dt.datetime(2005,9,8))]

VIPN2, DVIPN2, VIPE1, DVIPE1, time = [], [], [], [], []
with open("isr-jro-data/jul20050907_avg_150km.001.txt", "r") as f:
    lines = f.readlines()
    header = lines[0].split(" ")
    header = list(filter(lambda x: x != "" and x!="\n", header))
    print(header)
    for l in lines[1:]:
        l = list(filter(lambda x: x != "" and x!="\n", l.split(" ")))
        time.append(dt.datetime(int(l[0]),int(l[1]),int(l[2]),int(l[3]),int(l[4]),int(l[5])))
        if l[11] != "missing": VIPN2.append(float(l[11]))
        else: VIPN2.append(np.nan)
        if l[12] != "missing": DVIPN2.append(float(l[12]))
        else: DVIPN2.append(np.nan)
        if l[13] != "missing": VIPE1.append(float(l[13]))
        else: VIPE1.append(np.nan)
        if l[14] != "missing": DVIPE1.append(float(l[14]))
        else: DVIPE1.append(np.nan)

fd = FetchData("wal", [dt.datetime(2005,9,7), dt.datetime(2005,9,8)])
beams, _ = fd.fetch_data(v_params=["elv", "v", "w_l", "gflg", "p_l", "slist", "v_e"])
rec = fd.convert_to_pandas(beams)

fig, axes = plt.subplots(nrows=3,ncols=1,dpi=150,figsize=(5,9))
ax = axes[0]
ax.semilogy(goes.time_tag, goes.xl, "ro", lw=0.8, label="SXR", ms=0.8)
ax.semilogy(goes.time_tag, goes["xs"], "bo", lw=0.8, label="HXR", ms=0.8)
ax.set_ylabel(r"Flux, $Wm^{-2}$")
ax.legend(loc=1)
ax.set_xlim(dt.datetime(2005,9,7,17), dt.datetime(2005,9,7,18))

ax = axes[1]
ax.set_ylabel(r"Velocity, $ms^{-1}$")
ax.plot(rec.time, rec.v, "ro", ms=0.5, alpha=0.7)
ax.set_xlim(dt.datetime(2005,9,7,17), dt.datetime(2005,9,7,18))
ax.set_ylim(-50, 100)

ax = axes[2]
ax.set_ylabel(r"Velocity, $ms^{-1}$")
ax.set_xlabel(r"Time, UT")
ax.plot(time, VIPE1, "ro", lw=0.8, label="VIPE1", ms=0.8)
ax.legend(loc=2)
ax.set_ylim(-100,-50)
ax = ax.twinx()
ax.plot(time, VIPN2, "bo", lw=0.8, label="VIPN1", ms=0.8)
ax.legend(loc=1)
ax.set_ylim(-10,60)
ax.set_xlim(dt.datetime(2005,9,7,17), dt.datetime(2005,9,7,18))

fig.autofmt_xdate()
fig.subplots_adjust(hspace=0.3, wspace=0.3)
fig.savefig("op/2005_event.png", bbox_inches="tight")
