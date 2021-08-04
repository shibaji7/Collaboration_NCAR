import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import netCDF4 as nc
import datetime as dt

import os
from pysolar.solar import get_altitude
from scipy.io import loadmat
import glob
from scipy.integrate import trapz
from scipy import signal

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

def _estimate_dop_delh_(x, y, freq, phi=0):
    dh = (np.max(x.height) - np.max(y.height)) * 1000.
    xf = (-2.*freq*1e6/3e8) * (dh/60.) * np.cos(np.deg2rad(phi))
    xd = 0.5 * xf * 3e8 / (freq * 1e6)
    return xd

def get_vdeta(d, rlim, freq):
    d = d[(d.height>=rlim[0]) & (d.height<rlim[1])]
    f = trapz(signal.resample(d.dop,INT_F))
    v = (0.5 * f * 3e8 / (freq * 1e6))
    return v

def create_distributions():
    x = pd.read_csv("/home/shibaji/Collaboration_NCAR/code_rt_sd/config/finescale_simulate.csv")
    x = x[(np.abs(x.vT)>30) & (np.abs(x.vT)<300)]
    vd, ve, vf = np.array(x.vD/x.vT), np.array(x.vE/x.vT), np.array((x.vF + x.vFh)/x.vT)
    vd = vd[vd>.1]
    vn, vh = np.array((x.vD + x.vE + x.vF)/x.vT), np.array(x.vFh/x.vT)

    fig = plt.figure(figsize=(8,8), dpi=120)
    ax = fig.add_subplot(221)
    ax.hist(vn, bins=np.arange(0,1,.01), color="r", alpha=0.5, density=True, 
            label=r"$\frac{V_{d\eta}}{V_T}[\mu=%.2f]$"%np.mean(vn), histtype="step")
    ax.axvline(np.mean(vn), ls="--", color="r", lw=0.6)
    ax.axvline(np.mean(vh), ls="--", color="b", lw=0.6)
    ax.hist(vh, bins=np.arange(0,1,.01), color="b", alpha=0.5, density=True, 
            label=r"$\frac{V_{dh}}{V_T}[\mu=%.2f]$"%np.mean(vh), histtype="step")
    ax.text(0.1,0.9, "(a.1)", horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
    ax.set_xlim(0,1)
    ax.set_ylim(0,20)
    ax.legend(loc=1, prop={"size": 8})
    ax.set_ylabel(r"Density $\left(\frac{V_x}{V_T}\right)$")

    ax = fig.add_subplot(222)
    ax.axvline(np.mean(vd), ls="--", color="r", lw=0.6)
    ax.hist(vd, bins=np.arange(0,1,.01), color="r", alpha=0.5, density=True, 
            label=r"$\frac{V_D}{V_T}[\mu=%.2f]$"%np.mean(vd), histtype="step")
    ax.axvline(np.mean(ve), ls="--", color="g", lw=0.6)
    ax.hist(ve, bins=np.arange(0,1,.01), color="g", alpha=0.5, density=True, 
            label=r"$\frac{V_E}{V_T}[\mu=%.2f]$"%np.mean(ve), histtype="step")
    ax.axvline(np.mean(vf), ls="--", color="b", lw=0.6)
    ax.hist(vf, bins=np.arange(0,1,.01), color="b", alpha=0.5, density=True, 
            label=r"$\frac{V_F}{V_T}[\mu=%.2f]$"%np.mean(vf), histtype="step")
    ax.set_ylim(0,50)
    ax.text(0.1,0.9, "(b.1)", horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
    ax.set_xlim(0,1)
    ax.legend(loc=1, prop={"size": 8})
    ax.set_xlabel(r"$\frac{V_x}{V_T}$")

    x = pd.read_csv("../op/finescale_simulate_total_A.csv")
    x = x[(np.abs(x.vT)>30) & (np.abs(x.vT)<300)]
    vd, ve, vf = np.array(x.vD/x.vT), np.array(x.vE/x.vT), np.array((x.vF + x.vFh)/x.vT)
    vn, vh = np.array((x.vD + x.vE + x.vF)/x.vT), np.array(x.vFh/x.vT)
    ax = fig.add_subplot(223)
    ax.hist(vn, bins=np.arange(0,1,.01), color="r", alpha=0.5, density=True,
            label=r"$\frac{V_{d\eta}}{V_T}[\mu=%.2f]$"%np.mean(vn), histtype="step")
    ax.axvline(np.mean(vn), ls="--", color="r", lw=0.6)
    ax.axvline(np.mean(vh), ls="--", color="b", lw=0.6)
    ax.hist(vh, bins=np.arange(0,1,.01), color="b", alpha=0.5, density=True,
            label=r"$\frac{V_{dh}}{V_T}[\mu=%.2f]$"%np.mean(vh), histtype="step")
    ax.text(0.1,0.9, "(a.2)", horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
    ax.set_xlim(0,1)
    ax.set_ylim(0,20)
    ax.legend(loc=1, prop={"size": 8})
    ax.set_ylabel(r"Density $\left(\frac{V_x}{V_T}\right)$")
    ax.set_xlabel(r"$\frac{V_x}{V_T}$")

    ax = fig.add_subplot(224)
    ax.axvline(np.mean(vd), ls="--", color="r", lw=0.6)
    ax.hist(vd, bins=np.arange(0,1,.01), color="r", alpha=0.5, density=True,
            label=r"$\frac{V_D}{V_T}[\mu=%.2f]$"%np.mean(vd), histtype="step")
    ax.axvline(np.mean(ve), ls="--", color="g", lw=0.6)
    ax.hist(ve, bins=np.arange(0,1,.01), color="g", alpha=0.5, density=True,
            label=r"$\frac{V_E}{V_T}[\mu=%.2f]$"%np.mean(ve), histtype="step")
    ax.axvline(np.mean(vf), ls="--", color="b", lw=0.6)
    ax.hist(vf, bins=np.arange(0,1,.01), color="b", alpha=0.5, density=True,
            label=r"$\frac{V_F}{V_T}[\mu=%.2f]$"%np.mean(vf), histtype="step")
    ax.set_ylim(0,50)
    ax.text(0.1,0.9, "(b.2)", horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
    ax.set_xlim(0,1)
    ax.legend(loc=1, prop={"size": 8})
    ax.set_xlabel(r"$\frac{V_x}{V_T}$")

    fig.savefig("figures/histogram.png", bbox_inches="tight")
    return

def limb_disk_flares():
    def to_sza(lat, lon, dn):
        d = dn.replace(tzinfo=dt.timezone.utc)
        sza = 90.-get_altitude(lats[_i], lons[_j], d)
        return

    def to_alt(H, Ro = 6356.766):
        H = H / 1.e5
        z = Ro * H/(Ro - H)
        return z

    pwd = "/home/shibaji/Collaboration_NCAR/code_rt_sd/data/synthetic/lqian_2010/"
    ds = nc.Dataset(pwd + "2005_fism2_daily_23.nc")
    ds = nc.Dataset(pwd + "2005_fism2_daily_and_SeptemberFlares_23.nc")

    limb = nc.Dataset(pwd + "syr05250_17_limb.nc") 
    print(limb.variables) 
    return

def sza_distribution():
    sza = np.loadtxt("../op/sza/sza.txt")
    _xx = np.loadtxt("../op/sza/xx.txt")
    _yy = np.loadtxt("../op/sza/yy.txt")
    vDx = np.loadtxt("../op/sza/vD.txt")

    fig = plt.figure(dpi=100, figsize=(8,5))
    import sys
    sys.path.append("/home/shibaji/Collaboration_NCAR/code_rt_sd/sd_cartopy/")
    import fov
    ax, fig, _to, _from = fov.get_globe(fig, 111, dt.datetime(2015,5,5,22,11))
    ax.text(0.01, 1.05, "Geographic Coordinates", horizontalalignment="left",
            verticalalignment="center", transform=ax.transAxes)
    ax.text(0.99, 1.05, "22:11 UT 5 May, 2015", horizontalalignment="right",
                        verticalalignment="center", transform=ax.transAxes)
    ax.text(-0.03, 0.6, "Propagation: NVIS, 12MHz", horizontalalignment="center",
                        verticalalignment="center", transform=ax.transAxes, rotation=90)
    #ax.grid_on()
    x, y = np.zeros_like(_xx), np.zeros_like(_yy)
    for _i in range(_xx.shape[0]):
        for _j in range(_xx.shape[1]):
            x[_i, _j], y[_i, _j] = _to.transform_point(_xx[_i, _j], _yy[_i, _j], _from)
    c = ax.pcolormesh(_xx, _yy, vDx, cmap="jet")
    cb = fig.colorbar(c, ax=ax, shrink=0.5)
    cs = ax.contour(_xx, _yy, sza)
    ax.clabel(cs, inline=True, fontsize=10, fmt=r"$%d^{o}$")
    cb.set_label(r"$v_D$, $ms^{-1}$")
    fig.savefig("figures/sza_dist.png", bbox_inches="tight")
    return

def intensity_distribution():
    def vals(x):
        expo = 1e-4 if x[0] == "X" else 1e-5
        mant = float(x[1:])
        return mant*expo
    d_reg, e_reg, f_reg = [60,90], [100,130], [140,300]
    dat = pd.read_csv("../op/radar_event_list.csv", parse_dates=["date"])
    intenisty = dat["type"].apply(vals)
    velocity, It  = [], []
    for d, r, _i_ in zip(dat.date.tolist(), dat.rad.tolist(), intenisty):
        freq = get_freq(d, r)
        vD, vE, vF, vFh, vT, SZA = [], [], [], [], [], []
        for bm in range(24):
            bearing_file = "../../data/op/%s/waccmx/%s/bm.%02d/bearing.mat"%(d.strftime("%Y.%m.%d.%H.%M"), r, bm)
            obj = loadmat(bearing_file)
            lat, lon = obj["lat"], obj["lon"]
            sza = np.mean(calculate_sza(d, lat, lon))
            dic = "../../data/op/{dn}/waccmx/{r}/bm.{bm}/".format(r=r,bm="%02d"%bm,dn=d.strftime("%Y.%m.%d.%H.%M"))
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
                        except: pass
                        pass
                    pass
                pass
            pass
        vT = np.array(vT)
        vT = vT[vT<=500]
        velocity.extend(vT)
        It.extend(len(vT)*[np.round(np.log10(_i_),3)])
        pass
    df = pd.DataFrame()
    df["intensity"], df["velocity"] = It, velocity
    fig = plt.figure(dpi=100, figsize=(5,5))
    ax = fig.add_subplot(111)
    boxplot = df.boxplot(by="intensity", ax=ax)
    ax.set_ylabel(r"Modeled Doppler Velocity, $ms^{-1}$")
    ax.set_xlabel(r"Peak Soft X-ray, $Wm^{-2}$")
    ax.set_title("")
    fig.savefig("figures/intensity.png", bbox_inches="tight")
    return


if __name__ == "__main__":
    #create_distributions()
    #limb_disk_flares()
    #sza_distribution()
    #intensity_distribution()
    os.system("rm -fv *.log")
    pass
