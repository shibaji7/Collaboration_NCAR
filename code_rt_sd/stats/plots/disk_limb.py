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

C = {
    "boltz" : 1.38066e-23,       # Boltzmann constant  in Jule K^-1
    "h"     : 6.626e-34  ,       # Planks constant  in ergs s
    "c"     : 2.9979e08  ,       # in m s^-1
    "avo"   : 6.023e23   ,       # avogadro's number
    "Re"    : 6371.e3    ,
    "amu"   : 1.6605e-27 ,
    "q_e"   : 1.602e-19  ,       # Electron charge in C
    "m_e"   : 9.109e-31  ,       # Electron mass in kg
    "g"     : 9.81       ,       # Gravitational acceleration on the surface of the Earth
    "eps0"  : 1e-9/(36*np.pi),
    "R"     : 8.31       ,       # J mol^-1 K^-1
    }

def plot_fism_enhancement():
    from scipy.io import readsav
    cent = readsav("../../data/synthetic/lqian_2010/FISM_60sec_2005250_v00_01_clv_cent.sav")
    limb = readsav("../../data/synthetic/lqian_2010/FISM_60sec_2005250_v00_01_clv_limb.sav")
    no = readsav("../../data/synthetic/lqian_2010/FISM_60sec_2005250_v00_01_no_flr.sav")
    secs = 17*60 + 37
    fig = plt.figure(dpi=180, figsize=(6,7))
    ax = fig.add_subplot(211)
    date = dt.datetime(int(str(cent["date"])[:4]), 1, 1) + dt.timedelta(int(str(cent["date"])[4:])-1)
    ax.text(0.01, 1.05, "Date:%s"%date.strftime("%Y-%m-%d"), ha="left", va="center", transform=ax.transAxes)
    ax.text(0.01, 0.9, "(a)", ha="left", va="center", transform=ax.transAxes)
    ax.semilogy(cent["fism_wv"], cent["fism_pred"][:, secs], "r", lw=1.2, label="Center Flare (peak at 17:37 UT)")
    ax.semilogy(limb["fism_wv"], limb["fism_pred"][:, secs], "b", lw=1.0, label="Limb Flare (peak at 17:37 UT)")
    ax.semilogy(no["fism_wv"], no["fism_pred"][:, secs], "k", lw=0.6, label="Pre-flare (17:20 UT)")
    ax.legend(loc=1)
    ax.set_ylim(1e-6, 1e-1)
    ax.set_xlim(0,190)
    ax.set_ylabel(r"Solar Irradiance, $Wm^2$")
    ax = fig.add_subplot(212)
    ax.text(0.01, 0.9, "(b)", ha="left", va="center", transform=ax.transAxes)
    ax.semilogy(cent["fism_wv"],100*(cent["fism_pred"][:, secs]-no["fism_pred"][:, secs])/no["fism_pred"][:, secs], "r", lw=1.2, label="Center Flare")
    ax.semilogy(limb["fism_wv"],100*(limb["fism_pred"][:, secs]-no["fism_pred"][:, secs])/no["fism_pred"][:, secs], "b", lw=1.0, label="Limb Flare")
    ax.legend(loc=1)
    ax.set_xlabel(r"Wavelength, $nm$")
    ax.set_ylabel("Percentage Increase")
    ax.set_ylim(1e0, 1e5)
    ax.set_xlim(0,190)
    fig.savefig("figures/limb_disk.png", bbox_inches="tight")
    return

def plot_2D_distribution(alt=200, zipd=True, date=dt.datetime(2003,1,1), Re=6.356766e6, dlat=2, dlon=4, loc=[10,-100]):
    import fov
    import matplotlib.colors as colors
    mdata = sza.from_pickle("data_DL.pkl")
    cent, limb = {}, {}
    cent0, limb0 = {}, {}
    time_index, alt = 7, 200.
    if len(mdata.keys()) == 0:
        ds = Dataset("../../data/synthetic/lqian_2010/syr05250_17_disk.nc")
        cent["times"] = [date + dt.timedelta(days=x-1) for x in ds.variables["time"][:]] 
        Z = ds.variables["Z"][:] * 1e-3
        h = 1e-2*(Re*Z/(Re-Z))
        cent["nex"], cent["alts"], cent["latx"], cent["lonx"] = sza._intp_(h[time_index, :, :, :], ds.variables["lat"][:], 
                ds.variables["lon"][:], ds.variables["NE"][time_index,:,:,:], hd=[50,350,1], dlat=dlat, dlon=dlon, 
                scale="log", kind="cubic", v=True)
        cent0["nex"], cent0["alts"], cent0["latx"], cent0["lonx"] = sza._intp_(h[time_index-1, :, :, :], ds.variables["lat"][:],
                ds.variables["lon"][:], ds.variables["NE"][time_index-1,:,:,:], hd=[50,350,1], dlat=dlat, dlon=dlon,
                scale="log", kind="cubic", v=True)
        
        ds = Dataset("../../data/synthetic/lqian_2010/syr05250_17_limb.nc")
        limb["times"] = [date + dt.timedelta(days=x-1) for x in ds.variables["time"][:]]
        Z = ds.variables["Z"][:] * 1e-3
        h = 1e-2*(Re*Z/(Re-Z))
        limb["nex"], limb["alts"], limb["latx"], limb["lonx"] = sza._intp_(h[time_index, :, :, :], ds.variables["lat"][:], 
                ds.variables["lon"][:], ds.variables["NE"][time_index,:,:,:], hd=[50,350,1], dlat=dlat, dlon=dlon, 
                scale="log", kind="cubic", v=True)
        limb0["nex"], limb0["alts"], limb0["latx"], limb0["lonx"] = sza._intp_(h[time_index-1, :, :, :], ds.variables["lat"][:],
                ds.variables["lon"][:], ds.variables["NE"][time_index-1,:,:,:], hd=[50,350,1], dlat=dlat, dlon=dlon,
                scale="log", kind="cubic", v=True)
        mdata["cent"] = cent
        mdata["cent0"] = cent0
        mdata["limb"] = limb
        mdata["limb0"] = limb0
        sza.to_pickle(mdata, "data_DL.pkl")
    else: cent, limb = mdata["cent"], mdata["limb"]
    if zipd:
        alt_index = np.argmin(np.abs(cent["alts"]-200.))
        time = cent["times"][time_index]
        _sza = sza.get_SZA_array(cent["latx"], cent["lonx"], time)
        fig = plt.figure(dpi=150, figsize=(15,15))
        ax, fig, _to, _from = fov.get_globe(fig, 121, date)
        Z = cent["nex"][alt_index,:,:]*1e6
        xx, yy = np.meshgrid(cent["lonx"],cent["latx"])
        c = ax.pcolormesh(xx, yy, Z, cmap="plasma", norm=colors.LogNorm(vmin=1e8, vmax=1e12))
        cb = fig.colorbar(c, ax=ax, shrink=0.1)
        cb.set_label(r"$N_e^{Disk}$, $m^{-3}$")
        cs = ax.contour(cent["lonx"],cent["latx"], _sza)
        ax.clabel(cs, inline=True, fontsize=8, fmt=r"$%d^{o}$")
        ax.scatter([loc[1]],[loc[0]],color="g",s=15,marker="D")
        ax.text(0.01, 1.05, "Geographic Coordinates", horizontalalignment="left",
                verticalalignment="center", transform=ax.transAxes)
        ax.text(0.99, 1.05, "Disk Flare", ha="right", va="center", transform=ax.transAxes)
        ax.set_xlim(np.min(cent["lonx"]), np.max(cent["lonx"]))
        ax.set_ylim(np.min(cent["latx"]), np.max(cent["latx"]))
        ax, fig, _to, _from = fov.get_globe(fig, 122, date)
        Z = limb["nex"][alt_index,:,:]*1e6
        xx, yy = np.meshgrid(limb["lonx"],limb["latx"])
        c = ax.pcolormesh(xx, yy, Z, cmap="plasma", norm=colors.LogNorm(vmin=1e8, vmax=1e12))
        cb = fig.colorbar(c, ax=ax, shrink=0.1)
        cb.set_label(r"$N_e^{Limb}$, $m^{-3}$")
        cs = ax.contour(cent["lonx"],cent["latx"], _sza)
        ax.clabel(cs, inline=True, fontsize=8, fmt=r"$%d^{o}$")
        ax.scatter([loc[1]],[loc[0]],color="g",s=15,marker="D")
        ax.text(0.99, 1.05, "Limb Flare", ha="right", va="center", transform=ax.transAxes)
        ax.text(0.01, 1.05, "Height=%d km"%alt, ha="left", va="center", transform=ax.transAxes)
        ax.text(0.01, -0.05, "Date:%s UT"%time.strftime("%Y-%m-%d %H:%M"), ha="left", va="center", transform=ax.transAxes)
        ax.set_xlim(np.min(limb["lonx"]), np.max(limb["lonx"]))
        ax.set_ylim(np.min(limb["latx"]), np.max(limb["latx"]))
        fig.savefig("figures/disk_limb_dist_1.png", bbox_inches="tight")
    return

def plot_1D(frq=6.5e6, loc=[10,-100], date=dt.datetime(2003,9,7,17,37)):
    import sza
    mdata = sza.from_pickle("data_DL.pkl")
    print(mdata.keys())
    cent, limb = mdata["cent"], mdata["limb"]
    cent0, limb0 = mdata["cent0"], mdata["limb0"]

    fig = plt.figure(dpi=150, figsize=(6,3))
    ax = fig.add_subplot(121)
    ax.text(0.01, 0.9, "(c)", ha="left", va="center", transform=ax.transAxes)
    i,j = sza._get_ij_(cent["latx"], cent["lonx"], loc[0], loc[1])
    ax.semilogx(cent["nex"][:,i,j]*1e6, cent["alts"], lw=0.8, color="r", label="Disk")
    ax.semilogx(limb["nex"][:,i,j]*1e6, limb["alts"], lw=0.8, color="b", label="Limb")
    c_indx = np.argmin(np.abs(sza.calculate_fp(cent["nex"][:,i,j])-frq))
    ax.axhline(cent["alts"][c_indx], color="r", lw=0.4)
    l_indx = np.argmin(np.abs(sza.calculate_fp(limb["nex"][:,i,j])-frq))
    ax.axhline(limb["alts"][l_indx], color="b", lw=0.4)
    ax.text(0.8,0.9, "$\chi=%.1f^o$"%sza.get_SZA(loc[0], loc[1], date), ha="center", va="center", transform=ax.transAxes)
    ax.set_ylim(90, 300)
    ax.set_xlim(1e9,1e12)
    ax.set_xlabel(r"$N_e$, $m^{-3}$")
    ax.set_ylabel(r"Height, km")
    ax = fig.add_subplot(122)
    ax.set_ylim(90, 300)
    ax.axhline(cent["alts"][c_indx], color="r", lw=0.4)
    ax.axhline(limb["alts"][l_indx], color="b", lw=0.4)
    ax.text(1.05,0.5,"$f_0=%.1f$ MHz"%(frq/1e6), ha="center", va="center", transform=ax.transAxes, rotation=90)
    n_diff = -1*(sza.calculate_eta(cent["nex"][:,i,j], cent["alts"], loc[0], loc[1], date, frq) -\
            sza.calculate_eta(cent0["nex"][:,i,j], cent0["alts"], loc[0], loc[1], date, frq))
    n_diff[c_indx:] = np.nan
    v_diff = 1e3*2*frq*n_diff/3e8
    ax.axhline(cent["alts"][np.nanargmax(v_diff)], color="r", lw=0.8, ls="--")
    ax.plot(v_diff, cent["alts"], lw=1, color="r")

    n_diff = -1*(sza.calculate_eta(limb["nex"][:,i,j], limb["alts"], loc[0], loc[1], date, frq) -\
            sza.calculate_eta(limb0["nex"][:,i,j], limb0["alts"], loc[0], loc[1], date, frq))
    n_diff[l_indx:] = np.nan
    v_diff = 1e3*2*frq*n_diff/3e8
    ax.axhline(cent["alts"][np.nanargmax(v_diff)], color="b", lw=0.8, ls="--")
    ax.text(0.01, 0.9, "(d)", ha="left", va="center", transform=ax.transAxes)
    ax.plot(v_diff, limb["alts"], lw=1, color="b")
    ax.set_xlim(0,20)
    ax.set_xlabel(r"Doppler Velocity, $V_x (\times 1000, ms^{-1})$")
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    fig.savefig("figures/disk_limb_dist_2.png", bbox_inches="tight")
    return


if __name__ == "__main__":
    plot_fism_enhancement()
    #plot_2D_distribution(alt=200, loc=[10,-100])
    plot_1D(loc=[10,-100])
    os.system("rm -rf __pycache__")
    pass
