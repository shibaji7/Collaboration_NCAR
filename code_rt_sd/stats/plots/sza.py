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


def smooth(x, window_len=51, window="hanning"):
    if x.ndim != 1: raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len: raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<3: return x
    if not window in ["flat", "hanning", "hamming", "bartlett", "blackman"]: raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    s = np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    if window == "flat": w = numpy.ones(window_len,"d")
    else: w = eval("np."+window+"(window_len)")
    y = np.convolve(w/w.sum(),s,mode="valid")
    d = window_len - 1
    y = y[int(d/2):-int(d/2)]
    return y

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


def calculate_fp(ne):
    fp = (1/(2*np.pi))*np.sqrt(ne*1e6*C["q_e"]**2/(C["m_e"]*C["eps0"]))
    return fp

def to_pickle(data, fname):
    import pickle
    with open(fname, "wb") as f:
        pickle.dump(data, f)
    return

def from_pickle(fname):
    import pickle
    data = {}
    if os.path.exists(fname):
        with open(fname, "rb") as f:
            data = pickle.load(f)
    return data


def extrap1d(x,y,kind="linear"):
    """ This method is used to extrapolate 1D paramteres """
    interpolator = interp1d(x,y,kind=kind)
    xs = interpolator.x
    ys = interpolator.y
    def pointwise(x):
        if x < xs[0]: return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]: return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else: return interpolator(x)
    def ufunclike(xs):
        return array(list(map(pointwise, array(xs))))
    return ufunclike


def _get_ij_(lats, lons, lat, lon):
    """ Get (lat, lon) index """
    _ij_ = (np.argmin(np.abs(lats-lat)), np.argmin(np.abs(lons-lon)))
    return _ij_

def _intp_heights_(h, hx, param, scale="linear", kind="cubic"):
    if scale == "linear": pnew = interp1d(h, param, kind=kind)(hx)
    if scale == "log": pnew = 10**interp1d(h, np.log10(param), kind=kind)(hx)
    return pnew

def _intp_latlon_(lats, lons, latx, lonx, param, scale="linear", kind="cubic"):
    if scale == "linear": pnew = interp2d(lats, lons, param.T, kind=kind)(latx, lonx).T
    if scale == "log": pnew = 10**interp2d(lats, lons, np.log10(param.T), kind=kind)(latx, lonx).T
    return pnew

def _intp_(h, lats, lons, param, hd=[50,300,1], dlat=0.5, dlon=1, scale="log", kind="cubic", v=False):
    lonx = np.arange(np.min(lons), np.max(lons), dlon)
    latx = np.arange(np.min(lats), np.max(lats), dlat)
    param_intm = np.zeros((len(h), len(latx), len(lonx)))
    h_intm = np.zeros((len(h), len(latx), len(lonx)))
    for k,_ in enumerate(h):
        param_intm[k, :, :] = _intp_latlon_(lats, lons, latx, lonx, param[k, :, :], scale, kind)
        h_intm[k, :, :] = _intp_latlon_(lats, lons, latx, lonx, h[k, :, :], scale, kind)
    if v: print("\tLatlon convt.")
    hx = np.arange(hd[0],hd[1],hd[2])
    pnew = np.zeros((len(hx), len(latx), len(lonx)))
    for i,_ in enumerate(latx):
        for j,_ in enumerate(lonx):
            pnew[:,i,j] = 10**(extrap1d(h_intm[:, i, j], np.log10(param_intm[:, i, j]))(hx))
    if v: print("\tHeight convt.")
    return pnew, hx, latx, lonx

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

def calculate_eta(ne, alts, lat, lon, dn, f):
    from nrlmsise00 import msise_model
    k = (2 * np.pi * f) / C["c"]
    nn, T = [], []
    for alt in alts:
        dens, t  = msise_model(dn, alt, lat, lon, 150, 150, 4)
        del dens[5]
        nn.append(np.sum(dens)*1e6)
        T.append(t[1])
    p = np.array(nn) * np.array(T) *C["boltz"]
    nu = (2.637e6/np.sqrt(T) + 4.945e5) * p * 2.5
    w = 2*np.pi*f
    X,Z = (ne *1e6 * C["q_e"]**2) / (C["eps0"] * C["m_e"] * w**2), nu / w
    x,jz = X,Z*1.j
    n = np.sqrt(1 - (x / (1-jz)))
    real = n.real
    return real

def create_sza_plots(bgc, flare, fbgc, pbgc, pflare, date, frq, sza_small=[-150, 25], sza_high=[-70,25], zipd=True):
    import fov
    import matplotlib.colors as colors
    if zipd:
        sza = get_SZA_array(bgc["latx"], bgc["lonx"], date)
        fig = plt.figure(dpi=150, figsize=(7,3.5))
        ax, fig, _to, _from = fov.get_globe(fig, 111, date)
        Z_bgc = pflare["nex"][60,:,:]*1e6
        Z_flare = flare["nex"][60,:,:]*1e6
        xx, yy = np.meshgrid(bgc["lonx"],bgc["latx"])
        c = ax.pcolormesh(xx, yy, Z_flare-Z_bgc, cmap="plasma", norm=colors.SymLogNorm(linthresh=1e1, linscale=1e1, vmin=1e8, vmax=1e11))
        cb = fig.colorbar(c, ax=ax, shrink=0.7)
        cb.set_label(r"$\Delta N_e = N_e^{22:11}-N_e^{22:10}$, $m^{-3}$", fontsize=8)
        ax.scatter([sza_small[0]],[sza_small[1]],color="b",s=20,marker="D")
        ax.scatter([sza_high[0]],[sza_high[1]],color="r",s=20,marker="D")
        cs = ax.contour(bgc["lonx"],bgc["latx"], sza)
        ax.clabel(cs, inline=True, fontsize=8, fmt=r"$%d^{o}$")
        ax.text(-0.02, 0.99, "Geographic Coordinates", horizontalalignment="center",
                verticalalignment="top", transform=ax.transAxes, rotation=90, fontdict={"size":8})
        ax.text(0.01,1.05, "(a) 2015-05-05 22:11", ha="left", va="center", transform=ax.transAxes, fontdict={"size":8})
        ax.text(0.99,1.06, "Flare Class = X2.7\nHeight = %d km"%flare["alts"][60], ha="right", va="center", transform=ax.transAxes
                , fontdict={"size":8})
        fig.subplots_adjust(hspace=0.1, wspace=0.1)
        fig.savefig("figures/sza_dist_1.png", bbox_inches="tight")
    
    fig = plt.figure(dpi=150, figsize=(6,3))
    ax = fig.add_subplot(121)
    i, j = _get_ij_(flare["latx"], flare["lonx"], sza_small[1], sza_small[0])
    z = (flare["nex"][:,i,j] - bgc["nex"][:,i,j])
    z = np.sign(z) * np.log10(np.abs(z))
    ax.plot(z, flare["alts"], lw=0.8, color="b", label="$\chi=%.1f^o$"%get_SZA(sza_small[1], sza_small[0], date))
    #ax.semilogx(pflare["nex"][:,i,j], flare["alts"], lw=0.8, color="b", ls="--")
    fps = calculate_fp(flare["nex"][:,i,j])
    f_indx = np.argmin(np.abs(fps-frq))
    ax.axhline(flare["alts"][f_indx], color="b", lw=0.8, ls="--")
    ax.set_ylim(60, 300)
    #ax.set_xlim(1e8,1e12)
    ax.set_ylabel(r"Height, km")
    ax = fig.add_subplot(122)
    n_diff = -1*(calculate_eta(flare["nex"][:,i,j], flare["alts"], sza_small[1], sza_small[0], date, frq) -\
            calculate_eta(fbgc["nex"][:,i,j], fbgc["alts"], sza_small[1], sza_small[0], date, frq))
    n_diff[f_indx:] = np.nan
    ax.axhline(flare["alts"][f_indx], color="b", lw=0.8, ls="--")
    ax.plot(n_diff, flare["alts"], lw=1, color="b")
    ax.set_ylim(100, 300)
    ax.set_xlim(-0.2,.2)
    ax.text(0.01,1.05,"2015-05-05 22:11", ha="left", va="center", transform=ax.transAxes)
    ax.text(1.05,0.5,"$f_0=%.1f$ MHz"%(frq/1e6), ha="center", va="center", transform=ax.transAxes, rotation=90)
    
    ax = fig.add_subplot(121)
    i, j = _get_ij_(flare["latx"], flare["lonx"], sza_high[0], sza_high[0])
    z = (flare["nex"][:,i,j] - bgc["nex"][:,i,j])
    z = np.sign(z) * np.log10(np.abs(z))
    ax.plot(smooth(z), flare["alts"], lw=0.8, color="r", label="$\chi=%.1f^o$"%get_SZA(sza_high[1], sza_high[0], date))
    #ax.semilogx(pflare["nex"][:,i,j], flare["alts"], lw=0.8, color="r", ls="--")
    #ax.plot(z, flare["alts"], lw=0.8, color="r", label="$\chi=%.1f^o$"%get_SZA(sza_high[1], sza_high[0], date))
    #ax.semilogx(flare["nex"][:,i,j]*1e6, flare["alts"], lw=0.8, color="b")
    #ax.semilogx(bgc["nex"][:,i,j]*1e6, bgc["alts"], lw=0.6, color="r")
    fps = calculate_fp(flare["nex"][:,i,j])
    f_indx = np.argmin(np.abs(fps-frq))
    ax.axhline(flare["alts"][f_indx], color="r", lw=0.8, ls="--")
    ax.set_ylim(60, 300)
    ax.legend(loc=2, fontsize=8)
    ax.set_xlabel(r"$\log_{10}(\Delta N_e$), $m^{-3}$")
    ax.set_ylabel(r"Height, km")
    ax = fig.add_subplot(122)
    n_diff = smooth(-1*(calculate_eta(flare["nex"][:,i,j], flare["alts"], sza_high[1], sza_high[0], date, frq) -\
            calculate_eta(bgc["nex"][:,i,j], bgc["alts"], sza_high[1], sza_high[0], date, frq)))
    n_diff[f_indx:] = np.nan
    ax.axhline(flare["alts"][f_indx], color="r", lw=0.8, ls="--")
    ax.plot(n_diff*1e3, flare["alts"], lw=1, color="r")
    ax.set_ylim(60, 300)
    ax.set_xlim(-0.2,.2)
    #ax.text(0.8,0.8, "$\chi=%.1f^o$"%get_SZA(sza_high[1], sza_high[0], date), ha="center", va="center", transform=ax.transAxes)
    ax.set_xlabel(r"Change in refractive index, $\delta\eta$")
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    fig.savefig("figures/sza_dist_2.png", bbox_inches="tight")
    return

dn = dt.datetime(2015,5,5,22,11)
files = glob.glob("../../data/op/2015.05.05.22.11/waccmx/*.nc.gz")
files.sort()
fo = 4e6
dlat, dlon = 2, 4
iri = False

mdata = from_pickle("data.pkl")
bgc, fbgc, flare, pbgc, pflare = {}, {}, {}, {}, {}
if len(mdata.keys()) == 0:
    for _i, f in enumerate(files):
        if _i == 0:
            print("File: ", f)
            os.system("gzip -d " + f)
            f  = f.replace(".gz","")
            ds = Dataset(f)
            os.system("gzip " + f)
            bgc = {}
            bgc["nex"], bgc["alts"], bgc["latx"], bgc["lonx"] = _intp_(ds.variables["ZGd"][0,:,:,:]*1e-3,
                    ds.variables["lat"][:], ds.variables["lon"][:], ds.variables["NEd"][0,:,:,:],
                    hd=[50,350,1], dlat=dlat, dlon=dlon, scale="log", kind="cubic", v=True)
            mdata["bgc"] = bgc
        if _i == 18:
            print("File: ", f)
            os.system("gzip -d " + f)
            f  = f.replace(".gz","")
            ds = Dataset(f)
            os.system("gzip " + f)
            fbgc["nex"], fbgc["alts"], fbgc["latx"], fbgc["lonx"] = _intp_(ds.variables["ZGd"][0,:,:,:]*1e-3, 
                    ds.variables["lat"][:], ds.variables["lon"][:], ds.variables["NEd"][0,:,:,:], 
                    hd=[50,350,1], dlat=dlat, dlon=dlon, scale="log", kind="cubic", v=True)
            flare["nex"], flare["alts"], flare["latx"], flare["lonx"] = _intp_(ds.variables["ZGf"][0,:,:,:]*1e-3,
                    ds.variables["lat"][:], ds.variables["lon"][:], ds.variables["NEf"][0,:,:,:], 
                    hd=[50,350,1], dlat=dlat, dlon=dlon, scale="log", kind="cubic", v=True)
            mdata["fbgc"] = fbgc
            mdata["flare"] = flare
        if _i == 17:
            print("File: ", f)
            os.system("gzip -d " + f)
            f  = f.replace(".gz","")
            ds = Dataset(f)
            os.system("gzip " + f)
            pbgc["nex"], pbgc["alts"], pbgc["latx"], pbgc["lonx"] = _intp_(ds.variables["ZGd"][0,:,:,:]*1e-3,
                    ds.variables["lat"][:], ds.variables["lon"][:], ds.variables["NEd"][0,:,:,:], 
                    hd=[50,350,1], dlat=dlat, dlon=dlon, scale="log", kind="cubic", v=True)
            pflare["nex"], pflare["alts"], pflare["latx"], pflare["lonx"] = _intp_(ds.variables["ZGf"][0,:,:,:]*1e-3,
                    ds.variables["lat"][:], ds.variables["lon"][:], ds.variables["NEf"][0,:,:,:],
                    hd=[50,350,1], dlat=dlat, dlon=dlon, scale="log", kind="cubic", v=True)
            mdata["pbgc"] = pbgc
            mdata["pflare"] = pflare
    to_pickle(mdata, "data.pkl")
else: bgc, flare, fbgc, pbgc, pflare = mdata["bgc"], mdata["flare"], mdata["fbgc"], mdata["pbgc"], mdata["pflare"]
create_sza_plots(bgc, flare, fbgc, pbgc, pflare, dt.datetime(2015,5,5,22,11), fo)

os.system("rm -rf *.log")
