import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.append("../sd")
import os
import glob
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy import array
import datetime as dt

from pysolar.solar import get_altitude
from nrlmsise00.dataset import msise_4d
from iri2016 import IRI

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

def calculate_eta(o, dn, f=3e6):
    ne, alts, lats, lons = o["nex"], o["alts"], o["latx"], o["lonx"]
    ds = msise_4d(dn, alts, lats, lons)
    nn = 1e6*(ds["He"] + ds["O"] + ds["N2"] + ds["O2"] + ds["Ar"] + ds["N"] + ds["H"] + ds["AnomO"])[0,:,:,:]
    T = ds["Talt"][0,:,:,:]
    p = nn * T *C["boltz"]
    nu = (2.637e6/np.sqrt(T) + 4.945e5) * p * 2.5
    w = 2*np.pi*f
    X,Z = (ne *1e6 * C["q_e"]**2) / (C["eps0"] * C["m_e"] * w**2), nu / w
    x,jz = X,Z*1.j
    n = np.sqrt(1 - (x / (1-jz)))
    real = n.real
    return real

def calculate_fp(ne):
    fp = (1/(2*np.pi))*np.sqrt(ne*1e6*C["q_e"]**2/(C["m_e"]*C["eps0"]))
    return fp

def get_R(fp, h, lats, lons, f=3e6):
    R = np.zeros((latN, lonN)) * np.nan
    for _i in range(latN):
        for _j in range(lonN):
            d = dn.replace(tzinfo=dt.timezone.utc)
            sza = get_altitude(lats[_i], lons[_j], d)
            R[_i, _j] = np.argmin(np.abs(fp[_i, _j]-f))
    return R

def intg(fDx, hx, lats, lons, R):
    fD = np.zeros((latN, lonN)) * np.nan
    for _i in range(latN):
        for _j in range(lonN):
            r = R[_i, _j]
            d = dn.replace(tzinfo=dt.timezone.utc)
            sza = 90.-get_altitude(lats[_i], lons[_j], d)
            if abs(sza) < 90.:
                fx = fDx[:int(r)+1, _i, _j]
                dh = hx[:int(r)+1]
                fD[_i, _j] = np.trapz(fx, dh)
            else: fD[_i, _j] = 0. 
    fD = np.ma.masked_array(fD, mask=np.isnan(fD))
    return fD

def to_interpolate(zz, SZA):
    zmax = np.max(zz)
    cosine = np.cos(np.deg2rad(SZA))**2
    cosine[np.abs(SZA) > 90] = 0.
    return cosine * zmax



dn = dt.datetime(2015,5,5,22,11)
files = glob.glob("../data/op/2015.05.05.22.11/waccmx/*.nc.gz")
files.sort()
fo = 12e6
dlat, dlon = 2, 4
iri = False

for _i, f in enumerate(files):
    if _i == 20:
        print("File: ", f)
        os.system("gzip -d " + f)
        f  = f.replace(".gz","")
        ds = Dataset(f)
        os.system("gzip " + f)
        bgc, flare = {}, {}
        if iri:
            pass
        else:
            bgc["nex"], bgc["alts"], bgc["latx"], bgc["lonx"] = _intp_(ds.variables["ZGd"][0,:,:,:]*1e-3, ds.variables["lat"][:], 
                ds.variables["lon"][:], ds.variables["NEd"][0,:,:,:], hd=[50,350,1], dlat=dlat, dlon=dlon, scale="log", kind="cubic", v=True)
            flare["nex"], flare["alts"], flare["latx"], flare["lonx"] = _intp_(ds.variables["ZGf"][0,:,:,:]*1e-3, ds.variables["lat"][:],
                ds.variables["lon"][:], ds.variables["NEf"][0,:,:,:], hd=[50,350,1], dlat=dlat, dlon=dlon, scale="log", kind="cubic", v=True)
        latN, lonN = len(flare["latx"]), len(flare["lonx"])
        bgc["fp"], flare["fp"] = calculate_fp(bgc["nex"]), calculate_fp(flare["nex"])
        R = get_R(flare["fp"], flare["alts"], flare["latx"], flare["lonx"], f=fo)
        bgc["eta"], flare["eta"] = calculate_eta(bgc, dn, f=fo), calculate_eta(flare, dn, f=fo)
        d_eta, d_ne = (flare["eta"] - bgc["eta"])/60., (flare["nex"] - bgc["nex"])*1e6/60.
        fDx = np.abs(2.*fo*d_eta/C["c"])
        fD = intg(fDx, flare["alts"], flare["latx"], flare["lonx"], R)
        vD = 2*C["c"]*fD
    else: continue
    fig, axes = plt.subplots(dpi=180, nrows=2, ncols=1, figsize=(10,12))
    ax = axes[0]
    _xx, _yy = np.meshgrid(flare["lonx"], flare["latx"])
    SZA = np.zeros_like(_xx)
    for _ix in range(len(flare["latx"])):
        for _jx in range(len(flare["lonx"])):
            d = dn.replace(tzinfo=dt.timezone.utc)
            SZA[_ix, _jx] = 90.-get_altitude(flare["latx"][_ix], flare["lonx"][_jx], d)
    #_fD = interp2d(_xx, _yy, fD, kind="cubic")
    cosine = np.cos(np.deg2rad(SZA))**2
    cosine[np.abs(SZA) > 90] = 0.
    vDx = to_interpolate(vD*.01, SZA) + (cosine*vD*.002)
    c = ax.pcolormesh(_xx, _yy, vDx, cmap="jet")
    cb = fig.colorbar(c, ax=ax, shrink=0.5)
    cb.set_label(r"$v_D$, $ms^{-1}$")
    ax.set_ylabel("Latitude")
    cs = ax.contour(_xx, _yy, SZA)
    ax.clabel(cs, inline=True, fontsize=10)
    ax = axes[1]
    c = ax.pcolormesh(_xx, _yy, flare["nex"].sum(axis=0)*1e6, cmap="jet")
    cb = fig.colorbar(c, ax=ax, shrink=0.5)
    cb.set_label(r"$n_ef$, $m^{-3}$")
    ax.set_ylabel("Latitude")
    np.savetxt("op/sza/sza.txt", SZA)
    np.savetxt("op/sza/xx.txt", _xx)
    np.savetxt("op/sza/yy.txt", _yy)
    np.savetxt("op/sza/vD.txt", vDx)
    ne = to_interpolate(vD*.01, SZA) + (cosine*vD*.002) 
    dne = to_interpolate(vD*.01, SZA) + (cosine*vD*.002)
    d_eta = to_interpolate(vD*.01, SZA) + (cosine*vD*.002)
    np.savetxt("op/sza/Ne.txt", )
    #ax = axes[2]
    #c = ax.pcolormesh(_xx, _yy, bgc["nex"].sum(axis=0)*1e6, cmap="jet")
    #cb = fig.colorbar(c, ax=ax, shrink=0.5)
    #cb.set_label(r"$n_eb$, $m^{-3}$")
    #ax.set_ylabel("Latitude")
    #ax = axes[3]
    #c = ax.pcolormesh(_xx, _yy, d_ne.sum(axis=0), cmap="jet")
    #cb = fig.colorbar(c, ax=ax, shrink=0.5)
    #cb.set_label(r"$dn_e$, $m^{-3}$")
    #ax.set_ylabel("Latitude")
    #ax = axes[4]
    #c = ax.pcolormesh(_xx, _yy, d_eta.sum(axis=0), cmap="jet")
    #cb = fig.colorbar(c, ax=ax, shrink=0.5)
    #cb.set_label(r"$d\eta$")
    #ax.set_ylabel("Latitude")
    #ax = axes[2]
    #c = ax.pcolormesh(_xx, _yy, SZA, cmap="jet")
    #cb = fig.colorbar(c, ax=ax, shrink=0.5)
    #cb.set_label(r"$SZA$, $^o$")
    #ax.set_ylabel("Latitude")
    #ax.set_xlabel("Longitude")
    fig.savefig("op/sza/%s_%d.png"%(dn.strftime("%Y.%m.%d.%H.%M"), _i), bbox_inches="tight")
    plt.close()

