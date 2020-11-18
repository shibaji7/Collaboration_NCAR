#!/usr/bin/env python

"""model_analysis.py: Analysis of e-Density and spectral irradiance test experiments python program for Doppler shift"""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import matplotlib
matplotlib.use("Agg")
import os
import glob
import numpy as np

import sys
sys.path.append("sd/")
import plotlib

case = 2

if case == 0:
    import pandas as pd
    import numpy as np
    from scipy import io
    from netCDF4 import Dataset
    import plotlib
    fname = "data/op/{dn}/waccmx/{rad}/bm.{bm}/ne.ti({ti}).f.mat"
    events = pd.read_csv("config/events.csv", parse_dates=["dn"])
    rads = pd.read_csv("config/radars.csv")
    eDensPC, eDensAC = {16: [], 17: [], 18: [], 19: [], 20: []}, {16: [], 17: [], 18: [], 19: [], 20: []}
    for d in events.dn.tolist():
        for r in rads.rad.tolist():
            for bm in range(24):
                for i in range(15,21):
                    f = fname.format(dn=d.strftime("%Y.%m.%d.%H.%M"), rad=r, bm="%02d"%bm, ti="%02d"%i)
                    if os.path.exists(f):
                        o = io.loadmat(f)["ne"]
                        ne = o[:,0].ravel()
                        print(ne.shape)
                        if i >= 16: 
                            eDensPC[i].append(100*((ne-b)/b))
                            eDensAC[i].append((ne-b))
                        b = ne
                        print(f)
                #break
            #break
        #break
    plotlib.plot_edens_versus_height(eDensPC, eDensAC, ylim=[50,350])


    irr = np.array([])
    for d in events.dn.tolist():
        year = d.year
        fname = "../FISM/%d_fism2_flare_waccmx_23_1min_fixed_c200828.nc"%year
        if os.path.exists(fname):
            time = np.array([d.strftime("%Y%m%d")]).astype(np.int64)
            print(time)
            nc = Dataset(fname)
            date = nc.variables["date"][:]
            ssi = nc.variables["ssi"][:]
            index = np.where(date==time)
            secs = nc.variables["datesec"][index]
            ssi = (ssi[index,:]).reshape((ssi[index,:].shape[1], ssi[index,:].shape[2]))
            wavelength = nc.variables["wavelength"][:]
            if len(irr) == 0: irr = ssi
            else: irr = np.append(ssi, irr, axis=0)
            print(irr.shape)
    plotlib.plot_ssi_versus_bins(irr, wavelength, ylim=[50,350])
if case == 1:
    ## Source activate aacgm_carto
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import sys
    sys.path.append("sd_cartopy/")
    from fov import *

    import datetime as dt
    from netCDF4 import Dataset
    from utils import InterpolateData
    fglob = "data/op/2015.05.05.22.11/waccmx/*.gz"
    files = glob.glob(fglob)
    files.sort()

    dates = [dt.datetime(2015,5,5,22), dt.datetime(2015,5,5,22,8), dt.datetime(2015,5,5,22,9), dt.datetime(2015,5,5,22,10), dt.datetime(2015,5,5,22,11)]
    dastrs = [x.strftime("%Y.%m.%d.%H.%M") for x in dates]
    fig = plt.figure(figsize=(12, 12), dpi=150)
    dsb = None
    for _i, dx in enumerate(dastrs):
        if _i>0:
            idx = 220+_i
            ax = fig.add_subplot(idx, projection="fovcarto",coords="geo", plot_date=dates[_i])
            ax.coastlines()
            ax.overlay_radar()
            ax.overlay_fov()
            ax.grid_on()
            ax.enum()
        for f in files:
            if dx in f:
                os.system("cp {f} .".format(f=f))
                f = f.split("/")[-1]
                os.system("gzip -d {f}".format(f=f))
                ds = Dataset(f.replace(".gz", ""))
                if dsb is None:
                    intp = InterpolateData()
                    dsb = ds
                    fx = plt.figure(figsize=(5,5), dpi=150)
                    axf = fx.add_subplot(111, projection="fovcarto",coords="geo", plot_date=dates[_i])
                    axf.coastlines()
                    axf.overlay_radar()
                    axf.overlay_fov()
                    axf.grid_on()
                    axf.enum()
                    px = intp._intrp_(dsb.variables["ZGf"][0,:,:,:]*1e-3, dsb.variables["lat"][:], dsb.variables["lon"][:],
                                                        dsb.variables["WIf"][0,:,:,:], hd=300)
                    axf.overlay_data(dsb.variables["lat"][:], dsb.variables["lon"][:], px,
                            tx=cartopy.crs.PlateCarree(), colorbar_label="$\omega_I, ms^{-1}$")
                    fx.savefig("data/WI.png", bbox_inches="tight")
                    plt.close()
                if _i>0:
                    intp = InterpolateData()
                    px = intp._intrp_(ds.variables["ZGf"][0,:,:,:]*1e-3, ds.variables["lat"][:], ds.variables["lon"][:],
                            ds.variables["WIf"][0,:,:,:]-dsb.variables["WIf"][0,:,:,:], hd=300)
                    if _i==4:
                        ax.overlay_data(ds.variables["lat"][:], ds.variables["lon"][:], px,
                        tx=cartopy.crs.PlateCarree(), p_max=4, p_min=-4, p_step=0.5, p_ub=8, p_lb=-8,
                        cmap=matplotlib.pyplot.get_cmap("Spectral"), add_colorbar=True, colorbar_label="$\Delta \omega_I, ms^{-1}$")
                    else:
                        ax.overlay_data(ds.variables["lat"][:], ds.variables["lon"][:], px,
                               tx=cartopy.crs.PlateCarree(), p_max=4, p_min=-4, p_step=0.5, p_ub=8, p_lb=-8,
                               cmap=matplotlib.pyplot.get_cmap("Spectral"), add_colorbar=False, colorbar_label="Velocity [m/s]")
                os.system("rm {f}".format(f=f.replace(".gz", "")))
                break
    fig.savefig("data/wi.png", bbox_inches="tight")
if case==2:
    plotlib.plot_ray_edens()
    plotlib.plot_ray_edens(time=12, diff=True)
    plotlib.plot_ray_edens(time=9, diff=False)
