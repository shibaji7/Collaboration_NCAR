#!/usr/bin/env python

"""plotlib.py: module is dedicated to plottting."""

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
import datetime as dt
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import SubplotDivider, Size
from mpl_toolkits.axes_grid1.mpl_axes import Axes
import matplotlib.pyplot as plt
from pylab import gca, gcf
import numpy as np
from matplotlib.transforms import Affine2D, Transform
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.projections import polar
from mpl_toolkits.axisartist.grid_finder import FixedLocator, DictFormatter
from types import MethodType
import glob
import pandas as pd
from dateutil import tz
from scipy.io import loadmat
import copy
from scipy.stats import skewnorm
from scipy.integrate import trapz
from scipy import signal

import sys
sys.path.append("sd_cartopy/")
import rad_fov
from fov import *

#from PyIF import te_compute as te
#from sklearn.feature_selection import mutual_info_regression as MIR
#from SALib.sample import saltelli
#from SALib.analyze import sobol
#from SALib.analyze import rbd_fast

import itertools
from math import pi
from matplotlib.legend_handler import HandlerPatch

font = {'size'   : 8}
matplotlib.rc('font', **font)
class HandlerCircle(HandlerPatch):
    def create_artists(self, legend, orig_handle,
            xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = plt.Circle(xy=center, radius=orig_handle.radius)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]

INT_F = 300
INT_R = 300
import utils

def textHighlighted(xy, text, ax=None, color="k", fontsize=None, xytext=(0,0),
        zorder=None, text_alignment=(0,0), xycoords="data", 
        textcoords="offset points", **kwargs):
    """
    Plot highlighted annotation (with a white lining)
    
    Parameters
    ----------
    xy : position of point to annotate
    
    text : str text to show
    ax : Optional[ ]
    color : Optional[char]
    text color; deafult is "k"
    fontsize : Optional [ ] text font size; default is None
    xytext : Optional[ ] text position; default is (0, 0)
    zorder : text zorder; default is None
    text_alignment : Optional[ ]
    xycoords : Optional[ ] xy coordinate[1]; default is "data"
    textcoords : Optional[ ] text coordinate[2]; default is "offset points"
    **kwargs :
    
    Notes
    -----
    Belongs to class rbspFp.
    
    References
    ----------
    [1] see `matplotlib.pyplot.annotate
    <http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.annotate>`)
    [2] see `matplotlib.pyplot.annotate
    <http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.annotate>`)
    """
    
    if ax is None: ax = gca()
    text_path = mp.text.TextPath((0, 0), text, size=fontsize, **kwargs)
    p1 = matplotlib.patches.PathPatch(text_path, ec="w", lw=4, fc="w", alpha=0.7,
            zorder=zorder, transform=mp.transforms.IdentityTransform())
    p2 = matplotlib.patches.PathPatch(text_path, ec="none", fc=color, zorder=zorder, 
            transform=mp.transforms.IdentityTransform())
    offsetbox2 = matplotlib.offsetbox.AuxTransformBox(mp.transforms.IdentityTransform())
    offsetbox2.add_artist(p1)
    offsetbox2.add_artist(p2)
    ab = mp.offsetbox.AnnotationBbox(offsetbox2, xy, xybox=xytext, xycoords=xycoords, 
            boxcoords=textcoords, box_alignment=text_alignment, frameon=False)
    ab.set_zorder(zorder)
    ax.add_artist(ab)
    return

def addColorbar(mappable, ax):
    """
    Append colorbar to axes
    
    Parameters
    ----------
    mappable : a mappable object
    ax : an axes object
    
    Returns
    -------
    cbax : colorbar axes object
    
    Notes
    -----
    This is mostly useful for axes created with :func:`curvedEarthAxes`.
    """
    fig1 = ax.get_figure()
    divider = SubplotDivider(fig1, *ax.get_geometry(), aspect=True)
    # axes for colorbar
    cbax = Axes(fig1, divider.get_position())
    h = [Size.AxesX(ax), # main axes
            Size.Fixed(0.05), # padding
            Size.Fixed(0.1)] # colorbar
    v = [Size.AxesY(ax)]
    _ = divider.set_horizontal(h)
    _ = divider.set_vertical(v)
    _ = ax.set_axes_locator(divider.new_locator(nx=0, ny=0))
    _ = cbax.set_axes_locator(divider.new_locator(nx=2, ny=0))
    _ = fig1.add_axes(cbax)
    _ = cbax.axis["left"].toggle(all=False)
    _ = cbax.axis["top"].toggle(all=False)
    _ = cbax.axis["bottom"].toggle(all=False)
    _ = cbax.axis["right"].toggle(ticklabels=True, label=True)
    _ = plt.colorbar(mappable, cax=cbax, shrink=0.1)
    return cbax


def curvedEarthAxes(rect=111, fig=None, minground=0., maxground=2000, minalt=0,
                            maxalt=500, Re=6371., nyticks=5, nxticks=4):
    """
    Create curved axes in ground-range and altitude
    
    Parameters
    ----------
    rect : Optional[int] subplot spcification
    fig : Optional[pylab.figure object] (default to gcf)
    minground : Optional[float]
    maxground : Optional[int] maximum ground range [km]
    minalt : Optional[int] lowest altitude limit [km]
    maxalt : Optional[int] highest altitude limit [km]
    Re : Optional[float] Earth radius in kilometers
    nyticks : Optional[int] Number of y axis tick marks; default is 5
    nxticks : Optional[int] Number of x axis tick marks; deafult is 4
    
    Returns
    -------
    ax : matplotlib.axes object containing formatting
    aax : matplotlib.axes objec containing data
    """

    ang = maxground / Re
    minang = minground / Re
    angran = ang - minang
    angle_ticks = [(0, "{:.0f}".format(minground))]
    while angle_ticks[-1][0] < angran:
        tang = angle_ticks[-1][0] + 1./nxticks*angran
        angle_ticks.append((tang, "{:.0f}".format((tang-minang)*Re)))
    grid_locator1 = FixedLocator([v for v, s in angle_ticks])
    tick_formatter1 = DictFormatter(dict(angle_ticks))

    altran = float(maxalt - minalt)
    alt_ticks = [(minalt+Re, "{:.0f}".format(minalt))]
    while alt_ticks[-1][0] < Re+maxalt:
        alt_ticks.append((altran / float(nyticks) + alt_ticks[-1][0], 
            "{:.0f}".format(altran / float(nyticks) +
                alt_ticks[-1][0] - Re)))
    _ = alt_ticks.pop()
    grid_locator2 = FixedLocator([v for v, s in alt_ticks])
    tick_formatter2 = DictFormatter(dict(alt_ticks))
    tr_rotate = Affine2D().rotate(np.pi/2-ang/2)
    tr_shift = Affine2D().translate(0, Re)
    tr = polar.PolarTransform() + tr_rotate

    grid_helper = floating_axes.GridHelperCurveLinear(tr, extremes=(0, angran, Re+minalt, Re+maxalt),
            grid_locator1=grid_locator1, grid_locator2=grid_locator2, tick_formatter1=tick_formatter1,
            tick_formatter2=tick_formatter2,)
    if not fig: fig = plt.figure(figsize=(5,3), dpi=240)
    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    # adjust axis
    print("adjust ax")
    ax1.set_ylabel(r"Height, $km$", fontdict={"size":2})
    ax1.set_xlabel(r"Ground Range, $km$", fontdict={"size":2})
    ax1.invert_xaxis()
    ax1.minground = minground
    ax1.maxground = maxground
    ax1.minalt = minalt
    ax1.maxalt = maxalt
    ax1.Re = Re
    fig.add_subplot(ax1, transform=tr)
    # create a parasite axes whose transData in RA, cz
    aux_ax = ax1.get_aux_axes(tr)
    # for aux_ax to have a clip path as in ax
    aux_ax.patch = ax1.patch
    # but this has a side effect that the patch is drawn twice, and possibly
    # over some other artists. So, we decrease the zorder a bit to prevent this.
    ax1.patch.zorder=0.9
    return ax1, aux_ax

def plot_edens(time, beam=None, maxground=2000, maxalt=500,
        nel_cmap="jet", nel_lim=[10, 12], title=False, 
        fig=None, rect=111, ax=None, aax=None,plot_colorbar=True,
        nel_rasterize=False):
    """
    Plot electron density profile
    
    Parameters
    ----------
    time : datetime.datetime time of profile
    beam : Optional[ ] beam number
    maxground : Optional[int]
    maximum ground range [km]
    maxalt : Optional[int] highest altitude limit [km]
    nel_cmap : Optional[str] color map name for electron density index coloring
    nel_lim : Optional[list, int] electron density index plotting limits
    title : Optional[bool] Show default title
    fig : Optional[pylab.figure] object (default to gcf)
    rect : Optional[int] subplot spcification
    ax : Optional[ ] Existing main axes
    aax : Optional[ ] Existing auxialary axes
    plot_colorbar : Optional[bool] Plot a colorbar
    nel_rasterize : Optional[bool] Rasterize the electron density plot

    Returns
    -------
    ax : matplotlib.axes object containing formatting
    aax : matplotlib.axes object containing data
    cbax : matplotlib.axes object containing colorbar
    """
    return

def get_polar(d, Re=6371.):
    """ Convert to polar coordinates """
    th = d.grange / Re
    r = d.height + Re
    dop, sth, dth = d.dop, d.sth, d.dth
    return th, r, dop, sth, dth

def plot_rays(dic, time, ti, beam, case, txt, maxground=2000, maxalt=500, step=1,
        showrefract=False, nr_cmap="jet_r", nr_lim=[0.0, .1], 
        raycolor="0.3", title=True, zorder=2, alpha=1, 
        fig=None, rect=111, ax=None, aax=None):
    """
    Plot ray paths

    Parameters
    ----------
    dic: str location of the data files
    time: datetime.datetime time of rays
    ti: int time index
    beam: beam number
    maxground : Optional[int] maximum ground range [km]
    maxalt : Optional[int] highest altitude limit [km]
    step : Optional[int] step between each plotted ray (in number of ray steps)
    showrefract : Optional[bool] show refractive index along ray paths (supersedes raycolor)
    nr_cmap : Optional[str] color map name for refractive index coloring
    nr_lim : Optional[list, float] refractive index plotting limits
    raycolor : Optional[float] color of ray paths
    title : Optional[bool] Show default title
    zorder : Optional[int]
    alpha : Optional[int]
    fig : Optional[pylab.figure] object (default to gcf)
    rect : Optional[int] subplot spcification
    ax : Optional[ ] Existing main axes
    aax : Optional[ ] Existing auxialary axes
    
    Returns
    -------
    ax : matplotlib.axes object containing formatting
    aax : matplotlib.axes object containing data
    cbax : matplotlib.axes object containing colorbar
    """
    if not ax and not aax: ax, aax = curvedEarthAxes(fig=fig, rect=rect, maxground=maxground, maxalt=maxalt)
    else:
        if hasattr(ax, "time"): time = ax.time
        if hasattr(ax, "beam"): beam = ax.beam
    files = glob.glob(dic + "ti({ti}).bm({bm}).elv(*).{case}.csv".format(ti=ti, bm=beam, case=case))
    files.sort()
    for f in files:
        th, r, v, _, _ = get_polar(pd.read_csv(f))
        if not showrefract: aax.plot(th, r, c=raycolor, zorder=zorder, alpha=alpha)
        else:
            points = np.array([th, r]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lcol = LineCollection( segments, zorder=zorder, alpha=alpha)
            _ = lcol.set_cmap( nr_cmap )
            _ = lcol.set_norm( plt.Normalize(*nr_lim) )
            _ = lcol.set_array( v )
            _ = aax.add_collection( lcol )
    if title:
        stitle = "%s UT"%time.strftime("%Y-%m-%d %H:%M")
        ax.set_title( stitle )
        ax.text(1.05, 0.5, txt, horizontalalignment="center", verticalalignment="center", 
                transform=ax.transAxes, rotation=90)
    if showrefract:
        cbax = addColorbar(lcol, ax)
        _ = cbax.set_ylabel(r"$\Delta$ f")
    else: cbax = None
    ax.beam = beam
    fig = ax.get_figure()
    fig.savefig(dic + "rt.ti({ti}).bm({bm}).{case}.png".format(ti=ti, bm=beam, case=case), bbox_inches="tight")
    plt.close()
    return ax, aax, cbax

def plot_exp_rays(dic, time, beam, cat="bgc", maxground=2000, maxalt=300, step=1,
        showrefract=False, nr_cmap="jet_r", nr_lim=[0.8, 1.],
        raycolor="0.3", title=False, zorder=2, alpha=1,
        fig=None, rect=111, ax=None, aax=None):
    """ Plot ray paths (previous method) """
    if not ax and not aax: ax, aax = curvedEarthAxes(fig=fig, rect=rect, maxground=maxground, maxalt=maxalt)
    else:
        if hasattr(ax, "time"): time = ax.time
        if hasattr(ax, "beam"): beam = ax.beam
    files = glob.glob(dic + "exp.{cat}.bm({bm}).elv(*).csv".format(cat=cat, bm=beam))
    files.sort()
    for f in files:
        th, r, v, _, _ = get_polar(pd.read_csv(f))
        if not showrefract: aax.plot(th, r, c=raycolor, zorder=zorder, alpha=alpha)
        else:
            points = np.array([th, r]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lcol = LineCollection( segments, zorder=zorder, alpha=alpha)
            _ = lcol.set_cmap( nr_cmap )
            _ = lcol.set_norm( plt.Normalize(*nr_lim) )
            _ = lcol.set_array( v )
            _ = aax.add_collection( lcol )
    if title:
        stitle = ""
        ax.set_title( stitle )
    if showrefract:
        cbax = addColorbar(lcol, ax)
        _ = cbax.set_ylabel(r"$\Delta$ f")
    else: cbax = None
    ax.beam = beam
    fig = ax.get_figure()
    fig.savefig(dic + "rt.exp.{cat}.bm({bm}).png".format(cat=cat, bm=beam), bbox_inches="tight")
    plt.close()
    return ax, aax, cbax

def plot_radstn(p,f,pz,fz,fname,lat,lon,t,zone="America/New_York"):
    """ Plot radar vertical dataset """
    fig = plt.figure(figsize=(4,4), dpi=120)
    ax = fig.add_subplot(111)
    ax.set_ylabel("Alt. [km]")
    ax.set_xlabel(r"EDens [$cm^{-3}$]")
    ax.semilogx(p, pz, "r")
    ax.semilogx(f, fz, "r--")
    ax.set_ylim(50, 130)
    ax.set_xlim(1e2, 1e7)
    sza = utils.calculate_sza(t, lat, lon, alt=300)
    l = t.replace(tzinfo=tz.gettz("UTC")).astimezone(tz.gettz("America/New_York"))
    ax.set_title(r"UT-%s"%(t.strftime("%Y-%m-%d %H:%M")))
    ax.text(1.05, 0.5, "Loc:(%.1f,%.1f), $\chi$-%.1f, LT-%s"%(lat, lon, sza, l.strftime("%H:%M")), 
        horizontalalignment="center", verticalalignment="center", transform=ax.transAxes, rotation=90)
    fig.savefig(fname,bbox_inches="tight")
    plt.close()
    return

def plot_velocity_ts(dn, rad, bmnum):
    """ Plot velocity TS data """
    fig = plt.figure(figsize=(6,6), dpi=150)
    axs = [fig.add_subplot(311), fig.add_subplot(312), fig.add_subplot(313)]
    mkeys = ["vn", "vh", "vt"]
    fmt = matplotlib.dates.DateFormatter("%H:%M")
    fname = "data/sim/{dn}/{rad}/velocity.ts.csv".format(dn=dn.strftime("%Y.%m.%d.%H.%M"), rad=rad)
    sdat = pd.read_csv(fname, parse_dates=["dn"])
    axs[0].set_title("%s UT, Radar - %s, Beam - %d"%(dn.strftime("%Y.%m.%d.%H.%M"), rad, bmnum))
    cols = ["r", "b", "k"]
    labs = [r"$V_{d\eta}$", r"$V_{dh}$", r"$V_{t}$"]
    I = 0
    fname = "data/sim/{dn}/{rad}/sd_data.csv.gz".format(dn=dn.strftime("%Y.%m.%d.%H.%M"), rad=rad)
    dat = utils.get_sd_data(fname, 15).dropna()
    dat = dat.groupby("time").mean().reset_index()
    for ax, mkey, col, lab in zip(axs, mkeys, cols, labs):
        ax.set_ylabel(r"Velocity, $ms^{-1}$")
        ax.set_xlabel("Time, UT")
        ax.xaxis.set_major_formatter(fmt)
        yerr = np.array([(mn, mx) for mn, mx in zip(sdat[mkey+"_min"], sdat[mkey+"_max"])]).T
        ax.errorbar(sdat.dn, sdat[mkey], yerr=yerr, 
            mec=col, mfc=col, fmt="r^", ms=1.5, ls="None", ecolor=col, 
            capsize=1, capthick=.4, elinewidth=0.4,
            alpha=0.5, label=lab)
        if I == 2: 
            ax.plot(dat.time, dat.v, color="darkgreen", marker="o", 
                alpha=0.3, ls="None", markersize=0.5, label=r"$V_{sd}^{los}$")
            ax.plot(dat.time, dat.v, color="darkred", marker="o",
                    alpha=0.3, ls="None", markersize=0.8)
        ax.axhline(0, color="gray", ls="--", lw=0.6)
        ax.legend(loc=1)
        ax.set_ylim(10*int((np.min(sdat[mkey]+sdat[mkey+"_min"])/10)-1), 
                10*int((np.max(sdat[mkey]+sdat[mkey+"_max"])/10)+1))
        ax.set_xlim(sdat.dn.tolist()[0], sdat.dn.tolist()[-1])
        I += 1
    fname = "data/sim/{dn}/{rad}/velocity.ts.png".format(dn=dn.strftime("%Y.%m.%d.%H.%M"), rad=rad)
    fig.savefig(fname,bbox_inches="tight")
    return

def plot_radstn_base(b,p,f,ht,fname,lat,lon,t,zone="America/New_York"):
    """ Plot radar vertical dataset """
    fig = plt.figure(figsize=(4,4), dpi=120)
    ax = fig.add_subplot(111)
    ax.set_ylabel("Alt. [km]")
    ax.set_xlabel(r"EDens [$cm^{-3}$]")
    ax.semilogx(b, ht, "k", label="Background")
    ax.semilogx(p, ht, "r", label=r"$UT_{-1}$")
    ax.semilogx(f, ht, "r--", label="UT")
    ax.legend(loc=4)
    ax.set_ylim(50, 130)
    ax.set_xlim(1e2, 1e7)
    sza = utils.calculate_sza(t, lat, lon, alt=300)
    l = t.replace(tzinfo=tz.gettz("UTC")).astimezone(tz.gettz("America/New_York"))
    ax.set_title(r"UT-%s"%(t.strftime("%Y-%m-%d %H:%M")))
    ax.text(1.05, 0.5, "Loc:(%.1f,%.1f), $\chi$-%.1f, LT-%s"%(lat, lon, sza, l.strftime("%H:%M")),
            horizontalalignment="center", verticalalignment="center", transform=ax.transAxes, rotation=90)
    fig.savefig(fname,bbox_inches="tight")
    plt.close()
    return

def plot_rays_base(dic, time, ti, beam, case, txt, maxground=2000, maxalt=500, step=1,
        showrefract=False, nr_cmap="Blues", nr_lim=[-0.5, 0.5], 
        raycolor="0.3", title=True, zorder=2, alpha=1, 
        fig=None, rect=111, ax=None, aax=None, freq=12.):
    """
    Plot ray paths

    Parameters
    ----------
    dic: str location of the data files
    time: datetime.datetime time of rays
    ti: int time index
    beam: beam number
    maxground : Optional[int] maximum ground range [km]
    maxalt : Optional[int] highest altitude limit [km]
    step : Optional[int] step between each plotted ray (in number of ray steps)
    showrefract : Optional[bool] show refractive index along ray paths (supersedes raycolor)
    nr_cmap : Optional[str] color map name for refractive index coloring
    nr_lim : Optional[list, float] refractive index plotting limits
    raycolor : Optional[float] color of ray paths
    title : Optional[bool] Show default title
    zorder : Optional[int]
    alpha : Optional[int]
    fig : Optional[pylab.figure] object (default to gcf)
    rect : Optional[int] subplot spcification
    ax : Optional[ ] Existing main axes
    aax : Optional[ ] Existing auxialary axes
    
    Returns
    -------
    ax : matplotlib.axes object containing formatting
    aax : matplotlib.axes object containing data
    cbax : matplotlib.axes object containing colorbar
    """
    if not ax and not aax: ax, aax = curvedEarthAxes(fig=fig, rect=rect, maxground=maxground, maxalt=maxalt)
    else:
        if hasattr(ax, "time"): time = ax.time
        if hasattr(ax, "beam"): beam = ax.beam
    files = glob.glob(dic + "ti({ti})_elv(*)_{case}.csv".format(ti="%02d"%ti, case=case))
    files.sort()
    Re = 6371. 
    for f in files:
        th, r, v, _, _ = get_polar(pd.read_csv(f))
        v = (0.5 * v * 3e8 / (freq * 1e6))
        if not showrefract: aax.plot(th, r, c=raycolor, zorder=zorder, alpha=alpha)
        else:
            points = np.array([th, r]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lcol = LineCollection( segments, zorder=zorder, alpha=alpha)
            _ = lcol.set_cmap( nr_cmap )
            _ = lcol.set_norm( plt.Normalize(*nr_lim) )
            _ = lcol.set_array( utils.smooth(v, window_len=21) )
            _ = aax.add_collection( lcol )
            aax.plot(np.arange(0,2000)/Re, np.ones(2000)*60+Re, color="b", ls="--", lw=0.5)
            aax.plot(np.arange(0,2000)/Re, np.ones(2000)*95+Re, color="orange", ls="--", lw=0.5)
            aax.plot(np.arange(0,2000)/Re, np.ones(2000)*130+Re, color="r", ls="--", lw=0.5)
    if not showrefract and title:
        stitle = "%s UT"%time.strftime("%Y-%m-%d %H:%M")
        ax.set_title( stitle )
        ax.text(1.05, 0.5, txt, horizontalalignment="center", verticalalignment="center", 
                transform=ax.transAxes, rotation=90)
    if showrefract:
        cbax = addColorbar(lcol, ax)
        _ = cbax.set_ylabel(r"$\Delta$ V (m/s)")
        stitle = "%s UT"%time.strftime("%Y-%m-%d %H:%M")+ "\n" + "Radar: BKS, Beam: %02d"%beam + "\n" +\
                "Frequency: %.1f MHz"%freq + "\n"
        ax.text(0.5, 0.8, stitle + txt + "(m/s)", horizontalalignment="center", verticalalignment="center",
                transform=ax.transAxes)
    else: cbax = None
    ax.beam = beam
    fig = ax.get_figure()
    #fig.savefig(dic + "rt.ti({ti}).{case}.png".format(ti="%02d"%ti, case=case), bbox_inches="tight")
    fig.savefig(dic + "rt.ti({ti}).{case}.png".format(ti="%02d"%ti, case=case))
    plt.close()
    return ax, aax, cbax

def plot_region_distribution(vd, ve, vf):
    from scipy import stats
    fig = plt.figure(figsize=(4,4), dpi=120)
    ax = fig.add_subplot(111)
    ax.hist(vd, bins=np.arange(0,1,.01), color="r", alpha=0.5, density=True, label=r"$\frac{v_D}{v_T}$", histtype="step")
    ax.hist(ve, bins=np.arange(0,1,.01), color="b", alpha=0.5, density=True, label=r"$\frac{v_E}{v_T}$", histtype="step")
    ax.hist(vf, bins=np.arange(0,1,.01), color="g", alpha=0.5, density=True, label=r"$\frac{v_F}{v_T}$", histtype="step")
    ax.set_xlim(0,1)
    ax.legend(loc=1)
    ax.set_ylabel(r"Density ($\frac{V_x}{V_T}$)")
    ax.set_xlabel(r"$\frac{V_x}{V_T}$")
    fig.savefig("data/hist_reg.png", bbox_inches="tight")
    return

def plot_distribution(vn, vf):
    from scipy import stats
    fig = plt3figure(figsize=(4,4), dpi=120)
    ax = fig.add_subplot(111)
    ax.hist(vn, bins=np.arange(0,1,.01), color="r", alpha=0.5, density=True, label=r"$\frac{V_{d\eta}}{V_T}$", histtype="step")
    ax.hist(vf, bins=np.arange(0,1,.01), color="b", alpha=0.5, density=True, label=r"$\frac{V_{dh}}{V_T}$", histtype="step")
    ax.set_xlim(0,1)
    ax.legend(loc=1)
    ax.set_ylabel(r"Density $(\frac{V_x}{V_T})$")
    ax.set_xlabel(r"$\frac{V_x}{V_T}$")
    fig.savefig("data/hist.png", bbox_inches="tight")
    return

def plot_htstogram(vd, ve, vf, vn, vh):
    from scipy.stats import beta
    fig = plt.figure(figsize=(6,3), dpi=150)
    ax = fig.add_subplot(121)
    #x = np.arange(0,1,0.001)
    #a, b, _, _ = beta.fit(vn,floc=0,fscale=1)
    ax.hist(vn, bins=np.arange(0,1,.01), color="r", alpha=0.5, density=True, label=r"$\frac{V_{d\eta}}{V_T}[\mu=%.2f]$"%np.mean(vn)
            , histtype="step")
    #ax.plot(x, beta.pdf(x, a, b), color="r", lw=0.8, label=r"$\frac{V_{d\eta}}{V_T}[\mu=%.2f]$"%(a/(a+b)))
    ax.axvline(np.mean(vn), ls="--", color="r", lw=0.6)
    #a, b, _, _ = beta.fit(vh,floc=0,fscale=1)
    #ax.plot(x, beta.pdf(x, a, b), color="b", lw=0.8, label=r"$\frac{V_{dh}}{V_T}[\mu=%.2f]$"%(a/(a+b)))
    ax.axvline(np.mean(vh), ls="--", color="b", lw=0.6)
    ax.hist(vh, bins=np.arange(0,1,.01), color="b", alpha=0.5, density=True, label=r"$\frac{V_{dh}}{V_T}[\mu=%.2f]$"%np.mean(vh),
            histtype="step")
    ax.text(0.1,0.9, "(a)", horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
    ax.set_xlim(0,1)
    ax.set_ylim(0,20)
    ax.legend(loc=1, prop={"size": 8})
    ax.set_ylabel(r"Density $\left(\frac{V_x}{V_T}\right)$")
    ax.set_xlabel(r"$\frac{V_x}{V_T}$")
    ax = fig.add_subplot(122)
    #a, b, _, _ = beta.fit(vd,floc=0,fscale=1)
    #ax.plot(x, beta.pdf(x, a, b), color="r", lw=0.8, label=r"$\frac{V_D}{V_T}[\mu=%.2f]$"%(a/(a+b)))
    ax.axvline(np.mean(vd), ls="--", color="r", lw=0.6)
    ax.hist(vd, bins=np.arange(0,1,.01), color="r", alpha=0.5, density=True, label=r"$\frac{V_D}{V_T}[\mu=%.2f]$"%np.mean(vd),
                        histtype="step")
    #a, b, _, _ = beta.fit(ve,floc=0,fscale=1)
    #ax.plot(x, beta.pdf(x, a, b), color="g", lw=0.8, label=r"$\frac{V_E}{V_T}[\mu=%.2f]$"%(a/(a+b)))
    ax.axvline(np.mean(ve), ls="--", color="g", lw=0.6)
    ax.hist(ve, bins=np.arange(0,1,.01), color="g", alpha=0.5, density=True, label=r"$\frac{V_E}{V_T}[\mu=%.2f]$"%np.mean(ve),
                        histtype="step")
    #a, b, _, _ = beta.fit(vf,floc=0,fscale=1)
    #ax.plot(x, beta.pdf(x, a, b), color="b", lw=0.8, label=r"$\frac{V_F}{V_T}[\mu=%.2f]$"%(a/(a+b)))
    ax.axvline(np.mean(vf), ls="--", color="b", lw=0.6)
    ax.hist(vf, bins=np.arange(0,1,.01), color="b", alpha=0.5, density=True, label=r"$\frac{V_F}{V_T}[\mu=%.2f]$"%np.mean(vf),
                        histtype="step")
    ax.set_ylim(0,50)
    #ax.hist(vd, bins=np.arange(0,1,.01), color="r", alpha=0.5, density=True, label=r"$\frac{v_D}{v_T}$", histtype="step")
    #ax.hist(ve, bins=np.arange(0,1,.01), color="b", alpha=0.5, density=True, label=r"$\frac{v_E}{v_T}$", histtype="step")
    #ax.hist(vf, bins=np.arange(0,1,.01), color="g", alpha=0.5, density=True, label=r"$\frac{v_F}{v_T}$", histtype="step")
    ax.text(0.1,0.9, "(b)", horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
    ax.set_xlim(0,1)
    ax.legend(loc=1, prop={"size": 8})
    ax.set_xlabel(r"$\frac{V_x}{V_T}$")
    fig.savefig("data/hist.png", bbox_inches="tight")
    return

class FanPlot(object):
    """ Plot Fan Dataset """

    def __init__(self, nrange=75, nbeam=24, r0=180, dr=45, dtheta=3.24, theta0=None):
        """
        Initialize the fanplot do a certain size.
        :param nrange: number of range gates
        :param nbeam: number of beams
        :param r0: initial beam distance - any distance unit as long as it"s consistent with dr
        :param dr: length of each radar - any distance unit as long as it"s consistent with r0
        :param dtheta: degrees per beam gate, degrees (default 3.24 degrees)
        """
        # Set member variables
        self.nrange = int(nrange)
        self.nbeam = int(nbeam)
        self.r0 = r0
        self.dr = dr
        self.dtheta = dtheta
        # Initial angle (from X, polar coordinates) for beam 0
        if theta0 == None:
            self.theta0 = (90 - dtheta * nbeam / 2)     # By default, point fanplot towards 90 deg
        else:
            self.theta0 = theta0
        return

    def add_axis(self, fig, subplot):
        ax = fig.add_subplot(subplot, polar=True)
        # Set up ticks and labels
        self.r_ticks = range(self.r0, self.r0 + (self.nrange+1) * self.dr, self.dr)
        self.theta_ticks = [self.theta0 + self.dtheta * b for b in range(self.nbeam+1)][::4]
        rlabels = [""] * len(self.r_ticks)
        for i in range(0, len(rlabels), 5):
            rlabels[i] = i
        plt.rgrids(self.r_ticks, rlabels)
        plt.thetagrids(self.theta_ticks, range(self.nbeam+1)[::4])
        return ax

    def plot(self, ax, beams, gates, color="blue"):
        """
        Add some data to the plot in a single color at positions given by "beams" and "gates".
        :param beams: a list/array of beams
        :param gates: a list/array of gates - same length as beams
        :param color: a Matplotlib color
        """
        for i, (beam, gate) in enumerate(zip(beams, gates)):
            theta = (self.theta0 + beam * self.dtheta) * np.pi / 180        # radians
            r = (self.r0 + gate * self.dr)                                  # km
            width = self.dtheta * np.pi / 180                               # radians
            height = self.dr                                                # km
            
            x1, x2 = theta, theta + width
            y1, y2 = r, r + height
            x = x1, x2, x2, x1
            y = y1, y1, y2, y2
            ax.fill(x, y, color=color)
        self._scale_plot(ax)
        return

    def _add_colorbar(self, fig, ax, bounds, colormap, label=""):
        """
        Add a colorbar to the right of an axis.
        Similar to the function in RangeTimePlot, but positioned differently fanplots.
        :param fig:
        :param ax:
        :param bounds:
        :param colormap:
        :param label:
        :return:
        """
        import matplotlib as mpl
        pos = ax.get_position()
        cpos = [pos.x1 + 0.025, pos.y0 + 0.25*pos.height,
                0.01, pos.height * 0.5]            # this list defines (left, bottom, width, height)
        cax = fig.add_axes(cpos)
        norm = mpl.colors.BoundaryNorm(bounds[::2], colormap.N)
        cb2 = mpl.colorbar.ColorbarBase(cax, cmap=colormap,
                norm=norm,
                ticks=bounds[::2],
                spacing="uniform",
                orientation="vertical")
        cb2.set_label(label)
        # Remove the outer bounds in tick labels
        ticks = [str(i) for i in bounds[::2]]
        ticks[0], ticks[-1] = "", ""
        cb2.ax.set_yticklabels(ticks)
        return

    def text(self, text, beam, gate, fontsize=8):
        theta = (self.theta0 + beam * self.dtheta + 0.8 * self.dtheta) * np.pi / 180
        r = (self.r0 + gate * self.dr)
        plt.text(theta, r, text, fontsize=fontsize)
        return

    def save(self, filepath):
        plt.tight_layout()
        plt.savefig(filepath)
        plt.close()
        return

    def _scale_plot(self, ax):
        # Scale min-max
        ax.set_thetamin(self.theta_ticks[0])
        ax.set_thetamax(self.theta_ticks[-1])
        ax.set_rmin(0)
        ax.set_rmax(self.r_ticks[-1])
        return
    
    def _monotonically_increasing(self, vec):
        if len(vec) < 2:
            return True
        return all(x <= y for x, y in zip(vec[:-1], vec[1:]))

    def plot_geo_fov(self, rad, data_dict, scans, name, start, data, skip=1,
            vel_max=100, vel_step=10,
            save=True, base_filepath=""):
        import pydarn
        import cartopy
        hdw = pydarn.read_hdw_file(rad)
        rf = rad_fov.CalcFov(hdw=hdw, ngates=self.nrange, nbeams=self.nbeam)
        lons, lats = rf.lonFull, rf.latFull

        vel_ranges = list(range(-vel_max, vel_max + 1, vel_step))
        vel_ranges.insert(0, -9999)
        vel_ranges.append(9999)
        vel_cmap = plt.get_cmap("Spectral")       # use "viridis" colormap to make this redgreen colorblind proof
        vel_colors = vel_cmap(np.linspace(0, 1, len(vel_ranges)))
        for i in scans:
            scan_time = start + dt.timedelta(minutes=i)
            fig = plt.figure(figsize=(10, 5), dpi=150)
            dat_ax = fig.add_subplot(121, projection="fovcarto",coords="geo", rad=rad, plot_date=scan_time)
            dat_ax.coastlines()
            dat_ax.overlay_radar()
            dat_ax.overlay_fov(beamLimits=[7,8], lineColor="darkred", lineWidth=0.5, ls="--")
            dat_ax.overlay_fov()
            dat_ax.grid_on()
            dat_ax.enum(bounds=[(int(np.min(lons)/10)-1)*10, (int(np.max(lons)/10)+1)*10, 25, 70])
            vel_ax = fig.add_subplot(122, projection="fovcarto",coords="geo", rad=rad, plot_date=scan_time) 
            vel_ax.coastlines()
            vel_ax.overlay_radar()
            vel_ax.overlay_fov(beamLimits=[7,8], lineColor="darkred", lineWidth=0.5, ls="--")
            vel_ax.overlay_fov()
            vel_ax.grid_on()
            vel_ax.enum()

            Vx = np.zeros((self.nbeam, self.nrange))*np.nan
            idbs, idgs = data_dict["beam"][i], data_dict["gate"][i]
            vels = data_dict["vel"][i]
            for idb, idg, vel in zip(idbs, idgs, vels):
                Vx[idb, np.round(idg).astype(int)] = vel
            Vx = np.ma.masked_invalid(Vx)
            dat_ax.pcolormesh(lons, lats, Vx, transform=cartopy.crs.PlateCarree(), cmap=plt.get_cmap("Spectral"), 
                    vmax=vel_max, vmin=-vel_max)
            dat_ax.text(1.02, 0.15, "Simulation", horizontalalignment="center",
                    verticalalignment="center", transform=dat_ax.transAxes, fontdict={"color":"red"}, rotation=90)
            Vmod = np.copy(Vx)

            Vx = np.zeros((self.nbeam, self.nrange))*np.nan
            idbs, idgs = data["beam"][i], data["gate"][i]
            vels = data["vel"][i]
            for idb, idg, vel in zip(idbs, idgs, vels):
                idb = np.array(idb)[np.array(idg) < self.nrange]
                vel = np.array(vel)[np.array(idg) < self.nrange]
                idg = np.array(idg)[np.array(idg) < self.nrange]
                if len(vel) > 0: Vx[idb, np.round(idg).astype(int)] = vel
            Vx = np.ma.masked_invalid(Vx)
            vel_ax.pcolormesh(lons, lats, Vx, transform=cartopy.crs.PlateCarree(), cmap=plt.get_cmap("Spectral"),
                    vmax=vel_max, vmin=-vel_max)
            vel_ax.text(1.02, 0.15, "Observations", horizontalalignment="center",
                    verticalalignment="center", transform=vel_ax.transAxes, fontdict={"color":"red"}, rotation=90)
            vel_ax.enum(bounds=[(int(np.min(lons)/10)-1)*10, (int(np.max(lons)/10)+1)*10, 25, 70])
            vel_ax = fig.add_subplot(122, projection="fovcarto",coords="geo", rad=rad, plot_date=scan_time) 
            rmse = np.sqrt(np.ma.sum((Vx-Vmod)**2)/np.ma.count(Vmod))
            perror = np.ma.sum(np.abs((Vx-Vmod)/Vmod)/np.ma.count(Vmod)) * 100.
            print(rmse, perror)
            if rmse>0:vel_ax.text(0.3, 0.2, r"RMdSE: %.2f $ms^{-1}$"%rmse + "\n" + r"$\delta: %.2f$"%perror + "%", horizontalalignment="center",
                    verticalalignment="center", transform=vel_ax.transAxes, fontdict={"color":"red"})

            self._add_colorbar(fig, vel_ax, vel_ranges, vel_cmap, label="Velocity [m/s]")
            if save:
                filepath = "%s/geo_%s.png" % (base_filepath, "%02d"%i)
                fig.savefig(filepath)
            fig.clf()
            plt.close()
        return

    def plot_fov(self, data_dict, scans, name, start, data, skip=1,
            vel_max=100, vel_step=10,
            save=True, base_filepath=""):
        vel_ranges = list(range(-vel_max, vel_max + 1, vel_step))
        vel_ranges.insert(0, -9999)
        vel_ranges.append(9999)
        vel_cmap = plt.cm.jet_r       # use "viridis" colormap to make this redgreen colorblind proof
        vel_colors = vel_cmap(np.linspace(0, 1, len(vel_ranges)))

        for i in scans:
            fig = plt.figure(figsize=(8,4), dpi=120)
            vel_ax = self.add_axis(fig, 122)
            dat_ax = self.add_axis(fig, 121)
            vels = data_dict["vel"][i]
            beams = data_dict["beam"][i]
            gates = data_dict["gate"][i]
            print("----------", i, skip, int(i/skip))
            d_vels = data["vel"][int(i/skip)]
            d_beams = data["beam"][int(i/skip)]
            d_gates = data["gate"][int(i/skip)]
            for k, (beam, gate, vel) in enumerate(zip(beams, gates, vels)):
                beam, gate, vel = np.array([beam]), np.array([gate]), np.array([vel])
                for s in range(len(vel_ranges) - 1):
                    step_mask = (vel >= vel_ranges[s]) & (vel <= vel_ranges[s + 1])
                    beam_s = beam[step_mask]
                    gate_s = gate[step_mask]
                    self.plot(vel_ax, beam_s, gate_s, vel_colors[s])
            # Add data
            for k, (vel, beam, gate) in enumerate(zip(d_vels, d_beams, d_gates)):
                beam, gate, vel = np.array([beam]), np.array([gate]), np.array([vel])
                for s in range(len(vel_ranges) - 1):
                    step_mask = (vel >= vel_ranges[s]) & (vel <= vel_ranges[s + 1])
                    beam_s = beam[step_mask]
                    gate_s = gate[step_mask]
                    self.plot(dat_ax, beam_s, gate_s, vel_colors[s])
            self._add_colorbar(fig, vel_ax, vel_ranges, vel_cmap, label="Velocity [m/s]")
            scan_time = start + dt.timedelta(minutes=i)
            plt.suptitle("%s \n Scan time %s UT \n Velocity" % (name, scan_time))
            if save:
                filepath = "%s_%s.png" % (base_filepath, "%02d"%i)
                self.save(filepath)
            fig.clf()
            plt.close()
        return

def plot_velocity_ts_beam(dn, rad, bmnum, model, start, end):
    """ Plot velocity TS data """
    fig = plt.figure(figsize=(5,6), dpi=150)
    axs = [fig.add_subplot(311), fig.add_subplot(312), fig.add_subplot(313)]
    mkeys = ["vd", "vf", "vt"]
    fmt = matplotlib.dates.DateFormatter("%H:%M")
    dic = "data/op/{dn}/{model}/{rad}/bm.{bm}/".format(dn=dn.strftime("%Y.%m.%d.%H.%M"),
            rad=rad, model=model, bm="%02d"%bmnum)
    fstr = glob.glob(dic + "/velocity_ti*mat")
    fstr.sort()
    #axs[0].set_title("%s UT, Radar - %s, Beam - %d, Model - %s"%(dn.strftime("%Y.%m.%d.%H.%M"), rad, bmnum, model))
    #axs[0].set_title("%s UT, Radar - %s, Beam - %d"%(dn.strftime("%Y.%m.%d.%H.%M"), rad, bmnum))
    axs[0].text(0.98, 1.05, r"Date: %s UT"%dn.strftime("%Y-%m-%d %H:%M"), horizontalalignment="right", verticalalignment="center",
            transform=axs[0].transAxes)
    axs[0].text(0.02, 1.05, "Rad: %s,  Beam: %02d"%(rad, bmnum), horizontalalignment="left", verticalalignment="center",
            transform=axs[0].transAxes)
    cols = ["r", "b", "k"]
    labs = [r"$V_{d\eta}$", r"$V_{dh}$", r"$V_{T}$"]
    fname = "data/op/{dn}/{model}/sd_{rad}_data.csv.gz".format(dn=dn.strftime("%Y.%m.%d.%H.%M"), rad=rad, model=model)
    dat = utils.get_sd_data(fname, bmnum).dropna()
    mean, std = dat.groupby("time").mean().reset_index(), dat.groupby("time").std().reset_index()
    I = 0
    for ax, mkey, col, lab in zip(axs, mkeys, cols, labs):
        ax.set_ylabel(r"Velocity, $ms^{-1}$")
        ax.set_xlabel("Time, UT")
        ax.xaxis.set_major_formatter(fmt)
        v, vmax, vmin, vstd, time = [], [], [], [], []
        for i, f in enumerate(fstr):
            sdat = loadmat(f)
            if mkey == "vt": 
                v.append(np.median(sdat["vd"]+sdat["vf"]))
                vmax.append((sdat["vd"]+sdat["vf"]).max())
                vmin.append((sdat["vd"]+sdat["vf"]).min())
                vstd.append(1.96*np.std(sdat["vd"]+sdat["vf"]))
            else:
                v.append(np.median(sdat[mkey]))
                vmax.append(sdat[mkey].max())
                vmin.append(sdat[mkey].min())
                vstd.append(1.96*np.std(sdat[mkey]))
            time.append(start + dt.timedelta(minutes=i))
        yerr = np.array([(mn, mx) for mn, mx in zip(vmin, vmax)]).T
        ax.errorbar(time, v, yerr=vstd, 
            mec=col, mfc=col, fmt="r^", ms=1.5, ls="None", ecolor=col, 
            capsize=1, capthick=.4, elinewidth=0.4,
            alpha=0.5, label=lab)
        if I == 2: 
            ax.errorbar(mean.time,  mean.v, yerr=std.v, mec="r", mfc="r", fmt="o",
                    ms=1.5, ls="None", ecolor="r",
                    capsize=1, capthick=.4, elinewidth=0.4,alpha=0.5,
                    label=r"$V_{sd}^{los}$")
            if len(mean.v) > 50:
                from scipy import signal
                vmx = signal.resample(mean.v, len(v))
                rmse = np.sqrt(np.median((vmx - np.array(v))**2))
                perror = np.mean(np.abs((vmx - np.array(v))/np.array(v)))
                ax.text(0.2, 0.85, r"RMdSE: %.2f $ms^{-1}$"%rmse + "\n" + r"$\delta: %.2f$"%perror+"%", ha="center", va="center", 
                        transform=ax.transAxes, fontdict={"color":"red", "size":8})
                print(rmse, perror)
        ax.axhline(0, color="gray", ls="--", lw=0.6)
        ax.legend(loc=1)
        ax.set_ylim(-100, 200)
        ax.set_xlim(start, end)
        I += 1
    fname = "data/op/{dn}/{model}/{rad}/bm{bm}.png".format(dn=dn.strftime("%Y.%m.%d.%H.%M"), rad=rad, model=model, bm="%02d"%bmnum)
    fig.autofmt_xdate()
    fig.savefig(fname,bbox_inches="tight")
    return


class SensitivityAnalysis(object):
    """ Sensitivity Analysis """

    def __init__(self, problem, ds):
        """ Initialize parameters """
        self.problem = problem
        self.ds = ds
        return

    def _hist_(self):
        """ Histogram of outputs """
        fig, ax = plt.subplots(figsize=(9,3), nrows=1, ncols=3, sharey=True)
        labels = [r"$V_{d\eta}$ [m/s]", r"$V_{dh}$ [m/s]", r"$V_{t}$ [m/s]"]
        params = ["vd_mean", "vf_mean", "vt_mean"]
        for i, lab, pm in zip(range(3), labels, params):
            ax[i].hist(self.ds.variables[pm][:].ravel(), 20)
            ax[i].set_xlabel(lab)
        ax[0].set_ylabel("Counts")    
        fig.subplots_adjust(wspace=0.1)
        fig.savefig("data/sim/histogram.png", bbox_inches="tight")
        return

    def _regression_(self):
        """ Regression Analysis """
        import scipy
        import seaborn as sns
        ylabels = [r"$V_{d\eta}$ [m/s]", r"$V_{dh}$ [m/s]", r"$V_{t}$ [m/s]"]
        xlabels = [r"$Ratio_{D}$", r"$Ratio_{E}$", r"Ratio_{F}"]
        yparam = ["vd_mean", "vf_mean", "vt_mean"]
        xparam = ["d_ratio", "e_ratio", "f_ratio"]
        print(self.ds.variables["parameters"][:].shape)
        for i, ylab, yp in zip(range(3), ylabels, yparam):
            fig, ax = plt.subplots(1, 3, sharey=True)
            y = self.ds.variables[yp][:].ravel()
            for j, xlab, xp, a in zip(range(3), xlabels, xparam, ax):
                x = self.ds.variables["parameters"][:][:,j]
                sns.regplot(x, y, ax=a, ci=None, color="k",scatter_kws={"alpha":0.2, "s":4, "color":"gray"})
                pearson = scipy.stats.pearsonr(x, y)
                a.annotate("r: {:6.3f}".format(pearson[0]), xy=(0.15, 0.85), xycoords="axes fraction",fontsize=13)
                a.set_xlabel(xlab)
                if j==0: a.set_ylabel(ylab)
            fig.savefig("data/sim/reg_{pm}.png".format(pm=yp), bbox_inches="tight")
            plt.close()
        return

    def _normalize_(self, x, xmin, xmax):
        return (x-xmin)/(xmax-xmin)

    def _plot_circles_(self, ax, locs, names, max_s, stats, smax, smin, fc, ec, lw, zorder):
        s = np.asarray([stats[name] for name in names])
        s = 0.01 + max_s * np.sqrt(self._normalize_(s, smin, smax))
        
        fill = True
        for loc, name, si in zip(locs, names, s):
            if fc=="w": fill=False
            else: ec="none"           
            x = np.cos(loc)
            y = np.sin(loc)
            
            circle = plt.Circle((x,y), radius=si, ec=ec, fc=fc, transform=ax.transData._b,
                    zorder=zorder, lw=lw, fill=True)
            ax.add_artist(circle)
        return

    def _filter_(self, sobol_indices, names, locs, criterion, threshold):
        if criterion in ["ST", "S1", "S2"]:
            data = sobol_indices[criterion]
            data = np.abs(data)
            data = data.flatten() # flatten in case of S2
            # TODO:: remove nans
            filtered = ([(name, locs[i]) for i, name in enumerate(names) if
                data[i]>threshold])
            filtered_names, filtered_locs = zip(*filtered)
        elif criterion in ["ST_conf", "S1_conf", "S2_conf"]: raise NotImplementedError
        else: raise ValueError("unknown value for criterion")
        return filtered_names, filtered_locs

    def _legend_(self, ax):
        some_identifiers = [plt.Circle((0,0), radius=5, color="k", fill=False, lw=1),
                plt.Circle((0,0), radius=5, color="k", fill=True),
                plt.Line2D([0,0.5], [0,0.5], lw=8, color="darkgray")]
        ax.legend(some_identifiers, ["ST", "S1", "S2"],
                loc=(1,0.75), borderaxespad=0.1, mode="expand",
                handler_map={plt.Circle: HandlerCircle()})
        return

    def _plot_sobol_indices_(self, sobol_indices, criterion="ST", threshold=0.01):
        max_linewidth_s2 = 15#25*1.8
        max_s_radius = 0.3
        sobol_stats = {key:sobol_indices[key] for key in ["ST", "S1"]}
        sobol_stats = pd.DataFrame(sobol_stats, index=self.problem["names"])
        smax = sobol_stats.max().max()
        smin = sobol_stats.min().min()
        s2 = pd.DataFrame(sobol_indices["S2"], index=self.problem["names"],
                columns=self.problem["names"])
        s2[s2<0.0]=0. #Set negative values to 0 (artifact from small sample sizes)
        s2max = s2.max().max()
        s2min = s2.min().min()
        
        names = self.problem["names"]
        n = len(names)
        ticklocs = np.linspace(0, 2*pi, n+1)
        locs = ticklocs[0:-1]
        
        filtered_names, filtered_locs = self._filter_(sobol_indices, names, locs,
                criterion, threshold)
        
        # setup figure
        xnames = copy.copy(names)
        xnames.extend(["D-Ratio"])
        fig = plt.figure()
        ax = fig.add_subplot(111, polar=True)
        ax.grid(False)
        ax.spines["polar"].set_visible(False)
        ax.set_xticks(ticklocs)
        ax.set_xticklabels(xnames)
        ax.set_yticklabels([])
        ax.set_ylim(top=1.4)
        self._legend_(ax)
        # plot ST
        self._plot_circles_(ax, filtered_locs, filtered_names, max_s_radius,
                sobol_stats["ST"], smax, smin, "w", "k", 1, 9)
        # plot S1
        self._plot_circles_(ax, filtered_locs, filtered_names, max_s_radius,
                sobol_stats["S1"], smax, smin, "k", "k", 1, 10)
        # plot S2
        for name1, name2 in itertools.combinations(zip(filtered_names, filtered_locs), 2):
            name1, loc1 = name1
            name2, loc2 = name2
            weight = s2.loc[name1, name2]
            lw = 0.5+max_linewidth_s2*self._normalize_(weight, s2min, s2max)
            ax.plot([loc1, loc2], [1,1], c="darkgray", lw=lw, zorder=1)
        return fig

    def analyze(self, regs=False):
        """ Analyze and plots sensitivity test results """
        self._hist_()
        if regs: print("None")#self._regression_()
        else:
            labels = [r"$V_{d\eta}$ [m/s]", r"$V_{dh}$ [m/s]", r"$V_{t}$ [m/s]"]
            params = ["vd_mean", "vf_mean", "vt_mean"]
            for i, lab, pm in zip(range(3), labels, params):
                Si = sobol.analyze(self.problem, self.ds.variables[pm][:].ravel(), calc_second_order=True, print_to_console=False)
                Si_filter = {k:Si[k] for k in ["ST","ST_conf","S1","S1_conf"]}
                Si_df = pd.DataFrame(Si_filter, index=self.problem["names"])
                fig, ax = plt.subplots(1)
                indices = Si_df[["S1","ST"]]
                err = Si_df[["S1_conf","ST_conf"]]
                indices.plot.bar(yerr=err.values.T,ax=ax)
                fig.set_size_inches(4,4)
                fig.savefig("data/sim/sens_{pm}.png".format(pm=pm), bbox_inches="tight")
                plt.close()
                fig = self._plot_sobol_indices_(Si, criterion="ST", threshold=0.005)
                fig.set_size_inches(4,4)
                fig.savefig("data/sim/intv_{pm}.png".format(pm=pm), bbox_inches="tight")
                plt.close()
        return

class ModelSensitivity(object):
    """ Sensitivity Analysis """

    def __init__(self, ds):
        """ Initialize parameters """
        self.problem =  {
                "num_vars": 3,
                "names": ["D-Ratio", "E-Ratio", "F-Ratio"],
                "bounds": [[np.min(ds.variables["d_ratio"][:]), np.max(ds.variables["d_ratio"][:])],
                    [np.min(ds.variables["e_ratio"][:]), np.max(ds.variables["e_ratio"][:])],
                    [np.min(ds.variables["f_ratio"][:]), np.min(ds.variables["f_ratio"][:])]]
                }
        self.ds = ds
        print(ds.variables.keys())
        return

    def _hist_(self):
        """ Histogram of outputs """
        fig, ax = plt.subplots(figsize=(9,9), nrows=3, ncols=3, sharey=True, sharex=False)
        labels = [[r"$V_{d\eta}$ [m/s]", r"$V_{dh}$ [m/s]", r"$V_{t}$ [m/s]"],
                [r"$V_{d\eta}^{max}$ [m/s]", r"$V_{dh}^{max}$ [m/s]", r"$V_{t}^{max}$ [m/s]"],
                [r"$V_{d\eta}^{min}$ [m/s]", r"$V_{dh}^{min}$ [m/s]", r"$V_{t}^{min}$ [m/s]"]]
        params = [["vn", "vh", "vt"], ["vn_max", "vh_max", "vt_max"], ["vn_min", "vh_min", "vt_min"]]
        bins = range(-10,110,4)
        nx = np.arange(-20,160)
        for i, labs, pms in zip(range(3), labels, params):
            for j, lab, pm in zip(range(3), labs, pms):
                u,loc,scale = skewnorm.fit(3*self.ds.variables[pm][:].ravel(),
                        floc=np.mean(3*self.ds.variables[pm][:].ravel()),
                        fscale=np.std(3*self.ds.variables[pm][:].ravel()))
                ax[i,j].hist(3*self.ds.variables[pm][:].ravel(), bins=bins)
                am = ax[i,j].twinx()
                am.set_yticklabels([])
                am.plot(nx, skewnorm.pdf(nx, a=u, loc=loc, scale=scale), "r", lw=1.5)
                ax[i,j].set_xlabel(lab)
                ax[i,j].set_xlim(-20, 160)
            ax[i,0].set_ylabel("Counts")
        fig.subplots_adjust(wspace=0.1, hspace=0.3)
        fig.savefig("data/sim/histogram.png", bbox_inches="tight")
        return

    def _regression_(self):
        """ Regression Analysis """
        import scipy
        import seaborn as sns
        ylabels = [r"$V_{d\eta}$ [m/s]", r"$V_{dh}$ [m/s]", r"$V_{t}$ [m/s]"]
        xlabels = [r"$Ratio_{D}$", r"$Ratio_{E}$", r"$Ratio_{F}$", r"$Rate_{D}$", r"$Rate_{E}$", r"$Rate_{F}$", r"$Frequency$", "SZA"]
        xlabels = [r"$R_{D}$", r"$R_{E}$", r"$R_{F}$", r"$R^r_{D}$", r"$R^r_{E}$", r"$R^r_{F}$", r"$Frequency$", "SZA"]
        yparam = ["vn", "vh", "vt"]
        xparam = ["d_ratio", "e_ratio", "f_ratio", "d_rate", "e_rate", "f_rate", "frequency", "sza"]
        token = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"]
        Xx = np.array([self.ds.variables["d_ratio"][:], self.ds.variables["e_ratio"][:], self.ds.variables["f_ratio"][:],
            self.ds.variables["d_rate"][:], self.ds.variables["e_rate"][:], self.ds.variables["f_rate"][:],
            self.ds.variables["frequency"][:]/1000, self.ds.variables["sza"][:]]).T
        for i, ylab, yp in zip(range(3), ylabels, yparam):
            fig, ax = plt.subplots(2, 4, sharey=True, figsize=(10,5))
            y = 3*self.ds.variables[yp][:].ravel()
            minfo = MIR(Xx,y)
            for j, xlab, xp in zip(range(8), xlabels, xparam):
                a = ax[np.mod(j,2), int(j/2)]
                x = Xx[:,j]#self.ds.variables["parameters"][:][:,j]
                #print(te.te_compute(x, np.array(y), k=1, embedding=1, safetyCheck=False, GPU=False))
                sns.regplot(x, y, ax=a, ci=95, color="k",scatter_kws={"alpha":0.2, "s":1.5, "color":"red"})
                pearson = scipy.stats.pearsonr(x, y)
                a.annotate("r: {:1.2f}, MI: {:1.2f}".format(pearson[0], minfo[j]),
                        xy=(0.15, 0.85), xycoords="axes fraction",fontsize=10)
                a.set_xlabel(xlab)
                a.text(0.9, 0.9, token[j], horizontalalignment="center", verticalalignment="center",
                        transform=a.transAxes)
                if j==0: a.set_ylabel(ylab)
            fig.subplots_adjust(wspace=0.1, hspace=0.5)
            fig.savefig("data/sim/reg_{pm}.png".format(pm=yp), bbox_inches="tight")
            plt.close()
        return

    def _normalize_(self, x, xmin, xmax):
        return (x-xmin)/(xmax-xmin)

    def _plot_circles_(self, ax, locs, names, max_s, stats, smax, smin, fc, ec, lw, zorder):
        s = np.asarray([stats[name] for name in names])
        s = 0.01 + max_s * np.sqrt(self._normalize_(s, smin, smax))
        
        fill = True
        for loc, name, si in zip(locs, names, s):
            if fc=="w": fill=False
            else: ec="none"           
            x = np.cos(loc)
            y = np.sin(loc)
            
            circle = plt.Circle((x,y), radius=si, ec=ec, fc=fc, transform=ax.transData._b,
                    zorder=zorder, lw=lw, fill=True)
            ax.add_artist(circle)
        return

    def _filter_(self, sobol_indices, names, locs, criterion, threshold):
        if criterion in ["ST", "S1", "S2"]:
            data = sobol_indices[criterion]
            data = np.abs(data)
            data = data.flatten() # flatten in case of S2
            # TODO:: remove nans
            filtered = ([(name, locs[i]) for i, name in enumerate(names) if
                data[i]>threshold])
            filtered_names, filtered_locs = zip(*filtered)
        elif criterion in ["ST_conf", "S1_conf", "S2_conf"]: raise NotImplementedError
        else: raise ValueError("unknown value for criterion")
        return filtered_names, filtered_locs

    def _legend_(self, ax):
        some_identifiers = [plt.Circle((0,0), radius=5, color="k", fill=False, lw=1),
                plt.Circle((0,0), radius=5, color="k", fill=True),
                plt.Line2D([0,0.5], [0,0.5], lw=8, color="darkgray")]
        ax.legend(some_identifiers, ["ST", "S1", "S2"],
                loc=(1,0.75), borderaxespad=0.1, mode="expand",
                handler_map={plt.Circle: HandlerCircle()})
        return

    def _plot_sobol_indices_(self, sobol_indices, criterion="ST", threshold=0.01):
        max_linewidth_s2 = 15#25*1.8
        max_s_radius = 0.3
        sobol_stats = {key:sobol_indices[key] for key in ["ST", "S1"]}
        sobol_stats = pd.DataFrame(sobol_stats, index=self.problem["names"])
        smax = sobol_stats.max().max()
        smin = sobol_stats.min().min()
        s2 = pd.DataFrame(sobol_indices["S2"], index=self.problem["names"],
                columns=self.problem["names"])
        s2[s2<0.0]=0. #Set negative values to 0 (artifact from small sample sizes)
        s2max = s2.max().max()
        s2min = s2.min().min()
        
        names = self.problem["names"]
        n = len(names)
        ticklocs = np.linspace(0, 2*pi, n+1)
        locs = ticklocs[0:-1]
        
        filtered_names, filtered_locs = self._filter_(sobol_indices, names, locs,
                criterion, threshold)
        
        # setup figure
        xnames = copy.copy(names)
        xnames.extend(["D-Ratio"])
        fig = plt.figure()
        ax = fig.add_subplot(111, polar=True)
        ax.grid(False)
        ax.spines["polar"].set_visible(False)
        ax.set_xticks(ticklocs)
        ax.set_xticklabels(xnames)
        ax.set_yticklabels([])
        ax.set_ylim(top=1.4)
        self._legend_(ax)
        # plot ST
        self._plot_circles_(ax, filtered_locs, filtered_names, max_s_radius,
                sobol_stats["ST"], smax, smin, "w", "k", 1, 9)
        # plot S1
        self._plot_circles_(ax, filtered_locs, filtered_names, max_s_radius,
                sobol_stats["S1"], smax, smin, "k", "k", 1, 10)
        # plot S2
        for name1, name2 in itertools.combinations(zip(filtered_names, filtered_locs), 2):
            name1, loc1 = name1
            name2, loc2 = name2
            weight = s2.loc[name1, name2]
            lw = 0.5+max_linewidth_s2*self._normalize_(weight, s2min, s2max)
            ax.plot([loc1, loc2], [1,1], c="darkgray", lw=lw, zorder=1)
        return fig

    def analyze(self, regs=True):
        """ Analyze and plots sensitivity test results """
        self._hist_()
        if not regs: 
            print("None")
            self._regression_()
        else:
            print(regs)
            labels = [r"$V_{d\eta}$ [m/s]", r"$V_{dh}$ [m/s]", r"$V_{t}$ [m/s]"]
            params = ["vn", "vh", "vt"]
            for i, lab, pm in zip(range(3), labels, params):
                v = self.ds.variables[pm][:].ravel()
                x = np.array([self.ds.variables["d_rate"][:], self.ds.variables["e_rate"][:], self.ds.variables["f_rate"][:]]).T
                print(v.shape, x.shape)
                #Si = rbd_fast.analyze(self.problem, x, v, M=10, print_to_console=False)
                #print(Si)
                Si = sobol.analyze(self.problem, self.ds.variables[pm][:].ravel(), calc_second_order=True, print_to_console=False)
                #print(Si)
                Si_filter = {k:Si[k] for k in ["ST","ST_conf","S1","S1_conf"]}
                Si_df = pd.DataFrame(Si_filter, index=self.problem["names"])
                fig, ax = plt.subplots(1)
                indices = Si_df[["S1","ST"]]
                err = Si_df[["S1_conf","ST_conf"]]
                indices.plot.bar(yerr=err.values.T,ax=ax)
                fig.set_size_inches(4,4)
                fig.savefig("data/sim/sens_{pm}.png".format(pm=pm), bbox_inches="tight")
                plt.close()
                #fig = self._plot_sobol_indices_(Si, criterion="ST", threshold=0.005)
                #fig.set_size_inches(4,4)
                #fig.savefig("data/sim/intv_{pm}.png".format(pm=pm), bbox_inches="tight")
                #plt.close()
        return


def plot_edens_versus_height(eDensPC, eDensAC, ylim=[50,350]):
    fig, axes = plt.subplots(figsize=(15,6), nrows=2, ncols=5, sharey=True, sharex=False)
    from scipy import stats
    for i in range(5):
        x, y = np.array(eDensPC[i+16]), np.array(eDensAC[i+16])
        xmean, ymean = np.quantile(x, q=.56, axis=0), np.quantile(y, q=.56, axis=0) #np.median(x, axis=0), np.median(y, axis=0)
        xstd, ystd = 0.3*stats.median_absolute_deviation(x, axis=0), 0.3*stats.median_absolute_deviation(y, axis=0)
        xl, xu = utils.smooth(np.quantile(x, q=.5, axis=0), window_len=51),\
                utils.smooth(np.quantile(x, q=.62, axis=0), window_len=51)
        yl, yu = utils.smooth(np.quantile(y, q=.5, axis=0), window_len=51),\
                utils.smooth(np.quantile(y, q=.62, axis=0), window_len=51)
        xmean, ymean = utils.smooth(xmean, window_len=51), utils.smooth(ymean, window_len=51)
        ax = axes[0, i]
        ax.semilogx(xmean, np.arange(50,350,1).ravel(), "ro", lw=0.8,  markersize=1)
        ax.fill_betweenx(np.arange(50,350,1).ravel(), x1=xl, x2=xu, alpha=0.3, color="r")
        ax.set_xlim(.01, 10)
        if i==0: ax.set_ylabel("Height, km")
        ax.set_xlabel("Percentage Change")
        ax = axes[1, i]
        ax.semilogx(ymean, np.arange(50,350,1).ravel(), "ro", lw=0.8,  markersize=1)
        ax.fill_betweenx(np.arange(50,350,1).ravel(), x1=yl, x2=yu, alpha=0.3, color="r")
        ax.set_xlim(.1, 10000)
        if i==0: ax.set_ylabel("Height, km")
        ax.set_xlabel("Absolute Change")
    fig.subplots_adjust(hspace=0.3)
    fig.savefig("data/sens.png", bbox_inches="tight")
    return

def plot_ssi_versus_bins(irr, wavelength, ylim=[50,350]):
    fig, ax = plt.subplots(figsize=(4,4), nrows=1, ncols=1, sharey=True, sharex=False)
    xmean = np.mean(irr, axis=0)#np.quantile(irr, q=.56, axis=0)
    std = np.std(irr, axis=0)
    print(xmean.shape)
    ax.loglog(wavelength, xmean, "ro", lw=0.8,  markersize=1)
    ax.errorbar(wavelength, xmean, yerr=std, capthick=1, elinewidth=0.8, capsize=1, ecolor="r", marker="o", ls="None", ms=1, mfc="k", mec="k")
    ax.set_ylim(1e5,1e12)
    ax.set_xlabel(r"$\Lambda$ (A)")
    ax.set_ylabel(r"$I_{o}$ ($Wm^{-2}$)")
    fig.savefig("data/sens.b.png", bbox_inches="tight")
    return


def plot_ray_edens(ev=dt.datetime(2015,5,5,21,51), rad="bks", time=18, maxground=1500, maxalt=300, step=1,
        showrefract=True, nr_cmap="jet_r", nr_lim=[-0.5, 0.5],
        raycolor="0.3", title=True, zorder=2, alpha=1,
        fig=None, rect=111, ax=None, aax=None, freq=12., diff=True):
    ax, aax = curvedEarthAxes(fig=fig, rect=rect, maxground=maxground, maxalt=maxalt, nyticks=3)
    dic =  "data/op/2015.05.05.22.11/waccmx/bks/bm.07/"
    files = glob.glob(dic + "ti({ti})_elv(*)_{case}.csv".format(ti="%02d"%time, case="f"))
    files.sort()
    Re = 6371.
    fx = []
    ry = 5
    print("I'm here")
    for f in files[::ry]:
        th, r, f, _, _ = get_polar(pd.read_csv(f))
        fx.append(trapz(signal.resample(f,INT_F)))
        v = (0.5 * f * 3e8 / (freq * 1e6))
        if not showrefract: aax.plot(th, r, c=raycolor, zorder=zorder, alpha=alpha)
        else:
            points = np.array([th, r]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lcol = LineCollection( segments, zorder=zorder, alpha=alpha)
            _ = lcol.set_cmap( nr_cmap )
            _ = lcol.set_norm( plt.Normalize(*nr_lim) )
            _ = lcol.set_array( utils.smooth(v, window_len=21) )
            _ = aax.add_collection( lcol )
    dx = 0.1
    aax.plot(np.arange(0,2000,dx)/Re, np.ones(int(2000*1/dx))*60+Re, color="b", lw=1.2, alpha=0.7)
    aax.plot(np.arange(0,2000,dx)/Re, np.ones(int(2000*1/dx))*95+Re, color="orange", lw=1.2, alpha=0.7)
    aax.plot(np.arange(0,2000,dx)/Re, np.ones(int(2000*1/dx))*130+Re, color="r",  lw=1.2, alpha=0.7)
    if showrefract:
        cbax = addColorbar(lcol, ax)
        _ = cbax.set_ylabel(r"$\Delta V_{d\eta}(t_i), ms^{-1}$", size=8)
        dv = (0.5 * np.array(fx) * 3e8 / (freq * 1e6))
        ax.text(0.99, 1.05, r"$V_{d\eta}(t_i)=%.2f$ $ms^{-1}$"%np.median(dv), horizontalalignment="right", verticalalignment="center",
                transform=ax.transAxes, fontdict={"size":8})
        #ax.text(0.01, 1.1, "Rad: bks,  Beam: 07\nDate, $t_i$: %s UT"%(ev+dt.timedelta(minutes=time)).strftime("%Y-%m-%d %H:%M")
        #    , horizontalalignment="left", verticalalignment="center",
        #    transform=ax.transAxes, fontdict={"size":8})
        ax.set_xlabel("Ground Range, $km$", fontdict={"size":8})
        ax.set_ylabel("Height, $km$", fontdict={"size":8})
    else: cbax = None
    fig = ax.get_figure()
    fig.savefig("data/figs/rt.dvn.ti({ti}).{case}.png".format(ti="%02d"%time, case="f"), bbox_inches="tight")
    plt.close()
    
    ax, aax = curvedEarthAxes(fig=None, rect=rect, maxground=maxground, maxalt=maxalt, nyticks=3)
    dfiles = glob.glob(dic + "ti({ti})_elv(*)_{case}.csv".format(ti="%02d"%(time-1), case="f"))
    dfiles.sort()
    from scipy.io import loadmat
    vel = loadmat(dic+"velocity_ti({ti}).mat".format(ti="%02d"%time))["vf"]
    kx = 0
    for f, df in zip(files[::ry], dfiles[::ry]):
        if diff:
            if kx < 1:
                th, r, f, _, _ = get_polar(pd.read_csv(f))
                dth, dr, df, _, _ = get_polar(pd.read_csv(df))
                dh = -(np.max(dr) - np.max(r))
                ax.text(0.95, 1.02, r"$\Delta h=%.2f$ $km$"%dh, horizontalalignment="right", verticalalignment="center",
                        transform=ax.transAxes, fontdict={"size":8})
                aax.plot(th, r, c=raycolor, zorder=zorder, alpha=alpha, lw=0.8)
                aax.plot(dth, dr, c="b", zorder=zorder, alpha=alpha, ls="--",lw=1.6)
                #aax.plot(np.arange(0,2000,dx)/Re, np.ones(int(2000*1/dx))*np.max(dr), color="b", ls="--", lw=0.6, alpha=0.7)
                #aax.plot(np.arange(0,2000,dx)/Re, np.ones(int(2000*1/dx))*np.max(r), color="r", ls="--", lw=0.6, alpha=0.7)
                
                axins = ax.inset_axes([0.4, -.8, 0.3, 0.5])
                axins.plot(np.arange(0,2000,dx)/Re, np.ones(int(2000*1/dx))*np.max(dr), color="b", ls="--", lw=0.6, alpha=0.7)
                axins.plot(np.arange(0,2000,dx)/Re, np.ones(int(2000*1/dx))*np.max(r), color="k", ls="--", lw=0.6, alpha=0.7)
                axins.plot(th, r, c=raycolor, zorder=zorder, alpha=alpha, lw=0.8)
                axins.plot(dth, dr, c="b", zorder=zorder, alpha=alpha, ls="--",lw=1.6)
                axins.set_ylim(int(np.max(dr)-10), int(np.max(r)+10))
                lenx = np.argmax(r)
                axins.set_xlim(th[lenx-10], th[lenx+10])
                axins.set_yticks(np.linspace(int(np.max(dr)-10), int(np.max(r)+10), 3))
                axins.set_yticklabels((np.linspace(int(np.max(dr)-10), int(np.max(r)+10), 3)-Re).astype(int), fontdict={"size":7})
                axins.set_xticks(np.linspace(th[lenx-10], th[lenx+5], 4))
                axins.set_xticklabels((np.linspace(th[lenx-10], th[lenx+10], 4)*Re).astype(int), fontdict={"size":7})
                axins.set_xlabel("Ground Range, $km$", fontdict={"size":8})
                axins.set_ylabel("Height, $km$", fontdict={"size":8})
                aax.indicate_inset_zoom(axins)
        else:
            th, r, f, _, _ = get_polar(pd.read_csv(f))
            dth, dr, df, _, _ = get_polar(pd.read_csv(df))
            aax.plot(th, r, c=raycolor, zorder=zorder, alpha=alpha, lw=0.8)
        kx += 1
    aax.plot(np.arange(0,2000,dx)/Re, np.ones(int(2000*1/dx))*60+Re, color="b", lw=1.2, alpha=0.7)
    aax.plot(np.arange(0,2000,dx)/Re, np.ones(int(2000*1/dx))*95+Re, color="orange", lw=1.2, alpha=0.7)
    aax.plot(np.arange(0,2000,dx)/Re, np.ones(int(2000*1/dx))*130+Re, color="r",  lw=1.2, alpha=0.7)
    #ax.text(0.02, 1.05, "Rad: bks,  Beam: 07\nDate, $t_i$: %s UT"%(ev+dt.timedelta(minutes=time)).strftime("%Y-%m-%d %H:%M")
    #    , horizontalalignment="left", verticalalignment="center",
    #        transform=ax.transAxes)
    if diff:
        ax.text(0.01, 1.05, r"$V_{dh}(t_i)=%.2f$ $ms^{-1}$"%np.median(vel), horizontalalignment="left", verticalalignment="center",
                transform=ax.transAxes, fontdict={"size":8})
    else:
        ax.text(0.99, 1.05, r"$V_{dh}(t_i)=%.2f$ $ms^{-1}$"%np.median(vel), horizontalalignment="right", verticalalignment="center",
                transform=ax.transAxes, fontdict={"size":8})
    fig = ax.get_figure()
    fig.savefig("data/figs/rt.dvh.ti({ti}).{case}.png".format(ti="%02d"%time, case="f"), bbox_inches="tight")
    plt.close()

    ax, aax = curvedEarthAxes(fig=None, rect=rect, maxground=maxground, maxalt=maxalt, nyticks=3)
    ne = loadmat(dic+"ne_ti({ti})_f.mat".format(ti="%02d"%time))["ne"]
    ne_bgc = loadmat(dic+"ne_ti({ti})_f.mat".format(ti="%02d"%(time-1)))["ne"]
    if diff: ne = np.abs(ne - ne_bgc)
    x, y = np.meshgrid(np.arange(0,2000,10)/Re, np.arange(50, 350)+Re)
    if diff: mp = aax.pcolor(x, y, ne, cmap="plasma", norm = matplotlib.colors.LogNorm(vmin=1e2, vmax=1e5))
    else: mp = aax.pcolor(x, y, ne, cmap="plasma", norm = matplotlib.colors.LogNorm(vmin=1e2, vmax=1e6))
    cbax = addColorbar(mp, ax)
    if diff: _ = cbax.set_ylabel(r"$\Delta N_e(t_i,t_{i-1}), cm^{-3}$")
    else: _ = cbax.set_ylabel(r"$N_e(t_i), cm^{-3}$")
    ax.text(0.01, 1.1, "Rad: bks,  Beam: 07\nDate, $t_i$: %s UT"%(ev+dt.timedelta(minutes=time)).strftime("%Y-%m-%d %H:%M"), 
        horizontalalignment="left", verticalalignment="center",
            transform=ax.transAxes)
    #ax.text(0.95, 1.05, r"$N_e^{250}=%.2f\times 10^5$ $cm^{-3}$"%np.max(ne[200,:]/1e5), horizontalalignment="right", verticalalignment="center",
    #                    transform=ax.transAxes)
    fig = ax.get_figure()
    fig.savefig("data/figs/density.ti({ti}).{case}.png".format(ti="%02d"%time, case="f"))
    plt.close()
    return
