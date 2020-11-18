#!/usr/bin/env python

"""plotutil.py: module is dedicated to plottting."""

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

from PyIF import te_compute as te
from sklearn.feature_selection import mutual_info_regression as MIR
from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.analyze import rbd_fast

import itertools
from math import pi
from matplotlib.legend_handler import HandlerPatch
class HandlerCircle(HandlerPatch):
    def create_artists(self, legend, orig_handle,
            xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = plt.Circle(xy=center, radius=orig_handle.radius)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]

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
            Size.Fixed(0.1), # padding
            Size.Fixed(0.2)] # colorbar
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
    _ = plt.colorbar(mappable, cax=cbax, shrink=0.7)
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
    if not fig: fig = plt.figure(figsize=(10,6))
    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    # adjust axis
    ax1.axis["left"].label.set_text(r"Alt. [km]")
    ax1.axis["bottom"].label.set_text(r"Ground range [km]")
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
        ax.set_ylabel("Velocity [m/s]")
        ax.set_xlabel("Time [UT]")
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
        showrefract=False, nr_cmap="jet_r", nr_lim=[-0.5, 0.5], 
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
    files = glob.glob(dic + "ti({ti}).elv(*).{case}.csv".format(ti="%02d"%ti, case=case))
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
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111)
    ax.hist(vd, bins=np.arange(0,1,.005), color="r", alpha=0.5, density=True, label=r"$\frac{v_D}{v_T}$", histtype="step")
    ax.hist(ve, bins=np.arange(0,1,.005), color="b", alpha=0.5, density=True, label=r"$\frac{v_E}{v_T}$", histtype="step")
    ax.hist(vf, bins=np.arange(0,1,.005), color="g", alpha=0.5, density=True, label=r"$\frac{v_F}{v_T}$", histtype="step")
    ax.set_xlim(0,1)
    ax.legend(loc=1)
    ax.set_ylabel("Density")
    ax.set_xlabel(r"$\frac{v_x}{v_T}$")
    fig.savefig("data/hist_reg.png", bbox_inches="tight")
    return

def plot_distribution(vn, vf):
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111)
    ax.hist(vn, bins=np.arange(0,1,.005), color="r", alpha=0.5, density=True, label=r"$\frac{v_{\eta}}{v_T}$", histtype="step")
    ax.hist(vf, bins=np.arange(0,1,.005), color="b", alpha=0.5, density=True, label=r"$\frac{v_{h}}{v_T}$", histtype="step")
    ax.set_xlim(0,1)
    ax.legend(loc=1)
    ax.set_ylabel("Density")
    ax.set_xlabel(r"$\frac{v_x}{v_T}$")
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
    fig = plt.figure(figsize=(6,6), dpi=150)
    axs = [fig.add_subplot(311), fig.add_subplot(312), fig.add_subplot(313)]
    mkeys = ["vd", "vf", "vt"]
    fmt = matplotlib.dates.DateFormatter("%H:%M")
    dic = "data/op/{dn}/{model}/{rad}/bm.{bm}/".format(dn=dn.strftime("%Y.%m.%d.%H.%M"),
            rad=rad, model=model, bm="%02d"%bmnum)
    fstr = glob.glob(dic + "/velocity.ti*mat")
    fstr.sort()
    axs[0].set_title("%s UT, Radar - %s, Beam - %d, Model - %s"%(dn.strftime("%Y.%m.%d.%H.%M"), rad, bmnum, model))
    cols = ["r", "b", "k"]
    labs = [r"$V_{d\eta}$", r"$V_{dh}$", r"$V_{t}$"]
    fname = "data/op/{dn}/{model}/sd_{rad}_data.csv.gz".format(dn=dn.strftime("%Y.%m.%d.%H.%M"), rad=rad, model=model)
    dat = utils.get_sd_data(fname, bmnum).dropna()
    dat = dat.groupby("time").mean().reset_index()
    I = 0
    for ax, mkey, col, lab in zip(axs, mkeys, cols, labs):
        ax.set_ylabel("Velocity [m/s]")
        ax.set_xlabel("Time [UT]")
        ax.xaxis.set_major_formatter(fmt)
        v, vmax, vmin, time = [], [], [], []
        for i, f in enumerate(fstr):
            sdat = loadmat(f)
            if mkey == "vt": 
                v.append(np.median(sdat["vd"]+sdat["vf"]))
                vmax.append((sdat["vd"]+sdat["vf"]).max())
                vmin.append((sdat["vd"]+sdat["vf"]).min())
            else:
                v.append(np.median(sdat[mkey]))
                vmax.append(sdat[mkey].max())
                vmin.append(sdat[mkey].min())
            time.append(start + dt.timedelta(minutes=i))
        yerr = np.array([(mn, mx) for mn, mx in zip(vmin, vmax)]).T
        ax.errorbar(time, v, yerr=yerr, 
            mec=col, mfc=col, fmt="r^", ms=1.5, ls="None", ecolor=col, 
            capsize=1, capthick=.4, elinewidth=0.4,
            alpha=0.5, label=lab)
        if I == 2: 
            ax.plot(dat.time, dat.v, color="darkred", marker="o",
                    alpha=0.7, ls="None", markersize=1.5, label=r"$V_{sd}^{los}$")
        ax.axhline(0, color="gray", ls="--", lw=0.6)
        ax.legend(loc=1)
        ax.set_ylim(-100, 200)
        ax.set_xlim(start, end)
        I += 1
    fname = "data/op/{dn}/{model}/{rad}/bm{bm}.png".format(dn=dn.strftime("%Y.%m.%d.%H.%M"), rad=rad, model=model, bm="%02d"%bmnum)
    fig.savefig(fname,bbox_inches="tight")
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
    ax.loglog(wavelength, xmean, "ro", lw=0.8,  markersize=1)
    ax.errorbar(wavelength, xmean, yerr=std, capthick=1, elinewidth=0.8, capsize=1, ecolor="r", marker="o", ls="None", ms=1, mfc="k", mec="k")
    ax.set_ylim(1e5,1e12)
    ax.set_xlabel(r"$\Lambda$ (A)")
    ax.set_ylabel(r"$I_{o}$ ($Wm^{-2}$)")
    fig.savefig("data/sens.b.png", bbox_inches="tight")
    return

def plot_rti(rad="bks", stime=dt.datetime(2015,5,5,21,50), etime=dt.datetime(2015,5,5,22,51), fs=None, bs=None):
    import os
    import utils
    from scipy.signal import medfilt2d
    fname = "data/op/2015.05.05.22.11/waccmx/sd_bks_data.csv.gz"
    os.system("gzip -d " + fname)
    df = pd.read_csv(fname.replace(".gz", ""), parse_dates=["time"])
    os.system("gzip " + fname.replace(".gz", ""))
    X, Y, Z = utils.get_gridded_parameters(df)
    fig, ax = plt.subplots(figsize=(5,2.5), nrows=1, ncols=1, sharey=True, sharex=False, dpi=150)
    fmt = matplotlib.dates.DateFormatter("%H:%M")
    ax.xaxis.set_major_formatter(fmt)
    ax.text(0.75, 1.05, r"%s| $f_o$=%.1f MHz"%(rad.upper(), np.median(df.tfreq)/1e3), horizontalalignment="center",
            verticalalignment="center", transform=ax.transAxes)
    if fs is not None: ax.axvline(fs, ls="--", color="b", lw=0.8)
    if bs is not None: ax.axvline(bs, ls="--", color="r", lw=0.8)
    ax.set_xlim(stime, etime)
    ax.set_ylim(0, 70)
    ax.set_xlabel("Time, UT")
    ax.set_ylabel("Range Gate")
    from matplotlib import colors as mpl_colors
    import matplotlib.cm as cm
    cmj = cm.get_cmap("jet")
    cmpr = matplotlib.cm.prism
    cmap = matplotlib.colors.ListedColormap([cmpr(.142), cmpr(.125),
        cmpr(.11), cmpr(.1),
        ".6", cmpr(.175),
        cmpr(.158), cmj(.32),
        cmj(.37)])
    bounds = np.round(np.linspace(-100, 100, 7))
    bounds[3] = -15.
    bounds = np.insert(bounds,4,15.)
    norm = mpl_colors.BoundaryNorm(bounds, cmap.N) 
    cmap.set_bad("w", alpha=0.0)
    pcoll = ax.pcolormesh(X.T, Y.T, utils.medfilt2D_weight(Z), lw=0.01, edgecolors="None", cmap=cmap, norm=norm)
    cb = fig.colorbar(pcoll,shrink=0.9)
    cb.set_label(r"Velocity, $ms^{-1}$")
    fig.autofmt_xdate()
    fig.savefig("data/rti.example.png", bbox_inches="tight")
    return


def plot_supermag():
    import os
    from scipy import signal
    import matplotlib.dates as mdates
    fname = "data/op/2015.05.05.22.11/waccmx/sd_bks_data.csv.gz"
    os.system("gzip -d " + fname)
    df = pd.read_csv(fname.replace(".gz", ""), parse_dates=["time"])
    os.system("gzip " + fname.replace(".gz", ""))
    fs=dt.datetime(2015,5,5,22,7,40)
    bs=dt.datetime(2015,5,5,22,9,30)
    df = df[(df.time>=fs) & (df.time<=bs)][["time", "v"]]
    df = df.groupby("time").v.agg(["median", "std"]).reset_index()
    dx = pd.read_csv("data/supermag.csv", parse_dates=["Date_UTC"])
    dx = dx[["Date_UTC", "IAGA", "GEOLON", "GEOLAT", "MAGON", "MAGLAT", "MLT", "MCOLAT", "SZA", "dbn_nez", "dbe_nez", "dbz_nez"]]
    stns = ["BOU", "TUC", "T25", "BSL", "M06", "VIC", "BSL", "M08", "C08",]
    fig, axes = plt.subplots(figsize=(6,12), nrows=3, ncols=3, sharey=True, sharex=True, dpi=150)
    print(set(dx.IAGA))
    j = 0
    for stn in set(dx.IAGA):
        if stn in stns:
            if len(stns) > 1: ax = axes[j]
            else: ax = axes
            du = dx[(dx.IAGA==stn)]
            keys, col, labs = ["dbn_nez", "dbe_nez", "dbz_nez"], ["r", "b", "k"], ["dN", "dE", "dZ"]
            for i, k in enumerate(keys):
                x = du[[k, "Date_UTC"]].set_index("Date_UTC").resample("2s").interpolate("linear").reset_index()
                y = signal.resample(np.array(df["median"]), len(x))
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%H-%M"))
                ax.plot(x.Date_UTC, x[k], col[i], ls="--", lw=0.8, label=labs[i])
                ax.legend(loc=4)
                ax.set_xlim(x.Date_UTC.tolist()[0], x.Date_UTC.tolist()[-1])
                ax.set_ylabel("nT")
                ax.set_xlabel("Time, UT")
                ax.set_title("Mag: %s (%.2f, %.2f, %.2f)"%(stn, 
                    du["GEOLON"].tolist()[0], du["GEOLAT"].tolist()[0], du["SZA"].tolist()[0]))
            j += 1
    fig.autofmt_xdate()
    fig.savefig("data/corr.example.png", bbox_inches="tight")
    return
