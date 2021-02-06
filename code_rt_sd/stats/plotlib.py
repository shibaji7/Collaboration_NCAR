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
import matplotlib.pyplot as plt
import numpy as np

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
    fig.savefig("op/hist_reg.png", bbox_inches="tight")
    return

def plot_distribution(vn, vf):
    from scipy import stats
    fig = plt.figure(figsize=(4,4), dpi=120)
    ax = fig.add_subplot(111)
    ax.hist(vn, bins=np.arange(0,1,.01), color="r", alpha=0.5, density=True, label=r"$\frac{V_{d\eta}}{V_T}$", histtype="step")
    ax.hist(vf, bins=np.arange(0,1,.01), color="b", alpha=0.5, density=True, label=r"$\frac{V_{dh}}{V_T}$", histtype="step")
    ax.set_xlim(0,1)
    ax.legend(loc=1)
    ax.set_ylabel(r"Density $(\frac{V_x}{V_T})$")
    ax.set_xlabel(r"$\frac{V_x}{V_T}$")
    fig.savefig("op/hist.png", bbox_inches="tight")
    return

def plot_htstogram(vd, ve, vf, vn, vh, type):
    from scipy.stats import beta
    fig = plt.figure(figsize=(6,3), dpi=120)
    ax = fig.add_subplot(121)
    bins = np.arange(0,1,.05)
    #x = np.arange(0,1,0.001)
    #a, b, _, _ = beta.fit(vn,floc=0,fscale=1)
    ax.hist(vn, bins=bins, color="r", alpha=0.5, density=True, label=r"$\frac{V_{d\eta}}{V_T}[\mu=%.2f]$"%np.mean(vn)
            , histtype="step")
    #ax.plot(x, beta.pdf(x, a, b), color="r", lw=0.8, label=r"$\frac{V_{d\eta}}{V_T}[\mu=%.2f]$"%(a/(a+b)))
    ax.axvline(np.mean(vn), ls="--", color="r", lw=0.6)
    #a, b, _, _ = beta.fit(vh,floc=0,fscale=1)
    #ax.plot(x, beta.pdf(x, a, b), color="b", lw=0.8, label=r"$\frac{V_{dh}}{V_T}[\mu=%.2f]$"%(a/(a+b)))
    ax.axvline(np.mean(vh), ls="--", color="b", lw=0.6)
    ax.hist(vh, bins=bins, color="b", alpha=0.5, density=True, label=r"$\frac{V_{dh}}{V_T}[\mu=%.2f]$"%np.mean(vh),
            histtype="step")
    ax.text(0.1,0.9, "(a)", horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
    ax.set_xlim(0,1)
    ax.set_ylim(0,15)
    ax.legend(loc=1, prop={"size": 8})
    ax.set_ylabel(r"Density $\left(\frac{V_x}{V_T}\right)$")
    ax.set_xlabel(r"$\frac{V_x}{V_T}$")
    ax = fig.add_subplot(122)
    #a, b, _, _ = beta.fit(vd,floc=0,fscale=1)
    #ax.plot(x, beta.pdf(x, a, b), color="r", lw=0.8, label=r"$\frac{V_D}{V_T}[\mu=%.2f]$"%(a/(a+b)))
    ax.axvline(np.mean(vd), ls="--", color="r", lw=0.6)
    ax.hist(vd, bins=bins, color="r", alpha=0.5, density=True, label=r"$\frac{V_D}{V_T}[\mu=%.2f]$"%np.mean(vd),
            histtype="step")
    #a, b, _, _ = beta.fit(ve,floc=0,fscale=1)
    #ax.plot(x, beta.pdf(x, a, b), color="g", lw=0.8, label=r"$\frac{V_E}{V_T}[\mu=%.2f]$"%(a/(a+b)))
    ax.axvline(np.mean(ve), ls="--", color="g", lw=0.6)
    ax.hist(ve, bins=bins, color="g", alpha=0.5, density=True, label=r"$\frac{V_E}{V_T}[\mu=%.2f]$"%np.mean(ve),
            histtype="step")
    
    #a, b, _, _ = beta.fit(vf,floc=0,fscale=1)
    #ax.plot(x, beta.pdf(x, a, b), color="b", lw=0.8, label=r"$\frac{V_F}{V_T}[\mu=%.2f]$"%(a/(a+b)))
    ax.axvline(np.mean(vf), ls="--", color="b", lw=0.6)
    ax.hist(vf, bins=bins, color="b", alpha=0.5, density=True, label=r"$\frac{V_F}{V_T}[\mu=%.2f]$"%np.mean(vf),
            histtype="step")
    ax.set_ylim(0,15)
    #ax.hist(vd, bins=np.arange(0,1,.01), color="r", alpha=0.5, density=True, label=r"$\frac{v_D}{v_T}$", histtype="step")
    #ax.hist(ve, bins=np.arange(0,1,.01), color="b", alpha=0.5, density=True, label=r"$\frac{v_E}{v_T}$", histtype="step")
    #ax.hist(vf, bins=np.arange(0,1,.01), color="g", alpha=0.5, density=True, label=r"$\frac{v_F}{v_T}$", histtype="step")
    ax.text(0.1,0.9, "(b)", horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
    ax.set_xlim(0,1)
    ax.legend(loc=1, prop={"size": 8})
    ax.set_xlabel(r"$\frac{V_x}{V_T}$")
    fig.savefig("op/hist_{kind}.png".format(kind=type), bbox_inches="tight")
    return