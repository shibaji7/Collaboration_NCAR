#!/usr/bin/env python

"""fov.py: Fields-of-View plots showing all the radar fields of view """

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
import datetime as dt


fonttext = {"family": "serif", "color":  "darkblue", "weight": "normal", "size": 12}
fonttextx = {"family": "serif", "color":  "b", "weight": "normal", "size": 12}
fontT = {"family": "serif", "color":  "k", "weight": "normal", "size": 6}
font = {"family": "serif", "color":  "black", "weight": "normal", "size": 8}
fontx = {"family": "serif", "weight": "normal", "size": 9}
fontL = {"family": "serif", "color":  "darkblue", "weight": "bold", "size": 8}

from matplotlib import font_manager
ticks_font = font_manager.FontProperties(family="serif", size=10, weight="normal")
matplotlib.rcParams["xtick.color"] = "k"
matplotlib.rcParams["ytick.color"] = "k"
matplotlib.rcParams["xtick.labelsize"] = 8
matplotlib.rcParams["ytick.labelsize"] = 8
matplotlib.rcParams["mathtext.default"] = "default"


fonttext["size"] = 10
maxgate=75
plot_date = dt.datetime(2015, 5, 5)
import sys
sys.path.append("sd_cartopy/")
from fov import *


fig = plt.figure()
ax = fig.add_subplot(projection="fovcarto",\
        coords="geo", plot_date=plot_date)
ax.coastlines()
ax.overlay_radar()
ax.overlay_fov(fovColor="r")
ax.grid_on()
ax.enum()
fig.savefig("data/fov.png", bbox_inches="tight")

import os
os.system("rm -rf *.log")
os.system("rm -rf __pycache__/")
os.system("rm -rf sd_cartopy/__pycache__/")
