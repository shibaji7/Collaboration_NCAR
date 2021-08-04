import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np
import pydarn
import rad_fov

from shapely.geometry import  MultiLineString, mapping, LineString, Polygon
from descartes.patch import PolygonPatch

def overlay_radar(rads, ax, _to, _from, color="k", zorder=2, marker="o", ms=2, font={"size":7}, north=True):
    times = -1 if north else 1    
    for rad in rads:
        hdw = pydarn.read_hdw_file(rad)
        lat, lon, _ = hdw.geographic
        tlat, tlon = lat+3*times, lon+3*times*np.sign(hdw.boresight*-1)
        x, y = _to.transform_point(lon, lat, _from)
        tx, ty = _to.transform_point(tlon, tlat, _from)
        ax.plot(x, y, color=color, zorder=zorder, marker=marker, ms=ms)
        ax.text(tx, ty, rad.upper(), ha="center", va="center", fontdict=font)
    return

def overlay_fov(rads, ax, _to, _from, maxGate=40, rangeLimits=None, beamLimits=None,
            model="IS", fov_dir="front", fovColor=None, fovAlpha=0.2,
            fovObj=None, zorder=2, lineColor="k", lineWidth=0.5, ls="-"):
    """ Overlay radar FoV """
    from numpy import transpose, ones, concatenate, vstack, shape
    for rad in rads:
        hdw = pydarn.read_hdw_file(rad)
        sgate = 0
        egate = hdw.gates if not maxGate else maxGate
        ebeam = hdw.beams
        if beamLimits is not None: sbeam, ebeam = beamLimits[0], beamLimits[1]
        else: sbeam = 0
        rfov = rad_fov.CalcFov(hdw=hdw, ngates=egate)
        x, y = np.zeros_like(rfov.lonFull), np.zeros_like(rfov.latFull)
        for _i in range(rfov.lonFull.shape[0]):
            for _j in range(rfov.lonFull.shape[1]):
                x[_i, _j], y[_i, _j] = _to.transform_point(rfov.lonFull[_i, _j], rfov.latFull[_i, _j], _from)
        contour_x = concatenate((x[sbeam, sgate:egate], x[sbeam:ebeam, egate],
                                 x[ebeam, egate:sgate:-1],
                                 x[ebeam:sbeam:-1, sgate]))
        contour_y = concatenate((y[sbeam, sgate:egate], y[sbeam:ebeam, egate],
                                 y[ebeam, egate:sgate:-1],
                                 y[ebeam:sbeam:-1, sgate]))
        ax.plot(contour_x, contour_y, color=lineColor, zorder=zorder, linewidth=lineWidth, ls=ls)
        if fovColor:
            contour = transpose(vstack((contour_x, contour_y)))
            polygon = Polygon(contour)
            patch = PolygonPatch(polygon, facecolor=fovColor, edgecolor=fovColor, alpha=fovAlpha, zorder=zorder)
            ax.add_patch(patch)
    return

def create_cartopy(prj):
    fig = plt.figure(dpi=120, figsize=(6,6))
    ax = fig.add_subplot(111, projection=prj)
    ax.add_feature(cartopy.feature.OCEAN, zorder=0, alpha=0.1)
    ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor="black", alpha=0.2, lw=0.3)
    ax.set_global()
    gl = ax.gridlines(linewidth=0.3, draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_left = False
    return fig, ax

def convert_to_map_lat_lon(xs, ys, _from, _to):
    lat, lon = [], []
    for x, y in zip(xs, ys):
        _lon, _lat = _to.transform_point(x, y, _from)
        lat.append(_lat)
        lon.append(_lon)
    return lat, lon

def create_plots(rads=["cvw","bks","wal","cve","fhe","fhw","kap","sas","pgr"], title="", pngfname="fov.png"):
    geodetic = ccrs.Geodetic()
    #orthographic = ccrs.Orthographic(-120,90)
    orthographic = ccrs.PlateCarree(central_longitude=-95)
    fig, ax = create_cartopy(orthographic)
    for rad, c in zip(rads, 6*["r"]+3*["b"]):
        overlay_fov([rad], ax, orthographic, geodetic, fovColor=c)
        overlay_radar([rad], ax, orthographic, geodetic, north=True)
    ax.set_title(title)
    ax.set_extent([-150, -50, 20, 70])
    ax.text(-0.02, 0.65, "Geographic Coordinates", horizontalalignment="center",
            verticalalignment="center", transform=ax.transAxes, rotation=90)
    fig.savefig(pngfname, bbox_inches="tight")
    return

if __name__ == "__main__":
    #create_dec4_eclipse()
    #create_june10_eclipse()
    create_plots()
    import os
    os.system("rm -rf *.log")
    os.system("rm -rf __pycache__")
    pass
