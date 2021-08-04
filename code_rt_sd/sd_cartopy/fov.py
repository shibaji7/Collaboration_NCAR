import matplotlib
matplotlib.use("Agg")

import pydarnio
import pydarn
import numpy as np

import cartopy
from cartopy.mpl.geoaxes import GeoAxes
import aacgmv2
import numpy
from shapely.geometry import  MultiLineString, mapping, LineString, Polygon
from descartes.patch import PolygonPatch
from matplotlib.projections import register_projection
import copy
import datetime as dt
import rad_fov

class FoVCarto(GeoAxes):
    name = "fovcarto"

    def __init__(self, *args, **kwargs):
        if "map_projection" in kwargs: self.map_projection = kwargs.pop("map_projection")
        else: self.map_projection = cartopy.crs.PlateCarree(central_longitude=0)
        if "rad" in kwargs: self.rad = kwargs.pop("rad")
        else: self.rad = "bks"
        if "plot_date" in kwargs: self.plot_date = kwargs.pop("plot_date")
        else: self.plot_date = None
        self.supported_coords = [ "geo" ]
        if "coords" in kwargs:
            self.coords = kwargs.pop("coords")
            if self.coords not in self.supported_coords:
                err_str = "coordinates not supported, choose from : "
                for _n,_sc in enumerate(self.supported_coords):
                    if _n + 1 != len(self.supported_coords): err_str += _sc + ", "
                    else: err_str += _sc
                raise TypeError(err_str)
        else: self.coords = "geo"
        super().__init__(map_projection=self.map_projection,*args, **kwargs)
        return

    def overaly_coast_lakes(self, resolution="50m", color="black", **kwargs):
        """  Overlay AACGM coastlines and lakes  """
        kwargs["edgecolor"] = color
        kwargs["facecolor"] = "none"
        # overaly coastlines
        feature = cartopy.feature.NaturalEarthFeature("physical", "coastline",
                resolution, **kwargs)
        self.add_feature( cartopy.feature.COASTLINE, **kwargs )
        self.add_feature( cartopy.feature.LAKES, **kwargs )
        return

    def coastlines(self,resolution="50m", color="black", **kwargs):
        # details!
        kwargs["edgecolor"] = color
        kwargs["linewidth"] = 0.8
        kwargs["facecolor"] = "none"
        feature = cartopy.feature.NaturalEarthFeature("physical", "coastline",
                resolution, **kwargs)
        return self.add_feature(feature, **kwargs)
    
    def add_feature(self, feature, **kwargs):
        if "edgecolor" not in kwargs: kwargs["edgecolor"] = "black"
        kwargs["facecolor"] = "none"
        super().add_feature(feature, **kwargs)
        return

    def grid_on(self, tx=cartopy.crs.PlateCarree(), draw_labels=True, linewidth=0.5, color="gray", alpha=0.5, linestyle="--"):
        """ Adding grids to map """
        import matplotlib.ticker as mticker
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
        gl = self.gridlines(crs=tx, draw_labels=draw_labels,
                linewidth=linewidth, color=color, alpha=alpha, linestyle=linestyle)
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 15))
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        return

    def overlay_radar(self, tx=cartopy.crs.PlateCarree(), zorder=2, markerColor="darkblue", markerSize=15, fontSize=10, font_color="k",
            xOffset=None, yOffset=-0.1, annotate=True):
        """ Adding the radar location """
        self.hdw = pydarn.read_hdw_file(self.rad)
        print(" Radar:", self.hdw.geographic.lon, self.hdw.geographic.lat)
        self.scatter([self.hdw.geographic.lon], [self.hdw.geographic.lat], s=markerSize, marker="D",
                color=markerColor, zorder=zorder, transform=tx)
        nearby_rad = [["adw", "kod", "cve", "fhe", "wal", "gbr", "pyk", "aze", "sys"],
                            ["ade", "ksr", "cvw", "fhw", "bks", "sch", "sto", "azw", "sye"]]
        if annotate:
            if self.rad in nearby_rad[0]: xOff, ha = 0.1 if not xOffset else xOffset, 0
            elif self.rad in nearby_rad[1]: xOff, ha = -0.1 if not xOffset else xOffset, 1
            else: xOff, ha = 0.0, 0.5
            lon, lat = self.hdw.geographic.lon, self.hdw.geographic.lat
            x, y = self.projection.transform_point(lon, lat, src_crs=tx)
        return

    def overlay_fov(self, tx=cartopy.crs.PlateCarree(), maxGate=75, rangeLimits=None, beamLimits=None,
            model="IS", fov_dir="front", fovColor=None, fovAlpha=0.2,
            fovObj=None, zorder=2, lineColor="k", lineWidth=1, ls="-"):
        """ Overlay radar FoV """
        lcolor = lineColor
        from numpy import transpose, ones, concatenate, vstack, shape
        self.hdw = pydarn.read_hdw_file(self.rad)
        sgate = 0
        egate = self.hdw.gates if not maxGate else maxGate
        ebeam = self.hdw.beams
        if beamLimits is not None: sbeam, ebeam = beamLimits[0], beamLimits[1]
        else: sbeam = 0
        self.rad_fov = rad_fov.CalcFov(hdw=self.hdw, ngates=egate)
        xyz = self.projection.transform_points(tx, self.rad_fov.lonFull, self.rad_fov.latFull)
        x, y = xyz[:, :, 0], xyz[:, :, 1]
        contour_x = concatenate((x[sbeam, sgate:egate], x[sbeam:ebeam, egate],
            x[ebeam, egate:sgate:-1],
            x[ebeam:sbeam:-1, sgate]))
        contour_y = concatenate((y[sbeam, sgate:egate], y[sbeam:ebeam, egate],
            y[ebeam, egate:sgate:-1],
            y[ebeam:sbeam:-1, sgate]))
        self.plot(contour_x, contour_y, color=lcolor, zorder=zorder, linewidth=lineWidth, ls=ls)
        if fovColor:
            contour = transpose(vstack((contour_x, contour_y)))
            polygon = Polygon(contour)
            patch = PolygonPatch(polygon, facecolor=fovColor, edgecolor=fovColor, alpha=fovAlpha, zorder=zorder)
            self.add_patch(patch)
        return

    def enum(self, bounds=[-120, -70, 25, 70]):
        if bounds is not None: self.set_extent(bounds)
        if self.coords == "geo": self.text(1.02, 0.75, "Geographic Coordinates", horizontalalignment="center",
                verticalalignment="center", transform=self.transAxes, rotation=90)
        self.text(0.1, 1.02, "Rad: "+self.rad, horizontalalignment="center",
                verticalalignment="center", transform=self.transAxes)
        if self.plot_date: self.text(0.75, 1.02, self.plot_date.strftime("%Y-%m-%d %H:%M") + " UT", horizontalalignment="center",
                verticalalignment="center", transform=self.transAxes)
        return

    def add_dn_terminator(self, **kwargs):
        """ Adding day night terminator """
        from cartopy.feature.nightshade import Nightshade
        if self.plot_date: 
            ns_feature = Nightshade(self.plot_date, alpha=0.2)
            super().add_feature(feature, **kwargs)
        return

    def overlay_radar_data(self, dat, p_name = "v", tx=cartopy.crs.PlateCarree(), 
            p_max=100, p_min=-100, p_step=10, p_ub=9999, p_lb=-9999,
            cmap=matplotlib.pyplot.get_cmap("Spectral"), add_colorbar=True, colorbar_label="Velocity [m/s]",
            **kwargs):
        """ 
        Adding radar data
        dat: dict()-> with "bmnum", "slist" and "v" or other parameter in list of list format
        """
        nbeam = np.max(np.max(dat["bmnum"])) + 1
        if self.maxGate: nrange = self.maxGate
        else: nrange = np.max(np.max(dat["slist"])) + 1
        hdw = pydarn.read_hdw_file(self.rad)
        rf = rad_fov.CalcFov(hdw=hdw, ngates=nrange, nbeams=nbeam)
        lons, lats = rf.lonFull, rf.latFull
        
        p_ranges = list(range(p_min, p_max + 1, p_step))
        p_ranges.insert(0, p_lb)
        p_ranges.append(p_ub)        
        
        Px = np.zeros((nbeam, nrange))*np.nan
        idbs, idgs = dat["bmnum"], dat["slist"]
        params = dat[p_name]
        for idb, idg, par in zip(idbs, idgs, params):
            idb = np.array(idb)[np.array(idg) < nrange]
            par = np.array(par)[np.array(idg) < nrange]
            idg = np.array(idg)[np.array(idg) < nrange]
            if len(par) > 0: Px[idb, np.round(idg).astype(int)] = par
        Px = np.ma.masked_invalid(Px)
        self.pcolormesh(lons, lats, Px, transform=tx, cmap=cmap,
                vmax=p_max, vmin=p_min, **kwargs)
        if add_colorbar: self._add_colorbar(p_ranges, cmap, label=colorbar_label)
        return
    
    def _add_colorbar(self, bounds, colormap, label="", mark=2):
        """ Add a colorbar to the right of an axis. """
        pos = self.get_position()
        cpos = [pos.x1 + 0.035, pos.y0 + 0.25*pos.height,
                0.01, pos.height * 0.5]            # this list defines (left, bottom, width, height)
        cax = matplotlib.pyplot.gcf().add_axes(cpos)
        norm = matplotlib.colors.BoundaryNorm(bounds[::mark], colormap.N)
        cb2 = matplotlib.colorbar.ColorbarBase(cax, cmap=colormap,
                norm=norm,
                ticks=bounds[::mark],
                spacing="uniform",
                orientation="vertical")
        cb2.set_label(label)
        # Remove the outer bounds in tick labels
        ticks = [str(i) for i in bounds[::mark]]
        ticks[0], ticks[-1] = "", ""
        cb2.ax.set_yticklabels(ticks)
        return
    
    def overlay_data(self, lats, lons, dat, tx=cartopy.crs.PlateCarree(), p_max=30, p_min=-30, p_step=5, p_ub=200, p_lb=-200,
            cmap=matplotlib.pyplot.get_cmap("Spectral"), add_colorbar=True, colorbar_label="Velocity [m/s]",
            **kwargs):
        """ Adding data """
        p_ranges = list(np.arange(p_min, p_max + 1, p_step))
        p_ranges.insert(0, p_lb)
        p_ranges.append(p_ub)
        dat = np.ma.masked_invalid(dat)
        self.pcolormesh(lons, lats, dat, transform=tx, cmap=cmap,
                vmax=p_max, vmin=p_min, **kwargs)
        if add_colorbar: self._add_colorbar(p_ranges, cmap, label=colorbar_label, mark=1)
        return

register_projection(FoVCarto)


def adding_plot(date, fname, coords="geo", rad="bks"):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(projection="fovcarto",\
            coords=coords, plot_date=date,  rad=rad)
    ax.coastlines()
    ax.overlay_radar()
    ax.overlay_fov()
    ax.grid_on()
    ax.enum()
    fig.savefig(fname, bbox_inches="tight")
    return ax, fig

def adding_plots(dates, fname, coords="geo", nrows=2, ncols=2):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(5*nrows, 5*ncols), dpi=150)
    share = None
    for i, date in enumerate(dates):
        idx = 100*nrows + 10*ncols + (1+i)
        if i == 0: 
            ax = fig.add_subplot(idx, projection="fovcarto",\
                    coords=coords, plot_date=date)
            share = ax
        else: ax = fig.add_subplot(idx, projection="fovcarto",\
                coords=coords, plot_date=date, sharex=share)
        ax.coastlines()
        ax.overlay_radar()
        ax.overlay_fov()
        ax.grid_on()
        ax.enum()
    fig.savefig(fname, bbox_inches="tight")
    return

def adding_radars(date, fname, coords="geo", rads=["bks", "wal", "fhw", "fhe"]):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(projection="fovcarto",\
            coords=coords, plot_date=date, map_projection=cartopy.crs.Orthographic(central_longitude=-95, central_latitude=90))
    ax.coastlines()
    for rad in rads:
        ax.rad = rad
        #ax.overlay_radar()
        #ax.overlay_fov()
    #ax.grid_on()
    #ax.enum(bounds=[-150, -40, 20, 90])
    fig.savefig(fname, bbox_inches="tight")
    return ax, fig

def get_globe(fig, num, date, coords="geo"):
    import matplotlib.pyplot as plt
    ax = fig.add_subplot(num, projection="fovcarto",\
            coords=coords, plot_date=date)
    ax.coastlines()
    return ax, fig, ax.map_projection, cartopy.crs.Geodetic()


if __name__ == "__main__":
    #plot_date = dt.datetime(2015,1,1,4)
    #adding_plot(plot_date, "carto_test.png", coords="geo")
    #adding_plots([dt.datetime(2015,5,5,22,0), dt.datetime(2015,5,5,22,0),
    #    dt.datetime(2015,5,5,22,0), dt.datetime(2015,5,5,22,0)], "carto_test.png", coords="geo", nrows=2, ncols=2)
    adding_radars(dt.datetime(2015,5,5,22,0), "fov.png", coords="geo")
    import os
    os.system("rm -rf *.log")
    os.system("rm -rf __pycache__/")
