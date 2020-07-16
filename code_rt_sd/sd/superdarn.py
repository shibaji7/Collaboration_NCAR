#!/usr/bin/env python

"""superdarn.py: module is dedicated to sd proc study."""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import os
import datetime as dt
import numpy as np
import pandas as pd
from geopy.distance import great_circle as GC
from scipy.interpolate import interp2d
from scipy.io import savemat, loadmat
import glob
from scipy.integrate import trapz
from scipy.stats import norm

from get_sd_data import FetchData
import utils
import plotlib


class SuperDARN(object):
    """ SuperDARN Model Estimate """

    def __init__(self, args):
        """ Initialze the parameters """
        for k in vars(args).keys():
            setattr(self, k, vars(args)[k])
        self.Re = 6371.
        utils.create_folder_structures(self.event, self.rad)
        self._estimate_bearing_()
        if hasattr(self, "save_radar") and self.save_radar: self._fetch_sd_()
        return
    
    def _estimate_bearing_(self):
        """ Estimate laitude and logitude bearings """
        fname = "data/sim/{dn}/{rad}/bearing.mat".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        m = {}
        lat, lon, bearing = utils.get_sd_radar(self.rad)
        p = (lat, lon)
        gc = GC(p, p)
        dist = np.linspace(0,self.mrange,self.nmrange)
        lats, lons = [], []
        for d in dist:
            x = gc.destination(p, bearing, distance=d)
            lats.append(x[0])
            lons.append(x[1])
        rinc = dist[1]-dist[0]
        m["dist"], m["lat"], m["lon"] = dist, np.array(lats), np.array(lons)
        m["olat"], m["olon"], m["rb"], m["num_range"], m["max_range"], m["range_inc"] = lat, lon, bearing, float(self.nmrange),\
                float(self.mrange), float(rinc)
        m["start_height"], m["height_inc"], m["num_heights"] = float(self.sheight), float(self.hinc),\
                float(len(np.arange(self.sheight,self.eheight,self.hinc)))
        m["ht"] = np.arange(self.sheight,self.eheight,self.hinc)
        m["freq"], m["tol"], m["nhops"] = float(self.frequency), float(1e-7), float(self.nhops)
        m["elev_s"], m["elev_i"], m["elev_e"] = float(self.selev), float(self.ielev), float(self.eelev)
        m["radius_earth"] = 6371.0
        m["d_ratio"], m["d_start"], m["d_end"], m["d_rtime"] = float(self.d_ratio), float(self.d_start),\
                float(self.d_end), float(self.d_rtime)
        m["f_ratio"], m["f_start"], m["f_end"], m["f_rtime"] = float(self.f_ratio), float(self.f_start),\
                float(self.f_end), float(self.f_rtime)
        m["e_ratio"] = float(self.e_ratio)
        savemat(fname, m)
        self.m, self.lat, self.lon, self.ht = m, m["lat"], m["lon"], m["ht"]
        return

    def _compute_doppler_(self):
        """ Compute doppler velocity """
        kconst, cconst, delt = 80.6, 3e8, self.d_rtime*60
        dic = "data/sim/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        grange, height = self.m["dist"], self.ht
        fis = glob.glob(dic + "exp.flare.bm({bm}).elv(*).csv".format(bm=self.bmnum))
        fis.sort()
        fni = "data/sim/{dn}/{rad}/exp.flare.bm({bm}).ne.mat".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), 
                rad=self.rad, bm=self.bmnum)
        fnj = "data/sim/{dn}/{rad}/exp.bgc.bm({bm}).ne.mat".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
                rad=self.rad, bm=self.bmnum)
        nei, nej = loadmat(fni)["iono_en_grid"]*1e6, loadmat(fnj)["iono_en_grid"]*1e6
        funci = interp2d(grange,height,np.log10(nei))
        funcj = interp2d(grange,height,np.log10(nej))
        for fi in fis:
            dop, dth = [], []
            elvi, rayi = float(fi.split(".")[-2][4:-1]), pd.read_csv(fi)
            gi, hi = np.array(rayi.grange), np.array(rayi.height)
            dth.append(0.)
            for k in range(len(gi[:-1])):
                dth.append(np.abs(hi[k]-hi[k+1]))
            dth = np.array(dth)*1000.
            for k, g, h in zip(range(len(gi)), gi, hi):
                if h > 50.:
                    dne = (10**funci(g,h) - 10**funcj(g, h))[0]
                    df = (kconst / (cconst * self.frequency*1e6)) * (dne / delt) * (dth[k] / np.cos(np.deg2rad(90.-elvi)))
                    if np.isnan(df):  df = 0.
                    dop.append(df)
                else: dop.append(0.)
            rayi["dop"] = dop
            rayi.to_csv(fi, header=True, index=False)
        return

    def _compute_velocity_(self):
        """ Compute Velocity """
        def _estimate_dop_delh_(x, y, phi=0):
            dh = np.abs(np.max(x.height) - np.max(y.height)) * 1000.
            xf = (2*self.frequency*1e6/3e8) * (dh/(self.f_rtime*60.)) * np.cos(np.deg2rad(phi))
            xd = 0.5 * xf * 3e8 / (self.frequency * 1e6)
            return xd
        dic = "data/sim/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        fis = glob.glob(dic + "exp.flare.bm({bm}).elv(*).csv".format(bm=self.bmnum))
        fjs = glob.glob(dic + "exp.bgc.bm({bm}).elv(*).csv".format(bm=self.bmnum))
        fis.sort()
        fjs.sort()
        vd, vf, itr = 0., 0., 0
        for fi, fj in zip(fis, fjs):
            elvi = float(fi.split(".")[-2][4:-1])
            d = pd.read_csv(fi)
            b = pd.read_csv(fj)
            f = np.abs(trapz(np.array(d.dop), np.array(d.grange)/np.cos(np.deg2rad(90-elvi))))
            vd += (0.5 * f * 3e8 / (self.frequency * 1e6)) 
            vf += _estimate_dop_delh_(d,b)
            itr += 1
        vd = vd / itr
        vf = vf / itr
        if self.verbose:
            print("\tDoppler velocity at {ev} -> Vd={vd} m/s, Vf={vf} m/s, Vt={vt} m/s".format(ev=
            self.event.strftime("%Y-%m-%d %H:%M"), vd=np.round(vd,1), vf=np.round(vf,1), vt=np.round(vf+vd,1)))
            print("\tParameter used - ")
            print("\t D Region e-density enhancement factor - ", self.d_ratio)
            print("\t D Region e-density rise time - ", self.d_rtime*60, " sec")
            print("\t F Region e-density enhancement factor - ", self.f_ratio)
            print("\t F Region e-density rise time - ", self.f_rtime*60, " sec")
        return vd, vf, vd+vf

    def _compute_(self):
        """ Compute RT using Pharlap """
        dic = "data/sim/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        fbgc = "data/sim/{dn}/{rad}/exp.bgc.bm({bm}).elv(<elv>).csv".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), 
                rad=self.rad, bm=self.bmnum)
        fflare = "data/sim/{dn}/{rad}/exp.flare.bm({bm}).elv(<elv>).csv".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
                rad=self.rad, bm=self.bmnum)
        cmd = "export DIR_MODELS_REF_DAT=/home/shibaji/Collaboration_NCAR/code_rt_sd/pharlap/pharlap_4.1.3/dat;\
                cd pharlap;\
                matlab -nodisplay -nodesktop -nosplash -nojvm -r \"UT=[{ut}];rad='{rad}';dic='{dic}';fbgc='{fbgc}';bm={bm};\
                fflare='{fflare}';rt_1D_sim;exit;\"".format(ut=self.event.strftime("%Y %m %d %H %S"), rad=self.rad,
                        dic=dic, bm=self.bmnum, fbgc=fbgc, fflare=fflare)
        os.system(cmd)
        return

    def _fetch_sd_(self):
        """ Fetch SuperDARN data and save to local """
        fname = "data/sim/{dn}/{rad}/sd_data.csv".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        if not os.path.exists(fname + ".gz"):
            fd = FetchData(self.rad, [self.event - dt.timedelta(minutes=8), self.event + dt.timedelta(minutes=8)])
            beams, _ = fd.fetch_data(v_params=["elv", "v", "w_l", "gflg", "p_l", "slist", "v_e"])
            self.rec = fd.convert_to_pandas(beams)
            self.rec.to_csv(fname, header=True, index=False)
            os.system("gzip " + fname)
        else:
            os.system("gzip -d " + fname+".gz")
            self.rec = pd.read_csv(fname, parse_dates=["time"])
            os.system("gzip " + fname)
        return

    def _exe_(self):
        """ Execute the RT model and save results"""
        print("\n Start simulation (using Pharlap) ...")
        dic = "data/sim/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        self._compute_()
        plotlib.plot_exp_rays(dic, self.event, self.bmnum, "bgc")
        plotlib.plot_exp_rays(dic, self.event, self.bmnum, "flare")
        if self.verbose: print("\n Processing Doppler.")
        self._compute_doppler_()
        vd, vf, vt = self._compute_velocity_()
        return vd, vf, vt
