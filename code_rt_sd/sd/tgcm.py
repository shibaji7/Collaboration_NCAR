#!/usr/bin/env python

"""tgcm.py: module is dedicated to download and process the TIME-GCM data."""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import os
from paramiko import SSHClient
from scp import SCPClient
import datetime as dt
import numpy as np
import pandas as pd
from netCDF4 import Dataset, date2num, num2date
from scipy.interpolate import interp1d, interp2d
from geopy.distance import great_circle as GC
from scipy.io import savemat, loadmat
from scipy import array
import glob
from scipy.integrate import trapz

from get_sd_data import FetchData
import utils
import plotlib


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

class TGCM(object):
    """ TIME-GCM data download from Cheyenne"""

    def __init__(self, args):
        """ Initialze the parameters """
        for k in vars(args).keys():
            setattr(self, k, vars(args)[k])
        self.Re = 6371.
        if self.tsim_start is None: self.tsim_start = 0
        if self.tsim_end is None: self.tsim_end = int((self.end-self.start).total_seconds()/60.)
        utils.create_folder_structures(self.event, self.rad)
        self.rtime = utils.get_rtime(self.event, th = self.threshold)
        self.con = False
        if hasattr(self, "clear") and self.clear: os.system("rm data/sim/{dn}/{rad}/*".format(
            dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad))
        return

    def _conn_(self):
        """ Create Conn """
        self.ssh = SSHClient()
        self.ssh.load_system_host_keys()
        self.ssh.connect(hostname="cheyenne.ucar.edu", port = 22, username="shibaji", password="0513-shibaji-cit")
        self.scp = SCPClient(self.ssh.get_transport())
        self.con = True
        return

    def _close_(self):
        """ Close connection """
        if self.con:
            self.scp.close()
            self.ssh.close()
        return

    def _download_(self, i):
        """ Download the files from """
        dfname = "data/sim/{dn}/tgcm.{du}.nc.gz".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), 
                du=(self.start + dt.timedelta(minutes=i)).strftime("%Y.%m.%d.%H.%M"))
        dn = self.start + dt.timedelta(minutes=i)
        if not os.path.exists(dfname):
            if not self.con: self._conn_()
            if self.verbose: print("\n Create SSH.SCP connection, post process, and fetch data.")
            self.ssh.exec_command("mkdir tmp/")
            self.scp.put("config/proc.py", "tmp/")
            if self.verbose: print(" Run -", "python tmp/proc.py -ev {ev} -dn {dn} -rn {rn}".format(
                ev=self.event.strftime("%Y.%m.%d.%H.%M"),
                dn=dn.strftime("%Y.%m.%d.%H.%M"), rn=self.sim_id))
            stdin, stdout, stderr = self.ssh.exec_command("source ncar_tgcm_waccm_proc/ncar/bin/activate"\
                    "\n python tmp/proc.py -ev {ev} -dn {dn} -rn {rn}".format(ev=self.event.strftime("%Y.%m.%d.%H.%M"),
                        dn=dn.strftime("%Y.%m.%d.%H.%M"), rn=self.sim_id), get_pty=True)
            if self.verbose:
                for line in iter(stdout.readline, ""):
                    print(line, end="")
            if self.verbose: print(" End run ")
            self.ssh.exec_command("deactivate")
            self.scp.get("tmp/tgcm.nc.gz", dfname)
            self.ssh.exec_command("rm -rf tmp/")
        return


    def _interpolate_(self, i, case="d"):
        """ Interpolate electron density """
        m = {}
        fname = "data/sim/{dn}/{rad}/ne.ti({ti}).bm({bm}).{case}.mat".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), 
                rad=self.rad, ti=i, bm=self.bmnum, case=case)
        if not os.path.exists(fname):
            dfname = "data/sim/{dn}/tgcm.{du}.nc.gz".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
                    du=(self.start + dt.timedelta(minutes=i)).strftime("%Y.%m.%d.%H.%M"))
            os.system("gzip -d " + dfname)
            self._nc_ = Dataset(dfname.replace(".gz",""))
            os.system("gzip " + dfname.replace(".gz",""))
            nex, alts, latx, lonx = _intp_(self._nc_.variables["ZG"+case][0,:,:,:]*1e-5,
                    self._nc_.variables["lat"][:], self._nc_.variables["lon"][:], self._nc_.variables["NE"+case][0,:,:,:], 
                    hd=[self.sheight,self.eheight,self.hinc], 
                    dlat=2, dlon=4, scale="log", kind="cubic", v=self.verbose)
            ne = np.zeros((len(self.ht), len(self.lat)))
            for ix, la, lo in zip(range(len(self.lat)), self.lat, self.lon):
                _ij_ = _get_ij_(latx, lonx, la, lo)
                ne[:, ix] = nex[:, _ij_[0], _ij_[1]]
            m["ne"] = ne
            savemat(fname, m)
            self._compute_(i, case)
            if case == "f": self._plot_radstn_(i)
        return

    def _estimate_bearing_(self):
        """ Estimate laitude and logitude bearings """
        fname = "data/sim/{dn}/{rad}/bearing.mat".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        m = {}
        lat, lon, bearing = utils.get_sd_radar(self.rad)
        p = (lat, lon)
        self.rlat, self.rlon = lat, lon
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
        savemat(fname, m)
        self.m, self.lat, self.lon, self.ht = m, m["lat"], m["lon"], m["ht"]
        return

    def _compute_doppler_(self, i):
        """ Compute doppler velocity """
        kconst, cconst, delt = 80.6, 3e8, self.rtime*60
        dic = "data/sim/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        grange, height = self.m["dist"], self.ht
        fnd = glob.glob(dic + "ti({ti}).bm({bm}).elv(*).d.csv".format(ti=i, bm=self.bmnum))
        fnf = glob.glob(dic + "ti({ti}).bm({bm}).elv(*).f.csv".format(ti=i, bm=self.bmnum))
        fnd.sort()
        fnf.sort()
        ned = loadmat("data/sim/{dn}/{rad}/ne.ti({ti}).bm({bm}).d.mat".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
            rad=self.rad, ti=i, bm=self.bmnum))["ne"]*1e6
        nef = loadmat("data/sim/{dn}/{rad}/ne.ti({ti}).bm({bm}).f.mat".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
            rad=self.rad, ti=i, bm=self.bmnum))["ne"]*1e6
        funcf = interp2d(grange,height,np.log10(nef))
        funcd = interp2d(grange,height,np.log10(ned))
        for fi in fnf:
            dop, dth = [], []
            elvi, rayi = float(fi.split(".")[-3][4:-1]), pd.read_csv(fi)
            gi, hi = np.array(rayi.grange), np.array(rayi.height)
            dth.append(0.)
            for k in range(len(gi[:-1])):
                dth.append(np.abs(hi[k]-hi[k+1]))
            dth = np.array(dth)*1000.
            for k, g, h in zip(range(len(gi)), gi, hi):
                if h > 50.:
                    dne = (10**funcf(g,h) - 10**funcd(g, h))[0]
                    df = (kconst / (cconst * self.frequency*1e6)) * (dne / delt) * (dth[k] / np.cos(np.deg2rad(90.-elvi)))
                    if np.isnan(df):  df = 0.
                    dop.append(df)
                else: dop.append(0.)
            rayi["dop"] = dop
            rayi.to_csv(fi, header=True, index=False)
        self._compute_velocity_(i)
        return

    def _compute_velocity_(self, i):
        """ Compute Velocity """
        def _estimate_dop_delh_(x, y, phi=0):
            dh = np.abs(np.max(x.height) - np.max(y.height)) * 1000.
            xf = (2*self.frequency*1e6/3e8) * (dh/(self.rtime*60.)) * np.cos(np.deg2rad(phi))
            xd = 0.5 * xf * 3e8 / (self.frequency * 1e6)
            return xd
        dic = "data/sim/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        fnd = glob.glob(dic + "ti({ti}).bm({bm}).elv(*).d.csv".format(ti=i, bm=self.bmnum))
        fnf = glob.glob(dic + "ti({ti}).bm({bm}).elv(*).f.csv".format(ti=i, bm=self.bmnum))
        fnd.sort()
        fnf.sort()
        vd, vf, itr = 0., 0., 0
        for fi, fj in zip(fnf, fnd):
            elvi = float(fi.split(".")[-3][4:-1])
            d = pd.read_csv(fi)
            b = pd.read_csv(fj)
            f = np.abs(trapz(np.array(d.dop), np.array(d.grange)/np.cos(np.deg2rad(90-elvi))))
            vd += (0.5 * f * 3e8 / (self.frequency * 1e6))
            vf += _estimate_dop_delh_(d,b)
            itr += 1
        self.vd = vd / itr
        self.vf = vf / itr
        if self.verbose:
            print("\tDoppler velocity at {ev} -> Vd={vd} m/s, Vf={vf} m/s, Vt={vt} m/s".format(ev=
                (self.start+dt.timedelta(minutes=i)).strftime("%Y-%m-%d %H:%M"), 
                    vd=np.round(self.vd,1), vf=np.round(self.vf,1), vt=np.round(self.vf+self.vd,1)))
        return

    def _compute_(self, i, case="d"):
        """ Compute RT using Pharlap """
        u = self.start + dt.timedelta(minutes=i)
        dic = "data/sim/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        fname = "data/sim/{dn}/{rad}/ti({ti}).bm({bm}).elv(<elv>).{case}.csv".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), 
                rad=self.rad, bm=self.bmnum, ti=i, case=case)
        cmd = "export DIR_MODELS_REF_DAT=/home/shibaji/Collaboration_NCAR/code_rt_sd/pharlap/pharlap_4.1.3/dat;\
                cd pharlap;\
                matlab -nodisplay -nodesktop -nosplash -nojvm -r \"UT=[{ut}];rad='{rad}';dic='{dic}';fname='{fname}';bm={bm};\
                ti={ti};cse='{case}';rt_1D;exit;\"".format(ut=u.strftime("%Y %m %d %H %S"),
                        rad=self.rad, dic=dic, bm=self.bmnum, ti=i, fname=fname, case=case)
        os.system(cmd)
        return

    def _fetch_sd_(self):
        """ Fetch SuperDARN data and save to local """
        fname = "data/sim/{dn}/{rad}/sd_data.csv".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        if not os.path.exists(fname + ".gz"):
            fd = FetchData(self.rad, [self.start, self.end])
            beams, _ = fd.fetch_data(v_params=["elv", "v", "w_l", "gflg", "p_l", "slist", "v_e"])
            self.rec = fd.convert_to_pandas(beams)
            self.rec.to_csv(fname, header=True, index=False)
            os.system("gzip " + fname)
        else:
            os.system("gzip -d " + fname+".gz")
            self.rec = pd.read_csv(fname, parse_dates=["time"])
            os.system("gzip " + fname)
        return

    def _plot_radstn_(self, i):
        """ Plot radar station """
        fname = "data/sim/{dn}/{rad}/tgcm({i}).png".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad, i=i)
        ij = _get_ij_(self._nc_.variables["lat"][:], self._nc_.variables["lon"][:], self.rlat, self.rlon)
        p = self._nc_.variables["NEd"][0,:,ij[0],ij[1]]
        f = self._nc_.variables["NEf"][0,:,ij[0],ij[1]]
        pz = self._nc_.variables["ZGd"][0,:,ij[0],ij[1]]*1e-5
        fz = self._nc_.variables["ZGf"][0,:,ij[0],ij[1]]*1e-5
        plotlib.plot_radstn(p,f,pz,fz,fname,self.rlat,self.rlon,self.start+dt.timedelta(minutes=i))
        return

    def _exe_(self):
        """ Execute the RT model and save results"""
        print("\n Start simulation (using Pharlap) ...")
        dic = "data/sim/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        self._nc_ = None
        self._estimate_bearing_()
        if hasattr(self, "save_radar") and self.save_radar: self._fetch_sd_()
        for i in range(self.tsim_start, self.tsim_end):
            self._download_(i)
            if self.verbose: print("\tProcess-", self.start + dt.timedelta(minutes=i))
            for c in ["d", "f"]:
                self._interpolate_(i, c)
                txt = ""
                if c == "f":
                    self._compute_doppler_(i)
                    txt = r"$V_{d\eta}$=%.1f, $V_{dh}$=%.1f"%(self.vd,self.vf)
                plotlib.plot_rays(dic, self.start + dt.timedelta(minutes=i), i, self.bmnum, c, txt)
        if self.verbose: print("\n Interpolation completed.")
        if self.verbose: print("\n Processing Doppler.")
        self._close_()
        if hasattr(self, "archive") and self.archive: self._arch_()
        return
    
    def _arch_(self):
        """ Archive the data and create movie """
        print("\n System archiving ...")
        arc =  "data/sim/{dn}/archive/{rad}/tgcm/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        os.system("rm -rf " + arc)
        os.system("mkdir -p " + arc)
        os.system("cp -r data/sim/{dn}/{rad}/* {arc}".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
            rad=self.rad, arc=arc))
        cmd = "ffmpeg -r 1 -i {arc}/tgcm\(%d\).png -c:v libx264 -vf \
                'scale=1420:-2,fps=3,format=yuv420p' {arc}/tgcm_edens_rad.mp4".format(arc=arc)
        os.system(cmd)
        cmd = "ffmpeg -r 1 -i {arc}/rt.ti\(%d\).bm\({bm}\).d.png -c:v libx264 -vf \
                'scale=1420:-2,fps=3,format=yuv420p' {arc}/tgcm_rt_daily.mp4".format(arc=arc, bm=self.bmnum)
        os.system(cmd)
        cmd = "ffmpeg -r 1 -i {arc}/rt.ti\(%d\).bm\({bm}\).f.png -c:v libx264 -vf \
                'scale=1420:-2,fps=3,format=yuv420p' {arc}/tgcm_rt_flare.mp4".format(arc=arc, bm=self.bmnum)
        os.system(cmd)
        print(" System archived!")
        return
