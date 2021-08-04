#!/usr/bin/env python

"""models.py: module is dedicated to download and process the WACCM-X / TIME-GCM data."""

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
from scipy import signal

from get_sd_data import FetchData
import utils
import plotlib
import sys
sys.path.append("sd_cartopy/")
import rad_fov

INT_F = 300
INT_R = 300
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

class Model(object):
    """ General model data analysis for WACCM-X / TIME-GCM"""

    def __init__(self, args):
        """ Initialze the parameters """
        for k in vars(args).keys():
            setattr(self, k, vars(args)[k])
        self.Re = 6371.
        self.files = {
                "base": "data/op/{dn}/{model}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), model=self.model),
                "run": "data/op/{dn}/{model}/{rad}/bm.{bm}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
                    rad=self.rad, model=self.model, bm="%02d"%self.bmnum),
                "eden.nc": "{du}.nc",
                "eden.gz": "data/op/{dn}/{model}/{du}.nc.gz",
                "eden.mat": "ne_ti({ti})_{case}.mat",
                }
        if self.tsim_start is None: self.tsim_start = 0
        if self.tsim_end is None: self.tsim_end = int((self.end-self.start).total_seconds()/60.)
        if not os.path.exists(self.files["run"]): os.system("mkdir -p " + self.files["run"])
        if hasattr(self, "clear") and self.clear: self._clear_()
        self.recs = []
        return

    def _clear_(self):
        """ Clear past estimations """
        os.system("rm -rf " + self.files["run"])
        os.system("mkdir -p " + self.files["run"]) 
        return

    def _interpolate_(self, i, case="d"):
        """ Interpolate electron density """
        m = {}
        fname = (self.files["run"] + self.files["eden.mat"]).format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), 
                model=self.model, rad=self.rad, ti="%02d"%i, bm=self.bmnum, case=case)
        if not os.path.exists(fname):
            frm = self.files["eden.gz"].format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
                    du=(self.start + dt.timedelta(minutes=i)).strftime("%Y.%m.%d.%H.%M"), model=self.model)
            to = self.files["run"] + self.files["eden.nc"].format(du=(self.start +\
                    dt.timedelta(minutes=i)).strftime("%Y.%m.%d.%H.%M"))
            os.system("cp " + frm + " " + to+".gz")
            os.system("gzip -d " + to+".gz")
            self._nc_ = Dataset(to)
            os.system("rm " + to)
            nex, alts, latx, lonx = _intp_(self._nc_.variables["ZG"+case][0,:,:,:]*1e-3,
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
        if case == "f":
            self._compute_doppler_(i)
            self._plot_radstn_(i)
            txt = r"$V_{d\eta}$=%.1f, $V_{dh}$=%.1f"%(self.vd,self.vf)
            plotlib.plot_rays_base(self.files["run"], self.start + dt.timedelta(minutes=i), i, self.bmnum, case, txt, showrefract=True, 
                    freq=self.frequency)
        return
    
    def _plot_radstn_(self, i):
        """ Plot radar station """
        fname = self.files["run"]+"eden_ti({i}).png".format(i="%02d"%i)
        p = loadmat(self.files["run"] + "ne_ti({ti})_d.mat".format(ti="%02d"%i))["ne"][:,0]
        f = loadmat(self.files["run"] + "ne_ti({ti})_f.mat".format(ti="%02d"%i))["ne"][:,0]
        b = loadmat(self.files["run"] + "ne_ti({ti})_d.mat".format(ti="%02d"%0))["ne"][:,0]
        plotlib.plot_radstn_base(b,p,f,self.ht,fname,self.rlat,self.rlon,self.start+dt.timedelta(minutes=i))
        return

    def _estimate_bearing_(self):
        """ Estimate laitude and logitude bearings """
        fname = self.files["run"] + "bearing.mat"
        m = {}
        lat, lon, bearing = utils.get_sd_radar(self.rad)
        # T0DO
        bearing += (self.bmnum - 12)*3.24
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

    def _compute_(self, i, case="d"):
        """ Compute RT using Pharlap """
        u = self.start + dt.timedelta(minutes=i)
        dic = self.files["run"]
        fname = dic + "ti({ti})_elv(<elv>)_{case}.csv".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), ti="%02d"%i, case=case)
        cmd = "export DIR_MODELS_REF_DAT=/home/shibaji/Collaboration_NCAR/code_rt_sd/pharlap/pharlap_4.1.3/dat;\
                cd pharlap;\
                matlab -nodisplay -nodesktop -nosplash -nojvm -r \"UT=[{ut}];rad='{rad}';dic='{dic}';fname='{fname}';bm={bm};\
                ti='{ti}';cse='{case}';rt_1D_mod;exit;\"".format(ut=u.strftime("%Y %m %d %H %S"),
                rad=self.rad, dic=dic, bm=self.bmnum, ti="%02d"%i, fname=fname, case=case)
        os.system(cmd)
        return

    def _compute_doppler_(self, i):
        """ Compute doppler velocity """
        kconst, cconst, delt = 80.6, 3e8, self.rtime*60
        dic = self.files["run"]
        grange, height = self.m["dist"], self.ht
        fnd = glob.glob(dic + "ti({ti})_elv(*)_d.csv".format(ti="%02d"%i))
        fnf = glob.glob(dic + "ti({ti})_elv(*)_f.csv".format(ti="%02d"%i))
        fnd.sort()
        fnf.sort()
        ned = loadmat(dic + "ne_ti({ti})_d.mat".format(ti="%02d"%i))["ne"]*1e6
        nef = loadmat(dic + "ne_ti({ti})_f.mat".format(ti="%02d"%i))["ne"]*1e6
        funcf = interp2d(grange,height,np.log10(nef))
        funcd = interp2d(grange,height,np.log10(ned))
        for fi in fnf:
            dop, dth, sth = [], [], []
            elvi, rayi = float(fi.split("_")[-2][4:-1]), pd.read_csv(fi)
            gi, hi = np.array(rayi.grange), np.array(rayi.height)
            dth.append(0.)
            sth.append(0.)
            for k in range(len(gi[:-1])):
                dth.append(np.abs(hi[k]-hi[k+1]))
                sth.append(np.sqrt((hi[k]-hi[k+1])**2+(gi[k]-gi[k+1])**2))
            dth, sth = np.array(dth)*1000., np.array(sth)*1000
            for k, g, h in zip(range(len(gi)), gi, hi):
                if h > 50.:
                    dne = (10**funcf(g,h) - 10**funcd(g, h))[0]
                    #df = (kconst / (cconst * self.frequency*1e6)) * (dne / delt) * (dth[k] / np.cos(np.deg2rad(90.-elvi)))
                    df = (kconst / (cconst * self.frequency*1e6)) * (dne / delt) * (dth[k])
                    if np.isnan(df):  df = 0.
                    dop.append(df)
                else: dop.append(0.)
            rayi["dop"] = np.array(dop) / np.cos(np.deg2rad(90.-elvi))
            rayi["sth"] = sth
            rayi["dth"] = dth
            rayi.to_csv(fi, header=True, index=False)
        self._compute_velocity_(i)
        return

    def _compute_velocity_(self, i):
        """ Compute Velocity """
        def _estimate_dop_delh_(x, y, phi=0):
            dh = (np.max(x.height) - np.max(y.height)) * 1000.
            xf = (-2.*self.frequency*1e6/3e8) * (dh/(self.rtime*60.)) * np.cos(np.deg2rad(phi))
            xd = 0.5 * xf * 3e8 / (self.frequency * 1e6)
            return xd
        elvrang = np.arange(self.selev_d, self.eelev_d+1, self.ielev_d)
        dic = self.files["run"]
        fnd = glob.glob(dic + "ti({ti})_elv(*)_d.csv".format(ti="%02d"%i))
        fnf = glob.glob(dic + "ti({ti})_elv(*)_f.csv".format(ti="%02d"%i))
        fnd.sort()
        fnf.sort()
        vd, vf, srng,itr = np.zeros(len(elvrang)), np.zeros(len(elvrang)), np.zeros(len(elvrang)), 0
        for fi, fj in zip(fnf, fnd):
            elvi = float(fi.split("_")[-2][4:-1])
            if elvi in elvrang:
                d = pd.read_csv(fi)
                b = pd.read_csv(fj)
                f = trapz(signal.resample(d.dop,INT_F))
                #f = trapz(np.array(d.dop), np.array(d.grange))
                #f = (trapz(np.array(d.dop), np.array(d.grange)/np.cos(np.deg2rad(90-elvi))))
                srng[itr] = trapz(signal.resample(d.sth,INT_R))
                vd[itr] = (0.5 * f * 3e8 / (self.frequency * 1e6))
                vf[itr] = _estimate_dop_delh_(d,b)
                itr += 1
        fn = dic + "velocity_ti({ti}).mat".format(ti="%02d"%i)
        m = {"vd": vd, "vf": vf, "srng": srng}
        savemat(fn, m)
        self.vd = np.mean(vd.max()+vd.min())
        self.vf = np.mean(vf.max()+vf.min())
        self.vt = np.mean((vf+vd).max()+(vf+vd).min())
        self.recs.append({"dn":self.start+dt.timedelta(minutes=i), "vn":self.vd, "vh": self.vf, "vn_max":np.max(vd), 
            "vn_min":np.min(vd), "vh_max":np.max(vf), "vh_min":np.min(vf), "vt":self.vt, "vt_max":np.max(vf+vd), "vt_min":np.min(vf+vd)})
        if self.verbose:
            print("\tDoppler velocity at {ev} -> Vd={vd}*pm*{vds} m/s, Vf={vf}*pm*{vfs} m/s, Vt={vt}*pm*{vts} m/s".format(ev=
                (self.start+dt.timedelta(minutes=i)).strftime("%Y-%m-%d %H:%M"), 
                    vd=np.round(self.vd,1), vf=np.round(self.vf,1), vt=np.round(self.vf+self.vd,1),
                    vds=np.round(np.std(vd),1), vfs=np.round(np.std(vf),1), vts=np.round(np.std(vf+vd),1)))
        return

    def _exe_(self):
        """ Execute data simulations """
        print("\n Start simulation (using Pharlap) ...")
        self._estimate_bearing_()
        for i in range(self.tsim_start, self.tsim_end):
            [self._interpolate_(i, case) for case in ["d", "f"]]
        return

class Base(object):
    """ Base model for generic functionalities """

    def __init__(self, args):
        """ Initialze the parameters """
        self.args = args
        for k in vars(args).keys():
            setattr(self, k, vars(args)[k])
        self.files = {
                "base": "data/op/{dn}/{model}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), model=self.model),
                "eden.gz": "data/op/{dn}/{model}/{du}.nc.gz",
                }
        self.con = False
        if not os.path.exists(self.files["base"]): os.system("mkdir -p " + self.files["base"])
        if self.tsim_start is None: self.tsim_start = 0
        if self.tsim_end is None: self.tsim_end = int((self.end-self.start).total_seconds()/60.)
        return

    def _conn_(self):
        """ Create Conn """
        self.ssh = SSHClient()
        self.ssh.load_system_host_keys()
        self.ssh.connect(hostname="cheyenne.ucar.edu", port = 22, username="shibaji", password="0515-shibaji-cit")
        self.scp = SCPClient(self.ssh.get_transport())
        self.con = True
        return

    def _close_(self):
        """ Close connection """
        if self.con:
            self.scp.close()
            self.ssh.close()
        return

    def _fetch_sd_(self):
        """ Fetch SuperDARN data and save to local """
        fname = "data/op/{dn}/{model}/sd_{rad}_data.csv".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad, model=self.model)
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
        
        fd = FetchData(self.rad, [self.start, self.end + dt.timedelta(minutes=2)])
        _, scans = fd.fetch_data(v_params=["elv", "v", "w_l", "gflg", "p_l", "slist", "v_e"], by="scan")
        raw_dic = {"vel":[], "beam":[], "gate":[]}
        for sc in scans:
            vel, beam, gate = [], [], [] 
            for b in sc.beams:
                l = len(b.slist)
                beam.append([b.bmnum]*l)
                vel.append(b.v)
                gate.append(b.slist)
            raw_dic["vel"].append(np.array(vel))
            raw_dic["beam"].append(np.array(beam))
            raw_dic["gate"].append(np.array(gate))
        setattr(self, "frequency", np.median(self.rec.tfreq)/1e3)
        self.raw_dic, self.skip = raw_dic, int(self.tsim_end/len(scans) + 1)
        return

    def _run_bmnum_(self):
        """ Parallel process all the beams """
        from joblib import Parallel, delayed
        import multiprocessing
        def run(k):
            args = self.args
            setattr(args, "bmnum", k)
            Model(args)._exe_()
            return
        num_cores = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(run)(k) for k in range(24))
        return

    def _download_(self, i):
        return

    def _plot_summary_(self):
        """ Plot summary data """
        for bm in range(24):
            plotlib.plot_velocity_ts_beam(self.event, self.rad, bm, self.model, self.start, self.end)
        path = self.files["base"] + self.rad + "/"
        def to_gate(srng, irng=180, drng=45):
            gt = (srng/1000 - 180)/45
            gt[-1] = gt[-2]
            return gt
        fanplot = plotlib.FanPlot()
        data_dic = {"vel":[], "beam":[], "gate":[]}
        for i in range(self.tsim_start, self.tsim_end):
            vel, beam, gate = [], [], []
            for bm in range(24):
                rdic = "data/op/{dn}/{model}/{rad}/bm.{bm}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
                        rad=self.rad, model=self.model, bm="%02d"%bm) + "velocity_ti({ti}).mat".format(ti="%02d"%i)
                m = loadmat(rdic)
                vt, srng, l = m["vd"][0] + m["vf"][0], m["srng"][0], len(m["srng"][0])
                gt = to_gate(srng)
                if np.min(gt) < np.max(gt):
                    vel.append(vt.tolist())
                    beam.append([bm]*l)
                    gate.append(gt.tolist())
            data_dic["vel"].append(np.array(vel))
            data_dic["beam"].append(np.array(beam))
            data_dic["gate"].append(np.array(gate))
        fanplot.plot_geo_fov(self.rad, data_dic, range(self.tsim_start, self.tsim_end), "", 
                self.start, self.raw_dic, skip=1, base_filepath=path)
        #fanplot.plot_fov(data_dic, range(self.tsim_start, self.tsim_end), "Radar- %s, Model- %s"%(self.rad, self.model),
        #        self.start, self.raw_dic, skip=self.skip, base_filepath=path)
        return

    def _to_latlon(self, _dict_):
        import pydarn
        hdw = pydarn.read_hdw_file(self.rad)
        rf = rad_fov.CalcFov(hdw=hdw, ngates=75)
        lons, lats = rf.lonFull, rf.latFull
        _dict_["lat"], _dict_["lon"] = [], []
        for bs, gs in zip(_dict_["beam"], _dict_["gate"]):
            lat, lon = [], []
            for b, g in zip(bs, gs):
                lon.append(lons[b, np.round(g).astype(int)].tolist())
                lat.append(lats[b, np.round(g).astype(int)].tolist())
            _dict_["lat"].append(np.array(lat))
            _dict_["lon"].append(np.array(lon))
        return _dict_

    def _exe_(self):
        """ Execute data simulations """
        try:
            self._fetch_sd_()
        except: print("Not able to download data!")
        if hasattr(self, "plot_summary") and self.plot_summary: self._plot_summary_()
        else:
            for i in range(self.tsim_start, self.tsim_end):
                self._download_(i)
            self._close_()
            #self._run_bmnum_()
        return

class WACCM(Base):
    """ WACCMX data download from Cheyenne"""

    def __init__(self, args):
        """ Initialze the parameters """
        super().__init__(args)
        return

    def _download_(self, i):
        """ Download the files from """
        dfname = self.files["eden.gz"].format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
                du=(self.start + dt.timedelta(minutes=i)).strftime("%Y.%m.%d.%H.%M"), model=self.model)
        dn = self.start + dt.timedelta(minutes=i)
        if not os.path.exists(dfname):
            if not self.con: self._conn_()
            if self.verbose: print("\n Create SSH.SCP connection, post process, and fetch data.")
            self.ssh.exec_command("mkdir tmp/")
            self.scp.put("config/wxproc.py", "tmp/")
            if self.verbose: print(" Run -", "python tmp/wxproc.py -ev {ev} -dn {dn}".format(
                ev=self.event.strftime("%Y.%m.%d.%H.%M"), dn=dn.strftime("%Y.%m.%d.%H.%M")))
            stdin, stdout, stderr = self.ssh.exec_command("source ncar_tgcm_waccm_proc/ncar/bin/activate"\
                    "\n python tmp/wxproc.py -ev {ev} -dn {dn}".format(ev=self.event.strftime("%Y.%m.%d.%H.%M"),
                        dn=dn.strftime("%Y.%m.%d.%H.%M")), get_pty=True)
            if self.verbose:
                for line in iter(stdout.readline, ""):
                    print(line, end="")
            if self.verbose: print(" End run ")
            self.ssh.exec_command("deactivate")
            self.scp.get("tmp/waccmx.nc.gz", dfname)
            self.ssh.exec_command("rm -rf tmp/")
        return

class TGCM(Base):
    """ TIME-GCM data download from Cheyenne"""

    def __init__(self, args):
        """ Initialze the parameters """
        super().__init__(args)
        return

    def _download_(self, i):
        """ Download the files from """
        dfname = self.files["eden.gz"].format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
                du=(self.start + dt.timedelta(minutes=i)).strftime("%Y.%m.%d.%H.%M"), model=self.model)
        dn = self.start + dt.timedelta(minutes=i)
        if not os.path.exists(dfname):
            if not self.con: self._conn_()
            if self.verbose: print("\n Create SSH.SCP connection, post process, and fetch data.")
            self.ssh.exec_command("mkdir tmp/")
            self.scp.put("config/proc.py", "tmp/")
            if self.verbose: print(" Run -", "python tmp/proc.py -ev {ev} -dn {dn}".format(
                ev=self.event.strftime("%Y.%m.%d.%H.%M"), dn=dn.strftime("%Y.%m.%d.%H.%M")))
            stdin, stdout, stderr = self.ssh.exec_command("source ncar_tgcm_waccm_proc/ncar/bin/activate"\
                    "\n python tmp/proc.py -ev {ev} -dn {dn}".format(ev=self.event.strftime("%Y.%m.%d.%H.%M"),
                        dn=dn.strftime("%Y.%m.%d.%H.%M")), get_pty=True)
            if self.verbose:
                for line in iter(stdout.readline, ""):
                    print(line, end="")
            if self.verbose: print(" End run ")
            self.ssh.exec_command("deactivate")
            self.scp.get("tmp/tgcm.nc.gz", dfname)
            self.ssh.exec_command("rm -rf tmp/")
        return
