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
import glob
from scipy.integrate import trapz

from get_sd_data import FetchData
import utils
import plotlib


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

def _intp_(h, lats, lons, param, hd=[50,350,1], dlat=0.5, dlon=1, scale="log", kind="cubic", v=False):
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
            pnew[:,i,j] = _intp_heights_(h_intm[:, i, j], hx, param_intm[:, i, j], scale, kind)
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
        self._download_()
        return

    def _download_(self):
        """ Download the files from """
        doy = self.event.timetuple().tm_yday
        fname = "/glade/work/shibaji/timegcm/cases/{dn}-{sim}.2/timegcm_trunk/"\
                "timegcm-ch/timegcm.s_{tm}.nc".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), 
                        sim=self.sim_id, tm="%03d"%doy)
        dfname = "data/sim/{dn}/timegcm.nc".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"))
        gzfile = dfname + ".gz"
        start, end = (self.start.hour*60) + self.start.minute, (self.end.hour*60) + (self.end.minute+1)
        if not os.path.exists(gzfile):
            if self.verbose: print("\n Create SSH.SCP connection, post process, and fetch data.")
            if self.verbose: print(" Fname - ", fname)
            ssh = SSHClient()
            ssh.load_system_host_keys()
            ssh.connect(hostname="cheyenne.ucar.edu", port = 22, username="shibaji", password="0513-shibaji-cit")
            ssh.exec_command("mkdir tmp/")
            scp = SCPClient(ssh.get_transport())
            scp.put("config/proc.py", "tmp/")
            if self.verbose: print(" Run -", "python tmp/proc.py -f {f} -s {s} -e {e}".format(f=fname, s=start, e=end))
            stdin, stdout, stderr = ssh.exec_command("source ncar_tgcm_waccm_proc/ncar/bin/activate"\
                    "\n python tmp/proc.py -f {f} -s {s} -e {e}".format(f=fname, s=start, e=end), get_pty=True)
            if self.verbose:
                for line in iter(stdout.readline, ""):
                    print(line, end="")
            ssh.exec_command("deactivate")
            scp.get("tmp/timegcm.nc", dfname)
            ssh.exec_command("rm -rf tmp/")
            scp.close()
            ssh.close()
            os.system("gzip " + dfname)
        return

    def _interpolate_(self, i):
        """ Interpolate electron density """
        m = {}
        if self._nc_ is None:
            dfname = "data/sim/{dn}/timegcm.nc".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"))
            gzfile = dfname + ".gz"
            os.system("gzip -d " + gzfile)
            self._nc_ = Dataset(dfname)
            self.tx = num2date(self._nc_.variables["time"][:], units=self._nc_.variables["time"].units)
            self.tx = np.array([x._to_real_datetime() for x in self.tx]).astype("datetime64[ns]")
            self.tx = [dt.datetime.utcfromtimestamp(x.astype(int) * 1e-9) - dt.timedelta(days=1) for x in self.tx]
            os.system("gzip " + dfname)
        fname = "data/sim/{dn}/{rad}/ne.ti({ti}).bm({bm}).mat".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), 
                rad=self.rad, ti=i, bm=self.bmnum)
        if not os.path.exists(fname):
            nex, alts, latx, lonx = _intp_(self._nc_.variables["ZG"][i,:,:,:]*1e-5,
                    self._nc_.variables["lat"][:], self._nc_.variables["lon"][:], self._nc_.variables["NE"][i,:,:,:], 
                    hd=[self.sheight,self.eheight,self.hinc], 
                    dlat=2, dlon=4, scale="log", kind="cubic", v=self.verbose)
            ne = np.zeros((len(self.ht), len(self.lat)))
            for ix, la, lo in zip(range(len(self.lat)), self.lat, self.lon):
                _ij_ = _get_ij_(latx, lonx, la, lo)
                ne[:, ix] = nex[:, _ij_[0], _ij_[1]]
            m["ne"] = ne
            savemat(fname, m)
            self._compute_(i)
            dic = "data/sim/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
            plotlib.plot_rays(dic, self.start + dt.timedelta(minutes=i), i, self.bmnum)
        else: ne = loadmat(fname)["ne"]
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

    def _compute_doppler_(self, i, j=0):
        """ Compute doppler velocity """
        if self.verbose: print("\tCompute Doppler for -", self.start + dt.timedelta(minutes=i), " against - ", 
                self.start + dt.timedelta(minutes=j))
        k, c, dth, delt = 8.06e-11, 3e8, 30*1000, (i-j)*60
        dic = "data/sim/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        grange, height = self.m["dist"], self.ht
        fis = glob.glob(dic + "ti({ti}).bm({bm}).elv(*).csv".format(ti=i, bm=self.bmnum))
        fis.sort()
        fni = "data/sim/{dn}/{rad}/ne.ti({ti}).bm({bm}).mat".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
                rad=self.rad, ti=i, bm=self.bmnum)
        fnj = "data/sim/{dn}/{rad}/ne.ti({ti}).bm({bm}).mat".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"),
                rad=self.rad, ti=j, bm=self.bmnum)
        nei, nej = loadmat(fni)["ne"]*1e6, loadmat(fnj)["ne"]*1e6
        for fi in fis:
            dop, dist = [], []
            elvi, rayi = float(fi.split(".")[-2][4:-1]), pd.read_csv(fi)
            gi, hi = np.array(rayi.grange), np.array(rayi.height)
            funci = interp2d(grange,height,nei,kind="cubic")
            funcj = interp2d(grange,height,nej,kind="cubic")
            dist.append(0.)
            for k in range(len(gi[:-1])):
                dist.append(np.sqrt((gi[k]-gi[k+1])**2+(hi[k]-hi[k+1])**2))
            for g, h in zip(gi, hi):
                if h > 50.:
                    dne = (funci(g,h) - funcj(g, h))[0]
                    df = (k / (c * self.frequency*1e6)) * (dne / delt) * (dth / np.cos(np.deg2rad(90.-elvi)))
                    dop.append(df)
                else: dop.append(0.)
            rayi["dop"] = dop
            rayi["dist"] = dist
            rayi.to_csv(fi, header=True, index=False)
        #self._compute_velocity_(i, j)
        return

    def _compute_velocity_(self, i, j=0):
        """ Compute Velocity """
        dic = "data/sim/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        fis = glob.glob(dic + "ti({ti}).bm({bm}).elv(*).csv".format(ti=i, bm=self.bmnum))
        fis.sort()
        f = 0.
        for fi in fis:
            d = pd.read_csv(fi)
            print(d.head(30))
            print(np.abs(trapz(np.array(d.dop), np.array(d.dist))))
            f += np.abs(trapz(np.array(d.dop), np.array(d.dist)))
            break
        vd = 0.5 * f * 3e8 / (12e6*len(fis))
        if self.verbose: print("\tDoppler velocity -", self.start + dt.timedelta(minutes=i), " against - ",
                self.start + dt.timedelta(minutes=j), " V.d = ", vd, "m/s")
        return

    def _compute_(self, i):
        """ Compute RT using Pharlap """
        u = self.start + dt.timedelta(minutes=i)
        dic = "data/sim/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        fname = "data/sim/{dn}/{rad}/ti({ti}).bm({bm}).elv(<elv>).csv".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), 
                rad=self.rad, bm=self.bmnum, ti=i)
        cmd = "export DIR_MODELS_REF_DAT=/home/shibaji/Collaboration_NCAR/code_rt_sd/pharlap/pharlap_4.1.3/dat;\
                cd pharlap;\
                matlab -nodisplay -nodesktop -nosplash -nojvm -r \"UT=[{ut}];rad='{rad}';dic='{dic}';fname='{fname}';bm={bm};\
                ti={ti};rt_1D;exit;\"".format(ut=u.strftime("%Y %m %d %H %S"),rad=self.rad, dic=dic, bm=self.bmnum, ti=i, fname=fname)
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

    def _plot_radstn_(self):
        """ Plot radar station """
        fname = "data/sim/{dn}/{rad}/tgcm.png".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        ij = _get_ij_(self._nc_.variables["lat"][:], self._nc_.variables["lon"][:], self.rlat, self.rlon)
        xi, xj = self.tx.index(self.start), self.tx.index(self.start + dt.timedelta(minutes=7))
        p = self._nc_.variables["NE"][xi,:,ij[0],ij[1]]
        f = self._nc_.variables["NE"][xj,:,ij[0],ij[1]]
        pz = self._nc_.variables["ZG"][xi,:,ij[0],ij[1]]*1e-5
        fz = self._nc_.variables["ZG"][xj,:,ij[0],ij[1]]*1e-5
        print(p.tolist())
        print(f.tolist())
        plotlib.plot_radstn(p,f,pz,fz,fname,self.rlat,self.rlon,self.event)
        return

    def _exe_(self):
        """ Execute the RT model and save results"""
        print("\n Start simulation (using Pharlap) ...")
        self._nc_ = None
        self._estimate_bearing_()
        if hasattr(self, "save_radar") and self.save_radar: self._fetch_sd_()
        for i in range(self.tsim_start, self.tsim_end):
            if self.verbose: print("\tProcess-", self.start + dt.timedelta(minutes=i))
            self._interpolate_(i)
        if self.verbose: print("\n Interpolation completed.")
        if self.verbose: print("\n Processing Doppler.")
        self._plot_radstn_()
        self._compute_doppler_(7)
        return
