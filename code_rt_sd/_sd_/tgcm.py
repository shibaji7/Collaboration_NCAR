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
from scipy.io import savemat

from get_sd_data import FetchData
import utils

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
        dfname = "_data_/_sim_/{dn}/timegcm.nc".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"))
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
            scp.put("_config_/proc.py", "tmp/")
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
            dfname = "_data_/_sim_/{dn}/timegcm.nc".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"))
            gzfile = dfname + ".gz"
            os.system("gzip -d " + gzfile)
            self._nc_ = Dataset(dfname)
            os.system("gzip " + dfname)
        nex, alts, latx, lonx = _intp_(self._nc_.variables["ZG"][i,:,:,:]*1e-5,
                self._nc_.variables["lat"][:], self._nc_.variables["lon"][:], self._nc_.variables["NE"][i,:,:,:], 
                hd=[self.sheight,self.eheight,self.hinc], 
                dlat=2, dlon=4, scale="log", kind="cubic", v=self.verbose)
        m["ne"], m["h"], m["la"], m["lo"] = nex, alts, latx, lonx
        fname = "_data_/_sim_/{dn}/{rad}/tmp.mat".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        savemat(fname, m)
        return

    def _estimate_bearing_(self):
        """ Estimate laitude and logitude bearings """
        fname = "_data_/_sim_/{dn}/{rad}/bearing.mat".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        m = {}
        if not os.path.exists(fname):
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
            m["freq"], m["tol"], m["nhops"] = float(self.frequency), float(1e-7), float(self.nhops)
            m["elev_s"], m["elev_i"], m["elev_e"] = float(self.selev), float(self.ielev), float(self.eelev)
            m["radius_earth"] = 6371.0
            savemat(fname, m)
        return

    def _compute_(self, i):
        """ Compute RT using Pharlap """
        u = self.start + dt.timedelta(minutes=i)
        dic = "_data_/_sim_/{dn}/{rad}/".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
        fname = "_data_/_sim_/{dn}/{rad}/ti({ti}).bm({bm}).elv(<elv>).csv".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), 
                rad=self.rad, bm=self.bmnum, ti=i)
        cmd = "export DIR_MODELS_REF_DAT=/home/shibaji/Collaboration_NCAR/code_rt_sd/_pharlap_/pharlap_4.1.3/dat;\
                cd _pharlap_;\
                matlab -nodisplay -nodesktop -nosplash -nojvm -r \"UT=[{ut}];rad='{rad}';dic='{dic}';fname='{fname}';bm='{bm}';\
                ti={ti};rt_1D;exit;\"".format(ut=u.strftime("%Y %m %d %H %S"),rad=self.rad, dic=dic, bm=self.bmnum, ti=i, fname=fname)
        os.system(cmd)
        os.remove(dic + "tmp.mat")
        return

    def _fetch_sd_(self):
        """ Fetch SuperDARN data and save to local """
        fname = "_data_/_sim_/{dn}/{rad}/sd_data.csv".format(dn=self.event.strftime("%Y.%m.%d.%H.%M"), rad=self.rad)
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

    def _exe_(self):
        """ Execute the RT model and save results"""
        print("\n Start simulation (using Pharlap) ...")
        self._nc_ = None
        self._estimate_bearing_()
        if hasattr(self, "save_radar") and self.save_radar: self._fetch_sd_()
        for i in range(self.tsim_start, self.tsim_end):
            if self.verbose: print("\tProcess-", self.start + dt.timedelta(minutes=i))
            self._interpolate_(i)
            self._compute_(i)
            break
        return
