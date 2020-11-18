#!/usr/bin/env python

"""get_sd_data.py: module is dedicated to fetch fitacf data from files."""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import numpy as np
import pandas as pd
import datetime as dt
import glob
import bz2
import pydarn

class Gate(object):
    """Class object to hold each range cell value"""

    def __init__(self, bm, i, params=["v", "w_l", "gflg", "p_l", "v_e"], gflg_type=-1):
        """
        initialize the parameters which will be stored
        bm: beam object
        i: index to store
        params: parameters to store
        """
        for p in params:
            setattr(self, p, getattr(bm, p)[i])
        if gflg_type >= 0: setattr(self, "gflg", getattr(bm, "gsflg")[gflg_type][i])
        return


class Beam(object):
    """Class to hold one beam object"""

    def __init__(self):
        """
        initialize the instance
        """
        return

    def set(self, time, d, s_params=["bmnum", "noise.sky", "tfreq", "scan", "nrang"], 
            v_params=["pwr0", "v", "w_l", "gflg", "p_l", "slist", "v_e"]):
        """
        Set all parameters
        time: datetime of beam
        d: data dict for other parameters
        s_param: other scalar params
        v_params: other list params
        """
        self.time = time
        for p in s_params:
            if p in d.keys(): setattr(self, p, d[p])
            else: setattr(self, p, None)
        for p in v_params:
            if p in d.keys(): setattr(self, p, d[p])
            else: setattr(self, p, [])
        self.gs_estimation()
        return

    def copy(self, bm):
        """
        Copy all parameters
        """
        for p in bm.__dict__.keys():
            setattr(self, p, getattr(bm, p))
        return

    def gs_estimation(self):
        """
        Estimate GS flag using different criterion
        Cases - 
                0. Sundeen et al. |v| + w/3 < 30 m/s
                1. Blanchard et al. |v| + 0.4w < 60 m/s
                2. Blanchard et al. [2009] |v| - 0.139w + 0.00113w^2 < 33.1 m/s
        """
        self.gsflg = {}
        self.gsflg_prob = {}
        if len(self.v) > 0 and len(self.w_l) > 0: 
            self.gsflg[0] = ((np.abs(self.v) + (self.w_l/3.)) < 30.).astype(int)
            self.gsflg_prob[0] = 1. / (1 + np.exp(np.abs(self.v) + (self.w_l/3.) - 30.))
        if len(self.v) > 0 and len(self.w_l) > 0: 
            self.gsflg[1] = ((np.abs(self.v) + (self.w_l*0.4)) < 60.).astype(int)
            self.gsflg_prob[1] = 1. / (1 + np.exp(np.abs(self.v) + (self.w_l*0.4) - 60.))
        if len(self.v) > 0 and len(self.w_l) > 0: 
            self.gsflg[2] = ((np.abs(self.v) - (0.139*self.w_l) + (0.00113*self.w_l**2)) < 33.1).astype(int)
            self.gsflg_prob[2] = 1. / (1 + np.exp(np.abs(self.v) + (0.139*self.w_l) + (0.00113*self.w_l**2) - 33.1))
        return


class Scan(object):
    """Class to hold one scan (multiple beams)"""

    def __init__(self, stime=None, etime=None, stype="normal"):
        """
        initialize the parameters which will be stored
        stime: start time of scan
        etime: end time of scan
        stype: scan type
        """
        self.stime = stime
        self.etime = etime
        self.stype = stype
        self.beams = []
        return

    def update_time(self, up=True):
        """
        Update stime and etime of the scan.
        up: Update average parameters if True
        """
        self.stime = self.beams[0].time
        self.etime = self.beams[-1].time
        if up: self._populate_avg_params()
        return

    def _populate_avg_params(self):
        """
        Polulate average parameetrs
        """
        f, nsky = [], []
        for b in self.beams:
            f.append(getattr(b, "tfreq"))
            nsky.append(getattr(b, "noise.sky"))
        self.f, self.nsky = np.mean(f), np.mean(nsky)
        return


class FetchData(object):
    """Class to fetch data from fitacf files for one radar for atleast a day"""

    def __init__(self, rad, date_range, files=None, verbose=True):
        """
        initialize the vars
        rad = radar code
        date_range = [ start_date, end_date ]
        files = List of files to load the data from
        e.x :   rad = "sas"
                date_range = [
                    datetime.datetime(2017,3,17),
                    datetime.datetime(2017,3,18),
                ]
        """
        self.rad = rad
        self.date_range = date_range
        self.files = files
        self.verbose = verbose
        if (rad is not None) and (date_range is not None) and (len(date_range) == 2):
            self._create_files()
        return

    def _create_files(self):
        """
        Create file names from date and radar code
        """
        if self.files is None: self.files = []
        reg_ex = "/sd-data/{year}/fitacf/{rad}/{date}.*.{rad}.fitacf.bz2"
        days = (self.date_range[1] - self.date_range[0]).days + 1
        ent = -1
        for d in range(-1,days):
            e = self.date_range[0] + dt.timedelta(days=d)
            fnames = glob.glob(reg_ex.format(year=e.year, rad=self.rad, date=e.strftime("%Y%m%d")))
            fnames.sort()
            for fname in fnames:
                tm = fname.split(".")[1]
                sc = fname.split(".")[2]
                dus = dt.datetime.strptime(fname.split(".")[0].split("/")[-1] + tm + sc, "%Y%m%d%H%M%S")
                due = dus + dt.timedelta(hours=2)
                if (ent == -1) and (dus <= self.date_range[0] <= due): ent = 0
                if ent == 0: self.files.append(fname)
                if (ent == 0) and (dus <= self.date_range[1] <= due): ent = -1
        return

    def _parse_data(self, data, s_params, v_params, by, scan_prop):
        """
        Parse data by data type
        data: list of data dict
        params: parameter list to fetch
        by: sort data by beam or scan
        scan_prop: provide scan properties if by='scan'
                        {"stype": type of scan, "dur": duration in min}
        """
        _b, _s = [], []
        if self.verbose: print("\n Started converting to beam data.")
        for d in data:
            time = dt.datetime(d["time.yr"], d["time.mo"], d["time.dy"], d["time.hr"], d["time.mt"], d["time.sc"], d["time.us"])
            if time >= self.date_range[0] and time <= self.date_range[1]:
                bm = Beam()
                bm.set(time, d, s_params,  v_params)
                _b.append(bm)
        if self.verbose: print("\n Converted to beam data.")
        if by == "scan":
            if self.verbose: print("\n Started converting to scan data.")
            scan, sc =  0, Scan(None, None, scan_prop["stype"])
            sc.beams.append(_b[0])
            for d in _b[1:]:
                if d.scan == 1:
                    sc.update_time()
                    _s.append(sc)
                    sc = Scan(None, None, scan_prop["stype"])
                    sc.beams.append(d)
                else: sc.beams.append(d)
            if self.verbose: print("\n Converted to scan data.")
        return _b, _s

    def convert_to_pandas(self, beams, s_params=["bmnum", "noise.sky", "tfreq", "scan", "nrang", "time"],
            v_params=["v", "w_l", "gflg", "p_l", "slist", "v_e"]):
        """
        Convert the beam data into dataframe
        """
        _o = dict(zip(s_params+v_params, ([] for _ in s_params+v_params)))
        for b in beams:
            l = len(getattr(b, "slist"))
            for p in v_params:
                _o[p].extend(getattr(b, p))
            for p in s_params:
                _o[p].extend([getattr(b, p)]*l)
        return pd.DataFrame.from_records(_o)

    def fetch_data(self, s_params=["bmnum", "noise.sky", "tfreq", "scan", "nrang"],
                        v_params=["pwr0", "v", "w_l", "gflg", "p_l", "slist", "v_e"],
                        by="beam", scan_prop={"dur": 1, "stype": "normal"}):
        """
        Fetch data from file list and return the dataset
        params: parameter list to fetch
        by: sort data by beam or scan
        scan_prop: provide scan properties if by='scan' 
                   {"stype": type of scan, "dur": duration in min}
        """
        data = []
        for f in self.files:
            with bz2.open(f) as fp:
                fs = fp.read()
            if self.verbose: print("Read file - ", f)
            reader = pydarn.SDarnRead(fs, True)
            records = reader.read_fitacf()
            data += records
        if by is not None: data = self._parse_data(data, s_params, v_params, by, scan_prop)
        return data

if __name__ == "__main__":
    fdata = FetchData( "sas", [dt.datetime(2015,3,17,3),
        dt.datetime(2015,3,17,3,20)] )
    fdata.fetch_data()
    fdata.fetch_data(by="scan", scan_prop={"dur": 2, "stype": "themis"})
    import os
    os.system("rm *.log")
    os.system("rm -rf __pycache__/")
