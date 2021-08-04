

import os
import datetime as dt
import numpy as np
import pandas as pd
import pytz
from netCDF4 import Dataset
from timezonefinder import TimezoneFinder
import aacgmv2
import traceback


tf = TimezoneFinder(in_memory=True)


def convert_geomag(x):
    mlat, mlon, mlt = aacgmv2.get_aacgm_coord(x["lat"], x["lon"], 150, x["date"])
    x["mlat"], x["mlon"], x["mlt"] = np.round( mlat, 2), np.round( mlon, 2), np.round( mlt, 2)
    return x

def get_isr_data(f="stats/isr-jro-data/jul20050907_avg_150km.001.txt"):
    dates, lat, lon, vipn2, dvipn2, vipe1, dvipe1 = [], [], [], [], [], [], []
    with open(f, "r") as f: lines = f.readlines()
    for l in lines[1:]:
        l = list(filter(None, l.replace("\n", "").replace("missing","NaN").split(" ")))
        d = dt.datetime(int(l[0]), int(l[1]), int(l[2]), int(l[3]), int(l[4]), int(l[5]))
        dates.append(d)
        lat.append(float(l[6]))
        lon.append(float(l[7]))
        vipn2.append(float(l[11]))
        dvipn2.append(float(l[12]))
        vipe1.append(float(l[13]))
        dvipe1.append(float(l[14]))
    u = pd.DataFrame()
    u["date"], u["lat"], u["lon"], u["vipn2"], u["dvipn2"], u["vipe1"], u["dvipe1"] = dates, lat, lon, vipn2, dvipn2, vipe1, dvipe1
    local_time_zone = tf.timezone_at(lng=u.lon.mean(), lat=u.lat.mean())
    timezone = pytz.timezone(local_time_zone)
    u.date = [timezone.localize(d).astimezone(pytz.utc) + dt.timedelta(hours=-5) for d in u.date]
    u.date = [d.replace(day=7) for d in u.date]
    u = u.apply(convert_geomag, axis=1)
    return u

def get_mlat_mlon_index(mlats, mlons, mlat, mlon):
    i, j = np.argmin(np.abs(mlats-mlat)), np.argmin(np.abs(mlons-mlon))
    return i, j

def get_field_data(dates, lats, lons, fname):
    if os.path.exists(fname): o = pd.read_csv(fname, parse_dates=["date"])
    else:
        pc, hc, ed1, ed2 = [], [], [], []
        for d in dates:
            print(" Date:", d)
            fname = "data/op/2005.09.07.17.40/waccmx/%s.nc.gz"%d.strftime("%Y.%m.%d.%H.%M")
            os.system("gzip -d " + fname)
            fname = fname.replace(".gz", "")
            try:
                ds = Dataset(fname)
                mlats, mlons = ds.variables["mlat"][:], ds.variables["mlon"][:]
                ED1, ED2 = ds.variables["ED1"][:], ds.variables["ED2"][:]
                PC, HC = ds.variables["PC"][:], ds.variables["HC"][:]
                for lat, lon in zip(lats, lons):
                    i, j = get_mlat_mlon_index(mlats, mlons, lat, lon)
                    pc.append(PC[0,i,j])
                    hc.append(HC[0,i,j])
                    ed1.append(ED1[0,i,j])
                    ed2.append(ED2[0,i,j])
            except: traceback.print_exc()
            os.system("gzip " + fname)
        o = pd.DataFrame()
        o["date"], o["pc"], o["hc"], o["ed1"], o["ed2"] = dates, pc, hc, ed1, ed2
        o["mlat"], o["mlon"] = lats, lons
        o.to_csv(fname, index=False, header=True)
    print(o.head())
    return o


case = "field-analysis"
if case == "run-cmd":
    cmd = "nohup python model_sim.py -r wal -ev 2005-09-07T17:40 -s 2005-09-07T17 -e 2005-09-07T18 > /dev/null 2>&1 &"
    os.system(cmd)
elif case == "run-plot":
    cmd = "nohup python model_sim.py -r wal -ev 2005-09-07T17:40 -s 2005-09-07T17 -e 2005-09-07T18 -ps > /dev/null 2>&1 &"
    os.system(cmd)
elif case == "field-analysis": 
    f = "stats/isr-jro-data/jul20050907_avg_150km.001.txt"
    dates = [dt.datetime(2005,9,7,17) + dt.timedelta(minutes=i) for i  in range(60)]
    u = get_isr_data(f)
    o = get_field_data(dates, np.array([u.mlat.tolist()[0]]), np.array([u.mlon.tolist()[0]]), "stats/isr-jro-data/jro_model.csv")
elif case == "ion-drift-analsys": pass
elif case == "conductance-analsys": pass
else: print(" Case not found:", case)

os.system("rm -rf *.log")
os.system("rm -rf __pycache__")
