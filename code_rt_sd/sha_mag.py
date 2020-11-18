#!/usr/bin/env python

"""sha_mag.py: SHA from magnetometer data python program for Sq"""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import pandas as pd
import numpy as np

def P(t, m, n):
    if n==0 and m==0: return 1
    if n==1 and m==0: return np.cos(t)
    if n==1 and m==1: return np.sin(t)
    if n==m and m>1: return np.sqrt((2*m-1)/(2*m)) * np.sin(t) * P(t, m-1, n-1)
    if n>m: return ( ((2*n-1)*np.cos(t)*P(t, m, n-1)) - (np.sqrt((n-1)**2-m**2)*P(t, m, n-2)) ) / (np.sqrt(n**2-m**2))
    else: return 0.

def dP(t, m, n):
    return ( (n*np.cos(t)*P(t, m, n)) - (np.sqrt(n**2-m**2) * P(t, m, n-1)) ) / np.sin(t)

Xhu, Xhl = lambda t, p, n, m: np.cos(m*p) * dP(t, m, n), lambda t, p, n, m: np.sin(m*p) * dP(t, m, n)
Yhu, Yhl = lambda t, p, n, m: m * np.sin(m*p) * P(t, m, n), lambda t, p, n, m: -m * np.cos(m*p) * P(t, m, n)
Zhu, Zhl = lambda t, p, n, m: np.cos(m*p) * P(t, m, n), lambda t, p, n, m: np.sin(m*p) * P(t, m, n)

def create_a1_mat(colat, lon, order=3):
    K, T = len(colat), np.sum([t+1 for t in range(1, order+1)])
    a1 = np.zeros((2*K, 2*T))
    for _k, t, p in zip(range(K), colat, lon):
        _t = 0
        for n in range(1, order+1):
            for m in range(n+1):
                a1[_k, _t] = Xhu(t,p,n,m)
                a1[_k, T+_t] = Xhl(t,p,n,m)
                a1[K+_k, _t] = Yhu(t,p,n,m)
                a1[K+_k, T+_t] = Yhl(t,p,n,m)
                _t += 1
    return a1

def create_a2_mat(colat, lon, order=3):
    K, T = len(colat), np.sum([t+1 for t in range(1, order+1)])
    a2 = np.zeros((K, 2*T))
    for _k, t, p in zip(range(K), colat, lon):
        _t = 0
        for n in range(1, order+1):
            for m in range(n):
                a2[_k, _t] = Zhu(t,p,n,m)
                a2[_k, T+_t] = Zhl(t,p,n,m)
    return a2

files = "data/supermag.csv"
stns = ["BOU", "TUC", "T25", "BSL", "M06", "VIC", "T20"]
dat = pd.read_csv(files, parse_dates=["Date_UTC"])[["Date_UTC", "IAGA", "GEOLON", "GEOLAT", "MAGON", "MAGLAT", "MLT", "MCOLAT", "SZA", "dbn_nez", "dbe_nez", "dbz_nez"]]
mins = int((dat.Date_UTC.tolist()[-1]-dat.Date_UTC.tolist()[0]).total_seconds()/60) + 1
dates = [dat.Date_UTC.tolist()[0] + dt.timedelta(minutes=i) for i in range(mins)]
O = 1

for d in dates:
    _o = dat[(dat.Date_UTC==d) & dat.IAGA.str.contains("|".join(stns))]
    colat, lon = np.deg2rad(90. - np.array(_o["GEOLAT"])), np.deg2rad(np.mod( (_o["GEOLON"] + 180), 360 ) - 180)
    dX, dY, dZ = np.array(_o["dbn_nez"]), np.array(_o["dbe_nez"]), np.array(_o["dbz_nez"])
    d1, d2 = np.reshape(np.append(dX,dY), (2*len(stns),1)), np.reshape(dZ, (len(stns), 1))
    A1, A2 = create_a1_mat(colat, lon, O), create_a2_mat(colat, lon, O)
    print(A1.shape, A1.T.shape, d1.shape)
    #M1 = np.linalg.inv(np.matmul(A1.T,A1))
    print(A1,np.linalg.det(np.matmul(A1.T,A1)))
    #M1, M2 = np.matmul(np.matmul(np.matmul(A1.T,A1), A1.T, d1)), np.matmul(np.matmul(np.matmul(A2.T,A2), A2.T, d2))
    break
