#!/usr/bin/env python

"""simulate.py: simulate python program for Doppler shift"""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
import sys
sys.path.append("sd/")
import datetime as dt
import argparse
from dateutil import parser as dparser

from tgcm import TGCM
from waccmx import WACCM

if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("-m", "--model", default="waccmx", help="Model name [tgcm/waccmx] (default tgcm)")
        parser.add_argument("-r", "--rad", default="bks", help="Radar code (default bks)")
        parser.add_argument("-ev", "--event", default=dt.datetime(2015,3,11,16,22), help="Start date (default 2015-3-11T16:22)",
                type=dparser.isoparse)
        parser.add_argument("-si", "--sim_id", default="L100", help="Simulation ID (default L100)")
        parser.add_argument("-b", "--bmnum", type=int, default=0, help="Beam num starting from 0 in center and +, - on right-left")
        parser.add_argument("-ts", "--tsim_start", type=int, default=None, help="Minutes to simulate start (0-...) (default None)")
        parser.add_argument("-te", "--tsim_end", type=int, default=None, help="Minutes to simulate end (1-...) (default None)")
        parser.add_argument("-s", "--start", default=dt.datetime(2015,3,11,16,15), help="Start date (default 2015-3-11T16:05)",
                type=dparser.isoparse)
        parser.add_argument("-e", "--end", default=dt.datetime(2015,3,11,16,24), help="End date (default 2015-3-11T16:25)",
                type=dparser.isoparse)
        parser.add_argument("-rd", "--save_radar", action="store_false", help="Save riometer data (default True)")
        parser.add_argument("-ps", "--plot_summary", action="store_true", help="Plot summary report (default False)")
        parser.add_argument("-sr", "--save_result", action="store_false", help="Save results (default True)")
        parser.add_argument("-c", "--clear", action="store_true", help="Clear pervious stored files (default False)")
        parser.add_argument("-v", "--verbose", action="store_false", help="Increase output verbosity (default True)")
        parser.add_argument("-a", "--archive", action="store_false", help="Archive and create movies (default True)")
        parser.add_argument("-fr", "--frequency", type=float, default=12., help="Frequency of oprrations in MHz (default 12 MHz)")
        parser.add_argument("-mr", "--mrange", type=float, default=1990., help="Max ground range in km (default 1000 km)")
        parser.add_argument("-nmr", "--nmrange", type=int, default=200, help="Number of steps in ground range (default 101)")
        parser.add_argument("-sh", "--sheight", type=float, default=50., help="Start height in km (default 50 km)")
        parser.add_argument("-nmh", "--eheight", type=float, default=350., help="End height in km (default 350 km)")
        parser.add_argument("-hinc", "--hinc", type=float, default=1., help="Step in height in km (default 1 km)")
        parser.add_argument("-es", "--selev", type=float, default=30., help="Start elevation angle deg")
        parser.add_argument("-ee", "--eelev", type=float, default=35., help="End elevation angle deg")
        parser.add_argument("-ei", "--ielev", type=float, default=1., help="Inc of elevation angle deg (default 20*)")
        parser.add_argument("-nhops", "--nhops", type=float, default=1, help="Number of hops (default 1)")
        parser.add_argument("-th", "--threshold", type=float, default=1.e-5, help="Threshold for risetime, (1.e-5)")
        args = parser.parse_args()
        if args.verbose:
            print("\n Parameter list for simulation ")
            for k in vars(args).keys():
                print("     " , k , "->" , str(vars(args)[k]))
        if args.model == "waccmx": WACCM(args)._exe_()
        if args.model == "tgcm": TGCM(args)._exe_()
        print("")
        if os.path.exists("sd/__pycache__/"):
            os.system("rm -rf sd/__pycache__/")
            os.system("rm -rf *.log")
