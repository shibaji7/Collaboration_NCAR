#!/usr/bin/env python

"""experiment.py: experiment python program for Doppler shift"""

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

from superdarn import SuperDARN

if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("-r", "--rad", default="bks", help="Radar code (default bks)")
        parser.add_argument("-ev", "--event", default=dt.datetime(2015,3,11,16,22), help="Start date (default 2015-3-11T16:22)",
                type=dparser.isoparse)
        parser.add_argument("-b", "--bmnum", type=int, default=0, help="Beam num starting from 0 in center and +, - on right-left")
        parser.add_argument("-dra", "--d_ratio", type=float, default=100., help="D region e-denisity peak times (100)")
        parser.add_argument("-ds", "--d_start", type=float, default=10., help="D region start (60.)")
        parser.add_argument("-de", "--d_end", type=float, default=35., help="D region end (80.)")
        parser.add_argument("-drt", "--d_rtime", type=float, default=2., help="D region e-denisity rise time (2.)")
        parser.add_argument("-era", "--e_ratio", type=float, default=1., help="E region e-denisity peak order (1.)")
        parser.add_argument("-fra", "--f_ratio", type=float, default=1.01, help="F region e-denisity peak order")
        parser.add_argument("-fs", "--f_start", type=float, default=130., help="F region start (180.)")
        parser.add_argument("-fe", "--f_end", type=float, default=190, help="F region end (240.)")
        parser.add_argument("-frt", "--f_rtime", type=float, default=0.5, help="F region e-denisity rise time (0.5)")
        parser.add_argument("-rd", "--save_radar", action="store_false", help="Save riometer data (default True)")
        parser.add_argument("-ps", "--plot_summary", action="store_true", help="Plot summary report (default False)")
        parser.add_argument("-sr", "--save_result", action="store_false", help="Save results (default True)")
        parser.add_argument("-c", "--clear", action="store_true", help="Clear pervious stored files (default False)")
        parser.add_argument("-v", "--verbose", action="store_false", help="Increase output verbosity (default True)")
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
        args = parser.parse_args()
        if args.verbose:
            print("\n Parameter list for simulation ")
            for k in vars(args).keys():
                print("     " , k , "->" , str(vars(args)[k]))
        SuperDARN(args)._exe_()
        print("")
        if os.path.exists("sd/__pycache__/"):
            os.system("rm -rf sd/__pycache__/")
            os.system("rm -rf *.log")
