import os
import sys
import datetime as dt
import pandas as pd

I = int(sys.argv[1])
events = pd.read_csv("config/events.csv", parse_dates=["dn"])
for i, event in events.iterrows():
    if i == I:
        evnt = event["dn"]
        start = evnt - dt.timedelta(minutes=20)
        end = evnt + dt.timedelta(minutes=41)
        if end.day != start.day: end = end.replace(minute=0)
        cmd = "nohup python model_sim.py -ev {ev} -s {s} -e {e} > /dev/null 2>&1 &"
        cmd = cmd.format(ev=evnt.strftime("%Y-%m-%dT%H:%M"), s=start.strftime("%Y-%m-%dT%H:%M"), e=end.strftime("%Y-%m-%dT%H:%M"))
        print(cmd)
        os.system(cmd)
