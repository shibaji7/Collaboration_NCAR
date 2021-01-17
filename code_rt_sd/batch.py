import os
import sys
import datetime as dt
import pandas as pd
import numpy as np

if len(sys.argv) > 1:
    ind = int(sys.argv[1])
    events = pd.read_csv("config/radar_event_list.csv", parse_dates=["date"])
    event = events.iloc[ind]
    evnt = event["date"]
    start = evnt - dt.timedelta(minutes=20)
    end = evnt + dt.timedelta(minutes=41)
    if end.day != start.day: end = end.replace(minute=0)
    cmd = "nohup python model_sim.py -r {r} -ev {ev} -s {s} -e {e} -ps > /dev/null 2>&1 &"
    cmd = cmd.format(r=event["rad"], ev=evnt.strftime("%Y-%m-%dT%H:%M"),
            s=start.strftime("%Y-%m-%dT%H:%M"), e=end.strftime("%Y-%m-%dT%H:%M"))
    print(cmd)
else:
    events = pd.read_csv("config/radar_event_list.csv", parse_dates=["date"])
    if os.path.exists("config/radar_event_list_proc.csv"): L = len(pd.read_csv("config/radar_event_list_proc.csv"))
    else: L = 0
    event = events.iloc[L]
    evnt = event["date"]
    start = evnt - dt.timedelta(minutes=20)
    end = evnt + dt.timedelta(minutes=41)
    if end.day != start.day: end = end.replace(minute=0)
    cmd = "nohup python model_sim.py -r {r} -ev {ev} -s {s} -e {e} > /dev/null 2>&1 &"
    cmd = cmd.format(r=event["rad"], ev=evnt.strftime("%Y-%m-%dT%H:%M"), 
            s=start.strftime("%Y-%m-%dT%H:%M"), e=end.strftime("%Y-%m-%dT%H:%M"))
    print(cmd)
    os.system(cmd)
    result = [{"date": evnt, "rad": event["rad"]}]
    o = pd.DataFrame.from_records(result)
    if not os.path.exists("config/radar_event_list_proc.csv"): o.to_csv("config/radar_event_list_proc.csv", header=True, index=False)
    else: o.to_csv("config/radar_event_list_proc.csv", mode="a", header=False, index=False)
