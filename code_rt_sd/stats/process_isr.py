import os

for i in range(28):
    print("python capture_stats.py -n %d"%i)
    os.system("python capture_stats.py -n %d"%i)
