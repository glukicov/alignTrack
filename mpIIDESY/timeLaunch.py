#!/usr/bin/env python

import sys 
import argparse
import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
import time
from scipy.optimize import curve_fit
import subprocess


#First run the programme for specified number of tracks, and get the running constant
tracksCut=(1000, 3000, 8000, 12000, 19000, 29000, 40000, 50000, 65000, 70000, 80000, 90000, 100800)

#move previous time-file
if ( os.path.isfile("Tracker_time.txt") ):
	os.rename("Tracker_time.txt", "Tracker_time.txt.BK")

#Run AlignTrack
for i in range(0, len(tracksCut)):
	subprocess.call(["./AlignTracker", "n" , str(tracksCut[i]) ])
	time.sleep(9.0)  

timeDelayCut=[]

with open("Tracker_time.txt") as f:
	for line in f:  #Line is a string
		number_str = line.split()
		timeDelayCut.append(float(number_str[0]))

slope, intercept = np.polyfit(tracksCut, timeDelayCut, 1)
abline_values = [slope * i + intercept for i in tracksCut]

plt.figure(2)
plt.scatter(tracksCut,timeDelayCut)
plt.title('Tracks vs Time - DCA Cut; slope= %.7f intercept= %.7f' %(slope, intercept) )
axes = plt.gca()
plt.plot(tracksCut, abline_values, 'b')
plt.ylabel("timeDelay [s]")
plt.xlabel("# Tracks [requested; accepted=~10%]")
print "Slope= ", slope, "intercept= ", intercept
if ( os.path.isfile("TimeConstants.txt") ):
	os.rename("TimeConstants.txt", "TimeConstants.txt.BK")

f = open('TimeConstants.txt', 'w')
f.write(str(slope))
f.write(" ")
f.write(str(intercept))  
f.close()  # you can omit in most cases as the destructor will call it
plt.savefig("TimeTrend.png")
print("File produced: TimeTrend.png")




