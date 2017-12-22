#!/usr/bin/env python

import sys 
import argparse
import glob
import os
import datetime
import time
import subprocess

tracksCut=(1000, 3000, 8000, 12000, 19000, 29000, 40000, 50000, 65000, 70000, 80000, 90000, 100000)
#tracksCut=(1000, 3000, 8000, 12000)

subprocess.call(["clear"])

print "Staring PEDE tests with DCA Cut"

if ( os.path.isfile("PEDE_Mis.txt") ):
	os.rename("PEDE_Mis.txt", "PEDE_Mis.txt.BK")

slope = 0.0 
intercept = 0.0
with open("TimeConstants.txt") as f:
	for line in f:  #Line is a string
		number_str = line.split()
		slope = float(number_str[0])
		intercept = float(number_str[1])

for i in range(0, len(tracksCut)):
	subprocess.call(["./AlignTracker", "n" , str(tracksCut[i]), "0.0", "0.0"])
	
	subprocess.call(["./pede", "Tracker_str.txt" ])
	
	subprocess.call(["./ConcatenatePEDE.py"])
	

subprocess.call(["./TrackerLaunch.py"])
print "PEDE Tests Complete. "

