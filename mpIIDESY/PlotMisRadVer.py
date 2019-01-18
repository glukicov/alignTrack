#!/usr/bin/python

####################################################################
# Plot Radial and Vertical Misalignment given the FHICL file 
#####################################################################
import argparse, sys, os

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-f', '--file', help='path')
args = parser.parse_args()
file = args.file

import matplotlib.pyplot as plt #for plotting 
import matplotlib.ticker as ticker
import numpy as np  # smart arrays 
import itertools # smart lines
import subprocess
import pandas as pd
import re

# CONSTANTS 
moduleN = 8 #number of movable detectors
stationN = 3 # S0, S12, S18 
globalN = 2 # X, Y
parN = stationN * moduleN

def getOffsets(f, name):
	offsets = [] #tmp storage buffer 
	for line in f:
	    if re.match(name, line):
	        copy = True
	        offsets=line 

	    else:
	        copy = False

	return offsets   
		

###Store the offsets from file
f=open(file, "r")
radial = (getOffsets(f, "services.Geometry.strawtracker.radialOffsetPerModule:"))
radial = radial.replace("services.Geometry.strawtracker.radialOffsetPerModule: [", " ") 
radial = radial.replace("]", "") 
radialOff = np.array([float(r) for r in radial.split(',')])
radialOff=radialOff[8:16] 
radialOff=radialOff*1e3 # mm -> um 
print("Radial:", radialOff)

f=open(file, "r")
vertical = (getOffsets(f, "services.Geometry.strawtracker.verticalOffsetPerModule:"))
vertical = vertical .replace("services.Geometry.strawtracker.verticalOffsetPerModule: [", "") 
vertical = vertical .replace("]", "") 
verticalOff = np.array([float(v) for v in vertical .split(',')])
verticalOff=verticalOff[8:16] 
verticalOff=verticalOff*1e3 # mm -> um 
print("Vertical:",  verticalOff)


###Get a nice label
extraLabel = file.replace("/gm2/app/users/glukicov/TrackerAlignment/gm2Dev_v9_15_00/srcs/gm2tracker/align/Systematics_new/", " ") # replace "/" with "_"
extraLabel = extraLabel.replace("/", "_") # replace "/" with "_"
extraLabel = extraLabel.replace("_RunTrackingDAQ.fcl", "") # replace "/" with "_"

print("Label:", extraLabel)


colours = ["green", "blue", "black", "orange", "purple"]
yMin = -250
yMax = 250
plt.subplot(211) # X 
plt.rcParams.update({'font.size': 10})
axes = plt.gca()
axes.set_xlim(0.5, 8.5)
axes.set_ylim(yMin, yMax)
plt.title("Misalignment: Radial "+str(extraLabel), fontsize=10)
plt.ylabel("Misalignment [um]", fontsize=10)
line = [[0.5,0.0], [8.5, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
for i_module in range(0, 8):
	plt.plot(i_module+1, radialOff[i_module], marker=".", color="red")
	line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')

#Legend (stats X)
avgMeanMis = sum(radialOff)/float(len(radialOff))
SDMis = np.std(radialOff)
textstr = '<Truth>=%s um\nSD Truth=%s um \n'%(int(round(avgMeanMis)), int(round(SDMis)))
plt.text(4.7, 80, textstr, fontsize=10, color="black")
plt.xlabel("Module", fontsize=10)

plt.subplot(212) # Y 
axes = plt.gca()
axes.set_xlim(0.5, 8.5)
axes.set_ylim(yMin, yMax)
plt.title("Misalignment: Vertical "+str(extraLabel), fontsize=10)
plt.ylabel("Misalignment [um]", fontsize=10)
line = [[0.5,0.0], [8.5, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
for i_module in range(0, 8):
	plt.plot(i_module+1, verticalOff[i_module], marker=".", color="red")
	line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	

#Legend (stats Y)
avgMeanMis = sum(verticalOff)/float(len(verticalOff))
SDMis = np.std(verticalOff)
textstr = '<Truth>=%s um\nSD Truth=%s um \n'%(int(round(avgMeanMis)), int(round(SDMis)))
plt.text(4.7, 80, textstr, fontsize=10, color="black")
plt.xlabel("Module", fontsize=10)

plt.tight_layout()

plt.savefig("Mis"+str(extraLabel)+".png", dpi=600)