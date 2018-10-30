#!/usr/bin/env python

#Plotter

from ROOT import *
import matplotlib.pyplot as plt #for plotting 
import numpy as np  # smart arrays 
import itertools # smart lines 
import argparse, sys
from math import log10, floor
from matplotlib.ticker import MaxNLocator
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
import subprocess
def round_sig(x, sig=2):
	return round(x, sig-int(floor(log10(abs(x))))-1)

selectDirNames = ( "-50 < <y> < -30", "-30 < <y> < -10", "-10 < <y> < 10" , "10 < <y> < 30")
wireGroups = ("S0-S15", "S16-S31")
colour = ("green", "purple", "orange")


f = TFile.Open('TrackerAlignment.root')
if f:
    print("is open")
else:
    print("Not found")


ticks=[]

yMin = -0.02
yMax = 0.02
plt.figure(1)
axes = plt.gca()
for i_wireGroup in range(0, len(wireGroups)):
	colourDirs=colour[i_wireGroup]
	wireGroupName = wireGroups[i_wireGroup]
	for i_dirNames in range(0, len(selectDirNames)):
		dirName = selectDirNames[i_dirNames]
		name = "TrackerAlignment/"+str(dirName)+"/Residual "+str(wireGroupName)+" "+str(dirName)
		print name
		t = f.Get(str(name))
		mean=t.GetMean()
		meanError=t.GetMeanError()
		tick=i_dirNames+1
		if (i_wireGroup==0):
			ticks.append(tick)
		
		plt.plot(tick, mean, marker="*", color=str(colourDirs))
		plt.errorbar(tick, mean, yerr=meanError, color=str(colourDirs))	

		line = [[i_dirNames+0.5,yMin], [i_dirNames+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'black')




axes.set_xlim(0.5, 4.5)
axes.set_ylim(yMin, yMax)
plt.xticks(ticks, selectDirNames)
plt.ylabel("Residual Mean [error = Error on the Mean]", fontsize=10)
plt.xlabel("Vertical Group [mm]", fontsize=10)
plt.savefig("Vertical_Group.png")

