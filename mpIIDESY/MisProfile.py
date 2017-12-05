#!/usr/bin/env python

import sys 
import argparse
import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
import time
import itertools
from scipy.optimize import curve_fit
import subprocess
from matplotlib.ticker import MaxNLocator

#First run the programme for specified number of tracks, and get the running constant
runs = (1.0, 1.1, 1.2, 2.0, 2.1, 2.2, 3.0, 3.1, 3.2, 4.0, 4.1, 4.2, 5.0, 5.1, 5.2)
chi1=(1.225475, 0.977459, 0.977459, 0.976604, 0.976587, 0.976587, 0.989997, 0.989964, 0.989964, 0.981541, 0.981720, 0.981720, 0.986597, 0.986591, 0.986591)

fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.set_xlim(runs[0]-0.2,runs[-1]+0.2)
ax1.set_ylim(0.974, 0.999)
ax1 .scatter(runs,chi1, color='red')
plt.plot([runs[0]-0.2, runs[-1]+0.2], [1, 1], color="green")

plt.title("Alignment Convergence with Run #")
axes = plt.gca()
plt.ylabel("$\chi^2 $")
plt.xlabel("Run #.Iteration #")
plt.savefig("MisProfile_1.png")


#First run the programme for specified number of tracks, and get the running constant
runs =  (1, 2, 3, 4, 5)
mis1Raw=(5.42, 2.10, 0.82, -4.095, 1.415)
mis2Raw=(-5.70, -2.73, 2.76, 1.72, -1.512)

mis1=[]
mis2=[]

for i in range(0, len(mis1Raw)):
	if (i==0):
		mis1.append(mis1Raw[0])
		mis2.append(mis2Raw[0])
	else:
		mis1.append( mis1[i-1] + mis1Raw[i] )
		mis2.append( mis2[i-1] + mis2Raw[i] )

print mis1
print mis2

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.set_xlim(runs[0]-0.2,runs[-1]+0.2)
ax1 .scatter(runs,mis1, color='red')
ax1 .scatter(runs,mis2, color='blue')	
plt.plot([runs[0]-0.2, runs[-1]+0.2], [0, 0], color="green")
plt.title("Alignment Convergence with Run #")
axes = plt.gca()
plt.ylabel("Misalignment [um]")
plt.xlabel("Run #")
plt.savefig("MisProfile_2.png")
print("File produced: MisProfile_1.png and MisProfile_2.png")





