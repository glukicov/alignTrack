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
runs = (1.0, 1.1, 1.2, 2.0, 2.1, 2.2)
chi1=(1.24057, 0.99261, 0.99261, 0.99168, 0.99158, 0.99158)

fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.set_xlim(runs[0]-0.2,runs[-1]+0.2)
ax1.set_ylim(0.991,0.993)
ax1 .scatter(runs,chi1, color='red')
plt.plot([runs[0]-0.2, runs[-1]+0.2], [1, 1], color="green")

plt.title("Alignment Convergence with Run #")
axes = plt.gca()
plt.ylabel("$\chi^2 $")
plt.xlabel("Run #.Iteration #")
plt.savefig("MisProfile_1.png")


#First run the programme for specified number of tracks, and get the running constant
runs =  (1, 2)
mis1=(9, 1)
mis2=(-8, 1)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.set_xlim(runs[0]-0.2,runs[-1]+0.2)
ax1 .scatter(runs,mis1, color='red')
ax1 .scatter(runs,mis2, color='blue')	
plt.plot([runs[0]-0.2, runs[-1]+0.2], [0, 0], color="green")

popt1, pcov1 = curve_fit(func, runs, mis1)
plt.plot(func(runs, *popt1))

plt.title("Alignment Convergence with Run #")
axes = plt.gca()
plt.ylabel("Misalignment [um]")
plt.xlabel("Run #")
plt.savefig("MisProfile_2.png")
print("File produced: MisProfile_1.png and MisProfile_2.png")





