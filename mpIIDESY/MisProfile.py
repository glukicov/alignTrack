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

def func(x, a, b, c):
	return a * np.exp(-b * x) + c

#First run the programme for specified number of tracks, and get the running constant
runs = (0, 1, 2)
mis1=(300, 9, 1)
mis2=(-300, -8, 1)

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
plt.savefig("MisProfile.png")
print("File produced: TimeTrend.png")




