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


x = [-6880.76, -6873.53, -6863.74, -6851.36, -6836.33, -6818.78, -6798.7,-6776.11  ,291.97    ,429.41    ,566.13    ,702.13    ,837.39    ,971.94    ,1105.79   ,1238.93 ] 

z= [-291.97, -429.41, -566.13, -702.13, -837.39, -971.94, -1105.79, -1238.93, -6880.76, -6873.53, -6863.74, -6851.36, -6836.33, -6818.78, -6798.7, -6776.11]


yMin = -7600
yMax = 7600
plt.figure(1)
axes = plt.gca()
for i in range(0, len(x)):
	plt.plot(z[i], x[i], marker="*", color="red", markersize=20)

ring = plt.Circle((0, 0), 7000, color='blue')
axes.add_artist(ring)

axes.set_xlim(yMin, yMax)
axes.set_ylim(yMin, yMax)
plt.xlabel("Z [mm]", fontsize=20)
plt.ylabel("X [mm]", fontsize=20)
plt.savefig("Ring.png")