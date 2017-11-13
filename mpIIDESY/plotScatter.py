#!/usr/bin/env python

import sys 
import argparse
import glob
import matplotlib.pyplot as plt
import numpy as np
import datetime
import time
from scipy.optimize import curve_fit


tracksNoCut=( 50, 250, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
timeDelayNoCut=( 1.0, 1.5, 1.8, 2.3, 3.0, 3.5, 4.5, 6.0, 7.0, 9.0, 10.0, 11.5, 12.5, 14.0 )

slope, intercept = np.polyfit(tracksNoCut, timeDelayNoCut, 1)
abline_values = [slope * i + intercept for i in tracksNoCut]

plt.figure(1)
plt.scatter(tracksNoCut,timeDelayNoCut)
plt.title('Tracks vs Time - NoCut; slope= %.7f intercept= %.7f' %(slope, intercept))
axes = plt.gca()
plt.plot(tracksNoCut, abline_values, 'b')
plt.ylabel("timeDelay [s]")
plt.xlabel("# Tracks")
plt.show()

tracksCut=( 500, 3000, 8000, 12000, 19000, 29000, 40000, 50000, 65000, 70000, 80000, 90000, 100800)
timeDelayCut=( 1.0, 1.5, 2.2, 2.8, 3.8, 5.2, 6.9, 8.0, 11.0, 11.8, 13.6, 14.5, 16.0)

slope, intercept = np.polyfit(tracksCut, timeDelayCut, 1)
abline_values = [slope * i + intercept for i in tracksCut]

plt.figure(2)
plt.scatter(tracksCut,timeDelayCut)
plt.title('Tracks vs Time - DCA Cut; slope= %.7f intercept= %.7f' %(slope, intercept) )
axes = plt.gca()
plt.plot(tracksCut, abline_values, 'b')
plt.ylabel("timeDelay [s]")
plt.xlabel("# Tracks [requested; accepted=~10%]")
plt.show()