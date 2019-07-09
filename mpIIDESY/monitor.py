#########################################################################
# Monitoring plot for Tracker Alignment: for many g-2 runs  
#
# Expects to have offsets written to FHICLs from running RobustTrackerLaunch 
#  on alignment data 
# 
# INPUT: Top Directory 
# Directory structure: e.g. MONITOR/15922/S12 
#   where dir=MONITOR (top directory), containing many runs 
#         15922 is one of many runs 
#         S12 folder has the OffsetsPerModuleS12.fcl
#         S17 folder has the OffsetsPerModuleS18.fcl 
#
# Created: 9 July 2019 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 9 July 2019 by Gleb
#######################################################################

import sys, os # print out to terminal 
from pathlib import Path # list dirs
import numpy as np # arrays 
import argparse # command line inputs sub
import matplotlib.pyplot as plt #for plotting  
import itertools # smart lines in plotting 
import RobustTrackerLaunch as RTL # importing a function; this code
# must have an "if main/def main" protection

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-d', '--dir', help='top directory', default="MONITOR", type=str) 
args = parser.parse_args()
top_dir=args.dir

#Define constants 
stations = ["S12", "S18"]
stationN = len(stations)
moduleN = 8 # per station 
fileName = "OffsetsPerModule"
fileExt = ".fcl"
FHICLServicePath = "services.Geometry.strawtracker." #full path of the tracking FHICL patch 
FHICLPatchName = ["strawModuleRShift", "strawModuleHShift"]
GlobalParNames = ["Radial", "Vertical"]
globalN=len(FHICLPatchName)

#Get list of runs as the subdirs (removing the top dir name by spitting)
p = Path(top_dir)
runs = [int(str(x).rsplit('/', 1)[-1]) for x in p.iterdir() if x.is_dir()]
runN = len(runs)
print("Monitoring", runN,"runs:",runs)

#Containers to fill with offsets per dof per station per run per module 
offsets = [[[[0 for i_module in range(0, moduleN)] for i_run in range(0, runN) ] for i_station in range(0, stationN)] for i_global in range(0, globalN) ]

#Get offsets from files per run for all modules in that run
for i_global in range(0, globalN):
    for i_station in range(0, stationN):
        for i_run in range(0, runN):
            path = top_dir+"/"+str(runs[i_run])+"/"+stations[i_station]
            filePath=path+"/"+fileName+stations[i_station]+fileExt
            file = open (filePath, "r")
            offset_input = RTL.getOffsets(file, FHICLServicePath+FHICLPatchName[i_global]+stations[i_station][1:3])
            offsets[i_global][i_station][i_run]=offset_input

#Containers to fill with offsets per dof per station per module per run
offsets_run = [[[[0 for i_run in range(0, runN)] for i_module in range(0, moduleN)] for i_station in range(0, stationN)] for i_global in range(0, globalN) ]

for i_global in range(0, globalN):
    for i_station in range(0, stationN):
        for i_module in range(0, moduleN):
            for i_run in range(0, runN):
                offsets_run[i_global][i_station][i_module][i_run]=offsets[i_global][i_station][i_run][i_module]

##### Plotting
colors = ["red", "green", "blue", "purple", "orange", "black", "brown", "grey"]
f = plt.figure(figsize=(7,int(globalN*2+1)))
yMax = [120, 120]
yMin = [-120, -220]
#Make subplot for each result 
i_total=1
for i_global in range(0, globalN):
    for i_station in range(0, stationN):
        plt.subplot( int( str(globalN)+str(stationN)+str(i_total)) )
        axes = plt.gca()
        plt.title(GlobalParNames[i_global]+" alignment in "+stations[i_station], fontsize=12)
        plt.ylabel(r"Alignment [$\mathrm{\mu}$m]")
        plt.xlabel("Run #", fontsize=12)
        plt.xticks(fontsize=10, rotation=0) 
        plt.yticks(fontsize=10, rotation=0)
        plt.minorticks_on()
        axes.tick_params(axis='x',which='minor',bottom=False)
        axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
        plt.tight_layout()
        #Plot data per module
        for i_module in range(0, moduleN):
            plt.plot(runs, offsets_run[i_global][i_station][i_module], marker=".", color=colors[i_module], linewidth=1, label="M"+str(i_module+1))
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 
        i_total+=1
plt.savefig("Monitoring.png", dpi=250)