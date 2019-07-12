#########################################################################
# Monitoring plot for Tracker Alignment: for many g-2 runs  
#
# Expects to have offsets written to FHICLs from running RobustTrackerLaunch 
#  on alignment data 
# 
# INPUT: Top Directory 
# Directory structure: e.g. MONITOR/15922/ 
#         where 15922 is one of many aligned runs for one of the two stations
#         e.g. for S12 folder has the OffsetsPerModuleS12.fcl
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
parser.add_argument('-s12', '--dirS12', help='top directory', type=str) 
parser.add_argument('-s18', '--dirS18', help='top directory', type=str) 
args = parser.parse_args()
dirS12=args.dirS12
dirS18=args.dirS18

#Define constants 
stationPath = (dirS12, dirS18)
stationN = len(stationPath)
stations=("S12", "S18")
moduleN = 8 # per station 
fileName = "OffsetsPerModule"
fileExt = ".fcl"
FHICLServicePath = "services.Geometry.strawtracker." #full path of the tracking FHICL patch 
FHICLPatchName = ["strawModuleRShift", "strawModuleHShift"]
GlobalParNames = ["Radial", "Vertical"]
globalN=len(FHICLPatchName)

###Define DS and run ranges as a dict
data_sets =  {
  "60h": [15921, 15992],
  "9D": [16355, 16539],
  "End Game": [16908, 17528],
  "High Kick" : [16110, 16256],
  "Low Kick" : [16669, 16714],
  "Serenity Now": [24376, 24462],
  "Lazarus": [24575, 24688] 
}

#Get list of runs as the subdirs (removing the top dir name by spitting)
p = Path(dirS12)
runs = [int(str(x).rsplit('/', 1)[-1]) for x in p.iterdir() if x.is_dir()]
runN = len(runs)
print("Monitoring", runN,"runs:",runs)

#Get unique data sets names 
unique_data_sets=set()
for run in runs:
    for data_set in data_sets:
        if (run > data_sets[data_set][0] and run < data_sets[data_set][1]):
            unique_data_sets.add(data_set)
print("In the following data sets:", unique_data_sets)   

#Containers to fill with offsets per dof per station per run per module 
offsets = [[[[0 for i_module in range(0, moduleN)] for i_run in range(0, runN) ] for i_station in range(0, stationN)] for i_global in range(0, globalN) ]

#Get offsets from files per run for all modules in that run
for i_global in range(0, globalN):
    for i_station in range(0, stationN):
        for i_run in range(0, runN):
            path = stationPath[i_station]+"/"+str(runs[i_run])+"/"
            filePath=path+"/"+fileName+stations[i_station]+fileExt
            try:
                file = open (filePath, "r")
                offset_input = RTL.getOffsets(file, FHICLServicePath+FHICLPatchName[i_global]+stations[i_station][1:3])
                offsets[i_global][i_station][i_run]=offset_input
            except Exception:
                continue

#Containers to fill with offsets per dof per station per module per run
offsets_run = [[[[0 for i_run in range(0, runN)] for i_module in range(0, moduleN)] for i_station in range(0, stationN)] for i_global in range(0, globalN) ]

for i_global in range(0, globalN):
    for i_station in range(0, stationN):
        for i_module in range(0, moduleN):
            for i_run in range(0, runN):
                offsets_run[i_global][i_station][i_module][i_run]=offsets[i_global][i_station][i_run][i_module]

print(offsets_run)
##### Plotting

#First plot: nominal values
#Second plot: relative to the first run in monitor
stateN=2
colors = ["red", "green", "blue", "purple", "orange", "black", "brown", "grey"]
yTitle = [r"Alignment [$\mathrm{\mu}$m]", r"$\Delta$ Alignment [$\mathrm{\mu}$m]"]
plotName = ["", "_relative"]

for i_state in range(0, stateN):
    f = plt.figure(figsize=(7,int(globalN*2+1)))
    #Make subplot for each result 
    i_total=1
    for i_global in range(0, globalN):
        for i_station in range(0, stationN):
            plt.subplot( int( str(globalN)+str(stationN)+str(i_total)) )
            axes = plt.gca()
            plt.title(GlobalParNames[i_global]+" alignment in "+stations[i_station], fontsize=10)
            plt.ylabel(yTitle[i_state], fontsize=10)
            plt.xlabel("Run #", fontsize=10)
            plt.xticks(fontsize=8, rotation=0) 
            plt.yticks(fontsize=8, rotation=0)
            plt.minorticks_on()
            axes.tick_params(axis='x',which='minor',bottom=False)
            axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
            plt.tight_layout()
            #Plot data per module
            for i_module in range(0, moduleN):
                if (i_state==0):
                    plt.plot(runs, offsets_run[i_global][i_station][i_module], marker=".", color=colors[i_module], linewidth=1, label="M"+str(i_module+1))
                    yMax=np.max(offsets_run[i_global][i_station])
                    yMin=np.min(offsets_run[i_global][i_station])
                if (i_state==1):
                    plt.plot(runs, offsets_run[i_global][i_station][i_module]-offsets_run[i_global][i_station][i_module][0], marker=".", color=colors[i_module], linewidth=1, label="M"+str(i_module+1))
                    yMax=10
                    yMin=-10
            axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 
            i_total+=1
            #Set y limits to add DS labels 
            axes.set_ylim(yMin*1.5, yMax*1.1)
            #Add dataset divisions:
            for data_set in data_sets:
                if data_set in unique_data_sets:
                    line = [[ data_sets[data_set][0], yMin*1.5], [ data_sets[data_set][0],yMax*1.1]]
                    plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green', linewidth=2, linestyle=":")
                    line = [[ data_sets[data_set][1], yMin*1.5], [ data_sets[data_set][1], yMax*1.1]]
                    plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green', linewidth=2, linestyle=":")
                    middle_range = (data_sets[data_set][1]+data_sets[data_set][0])/2
                    plt.text(middle_range, yMin*1.36, data_set, fontsize=8, color='green', style='italic')
                    axes.annotate("",xy=(data_sets[data_set][0], yMin*1.4), xytext=(middle_range+middle_range*0.001, yMin*1.4), arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color="green") )
                    axes.annotate("",xy=(data_sets[data_set][1], yMin*1.4), xytext=(middle_range, yMin*1.4), arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color="green") )
    plt.savefig("Monitoring"+plotName[i_state]+".png", dpi=250)