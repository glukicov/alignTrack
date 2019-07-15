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

import sys, os, re # print out to terminal 
from glob import glob
import numpy as np # arrays 
import argparse # command line inputs sub
import matplotlib.pyplot as plt #for plotting  
import itertools # smart lines in plotting 
import RobustTrackerLaunch as RTL # importing a function; this code
import datetime
# must have an "if main/def main" protection

### HELPER FUNCTIONS ######
## Read track number from the PEDE log file 
def getTracks(f, name):
    tracks = [] #tmp storage buffer 
    for line in f:
        if re.match(name, line):
            tracks=line 
        else:
            pass
    tracks = tracks.replace("NREC =", "") 
    tracks = tracks.replace("= number of records", "")
    return int(tracks)

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

### Duplicate all cout into a log file
class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("Monitoring.log", "w+")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        self.terminal.flush()
        self.log.flush()

sys.stdout = Logger()


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

#Get list of runs as the subdirs (removing the top dir name by splitting)
station_runs=[[] ,[]] # all runs 
station_runs_bool=[[], []] # runs with results (use as boolean mask)
runs=set()  # final runs for both stations 

print("Summary monitoring plots for S12 and S18:", dirS12, dirS18)
print("Starting on:", datetime.datetime.now())

#quickly check if we get the .fcl results file for that run, and remove run otherwise 
for i_station, i_path in enumerate(stationPath): 
    all_runs_station = [ int( os.path.join(i_path, o).rsplit('/', 1)[-1] ) for o in os.listdir(i_path) if os.path.isdir(os.path.join(i_path,o))]
    print("Found", len(all_runs_station),"runs:",all_runs_station,"for station", stations[i_station])
    for i_run in all_runs_station:
            path = i_path+"/"+str(i_run)+"/"
            filePath=path+"/"+fileName+stations[i_station]+fileExt
            #append to run array if we have result for that run 
            if (os.path.exists(filePath)): 
                station_runs[i_station].append(i_run)
            else:
                print("No alignment data found for run", i_run, stations[i_station])

# Log the failed and non-failed runs 
for i_station in range(0, stationN):
    print(stations[i_station],":", len(station_runs[i_station]), "runs have alignment data: ", station_runs[i_station])

#Finally form the set of runs for plotting (ok in both stations)
runs = set(station_runs[0]) & set(station_runs[1])
runN=len(runs)
print("Using the following",runN,"runs for plotting (ok in both stations): ", runs)

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
trackN=[ [], [] ] #Get track number per run 
for i_global in range(0, globalN):
    for i_station in range(0, stationN):
        for i_run, run in enumerate(runs):
            path = stationPath[i_station]+"/"+str(run)+"/"
            filePath=path+"/"+fileName+stations[i_station]+fileExt
            file = open (filePath, "r")
            offset_input = RTL.getOffsets(file, FHICLServicePath+FHICLPatchName[i_global]+stations[i_station][1:3])
            offsets[i_global][i_station][i_run]=offset_input
            #append track numbers just once per global 
            if (i_global == 0):
                filePath=path+"/"+"millepede.log"
                file = open (filePath, "r")
                trackN[i_station].append(getTracks(file, "\s+NREC\s+=\s+"))

#Containers to fill with offsets per dof per station per module per run
offsets_run = [[[[0 for i_run in range(0, runN)] for i_module in range(0, moduleN)] for i_station in range(0, stationN)] for i_global in range(0, globalN) ]

for i_global in range(0, globalN):
    for i_station in range(0, stationN):
        for i_module in range(0, moduleN):
            for i_run in range(0, runN):
                offsets_run[i_global][i_station][i_module][i_run]=offsets[i_global][i_station][i_run][i_module]

##### Plotting
#First plot: nominal values
#Second plot: relative to the first run in monitor
#Third plot: track number per run 
stateN=3
colors = ["red", "green", "blue", "purple", "orange", "black", "brown", "grey"]
yTitle = [r"Alignment [$\mathrm{\mu}$m]", r"$\Delta$ Alignment [$\mathrm{\mu}$m]"]
plotName = ["", "_relative", "_tracks"]
colorsS=["red", "blue"]

for i_state in range(0, stateN):
    if (i_state == 0 or i_state == 1):
        f = plt.figure(figsize=(7,int(globalN*2+1)))
    if(i_state == 2):
         f = plt.figure(figsize=(7,7))
         globalN=1 # only 1 dimension for track N
         plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #Make subplot for each result 
    i_total=1
    for i_global in range(0, globalN):
        for i_station in range(0, stationN):
            if (i_state == 0 or i_state == 1 ):
                plt.subplot( int( str(globalN)+str(stationN)+str(i_total)) )
                plt.title(GlobalParNames[i_global]+" alignment in "+stations[i_station], fontsize=10)
                plt.ylabel(yTitle[i_state], fontsize=10)
            if (i_state == 2):
                plt.title("Tracks used in alignment", fontsize=10)
                plt.ylabel("Tracks", fontsize=10)
            axes = plt.gca()
            plt.xlabel("Run #", fontsize=10)
            plt.xticks(fontsize=8, rotation=0) 
            plt.yticks(fontsize=8, rotation=0)
            plt.minorticks_on()
            axes.tick_params(axis='x',which='minor',bottom=False)
            axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
            plt.tight_layout()
            #Plot data per module
            if (i_state == 0 or i_state == 1 ):
                for i_module in range(0, moduleN):
                    if (i_state==0):
                        plt.plot(list(runs), offsets_run[i_global][i_station][i_module], marker=".", color=colors[i_module], linewidth=1, label="M"+str(i_module+1))
                        yMax=np.max(offsets_run[i_global][i_station])
                        yMin=np.min(offsets_run[i_global][i_station])
                    if (i_state==1):
                        plt.plot(list(runs), offsets_run[i_global][i_station][i_module]-offsets_run[i_global][i_station][i_module][0], marker=".", color=colors[i_module], linewidth=1, label="M"+str(i_module+1))
                        yMax=10
                        yMin=-10
                #Set y limits to add DS labels 
                axes.set_ylim(yMin*1.5, yMax*1.1)
            if (i_state == 2):
                plt.plot(list(runs), trackN[i_station], marker=".", color=colorsS[i_station], linewidth=1, label=stations[i_station] )
            axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 
            i_total+=1
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



