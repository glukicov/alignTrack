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
import datetime # print current time 
from matplotlib.gridspec import GridSpec # for broken axis 
from brokenaxes import brokenaxes # brake x-axis for DS 
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

def make_plot(colors, colorsS, plotName, font_size, plot_dpi, y_Max, n_rows, globalN, data_sets, share_x, main_title):
    fig, axes = plt.subplots(nrows=n_rows, sharex=share_x, figsize=(10, 8)) # n_rows subplots, with shared axis if passed 
    if(globalN == 1):
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) #e-5 notation 
        axes=[axes] # make an array to use indexing 
    #Make subplot for each result 
    i_total=0 # count axis 
    for i_global in range(0, globalN):  
        for i_station in range(0, stationN):
            if (globalN == 2): # if alignment 
                axes[i_total].set_ylabel(GlobalParNames[i_global]+" in "+stations[i_station]+r" [$\mathrm{\mu}$m]", fontsize=font_size-1)
            if (globalN == 1):
                plt.ylabel("Tracks", fontsize=font_size)
            plt.suptitle(main_title, fontsize=font_size) #title for all subplots 
            plt.xlabel("Run #", fontsize=font_size)
            plt.xticks(fontsize=font_size, rotation=0) 
            plt.minorticks_on()
            axes[i_total].tick_params(axis='x',which='minor',bottom=False, labelsize=font_size)
            axes[i_total].tick_params(axis='y', which='both', left=True, right=True, direction='inout', labelsize=font_size)
            plt.tight_layout()
            fig.subplots_adjust(top=0.95) # leave space for the title of alll 
            if(share_x):
                plt.subplots_adjust(wspace=0, hspace=0) # merge x-axis 
            #Plot data per module
            yMin, yMax = None, None # scope
            if (globalN == 2):
                for i_module in range(0, moduleN):
                    #absolute alignment 
                    if (plotName==""):
                        axes[i_total].plot(list(runs), offsets_run[i_global][i_station][i_module], marker=".", color=colors[i_module], linewidth=1, label="M"+str(i_module+1))
                        yMax=np.max(offsets_run[i_global][i_station])
                        yMin=np.min(offsets_run[i_global][i_station])
                    if (plotName=="_relative"):
                        axes[i_total].plot(list(runs), offsets_run[i_global][i_station][i_module]-np.mean(offsets_run[i_global][i_station][i_module]), marker=".", color=colors[i_module], linewidth=0, label="M"+str(i_module+1))
                        yMax=np.max(offsets_run[i_global][i_station])
                        yMin=np.min(offsets_run[i_global][i_station])
                #Set y limits to add DS labels
                axes[i_total].set_ylim(yMin*1.5, yMax*1.1)
            if (globalN == 1):
                plt.plot(list(runs), trackN[i_station], marker=".", color=colorsS[i_station], linewidth=1, label=stations[i_station] )
                if (i_station==0):
                    yMax=np.max(trackN[i_station])
                    yMin=np.min(trackN[i_station])*0.3 

            #add legend at the end 
            axes[i_total].legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 9}) # outside (R) of the plot 
            #Add dataset divisions:
            if (globalN==1 and i_station==1):
                        pass # do just one for tracks plot 
            #draw DS lines
            else:
                for data_set in data_sets:
                    if data_set in unique_data_sets:
                        line = [[ data_sets[data_set][0], yMin*1.5], [ data_sets[data_set][0],yMax*1.1]]
                        axes[i_total].plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green', linewidth=2, linestyle=":")
                        line = [[ data_sets[data_set][1], yMin*1.5], [ data_sets[data_set][1], yMax*1.1]]
                        axes[i_total].plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green', linewidth=2, linestyle=":")
                        middle_range = (data_sets[data_set][1]+data_sets[data_set][0])/2
                        axes[i_total].text(middle_range, yMin*1.36, data_set, fontsize=font_size, color='green', style='italic')
                        axes[i_total].annotate("",xy=(data_sets[data_set][0], yMin*1.4), xytext=(middle_range+middle_range*0.001, yMin*1.4), arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color="green") )
                        axes[i_total].annotate("",xy=(data_sets[data_set][1], yMin*1.4), xytext=(middle_range, yMin*1.4), arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color="green") ) 
            if (globalN==2):
                i_total+=1 # increment axes counter for alignment plots 
    fig.align_ylabels(axes[:])
    plt.savefig("Monitoring"+plotName+".png", dpi=plot_dpi)

#arguments 
parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-s12', '--dirS12', help='top directory', type=str) 
parser.add_argument('-s18', '--dirS18', help='top directory', type=str) 
parser.add_argument('-trackCut', '--trackCut', help='top directory', type=int, default=300000) 
args = parser.parse_args()
dirS12=args.dirS12
dirS18=args.dirS18
trackCut=args.trackCut

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

#Welcome message -> Logger 
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
runs=sorted(runs)
runN=len(runs)
print("Using the following",runN,"runs for plotting (ok in both stations): ", runs)

#Get unique data sets names 
unique_data_sets=set()
for run in runs:
    for data_set in data_sets:
        if (run > data_sets[data_set][0] and run < data_sets[data_set][1]):
            unique_data_sets.add(data_set)
print("In the following data sets:", unique_data_sets) 

selected_runs=[]
#Get and cut on the number of tracks 
trackN=[ [], [] ] #Get track number per run
for i_run, run in enumerate(runs):
    path_s12 = stationPath[0]+"/"+str(run)+"/"+"millepede.log"
    path_s18 = stationPath[1]+"/"+str(run)+"/"+"millepede.log"
    file_s12 = open (path_s12, "r")
    file_s18 = open (path_s18, "r")
    tracks = [getTracks(file_s12, "\s+NREC\s+=\s+"), getTracks(file_s18, "\s+NREC\s+=\s+")]
    if (tracks[0] > int(trackCut) and tracks[1] > int(trackCut)):
        trackN[0].append(tracks[0])
        trackN[1].append(tracks[1])
        selected_runs.append(run)

runs=selected_runs
runN=len(runs)
print("Using the following",runN,"runs for plotting (after the track cut of",trackCut,"tracks): ", runs)

#Containers to fill with offsets per dof per station per run per module 
offsets = [[[[0 for i_module in range(0, moduleN)] for i_run in range(0, runN) ] for i_station in range(0, stationN)] for i_global in range(0, globalN) ]

#Get offsets from files per run for all modules in that run
for i_global in range(0, globalN):
    for i_station in range(0, stationN):
        for i_run, run in enumerate(runs):
            path = stationPath[i_station]+"/"+str(run)+"/"
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
#First plot: nominal values
#Second plot: relative to the first run in monitor
#Third plot: track number per run 

#Define plotting constants 
stateN=3 # 3 plots 
colors = ["red", "green", "blue", "purple", "orange", "black", "brown", "grey"]
colorsS=["red", "blue"]
font_size=12
plot_dpi=250 
y_Max = 12 # for relative plot 
share_x=True

#Pass as a function of i_state:
n_rows_array = [4, 4, 1]
global_n = [2, 2, 1]
plotName = ["", "_relative", "_tracks"]
main_title = ["Alignment per run", r"Alignment-$\langle$Alignment$\rangle$ per run", "Tracks used in alignment"]

#Produce the 3 plots: absolute, relative and trackN
for i_state in range(stateN):
   make_plot(colors, colorsS, plotName[i_state], font_size, plot_dpi, y_Max, n_rows_array[i_state], global_n[i_state], data_sets, share_x, main_title[i_state])

print("Monitoring plots are ready!")


fig = plt.figure(figsize=(5,5))
sps1, sps2 = GridSpec(2,1)
bax = brokenaxes(xlims=((.1, .3),(.7, .8)), subplot_spec=sps1)
x = np.linspace(0, 1, 100)
bax.plot(x, np.sin(x*30), ls=':', color='m')
x = np.random.poisson(3, 1000)
bax = brokenaxes(xlims=((0, 2.5), (3, 6)), subplot_spec=sps2)
bax.hist(x, histtype='bar')
plt.savefig("Broken.png", dpi=plot_dpi)