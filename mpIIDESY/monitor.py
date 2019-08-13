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
import subprocess, shlex 
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
from collections import OrderedDict
# must have an "if main/def main" protection

#arguments 
parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-s12', '--dirS12', help='top directory', type=str) 
parser.add_argument('-s18', '--dirS18', help='top directory', type=str) 
parser.add_argument('-trackCut', '--trackCut', help='top directory', type=int, default=300000)
parser.add_argument('-pedeCut', '--pedeCut', help='top directory', type=bool, default=True)
args = parser.parse_args()
dirS12=args.dirS12
dirS18=args.dirS18
trackCut=args.trackCut
pedeCut=args.pedeCut

#Define constants 
stationPath = (dirS12, dirS18)
stationN = len(stationPath)
stations=("S12", "S18")
tracker=("Tracker 1", "Tracker 2")
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
data_sets = OrderedDict() # keep the same order for Python2
data_sets["60h"] = [15921, 15992] 
data_sets["9D"] = [16355, 16539] 
data_sets["End Game"] = [16908, 17528] 
data_sets["High Kick" ] = [16110, 16256] 
data_sets["Low Kick" ] = [16669, 16714] 
data_sets["Serenity Now"] = [24376, 24462] 
data_sets["Lazarus"] = [24575, 24688]
data_sets["Sybil"] = [26053, 26384] 

#Get list of runs as the subdirs (removing the top dir name by splitting)
station_runs=[[] ,[]] # all runs 
station_runs_bool=[[], []] # runs with results (use as boolean mask)
runs=set()  # final runs for both stations 

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

def make_plot(colors, colorsS, plotName, font_size, plot_dpi, y_Max, n_rows, globalN, data_sets, main_title, broken=False, plot_show=True):
    
    if (broken is False):
        fig, axes = plt.subplots(nrows=n_rows, figsize=(10, 8)) # n_rows subplots, with shared axis if passed 
        fig.subplots_adjust(top=0.95) # leave space for the title of alll 
        plt.subplots_adjust(hspace=0.4) # move x-axis closer
        plt.tight_layout()
        font_size_broken = font_size-4
        
    if (broken is True):
        fig = plt.figure(figsize=(9, 8))
        sps1, sps2, sps3, sps4 = GridSpec(n_rows, 1)
        bax1, bax2, bax3, bax4 = brokenaxes(xlims=limits, subplot_spec=sps1), brokenaxes(xlims=limits, subplot_spec=sps2), brokenaxes(xlims=limits, subplot_spec=sps3), brokenaxes(xlims=limits, subplot_spec=sps4)
        axes=[bax1, bax2, bax3, bax4]
        font_size_broken = font_size-4
        # plt.subplots_adjust(hspace=0.3) # move x-axis closer

    #Plot decorations 
    plt.suptitle(main_title, fontsize=font_size, weight='bold') #title for all subplots
    plt.xticks(fontsize=font_size, rotation=0) 
    plt.minorticks_on()
    # plt.xlabel("Run #", fontsize=font_size)
    fig.align_ylabels()
    
    if(globalN == 1):
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) #e-5 notation 
        axes=[axes] # make an array to use indexing 
    #Make subplot for each result 
    i_total=0 # count axis
    for i_global in range(0, globalN):  
        for i_station in range(0, stationN):
            if (globalN == 2): # if alignment 
                axes[i_total].set_ylabel(GlobalParNames[i_global]+" in "+tracker[i_station]+r" [$\mathrm{\mu}$m]", fontsize=font_size_broken, labelpad=2)
                axes[i_total].set_xlabel("Run #", fontsize=font_size, labelpad=20)
            if (globalN == 1):
                plt.ylabel("Tracks", fontsize=font_size)

            #set appprox. y-error:
            if (i_total==0 or i_total==1):
                y_error=1
            if (i_total==2 or i_total==3):
                y_error=2
            
            #Axis decorations 
            axes[i_total].tick_params(axis='x',which='minor',bottom=False, labelsize=font_size)
            axes[i_total].tick_params(axis='y', which='both', left=True, right=True, direction='inout', labelsize=font_size)
            
            #Plot data per module
            yMin, yMax = None, None # scope
            if (globalN == 2):
                for i_module in range(0, moduleN):
                    #absolute alignment 
                    if (plotName==""):
                        axes[i_total].plot(list(selected_runs), offsets_run[i_global][i_station][i_module], marker=".", color=colors[i_module], linewidth=1, label="M"+str(i_module+1))
                        yMax=np.max(offsets_run[i_global][i_station])
                        yMin=np.min(offsets_run[i_global][i_station])
                    if (plotName=="_relative"):
                        axes[i_total].plot(list(selected_runs), offsets_run[i_global][i_station][i_module]-np.mean(offsets_run[i_global][i_station][i_module]), marker=".", color=colors[i_module], linewidth=0)
                        axes[i_total].errorbar(list(selected_runs), offsets_run[i_global][i_station][i_module]-np.mean(offsets_run[i_global][i_station][i_module]), yerr=y_error, marker=".", color=colors[i_module], linewidth=0, elinewidth=1, label="M"+str(i_module+1))
                        yMax=np.max(offsets_run[i_global][i_station])
                        yMin=np.min(offsets_run[i_global][i_station])
            if (globalN == 1):
                plt.plot(list(selected_runs), trackN[i_station], marker=".", color=colorsS[i_station], linewidth=1, label=stations[i_station] )
                if (i_station==0):
                    yMax=np.max(trackN[i_station])
                    yMin=np.min(trackN[i_station])*0.3 

            axes[i_total].legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 9}) # outside (R) of the plot
            
            #if a smaller range is passed - use it for the relative plot 
            if (y_Max < yMax and plotName=="_relative"):
                axes[i_total].set_ylim(-y_Max, y_Max)
            else:
                #Set y limits to add DS labels
                axes[i_total].set_ylim(yMin*1.5, yMax*1.1)
            
            #Add dataset divisions:
            if (globalN==1 and i_station==1):
                        pass # do just one for tracks plot 
            #draw DS lines
            else:
                for data_set in unique_data_sets:
                    line = [[ data_sets[data_set][0], yMin*1.5], [ data_sets[data_set][0],yMax*1.1]]
                    axes[i_total].plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green', linewidth=2, linestyle=":")
                    line = [[ data_sets[data_set][1], yMin*1.5], [ data_sets[data_set][1], yMax*1.1]]
                    axes[i_total].plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green', linewidth=2, linestyle=":")
                    middle_range = (data_sets[data_set][1]+data_sets[data_set][0])/2
                    # axes[i_total].text(middle_range, yMin*1.36, data_set, fontsize=font_size, color='green', style='italic')
                    axes[i_total].annotate("",xy=(data_sets[data_set][0], yMin*1.4), xytext=(middle_range+middle_range*0.001, yMin*1.4), arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color="green") )
                    axes[i_total].annotate("",xy=(data_sets[data_set][1], yMin*1.4), xytext=(middle_range, yMin*1.4), arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color="green") ) 
            if (globalN==2):
                i_total+=1 # increment axes counter for alignment plots 
    plt.savefig("Monitoring"+plotName+".png", dpi=plot_dpi)
    if (plot_show):
        plt.show()
        
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


selected_runs=[]
#Get and cut on the number of tracks 
trackN=[ [], [] ] #Get track number per run
for i_run, run in enumerate(runs):
    path_s12 = stationPath[0]+"/"+str(run)+"/"
    path_s18 = stationPath[1]+"/"+str(run)+"/"
    #print(run)
    file_s12 = open (path_s12+"millepede.end", "r")
    file_s18 = open (path_s18+"millepede.end", "r")
    status = [int(file_s12.read(8)), int(file_s18.read(8))]
    # check status is correct - abort and re-run pede otherwise 
    # 1 = ok
    # 3 = bad matrix status, removed later 
    if  ( status[0]==-1 or status[1]==-1 ):
        print(run, "has incorrect status, check millepede.end and re-run")
        print("S12",status[0])
        print("S18",status[1])
        sys.exit()

        # input("Clear bad data for re-run?")
        # if (status[0]==-1):
        #     os.chdir(path_s12)
        #     subprocess.call(["rm", "millepede.his"])
        #     subprocess.call(["rm", "millepede.end"])
        #     subprocess.call(["rm", "OffsetsPerModuleS12.fcl"])
        #     os.chdir("../../")
        # if (status[1]==-1):
        #     os.chdir(path_s18)
        #     subprocess.call(["rm", "millepede.his"])
        #     subprocess.call(["rm", "millepede.end"])
        #     subprocess.call(["rm", "OffsetsPerModuleS18.fcl"])
        #     os.chdir("../../")        

    # otherwise, use all '1' status  
    file_s12 = open (path_s12+"millepede.log", "r")
    file_s18 = open (path_s18+"millepede.log", "r")
    tracks = [getTracks(file_s12, "\s+NREC\s+=\s+"), getTracks(file_s18, "\s+NREC\s+=\s+")]
    if ( (tracks[0] > int(trackCut) and tracks[1] > int(trackCut)) and ( status[0]!=3 and status[1]!=3 )  ):
        trackN[0].append(tracks[0])
        trackN[1].append(tracks[1])
        selected_runs.append(run)

runN=len(selected_runs)
print("Using the following",runN,"runs for plotting (after the track cut of",trackCut,"tracks): ", selected_runs)

#Get unique data sets names 
unique_data_sets=[]
print(data_sets)
for data_set in data_sets:
    for run in selected_runs:
        if (run > data_sets[data_set][0] and run < data_sets[data_set][1]):
            print("Found unique", run, data_set)
            unique_data_sets.append(data_set)
            break # do once only 
data_sets_number=len(unique_data_sets)
print("In the following",data_sets_number,"data sets:", unique_data_sets) 

limits = [] # limit for broken axis plot 
for data_set in unique_data_sets:
    limits.append( [data_sets[data_set][0]-2, data_sets[data_set][1]+2] )

#Containers to fill with offsets per dof per station per run per module 
offsets = [[[[0 for i_module in range(0, moduleN)] for i_run in range(0, runN) ] for i_station in range(0, stationN)] for i_global in range(0, globalN) ]

#Get offsets from files per run for all modules in that run
for i_global in range(0, globalN):
    for i_station in range(0, stationN):
        for i_run, run in enumerate(selected_runs):
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

#Pass as a function of i_state:
n_rows_array = [4, 4, 1]
global_n = [2, 2, 1]
plotName = ["", "_relative", "_tracks"]
main_title = ["Alignment per run", r"Alignment per run-$\langle$Alignment$\rangle$ ", "Tracks used in alignment"]

#Produce the 3 plots: absolute, relative and trackN
# for i_state in range(stateN):
   # make_plot(colors, colorsS, plotName[i_state], font_size, plot_dpi, y_Max, n_rows_array[i_state], global_n[i_state], data_sets, main_title[i_state], broken=False)

#Produce 2 broken axis plots for alignment
for i_state in range(stateN-1):
    make_plot(colors, colorsS, plotName[i_state], font_size, plot_dpi, y_Max, n_rows_array[i_state], global_n[i_state], data_sets, main_title[i_state], broken=False, plot_show=True)

print("Monitoring plots are ready!")