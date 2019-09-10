### Gleb 
# Looking into vertical beam from trackers 
from os import listdir
from os.path import isfile, join
import argparse, sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from ROOT import TFile

#Define constants
stations = ("12", "18")
colors = ("red", "blue")

parser = argparse.ArgumentParser()
parser.add_argument("--path", type=str) # offset scale (SD) [um]
args = parser.parse_args()
path=args.path

#get files 
files = [file for file in listdir(path) if isfile(join(path, file))]
files.sort()
print("Found", len(files),"files:", files)

# open and keep all in scope 
tfiles = []
for i_file in files:
    tfiles.append(TFile.Open(path+"/"+i_file))

histArray=[] #1D vertical beam 

# loop over files and stations
run_numbers = []
track_number = []
vertical = [ [], []]
vertical_error = [ [], []]
for i_file, file in enumerate(tfiles):
    run_numbers.append(int(file.Get("TrackSummary/RunInfo/TracksPerEvent").GetMean()))
    track_number.append(file.Get("Extrapolation/vertices/station12/h_verticalPos").GetEntries())
    for i_station, station in enumerate(stations):
        # histo_1D_noCut=file.Get("Extrapolation/vertices/station"+str(station)+"/h_verticalPos")
        histo_2D = file.Get("Extrapolation/vertices/station"+str(station)+"/h_verticalPos_vs_time")
        first_bin = histo_2D.GetXaxis().FindBin(30.0) # 30 us 
        tmpNameTH1 = "tmpNameTH1_"+str(i_file)+str(i_station) # assign a new "name pointer" to the TH1 object
        histo_1D = histo_2D.ProjectionY(tmpNameTH1, first_bin, -1) 
        histArray.append(histo_1D) # keep TH1D in scope if want them later 
        vertical[i_station].append(histo_1D.GetMean())
        vertical_error[i_station].append(histo_1D.GetMeanError())

print("Tracks #",track_number)

### Plotting 
fig, ax = plt.subplots()
x_ticks=x =np.char.replace(files, '.root', '')

# S12 and S18 
# for i_station, station in enumerate(stations):
#     ax.scatter(x_ticks, vertical[i_station], label="S"+station, color=colors[i_station])
#     plt.errorbar(x_ticks, vertical[i_station], yerr=vertical_error[i_station],  color=colors[i_station], markersize=12, elinewidth=1, capsize=2, linewidth=0)

# S18 - S12
# y_title=r"$\Delta$ (S18-S12) Vertical [mm]"
# label = r"$\Delta$ (S18-S12)"
# color="green"
# value = np.array(vertical[1])-np.array(vertical[0])
# error = np.sqrt( np.array(vertical_error[1])**2 + np.array(vertical_error[0])**2 )

# # <S18S12>
# y_title=r"$\langle$S18S12$\rangle$ Vertical [mm]"
# label = r"$\langle$S18S12$\rangle$"
# color="green"
# value =( np.array(vertical[1]) + np.array(vertical[0]) ) / 2 
# error = np.sqrt( np.array(vertical_error[1])**2 + np.array(vertical_error[0])**2 )

ax.scatter(x_ticks, value, label=label, color=color)
plt.errorbar(x_ticks, value, yerr=error, color=color, markersize=12, elinewidth=1, capsize=2, linewidth=0)

#massaging plot 
ax.set_ylabel(y_title, fontsize=14, fontweight='bold')
ax.set_xlabel("Setting", fontsize=14, fontweight='bold')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 12}) # outside (R) of the plot 
plt.tight_layout()
plt.savefig("vertical_run.png", dpi=100)