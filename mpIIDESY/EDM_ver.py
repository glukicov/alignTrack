### Gleb 
# Looking into vertical beam from trackers 


from os import listdir
from os.path import isfile, join
import argparse, sys
import matplotlib.pyplot as plt
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

# loop over files and stations
run_numbers = []
track_number = []
vertical = [ [], []]
vertical_error = [ [], []]
for file in tfiles:
    run_numbers.append(int(file.Get("TrackSummary/RunInfo/TracksPerEvent").GetMean()))
    track_number.append(file.Get("Extrapolation/vertices/station12/h_verticalPos").GetEntries())
    for i_station, station in enumerate(stations):
        vertical[i_station].append(file.Get("Extrapolation/vertices/station"+str(station)+"/h_verticalPos").GetMean())
        vertical_error[i_station].append(file.Get("Extrapolation/vertices/station"+str(station)+"/h_verticalPos").GetMeanError())


print("Run numbers:", run_numbers)
print("Tracks #",track_number)

### Plotting 
fig, ax = plt.subplots()
for i_station, station in enumerate(stations):
    ax.scatter(run_numbers, vertical[i_station], label="S"+station, color=colors[i_station])
    plt.errorbar(run_numbers, vertical[i_station], yerr=vertical_error[i_station],  color=colors[i_station], markersize=12, elinewidth=1, linewidth=0)

#massaging plot 
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.set_ylabel("Vertical Beam Position [mm]", fontsize=14, fontweight='bold')
ax.set_xlabel("Run #", fontsize=14, fontweight='bold')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 12}) # outside (R) of the plot 
plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)

plt.savefig("vertical_run.png", dpi=200)






