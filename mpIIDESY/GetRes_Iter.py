####################################################################
# Input: TrackerAlignment.root analysis-level TFile after alignment
# Output: Steerable by args 
#
# e.g. python GetRes_Iter.py -f1/Iter1/TrackerAlignment.root -f2/Iter3/TrackerAlignment.root
# Will produce UV Residuals and Residual SD plots per module in a station 
# 
# Created: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 26 November 2018 by Gleb
#####################################################################

from ROOT import TFile, TStyle, TCanvas, gStyle, TF1
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

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-m', '--moduleN', help='mode', default=-1) # # of removed module from tracking (if applicable)
parser.add_argument('-f1', '--fileN1', help='input ROOT file', default="TrackerAlignment.root")
parser.add_argument('-f2', '--fileN2', help='input ROOT file', default="TrackerAlignment.root")
parser.add_argument('-eL', '--eL', help='extra label', default="", type=str)
args = parser.parse_args()

eL=str(args.eL)

NModules=8
NLayers=4 # per modules
tracker=("Tracker 1", "Tracker 2")
markerShape=["8", "*"]
NTotalLayers=32
LayerNames = ["U0", "U1", "V0", "V1"]
moduleNames=np.arange(1, NModules+1) #1-8

fileName1 = str(args.fileN1)
fileName2 = str(args.fileN2)

#get just the filename...
fileName1_short = fileName1.rsplit('/', 1)[-1]
fileName2_short = fileName1.rsplit('/', 1)[-1]

regime = None 
if (fileName1_short == "TrackerAlignment.root" and fileName2_short== "TrackerAlignment.root"):
    print("Plotting Residuals from Alignment Tracks!")
    regime="align"
else:
    print("Not expected input file name!")
    print("expected both files named:TrackerAlignment.root")
    print("Please provide TrackerAlingment analysis-level plots.")
    print("Run alignment FHICL with 'monitor : false' and output gm2tracker_reco.root file")
    print("RunAlignmentPlots.fcl on that file, and use the result (TrackerAlignment.root) it in this script")
    sys.exit()

f1 = TFile.Open(fileName1)
f2 = TFile.Open(fileName2)
if f1 and f2:
    print(str(fileName1)+" and "+ str(fileName2)+" are open")
else:
     print(str(fileName1)+" or "+ str(fileName2)+" not found")

if (regime=="align"): # only alignment data has Pz/P and implicit station number 
    label_mean_1 = f1.Get("TrackerAlignment/Hits/Labels").GetMean()
    label_mean_2 = f2.Get("TrackerAlignment/Hits/Labels").GetMean()
    #print("Mean label 1:", round(label_mean_1))
    #print("Mean label 2:", round(label_mean_2))
    label_mean = (label_mean_1 + label_mean_2)/2
    if(label_mean < 1280 and label_mean > 1210):
        stationN = "S12"
        i_station=0
    elif(label_mean < 1880 and label_mean > 1810):
        stationN = "S18"
        i_station=1
    elif(label_mean < 1080 and label_mean > 1010):
        stationN = "S0"
    else:
        print("Mean mismatch..")

####### LAYERS ##############
#-------LayerResiudals_Zoom----------
i_totalLayer=0
means=[]
MeanErrors=[]
yMin = -60
yMax = 60
plt.figure(71)
axes = plt.gca()
for i in range(0, len(moduleNames)):
    i_module=moduleNames[i]
    if (regime=="align"):
        name = "TrackerAlignment/Modules/Residuals UV Module " + str(i_module)
    t1 = f1.Get(str(name))
    t2 = f2.Get(str(name))
    mean1 = t1.GetMean()
    mean2 = t2.GetMean()
    meanError1 = t1.GetMeanError()
    meanError2 = t2.GetMeanError()
    plt.errorbar(i_module, mean1*1e3, yerr=meanError1*1e3, color="purple", markersize=15, elinewidth=3) 
    plt.errorbar(i_module, mean2*1e3, yerr=meanError2*1e3, color="green", markersize=15, elinewidth=3) 
    plt.plot(i_module, mean1*1e3, marker=markerShape[0], color="purple", markersize=15, mew=1, linewidth=0, label="Misaligned" if i == 0 else "")
    plt.plot(i_module, mean2*1e3, marker=markerShape[1], color="green", markersize=15, mew=1, linewidth=0, label="Aligned" if i == 0 else "")
    line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey', linewidth=1)
    i_totalLayer+=1

line = [[0.5,0.0], [NModules+0.5, 0.0]]
plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'black', linewidth=1)
# axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 12}) # outside (R) of the plot 
axes.legend(title=stationN+":", title_fontsize=15, loc='upper left', prop={'size': 12})
axes.set_xlim(0.5, NModules+0.5)
axes.set_ylim(yMin, yMax)
# plt.title("UV Residuals in "+stationN+" "+eL, fontsize=18)
# plt.title("UV Residuals in "+tracker[i_station]+" "+eL, fontsize=18)
plt.ylabel(r"Residual Mean [$\mathrm{\mu}$m]", fontsize=18)
axes.tick_params(axis="y", labelsize=11)
axes.tick_params(axis="x", labelsize=11)
plt.xlabel("Module", fontsize=18)
plt.tight_layout()
plt.savefig("Residuals_L_Zoom"+str(stationN)+".png", dpi=600)


print("ROOT File analysed!")