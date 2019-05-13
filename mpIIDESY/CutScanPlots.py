#####
# Loop over cuts scans and create a FoM plot 
####

import argparse, sys, glob, os # glob is awesome! allows regex in process calls 
import subprocess, shlex 
import re 
# from pathlib import Path #check files
import pandas as pd # get data frame from text file 
import numpy as np # arrays 
import numpy.polynomial.polynomial as poly # linear fit 
import matplotlib.pyplot as plt #for plotting  
import itertools # smart lines in plotting 
from matplotlib.ticker import MaxNLocator

#Define constants (Capital + cammelCase)
ModuleN = 8 # per station 
ModuleArray=np.arange(1, ModuleN+1) #(1, 2,...,8) for plotting  
GlobalParNames = ["Radial", "Vertical", "Track number"] #only ever going to have 5 pars.
units = [r" [$\mathrm{\mu m}$]", r" [$\mathrm{\mu m}$]", " ", " [mrad]", " [mrad]"]

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-scan", "--scan") # scan study 
args = parser.parse_args()

scan=str(args.scan)

#Define tracker constants 
colors=["red", "blue"]
stationName=["S12", "S18"]
# colors=["green", "red", "blue"]
# stationName=["S0", "S12", "S18"]
stationN=len(stationName)

sigma = [2, 10] # um 

if (scan == "DCA"):
    cutScans = [0.3, 0.4, 0.5, 0.6, 0.7]
    defaultValue = 0.5
    cutName = "cutDCA"
    scan_units = " [mm]"
   
elif (scan == "P"):
    cutScans = [0, 800, 1000, 1300, 1500, 1700, 2000]
    # cutScans = [0, 800]
    defaultValue = 0.0
    cutName = "PCut"
    scan_units = " [MeV]"

elif (scan == "start-time"):
    cutScans = [0, 10000, 20000, 30000, 35000]
    defaultValue = 0.0
    cutName = "timeCut"
    scan_units = " [ns]"

elif (scan == "pValue"):
    cutScans = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.007, 0.01, 0.03, 0.05]
    defaultValue = 0.005
    cutName = "pValueCut"
    scan_units = " "

elif (scan == "PzP"):
    cutScans = [0.0, 0.9, 0.93, 0.95, 0.98]
    defaultValue = 0.93
    cutName = "pzCut"
    scan_units = " "

elif (scan == "Hits-per-track"):
    cutScans = [5, 7, 8, 9, 10, 11]
    defaultValue = 9
    cutName = "minHits"
    scan_units = " "

elif (scan == "Iteration"):
    cutScans = ["defDig_defTrack_Iter1", "defDig_defTrack_Iter1", "defDig_defTrack_Iter3", "defDig_defTrack_Iter4"]
    defaultValue = 9
    cutName = "minHits"
    scan_units = " "

elif (scan == "t0"):
    cutScans = [30, 32, 34, 35, 37, 38]
    scan_units = " [ns]"

else:
    print("Incorrect scan specified!")


print("Starting ",str(int(len(cutScans)*stationN))," scan reads")
 
 # Chi2 
metricX=[]
metricY=[]
trackN=[]
yLabel = "" # set inside the loop 

for i_total, i_cut in  enumerate(cutScans):
    for i_station in range(0, stationN):
    
        os.chdir(str(i_cut)+"/"+str(stationName[i_station])+"/")
        subprocess.call("pwd") # status printout 
        fileName = "metric"+str(stationName[i_station])+".txt"
        metricFile = open(fileName, "r")
        numbers = [float(x) for x in next(metricFile).split()] # read first line
        
        # # residual= meanAbsReco-meanAbsTruth
        metricX.append( (numbers[0]) )
        metricY.append( (numbers[1]) )
        yLabel = "|<Alignment>| [um]"

        # # Chi2 = residual^2 / sigma^2 
        # metricX.append( (numbers[0])**2/sigma[0]**2 )
        # metricY.append( (numbers[1])**2/sigma[1]**2 )
        # yLabel = r"$\eta^{2}_{\mathrm{ndf}}$"


        metricFile = open("trackN.txt", "r")
        numbers = [float(x) for x in next(metricFile).split()] # read first line
        trackN.append(numbers[0])

        os.chdir("../../") # go back to top dir 

metric=(metricX, metricY, trackN)
print(metric)

x_length = len(cutScans)
print("Total of", x_length, "scans per station")

# cutScans = [1, 2, 3 , 4]

globalN = len(metric)
# Plotting "constants"
f = plt.figure(figsize=(7,int(globalN*2+1)))

for i_global in range(0, globalN): 

    plt.subplot(int( str(globalN)+"1"+str(int(i_global+1)) )) 
    axes = plt.gca()
    # axes.set_xlim(cutScans[0]-cutScans[1]/10, cutScans[-1]*1.2)
    axes.set_xlim(29, 39)
    # axes.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.title(GlobalParNames[i_global]+" alignment performance for "+ str(scan) +" scan", fontsize=12)
    if(i_global == 2):
        plt.title(GlobalParNames[i_global]+" for "+ str(scan) +" scan", fontsize=12)
    plt.ylabel(yLabel)
    if(i_global == 2):
        plt.ylabel("trackN")
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel(str(scan)+str(scan_units), fontsize=12)
    plt.xticks(fontsize=10, rotation=0) 
    plt.yticks(fontsize=10, rotation=0)
    plt.minorticks_on()
    axes.tick_params(axis='x',which='minor',bottom=False)
    axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
    plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)
    if(i_global != 2):
        #Plot the 0th line 
        line = [[cutScans[0]-cutScans[1]/10, 0.0], [cutScans[-1]*1.2, 0.0]]
        #plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'black', linewidth=0.5)
        #Plot the 1.0 line 
        line = [[cutScans[0]-cutScans[1]/10, 1.0], [cutScans[-1]*1.2, 1.0]]
        #plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'black', linewidth=0.5, linestyle=":")
   
    #Plot data 
    for i_station in range(0, stationN):
        plt.plot(cutScans, metric[i_global][int(i_station)::stationN] , marker="+", color=colors[i_station], label=stationName[i_station])

    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 

plt.savefig("ScanResults_"+str(scan)+".png", dpi=250)