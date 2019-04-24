#####
# Loop over cuts scans and create a FoM plot 
####

import argparse, sys, glob, os # glob is awesome! allows regex in process calls 
import subprocess, shlex 
import re 
from pathlib import Path #check files
import pandas as pd # get data frame from text file 
import numpy as np # arrays 
import numpy.polynomial.polynomial as poly # linear fit 
import matplotlib.pyplot as plt #for plotting  
import itertools # smart lines in plotting 

#Define constants (Capital + cammelCase)
ModuleN = 8 # per station 
globalN = 2 # X, y
ModuleArray=np.arange(1, ModuleN+1) #(1, 2,...,8) for plotting  
GlobalParNames = ["Radial", "Vertical", 'Φ', 'ψ', 'θ'] #only ever going to have 5 pars.
units = [r" [$\mathrm{\mu m}$]", r" [$\mathrm{\mu m}$]", " [mrad]", " [mrad]", " [mrad]"]

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-scan", "--scan") # scan study 
args = parser.parse_args()

scan=str(args.scan)

#Define tracker constants 
# stationName=[S12", "S18"]
stationName=["S0", "S12", "S18"]
stationN=len(stationName)

sigma = 10 # um 

if (scan == "dca"):
    cutScans = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5]
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
    cutScans = [0, 10000, 20000, 25000, 30000, 35000, 40000]
    defaultValue = 0.0
    cutName = "timeCut"
    scan_units = " [ns]"

elif (scan == "p-value"):
    cutScans = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05]
    defaultValue = 0.005
    cutName = "pValueCut"
    scan_units = " "

elif (scan == "Pz/P"):
    cutScans = [0.0, 0.7, 0.8, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0]
    defaultValue = 0.93
    cutName = "pzCut"
    scan_units = " "

elif (scan == "Hits-per-track"):
    cutScans = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    defaultValue = 9
    pzCut = "minHits"
    scan_units = " "

else:
    print("Incorrect scan specified!")


print("Starting ",str(int(len(cutScans)*stationN))," scan reads")
 
 # Chi2 
metricX=[]
metricY=[]

for i_total, i_cut in  enumerate(cutScans):
    for i_station in range(0, stationN):
    
        os.chdir(str(i_cut)+"/"+str(stationName[i_station])+"/")
        subprocess.call("pwd") # status printout 
        fileName = "metric"+str(stationName[i_station])+".txt"
        metricFile = open(fileName, "r")
        numbers = [int(x) for x in next(metricFile).split()] # read first line
        # Chi2 = residual^2 / sigma^2 
        metricX.append( int(numbers[0])**2/sigma**2 )
        metricY.append( int(numbers[1])**2/sigma**2 )

        os.chdir("../../") # go back to top dir 

metric=(metricX, metricY)
print(metric)

x_length = len(cutScans)
print("Total of", x_length, "scans per station")

colors=["green", "red", "blue"]

# Plotting "constants"
f = plt.figure(figsize=(7,int(globalN*2+1)))

for i_global in range(0, globalN): 

    yMax = [1.2, 1.2]
    yMin = [-0.1 ,-0.1] 
    # if(max(metric[i_global]) > 0):
    #     yMax = [1.25*max(metric[i_global]), 1.25*max(metric[i_global])] 
    # else:
    #     yMax = [0.75*max(metric[i_global]), 0.75*max(metric[i_global])] 

    # if(min(metric[i_global]) < 0):
    #     yMin = [1.25*min(metric[i_global]), 1.25*min(metric[i_global])]
    # else:
    #     yMin = [0.75*min(metric[i_global]), 0.75*min(metric[i_global])]

    plt.subplot(int( str(globalN)+"1"+str(int(i_global+1)) )) 
    axes = plt.gca()
    axes.set_xlim(cutScans[0]-cutScans[1]/10, cutScans[-1]*1.2)
    axes.set_ylim(yMin[i_global], yMax[i_global])
    plt.title(GlobalParNames[i_global]+" alignment performance for "+ str(scan) +" scan", fontsize=12)
    plt.ylabel(r"$\chi^{2}_{\mathrm{ndf}}$"+ units[i_global])
    plt.xlabel(str(scan)+str(scan_units), fontsize=12)
    plt.xticks(fontsize=10, rotation=0) 
    plt.yticks(fontsize=10, rotation=0)
    plt.minorticks_on()
    axes.tick_params(axis='x',which='minor',bottom=False)
    axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
    plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)
    #Plot the 0th line 
    line = [[cutScans[0]-cutScans[1]/10,0.0], [cutScans[-1]*1.2, 0.0]]
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'black')
    #Plot the 0th line 
    line = [[cutScans[0]-cutScans[1]/10,1.0], [cutScans[-1]*1.2, 1.0]]
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'black', linestyle="--")

    #Plot data 
    for i_station in range(0, stationN):
        plt.plot(cutScans, metric[i_global][int(i_station)::3] , marker="+", color=colors[i_station], label=stationName[i_station])

    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 

plt.savefig("ScanResults_"+str(scan)+".png", dpi=250)



