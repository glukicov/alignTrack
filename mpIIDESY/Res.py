############
# A single script that runs residual plotter, PEDE, and alignment plotter
# 
# Assumed to be run from a folder with :
# 1) TrackerAlignemnt.root (alignment ana file from RunAlignmentPlots.fcl)
# 2) PEDE data file from Alignment producer (e.g. trackerAlignment_S12.fcl)
#
# Created: 28 March 2019 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 28 March 2019 by Gleb
##############
import subprocess
import argparse, sys, glob, os # glob is awesome! allows regex in process calls 
import numpy.polynomial.polynomial as poly # linear fit 
import matplotlib.pyplot as plt #for plotting 
import pandas as pd # get data frame from text file 
import numpy as np # arrays  
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mticker
from ROOT import TFile, TStyle, TCanvas, gStyle, TF1
import itertools # smart lines in plotting 

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-scan", "--scan", default="None") # scan study 
parser.add_argument("-mode", "--mode", default="None") # scan study 
args = parser.parse_args()
scan=str(args.scan)
mode=str(args.mode)


if (scan == "t0"):
    cutScans = ["29","30", "31", "32", "33", "34", "35", "36", "37", "38", "39"]
    stationName=["S12", "S18"]
    colors=["red", "blue"]
    nominal = "34"

if (scan == "t0_sim"):
    cutScans = ["24.6", "25.1", "25.6", "26.1"]
    stationName=["S0"]
    colors=["green"]
    nominal = "24.6"

stationN=len(stationName)

if (mode == "res"):
    print("Making residual plots")    
        
    for i_total, i_cut in  enumerate(cutScans):
        print("Calling residual plotter in")
        os.chdir(str(i_cut))
        subprocess.call("pwd") # status printout 
        for i_station in range(0, stationN):
            subprocess.call(["python3", "/Users/gleb/software/alignTrack/mpIIDESY/GetRes.py", "-f", "gm2tracker_ana.root", "-s", stationName[i_station], "-eL", scan+": "+str(i_cut)+" ns"])
        os.chdir("../") 


if (mode == "Summary"):
    print("Making residual summary plots")    
   
    tracks = [[0 for i_iter in range(0, len(cutScans)) ] for i_station in range(stationN)]
    maxdU = [[0 for i_iter in range(0, len(cutScans))]  for i_station in range(stationN)]
    maxdV = [[0 for i_iter in range(0, len(cutScans)) ]  for i_station in range(stationN)]
    meandU = [[0 for i_iter in range(0, len(cutScans))]  for i_station in range(stationN)]
    meandV = [[0 for i_iter in range(0, len(cutScans)) ]  for i_station in range(stationN)]
    resolution = [[0 for i_iter in range(0, len(cutScans)) ]  for i_station in range(stationN)]

    for i_iter in range(0, len(cutScans)):
            for i_station in range(0, stationN):
                print("Station:", stationName[i_station], "iter:", cutScans[i_iter])
                
                fileName = "gm2tracker_ana.root" 
                f = TFile.Open(cutScans[i_iter]+"/"+fileName)
                tracks[i_station][i_iter]=(f.Get("TrackSummary"+str(stationName[i_station])+"/FitResults/pValues").GetEntries())
                metricFile=open(cutScans[i_iter]+"/"+"metric_"+stationName[i_station]+".txt", "r")
                numbers = [float(x) for x in next(metricFile).split()] # read first line
                maxdU[i_station][i_iter]=( (numbers[0]) )
                maxdV[i_station][i_iter]=( (numbers[1]) )
                meandU[i_station][i_iter]=( (numbers[2]) )
                meandV[i_station][i_iter]=( (numbers[3]) )
                resolution[i_station][i_iter]=( (numbers[4]) )
    
    print("tracks", tracks)
    #scale tracks by the nominal
    print("Nominal t0: ", nominal)
    if nominal in cutScans:
        nominal_index = cutScans.index(nominal)
        print("nominal index: ", nominal_index)
        for i_station in range(0, stationN):
            print("nominal track#: ", tracks[i_station][nominal_index], "for ", stationName[i_station])
            tracks[i_station]=np.array(tracks[i_station]/np.array(tracks[i_station][nominal_index])) * 100

    
    metric_1 = (tracks)
    metric_2 = (resolution)
    metric_3 = (maxdU, maxdV, meandU, meandV)

    metric=(metric_1, metric_2, metric_3)
    UVlabel=[r"$\langle \Delta U \rangle$", r"$\langle \Delta V \rangle$"]
    #print(metric)

    x_length = len(cutScans)
    print("Total of", x_length, "scans per station")
   
    globalN = len(metric)

    GlobalParNames = ("Tracks [%]", r"$\sigma_{UV}$", "U or V <separation> [um]")
    
    # Plotting "constants"
    f = plt.figure(figsize=(7,int(globalN*3+1)))

    for i_global in range(0, globalN): 

        plt.subplot(int( str(globalN)+"1"+str(int(i_global+1)) )) 
        axes = plt.gca()
        plt.xlabel(str(scan)+ " [ns]", fontsize=12)
        plt.ylabel(GlobalParNames[i_global], fontsize=12)
        plt.xticks(fontsize=12, rotation=0) 
        plt.yticks(fontsize=12, rotation=0)
        plt.minorticks_on()
        axes.tick_params(axis='x',which='minor',bottom=False)
        axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
        plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)
       
        #Plot data 
        for i_station in range(0, stationN):
            if (i_global == 0):
                plt.plot(cutScans, metric[i_global][i_station] , marker="+", color=colors[i_station], label=stationName[i_station])
            if (i_global == 1):
                plt.plot(cutScans, metric[i_global][i_station] , marker="+", color=colors[i_station], label=stationName[i_station])
            if (i_global == 2):
                plt.plot(cutScans,np.abs(metric[i_global][2][i_station]) , marker="+", linestyle=":", color=colors[i_station], label=stationName[i_station]+" "+UVlabel[0])
                plt.plot(cutScans,np.abs(metric[i_global][3][i_station]) , marker="+", linestyle="-.", color=colors[i_station], label=stationName[i_station]+" "+UVlabel[1])

        # axes.legend(loc='top centre', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 
        if (i_global == 1):
            axes.legend(loc='upper center', fontsize=14) # outside (R) of the plot 
        if (i_global == 0 or i_global==2):
            axes.legend(loc='lower center', fontsize=14) # outside (R) of the plot 

    plt.savefig("SummaryResiduals_"+str(scan)+".png", dpi=250)



