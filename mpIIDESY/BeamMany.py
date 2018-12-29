################
# Script that produces the final FoM 
#
# Gleb Lukicov 27 December 2018 
############

import argparse, sys, os 
import subprocess
from ROOT import TFile, TStyle, TCanvas, gStyle, TF1, TPaveStats, gPad
from scipy import stats
import matplotlib.pyplot as plt 
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker
import numpy as np

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-mis", "--mis") # offset scale (SD) [um]
parser.add_argument("-trialN", "--trialN") # number of iterations
args = parser.parse_args()

# CONSTANTS 
trialN = int(args.trialN) 
misalignment = int(args.mis)
station12Path = "Extrapolation/vertices/station12/pValue>0.005_and_noVolumesHit/"
station18Path = "Extrapolation/vertices/station18/pValue>0.005_and_noVolumesHit/"
plotNames=["h_radialPos", "h_verticalPos"]


#Global containers
S12_rad=[]
S12_ver=[]
S18_rad=[]
S18_ver=[]
S12_rad_error=[]
S12_ver_error=[]
S18_rad_error=[]
S18_ver_error=[]
cases=[]

#First take the reference nominal point "0"
scr = "1/trackRecoPlots.root"
# scrFile = "/gm2/app/users/glukicov/TrackerAlignment/gm2Dev_v9_14_00/srcs/gm2tracker/align/Systematics/nominalAlign/trackRecoPlots.root"
scrFile = TFile.Open(scr)

#open the histos
s12_rad = scrFile.Get(str(station12Path)+plotNames[0])
s18_rad = scrFile.Get(str(station18Path)+plotNames[0])
s12_ver = scrFile.Get(str(station12Path)+plotNames[1])
s18_ver = scrFile.Get(str(station18Path)+plotNames[1])

S12_rad.append(s12_rad.GetMean())
S12_ver.append(s12_ver.GetMean())
S18_rad.append(s18_rad.GetMean())
S18_ver.append(s18_ver.GetMean())
S12_rad_error.append(s12_rad.GetMeanError())
S12_ver_error.append(s12_ver.GetMeanError())
S18_rad_error.append(s18_rad.GetMeanError())
S18_ver_error.append(s18_ver.GetMeanError())
cases.append(int(0))


#Fil containers in a loop 
for i_trial in range(0, trialN):

    scr = str(i_trial+1)+"/trackRecoPlots.root"
    #scrFile = "/pnfs/GM2/scratch/users/glukicov/Systematics_Plots/"+str(misalignment)+"/"+str(i_trial+1)+"/trackRecoPlots.root"
    scrFile = TFile.Open(scr)

    #open the histos
    s12_rad = scrFile.Get(str(station12Path)+plotNames[0])
    s18_rad = scrFile.Get(str(station18Path)+plotNames[0])
    s12_ver = scrFile.Get(str(station12Path)+plotNames[1])
    s18_ver = scrFile.Get(str(station18Path)+plotNames[1])

    S12_rad.append(s12_rad.GetMean())
    S12_ver.append(s12_ver.GetMean())
    S18_rad.append(s18_rad.GetMean())
    S18_ver.append(s18_ver.GetMean())
    S12_rad_error.append(s12_rad.GetMeanError())
    S12_ver_error.append(s12_ver.GetMeanError())
    S18_rad_error.append(s18_rad.GetMeanError())
    S18_ver_error.append(s18_ver.GetMeanError())
    cases.append(int(i_trial+1))


#Fill Plots
plt.figure(1) 

plt.subplot(211) # RADIAL 
axes=plt.gca()
axes.xaxis.set_major_locator(MaxNLocator(integer=True)) 
plt.plot(cases, S12_rad,  marker="+", color="purple")
plt.errorbar(cases, S12_rad, yerr=S12_rad_error,  color="purple", markersize=12, elinewidth=2, label="S12")
plt.plot(cases, S18_rad,  marker="+", color="orange")
plt.errorbar(cases, S18_rad, yerr=S18_rad_error,  color="orange", markersize=12, elinewidth=2, label="S18")
plt.legend(loc='center')
plt.title("Radial Position misalignment SD of "+str(misalignment), fontsize=10)
plt.xlabel("Misalignment Case", fontsize=10)
plt.ylabel("Radial Position [mm]", fontsize=10)

plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)

plt.subplot(212) # VERTICAL 
axes=plt.gca()
axes.xaxis.set_major_locator(MaxNLocator(integer=True)) 
plt.plot(cases, S12_ver,  marker="+", color="purple")
plt.errorbar(cases, S12_ver, yerr=S12_ver_error,  color="purple", markersize=12, elinewidth=2, label="S12")
plt.plot(cases, S18_ver,  marker="+", color="orange")
plt.errorbar(cases, S18_ver, yerr=S18_ver_error,  color="orange", markersize=12, elinewidth=2, label="S18")
plt.legend(loc='center')
plt.title("Vertical Position misalignment SD of "+str(misalignment), fontsize=10)
plt.xlabel("Misalignment Case", fontsize=10)
plt.ylabel("Vertical Position [mm]", fontsize=10)

plt.subplots_adjust(bottom=0.1)
plt.subplots_adjust(hspace=0.43)

plt.savefig(str(misalignment)+".png", dpi=600)

print(str(trialN)+" cases analysed!")