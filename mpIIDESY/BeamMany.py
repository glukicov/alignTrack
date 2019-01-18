################
# Script that produces the final FoM 
#
# Gleb Lukicov 27 December 2018 
############

import argparse, sys, os
import os.path
import subprocess
from ROOT import TFile, TStyle, TCanvas, gStyle, TF1, TPaveStats, gPad
from scipy import stats
import matplotlib.pyplot as plt 
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats
import itertools

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-mis", "--mis") # offset scale (SD) [um]
parser.add_argument("-trialN", "--trialN") # number of iterations
parser.add_argument("-mode", "--mode") # offset scale (SD) [um]
parser.add_argument("-dof", "--dof") # radial, vertical or both 
parser.add_argument("-shift", "--shift") # offset scale (mean) [um]
args = parser.parse_args()

# CONSTANTS 
trialN = int(args.trialN) 
misalignment = int(args.mis)
shift = int(args.shift) 
mode = str(args.mode)
dof=str(args.dof)


dirName=str(shift)+"_Off_"+str(misalignment)+"_Mis_"+str(mode)+"_"+str(dof)


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
RMS_S12_rad=[]
RMS_S12_ver=[]
RMS_S18_rad=[]
RMS_S18_ver=[]
RMS_S12_rad_error=[]
RMS_S12_ver_error=[]
RMS_S18_rad_error=[]
RMS_S18_ver_error=[]
cases=[]

#First take the reference nominal point "0"
scr = "trackRecoPlots.root"
# scr = "/gm2/app/users/glukicov/TrackerAlignment/gm2Dev_v9_14_00/srcs/gm2tracker/align/Systematics/nominalAlign/trackRecoPlots.root"
scrFile = TFile.Open(scr)

#open the histos
s12_rad = scrFile.Get(str(station12Path)+plotNames[0])
s18_rad = scrFile.Get(str(station18Path)+plotNames[0])
s12_ver = scrFile.Get(str(station12Path)+plotNames[1])
s18_ver = scrFile.Get(str(station18Path)+plotNames[1])

s12_rad.GetXaxis().SetRangeUser(-70, 70) # cutting out tails 
s18_rad.GetXaxis().SetRangeUser(-70, 70)
s12_ver.GetXaxis().SetRangeUser(-70, 70)
s18_ver.GetXaxis().SetRangeUser(-70, 70)

S12_rad.append(s12_rad.GetMean())
S12_ver.append(s12_ver.GetMean())
S18_rad.append(s18_rad.GetMean())
S18_ver.append(s18_ver.GetMean())
S12_rad_error.append(s12_rad.GetMeanError())
S12_ver_error.append(s12_ver.GetMeanError())
S18_rad_error.append(s18_rad.GetMeanError())
S18_ver_error.append(s18_ver.GetMeanError())
RMS_S12_rad.append(s12_rad.GetRMS())
RMS_S12_ver.append(s12_ver.GetRMS())
RMS_S18_rad.append(s18_rad.GetRMS())
RMS_S18_ver.append(s18_ver.GetRMS())
RMS_S12_rad_error.append(s12_rad.GetRMSError())
RMS_S12_ver_error.append(s12_ver.GetRMSError())
RMS_S18_rad_error.append(s18_rad.GetRMSError())
RMS_S18_ver_error.append(s18_ver.GetRMSError())
cases.append(int(0))


#Fil containers in a loop 
for i_trial in range(0, trialN):

    #scr = str(i_trial+1)+"/trackRecoPlots.root"
    scr = "/pnfs/GM2/scratch/users/glukicov/Systematics_Plots_new/"+str(dirName)+"/"+str(i_trial+1)+"/trackRecoPlots.root"

    #Check file exists
    scrFilePresent=os.path.isfile(str(scr))
    #print(scrFilePresent)  
      
    if(scrFilePresent):
        scrFile = TFile.Open(scr)
        #open the histos
        s12_rad = scrFile.Get(str(station12Path)+plotNames[0])
        s18_rad = scrFile.Get(str(station18Path)+plotNames[0])
        s12_ver = scrFile.Get(str(station12Path)+plotNames[1])
        s18_ver = scrFile.Get(str(station18Path)+plotNames[1])

        s12_rad.GetXaxis().SetRangeUser(-70, 70) # cutting out tails 
        s18_rad.GetXaxis().SetRangeUser(-70, 70)
        s12_ver.GetXaxis().SetRangeUser(-70, 70)
        s18_ver.GetXaxis().SetRangeUser(-70, 70) 

        S12_rad.append(s12_rad.GetMean())
        S12_ver.append(s12_ver.GetMean())
        S18_rad.append(s18_rad.GetMean())
        S18_ver.append(s18_ver.GetMean())
        S12_rad_error.append(s12_rad.GetMeanError())
        S12_ver_error.append(s12_ver.GetMeanError())
        S18_rad_error.append(s18_rad.GetMeanError())
        S18_ver_error.append(s18_ver.GetMeanError())
        RMS_S12_rad.append(s12_rad.GetRMS())
        RMS_S12_ver.append(s12_ver.GetRMS())
        RMS_S18_rad.append(s18_rad.GetRMS())
        RMS_S18_ver.append(s18_ver.GetRMS())
        RMS_S12_rad_error.append(s12_rad.GetRMSError())
        RMS_S12_ver_error.append(s12_ver.GetRMSError())
        RMS_S18_rad_error.append(s18_rad.GetRMSError())
        RMS_S18_ver_error.append(s18_ver.GetRMSError())
        cases.append(int(i_trial+1)) #already have the 1st element from nominal 
    
    #if doesn't exist, move on to the next case 
    else:
        continue

#combine in stations for ease of plotting 
S_rad = (S12_rad, S18_rad)
S_rad_error = (S12_rad_error, S18_rad_error)
S_ver = (S12_ver, S18_ver)
S_ver_error = (S12_ver_error, S18_ver_error)
RMS_S_rad = (RMS_S12_rad, RMS_S18_rad)
RMS_S_rad_error = (RMS_S12_rad_error, RMS_S18_rad_error)
RMS_S_ver = (RMS_S12_ver, RMS_S18_ver)
RMS_S_ver_error = (RMS_S12_ver_error, RMS_S18_ver_error)
label=("S12", "S18")
color=("purple", "orange")

#New containers for the mean values 
S_rad_mean = []
S_rad_error_mean = []
S_ver_mean = []
S_ver_error_mean = []

#Fill Plots MEAN 
plt.figure(1) 

plt.subplot(211) # RADIAL 
axes=plt.gca()
axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

for i in range(0, 2):
    avgMean=np.mean(S_rad[i][1:]) # compute mean for offset cases only
    avgMeanError=stats.sem(S_rad[i][1:])
    S_rad_mean.append(avgMean)
    S_rad_error_mean.append(avgMeanError)
    line = [[0,S_rad[i][0]], [trialN+1, S_rad[i][0]]]  # put line at mean position 
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = str(color[i]), linestyle="--") # plot the line 
    # line = [[1,avgMean], [trialN+1, avgMean]]  # put line at mean position 
    # plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="--") # plot the line 
    # plt.text(trialN-3, avgMean+0.02, str(round(avgMean, 3))+r" $\pm$ "+ str(round(avgMeanError,3)), fontsize=9) # plot the mean value above line 
    plt.text(trialN-2, S_rad[i][0]+0.02, str(round(S_rad[i][0], 3))+r" $\pm$ "+str(round(S_rad_error[i][0],3)), fontsize=9) # plot the reference value 
    plt.scatter(cases, S_rad[i],  marker="+", color=str(color[i])) # plot all cases and reference points 
    plt.errorbar(cases, S_rad[i], yerr=S_rad_error[i],  color=str(color[i]), markersize=14, elinewidth=3, linewidth=0, label=str(label[i])) # add error bar

plt.legend(loc='center')
plt.title("Radial Position case " + str(dirName) , fontsize=10)
plt.xlabel("Sample #", fontsize=10)
plt.ylabel("Mean Radial Position [mm]", fontsize=10)



plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)

plt.subplot(212) # VERTICAL 
axes=plt.gca()
axes.xaxis.set_major_locator(MaxNLocator(integer=True))
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')

for i in range(0, 2):
    avgMean=np.mean(S_ver[i][1:]) # compute mean for offset cases only 
    avgMeanError=stats.sem(S_ver[i][1:])
    S_ver_mean.append(avgMean)
    S_ver_error_mean.append(avgMeanError)
    line = [[0,S_ver[i][0]], [trialN+1, S_ver[i][0]]]
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = str(color[i]), linestyle="--")
    # line = [[1,avgMean], [trialN+1, avgMean]]
    # plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="--")
    #plt.text(trialN-3, avgMean+0.02, str(round(avgMean,3))+r" $\pm$ "+str(round(avgMeanError,3)), fontsize=9)
    plt.text(trialN-2, S_ver[i][0]+0.02, str(round(S_ver[i][0], 3))+r" $\pm$ "+str(round(S_ver_error[i][0],3)), fontsize=9)
    plt.scatter(cases, S_ver[i],  marker="+", color=str(color[i]))
    plt.errorbar(cases, S_ver[i], yerr=S_ver_error[i],  color=str(color[i]), markersize=14, elinewidth=3, linewidth=0, label=str(label[i]))

plt.legend(loc='center')
plt.title("Vertical Position case "+str(dirName), fontsize=10)
plt.xlabel("Sample #", fontsize=10)
plt.ylabel("Mean Vertical Position [mm]", fontsize=10)

plt.subplots_adjust(bottom=0.1)
plt.subplots_adjust(hspace=0.43)
plt.subplots_adjust(right=0.9)

plt.savefig("Mean_"+str(dirName)+".png", dpi=600)

print(str(trialN)+" cases analysed!")

print("For the case of misalignment with SD of "+str(misalignment)+" um.")
#Print final results
print("Radial Shifts:")
for i in range(0, 2):
    # print("Mean Radial Position [mm] "+str(label[i])+ ": ", str(round(S_rad_mean[i], 3))+r" $\pm$ "+ str(round(S_rad_error_mean[i],3)) )
    # print("Reference Radial Position [mm] "+str(label[i])+ ": ", str(round(S_rad[i][0], 3))+r" $\pm$ "+str(round(S_rad_error[i][0],3)) )
    value = S_rad_mean[i] - S_rad[i][0] 
    error = np.sqrt( S_rad_error_mean[i]**2 +  S_rad_error[i][0]**2 )
    # print(r"$\Delta$ Radial Position [mm] "+str(label[i])+ ": ", str(round(value, 3))+r" $\pm$ "+str(round(error,3)) )
    print(r"$\Delta$ Radial Position [um] "+str(label[i])+ ": ", str(round(value*1e3, 3))+" +/- "+str(round(error*1e3,3)) )

    print("Vertical Shifts:")
for i in range(0, 2):
    # print("Mean Vertical Position [mm] "+str(label[i])+ ": ", str(round(S_ver_mean[i], 3))+r" $\pm$ "+ str(round(S_ver_error_mean[i],3)) )
    # print("Reference Vertical Position [mm] "+str(label[i])+ ": ", str(round(S_ver[i][0], 3))+r" $\pm$ "+str(round(S_ver_error[i][0],3)) )
    value = S_ver_mean[i] - S_ver[i][0] 
    error = np.sqrt( S_ver_error_mean[i]**2 +  S_ver_error[i][0]**2 )
    # print(r"$\Delta$ Vertical Position [mm] "+str(label[i])+ ": ", str(round(value, 3))+r" $\pm$ "+str(round(error,3)) )
    print(r"$\Delta$ Vertical Position [um] "+str(label[i])+ ": ", str(round(value*1e3, 3))+" +/- "+str(round(error*1e3,3)) )

#New containers for the mean values 
RMS_S_rad_mean = []
RMS_S_rad_error_mean = []
RMS_S_ver_mean = []
RMS_S_ver_error_mean = []

#Fill Plots SD 
plt.figure(2) 

plt.subplot(211) # RADIAL 
axes=plt.gca()
axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

for i in range(0, 2):
    avgMean=np.mean(RMS_S_rad[i][1:]) # compute mean for offset cases only
    avgMeanError=stats.sem(RMS_S_rad[i][1:])
    RMS_S_rad_mean.append(avgMean)
    RMS_S_rad_error_mean.append(avgMeanError)
    line = [[1,RMS_S_rad[i][0]], [trialN+1, RMS_S_rad[i][0]]]  # put line at mean position 
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="-") # plot the line 
    #line = [[1,avgMean], [trialN+1, avgMean]]  # put line at mean position 
    #plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="--") # plot the line 
    #plt.text(trialN-3, avgMean+0.02, str(round(avgMean, 3))+r" $\pm$ "+ str(round(avgMeanError,3)), fontsize=9) # plot the mean value above line 
    plt.text(0, RMS_S_rad[i][0]+0.02, str(round(RMS_S_rad[i][0], 3))+r" $\pm$ "+str(round(RMS_S_rad_error[i][0],3)), fontsize=9) # plot the reference value 
    plt.scatter(cases, RMS_S_rad[i],  marker="+", color=str(color[i])) # plot all cases and reference points 
    plt.errorbar(cases, RMS_S_rad[i], yerr=RMS_S_rad_error[i],  color=str(color[i]), markersize=12, elinewidth=2, label=str(label[i])) # add error bar

plt.legend(loc='center')
plt.title("Radial Position misalignment SD of "+str(misalignment), fontsize=10)
plt.xlabel("Sample #", fontsize=10)
plt.ylabel("SD Radial Position [mm]", fontsize=10)

plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)

plt.subplot(212) # VERTICAL 
axes=plt.gca()
axes.xaxis.set_major_locator(MaxNLocator(integer=True))
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')

for i in range(0, 2):
    avgMean=np.mean(RMS_S_ver[i][1:]) # compute mean for offset cases only 
    avgMeanError=stats.sem(RMS_S_ver[i][1:])
    RMS_S_ver_mean.append(avgMean)
    RMS_S_ver_error_mean.append(avgMeanError)
    line = [[1,RMS_S_ver[i][0]], [trialN+1, RMS_S_ver[i][0]]]
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="-")
    #line = [[1,avgMean], [trialN+1, avgMean]]
    #plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="--")
    #plt.text(trialN-3, avgMean+0.02, str(round(avgMean,3))+r" $\pm$ "+str(round(avgMeanError,3)), fontsize=9)
    plt.text(0, RMS_S_ver[i][0]+0.02, str(round(RMS_S_ver[i][0], 3))+r" $\pm$ "+str(round(RMS_S_ver_error[i][0],3)), fontsize=9)
    plt.scatter(cases, RMS_S_ver[i],  marker="+", color=str(color[i]))
    plt.errorbar(cases, RMS_S_ver[i], yerr=RMS_S_ver_error[i],  color=str(color[i]), markersize=12, elinewidth=2, label=str(label[i]))

plt.legend(loc='center')
plt.title("Vertical Position misalignment SD of "+str(misalignment), fontsize=10)
plt.xlabel("Sample #", fontsize=10)
plt.ylabel("SD Vertical Position [mm]", fontsize=10)

plt.subplots_adjust(bottom=0.1)
plt.subplots_adjust(hspace=0.43)

plt.savefig("SD_"+str(dirName)+".png", dpi=600)