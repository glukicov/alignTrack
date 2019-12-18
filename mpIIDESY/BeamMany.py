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
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats
import itertools
import re
import pandas as pd
from ROOT import TH1F, TF1, TCanvas
#sys.path.append(os.path.abspath(str(os.environ['MRB_SOURCE'])+"/gm2tracker/align/macros"))
# sys.path.append(os.path.abspath("~"))
# from rootlogon import SetMyStyle, myStyle
# SetMyStyle()

def getOffsets(f, name):
    offsets = [] #tmp storage buffer 
    for line in f:
        if re.match(name, line):
            copy = True
            offsets=line 

        else:
            copy = False

    return offsets   


parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-mis", "--mis") # offset scale (SD) [um]
parser.add_argument("-shift", "--shift") # offset mean shift (mean) [um]
parser.add_argument("-sign", "--sign") # offset sign
parser.add_argument("-trialN", "--trialN", nargs='?', type=int, default=25) # number of iterations
parser.add_argument("-mode", "--mode", nargs='?', type=str, default="U") # offset scale (SD) [um]
parser.add_argument("-dof", "--dof", nargs='?', type=str, default="B") # radial, vertical or both 
args = parser.parse_args()

# CONSTANTS 
trialN = int(args.trialN) 
mode = str(args.mode)
dof=str(args.dof)
misalignment = int(args.mis)
shift = int(args.shift) 
sign=str(args.sign)

dirName=str(sign)+str(shift)+"_Off_"+str(misalignment)+"_Mis_"+str(mode)+"_"+str(dof)


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
scr = "/Users/gleb/software/alignTrack/mpIIDESY/AlignmentData/BK/Systematics/nominalAlign/trackRecoPlots.root"
# scr = str(os.environ['MRB_SOURCE'])+"/gm2tracker/align/Systematics/nominalAlign/trackRecoPlots.root"
scrFile = TFile.Open(scr)

#open the histos
s12_rad = scrFile.Get(str(station12Path)+plotNames[0])
s18_rad = scrFile.Get(str(station18Path)+plotNames[0])
s12_ver = scrFile.Get(str(station12Path)+plotNames[1])
s18_ver = scrFile.Get(str(station18Path)+plotNames[1])

s12_rad.GetXaxis().SetRangeUser(-50, 50) # cutting out tails 
s18_rad.GetXaxis().SetRangeUser(-50, 50)
s12_ver.GetXaxis().SetRangeUser(-30, 30)
s18_ver.GetXaxis().SetRangeUser(-30, 30)


# ####### FITTING ###########

# s12_ver.GetXaxis().SetRangeUser(-30, 30)
# s18_ver.GetXaxis().SetRangeUser(-30, 30)

# gF = TF1("gF", "gaus", -30, 30)

# s12_ver.Fit(gF, "QR")
# S12_ver.append(gF.GetParameter(1))
# S12_ver_error.append(gF.GetParError(1))
# RMS_S12_ver.append(gF.GetParameter(2))
# RMS_S12_ver_error.append(gF.GetParError(2))

# s18_ver.Fit(gF, "QR")
# S18_ver.append(gF.GetParameter(1))
# S18_ver_error.append(gF.GetParError(1))
# RMS_S18_ver.append(gF.GetParameter(2))
# RMS_S18_ver_error.append(gF.GetParError(2))

####### FITTING ###########
S12_ver.append(s12_ver.GetMean())
S18_ver.append(s18_ver.GetMean())
S12_ver_error.append(s12_ver.GetMeanError())
S18_ver_error.append(s18_ver.GetMeanError())
RMS_S12_ver.append(s12_ver.GetRMS())
RMS_S18_ver.append(s18_ver.GetRMS())
RMS_S12_ver_error.append(s12_ver.GetRMSError())
RMS_S18_ver_error.append(s18_ver.GetRMSError())

S12_rad.append(s12_rad.GetMean())
S18_rad.append(s18_rad.GetMean())
S12_rad_error.append(s12_rad.GetMeanError())
S18_rad_error.append(s18_rad.GetMeanError())
RMS_S12_rad.append(s12_rad.GetRMS())
RMS_S18_rad.append(s18_rad.GetRMS())
RMS_S12_rad_error.append(s12_rad.GetRMSError())
RMS_S18_rad_error.append(s18_rad.GetRMSError())


cases.append(int(0))


Mis_SD_rad=[]
Mis_Mean_rad=[]
dS_all_rad=[]
dM_all_rad=[]
dS_all_ver=[]
dM_all_ver=[]

#Fil containers in a loop 
for i_trial in range(0, trialN):

    print(i_trial)
    # if (i_trial==5 or i_trial==11 or i_trial==12):
    #     continue

    # scr = "/Users/gleb/software/alignTrack/mpIIDESY/AlignmentData/AUE/BK/Systematics/"+str(dirName)+"/"+str(i_trial+1)+"/trackRecoPlots.root"
    scr = /Users/gleb/software/alignTrack/mpIIDESY/AlignmentData/AUE/BK/P0_Off_100_Mis_U_B

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

        s12_rad.GetXaxis().SetRangeUser(-50, 50) # cutting out tails 
        s18_rad.GetXaxis().SetRangeUser(-50, 50)
        s12_ver.GetXaxis().SetRangeUser(-30, 30)
        s18_ver.GetXaxis().SetRangeUser(-30, 30) 


        # ####### FITTING ###########

        # s12_ver.GetXaxis().SetRangeUser(-30, 30)
        # s18_ver.GetXaxis().SetRangeUser(-30, 30)

        # gF = TF1("gF", "gaus", -30, 30)

        # s12_ver.Fit(gF, "QR")
        # S12_ver.append(gF.GetParameter(1))
        # S12_ver_error.append(gF.GetParError(1))
        # RMS_S12_ver.append(gF.GetParameter(2))
        # RMS_S12_ver_error.append(gF.GetParError(2))

        # s18_ver.Fit(gF, "QR")
        # S18_ver.append(gF.GetParameter(1))
        # S18_ver_error.append(gF.GetParError(1))
        # RMS_S18_ver.append(gF.GetParameter(2))
        # RMS_S18_ver_error.append(gF.GetParError(2))

        ####### FITTING ###########
        S12_ver.append(s12_ver.GetMean())
        S18_ver.append(s18_ver.GetMean())
        S12_ver_error.append(s12_ver.GetMeanError())
        S18_ver_error.append(s18_ver.GetMeanError())
        RMS_S12_ver.append(s12_ver.GetRMS())
        RMS_S18_ver.append(s18_ver.GetRMS())
        RMS_S12_ver_error.append(s12_ver.GetRMSError())
        RMS_S18_ver_error.append(s18_ver.GetRMSError())

        S12_rad.append(s12_rad.GetMean())
        S18_rad.append(s18_rad.GetMean())
        S12_rad_error.append(s12_rad.GetMeanError())
        S18_rad_error.append(s18_rad.GetMeanError())
        RMS_S12_rad.append(s12_rad.GetRMS())
        RMS_S18_rad.append(s18_rad.GetRMS())
        RMS_S12_rad_error.append(s12_rad.GetRMSError())
        RMS_S18_rad_error.append(s18_rad.GetRMSError())


        cases.append(int(i_trial+1)) #already have the 1st element from nominal

        '''
        #file = "Systematics_reco/"+str(dirName)+"/"+str(i_trial+1)+"/RunTrackingDAQ.fcl"
        file = str(os.environ['MRB_SOURCE'])+"/gm2tracker/align/Systematics_reco/"+str(dirName)+"/"+str(i_trial+1)+"/RunTrackingDAQ.fcl"

        f=open(file, "r")
        radial = (getOffsets(f, "services.Geometry.strawtracker.radialOffsetPerModule:"))
        radial = radial.replace("services.Geometry.strawtracker.radialOffsetPerModule: [", " ") 
        radial = radial.replace("]", "") 
        radialOff = np.array([float(r) for r in radial.split(',')])
        radialOff=radialOff[8:16] 
        radialOff=radialOff*1e3 # mm -> um 

        Mis_SD_rad.append(int(np.std(radialOff))) # to the nearest um 
        Mis_Mean_rad.append(int(np.mean(radialOff)))
        '''

        

    
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
S_rad_max = []
S_rad_error_max = []

S_ver_mean = []
S_ver_error_mean = []
S_ver_max = []
S_ver_error_max = []


#Fill Plots MEAN 
plt.figure(1) 

plt.subplot(211) # RADIAL 
axes=plt.gca()
axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

for i in range(0, 2):
    line = [[0,S_rad[i][0]], [trialN+1, S_rad[i][0]]]  # put line at mean position 
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = str(color[i]), linestyle="--") # plot the line 
    plt.text(0, S_rad[i][0]+0.02, str(round(S_rad[i][0], 3))+r" $\pm$ "+str(round(S_rad_error[i][0],3)), fontsize=9) # plot the reference value 
    
    # for i_trial in range(0, trialN+1):
        # plt.text(i_trial, S_rad[i][i_trial]+0.02, str(round(S_rad[i][i_trial], 3))+r" $\pm$ "+str(round(S_rad_error[i][i_trial],3)), fontsize=9) # plot the reference value 
    
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
axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')

for i in range(0, 2):
    line = [[0,S_ver[i][0]], [trialN+1, S_ver[i][0]]]
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = str(color[i]), linestyle="--")
    plt.text(0, S_ver[i][0]+0.02, str(round(S_ver[i][0], 3))+r" $\pm$ "+str(round(S_ver_error[i][0],3)), fontsize=9)
    
    # for i_trial in range(0, trialN+1):
    #     plt.text(i_trial, S_ver[i][i_trial]+0.02, str(round(S_ver[i][i_trial], 3))+r" $\pm$ "+str(round(S_ver_error[i][i_trial],3)), fontsize=9)
    
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


# '''

#New containers for the mean values 
RMS_S_rad_mean = []
RMS_S_rad_error_mean = []
RMS_S_rad_max = []
RMS_S_rad_error_max = []

RMS_S_ver_mean = []
RMS_S_ver_error_mean = []
RMS_S_ver_max = []
RMS_S_ver_error_max = []



#Fill Plots SD 
plt.figure(2) 

plt.subplot(211) # RADIAL 
axes=plt.gca()
axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

for i in range(0, 2):
    line = [[0,RMS_S_rad[i][0]], [trialN+1, RMS_S_rad[i][0]]]  # put line at mean position 
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = str(color[i]), linestyle="--") 
    # plt.text(0, RMS_S_rad[i][0]+0.02, str(round(RMS_S_rad[i][0], 3))+r" $\pm$ "+str(round(RMS_S_rad_error[i][0],3)), fontsize=9) # plot the reference value 
    # plt.text(trialN, RMS_S_rad[i][trialN]+0.02, str(round(RMS_S_rad[i][trialN], 3))+r" $\pm$ "+str(round(RMS_S_rad_error[i][trialN],3)), fontsize=9) # plot the reference value 
    plt.scatter(cases, RMS_S_rad[i],  marker="+", color=str(color[i])) # plot all cases and reference points 
    plt.errorbar(cases, RMS_S_rad[i], yerr=RMS_S_rad_error[i],  color=str(color[i]), markersize=12, elinewidth=2, linewidth=0, label=str(label[i])) # add error bar

plt.legend(loc='center')
plt.title("Radial Position case " + str(dirName) , fontsize=10)
plt.xlabel("Sample #", fontsize=10)
plt.ylabel("SD Radial Position [mm]", fontsize=10)

plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)

plt.subplot(212) # VERTICAL 
axes=plt.gca()
axes.xaxis.set_major_locator(MaxNLocator(integer=True))
axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')

for i in range(0, 2):
    line = [[0,RMS_S_ver[i][0]], [trialN+1, RMS_S_ver[i][0]]]
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = str(color[i]), linestyle="--")
    # plt.text(0, RMS_S_ver[i][0]+0.02, str(round(RMS_S_ver[i][0], 3))+r" $\pm$ "+str(round(RMS_S_ver_error[i][0],3)), fontsize=9)
    # plt.text(trialN, RMS_S_ver[i][trialN]+0.02, str(round(RMS_S_ver[i][trialN], 3))+r" $\pm$ "+str(round(RMS_S_ver_error[i][trialN],3)), fontsize=9)
    plt.scatter(cases, RMS_S_ver[i],  marker="+", color=str(color[i]))
    plt.errorbar(cases, RMS_S_ver[i], yerr=RMS_S_ver_error[i],  color=str(color[i]), markersize=12, elinewidth=2, linewidth=0, label=str(label[i]))

plt.legend(loc='center')
plt.title("Vertical Position case "+str(dirName), fontsize=10)
plt.xlabel("Sample #", fontsize=10)
plt.ylabel("SD Vertical Position [mm]", fontsize=10)

plt.subplots_adjust(bottom=0.1)
plt.subplots_adjust(hspace=0.43)

plt.savefig("SD_"+str(dirName)+".png", dpi=600)

# '''

# # Plot for covariance
# plt.figure(3) 
# plt.suptitle("Correlation plots for "+str(dirName) , fontsize=10)
# plt.subplot(221) # Mean Beam vs Mean Mis  
# axes=plt.gca()
# axes.set_xlim(min(Mis_Mean_rad)*0.9, max(Mis_Mean_rad)*1.1)
# # axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
# axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
# axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

# for i in range(0, 2):
#     dM = np.array(S_rad[i][1:])-S_rad[i][0] # take away the nominal per station 
#     dM_all_rad.append(dM)
#     dM_ver=np.array(S_ver[i][1:])-S_ver[i][0]
#     dM_all_ver.append(dM_ver)
#     avgMean=np.mean(abs(dM)) # compute mean for offset cases only
#     avgMeanError=stats.sem(abs(dM))
#     S_rad_mean.append(avgMean)
#     S_rad_error_mean.append(avgMeanError)
#     avgMax=np.max(abs(dM)) # compute mean for offset cases only
#     S_rad_max.append(avgMax)
#     S_rad_error_max.append(S_rad_error[i][np.argmax(dM)])
#     # print("dM=",dM, label[i])
#     plt.scatter(Mis_Mean_rad, dM,  marker="+", color=str(color[i])) # plot all cases and reference points 
#     plt.errorbar(Mis_Mean_rad, dM, yerr=S_rad_error[i][1:],  color=str(color[i]), markersize=14, elinewidth=3, linewidth=0, label=str(label[i])) # add error bar
#     cor = np.corrcoef(Mis_Mean_rad, dM)
#     plt.plot([], [], ' ', label=str(label[i])+" corr: "+str(round(cor[0][1],3)))
#     plt.plot(np.unique(Mis_Mean_rad), np.poly1d(np.polyfit(Mis_Mean_rad, dM, 1))(np.unique(Mis_Mean_rad)), color=str(color[i]))

# plt.legend(loc='center')
# #plt.title("d(Radial Position) vs. <Misalignment> " + str(dirName) , fontsize=10)
# plt.xlabel("<Misalignment> [um]", fontsize=10)
# plt.ylabel("d(Radial Position) [mm]", fontsize=10)

# # plt.subplot(222) # Mean Beam vs SD Mis  
# # axes=plt.gca()
# # axes.set_xlim(min(Mis_SD_rad)*0.9, max(Mis_SD_rad)*1.1)
# # # axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
# # axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
# # axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

# # for i in range(0, 2):
# #     dM = np.array(S_rad[i][1:])-S_rad[i][0] # take away the nominal per station 
# #     plt.scatter(Mis_SD_rad, dM,  marker="+", color=str(color[i])) # plot all cases and reference points 
# #     plt.errorbar(Mis_SD_rad, dM, yerr=S_rad_error[i][1:],  color=str(color[i]), markersize=14, elinewidth=3, linewidth=0, label=str(label[i])) # add error bar
# #     cor = np.corrcoef(Mis_SD_rad, dM)
# #     plt.plot([], [], ' ', label=str(label[i])+" corr: "+str(round(cor[0][1],3)))
# #     plt.plot(np.unique(Mis_SD_rad), np.poly1d(np.polyfit(Mis_SD_rad, dM, 1))(np.unique(Mis_SD_rad)), color=str(color[i]))

# # plt.legend(loc='center')
# # #plt.title("d(Radial Position) vs. SD(Misalignment) " + str(dirName) , fontsize=10)
# # plt.xlabel("SD(Misalignment) [um]", fontsize=10)
# # plt.ylabel("d(Radial Position) [mm]", fontsize=10)


# plt.subplot(223) # SD Beam vs Mean Mis  
# axes=plt.gca()
# axes.set_xlim(min(Mis_Mean_rad)*0.9, max(Mis_Mean_rad)*1.1)
# # axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
# axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
# axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

# for i in range(0, 2):
#     dS = np.array(RMS_S_rad[i][1:])-RMS_S_rad[i][0] # take away the nominal per station 
#     dS_all_rad.append(dS)
#     dS_ver = np.array(RMS_S_ver[i][1:])-RMS_S_ver[i][0] # take away the nominal per station 
#     dS_all_ver.append(dS_ver)
#     avgMean=np.mean(abs(dS)) # compute mean for offset cases only
#     avgMeanError=stats.sem(abs(dS))
#     RMS_S_rad_mean.append(avgMean)
#     RMS_S_rad_error_mean.append(avgMeanError)
#     avgMax=np.max(abs(dS)) # compute mean for offset cases only
#     RMS_S_rad_max.append(avgMax)
#     RMS_S_rad_error_max.append(RMS_S_rad_error[i][np.argmax(dS)])
#     plt.scatter(Mis_Mean_rad, dS,  marker="+", color=str(color[i])) # plot all cases and reference points 
#     plt.errorbar(Mis_Mean_rad, dS, yerr=RMS_S_rad_error[i][1:],  color=str(color[i]), markersize=14, elinewidth=3, linewidth=0, label=str(label[i])) # add error bar
#     cor = np.corrcoef(Mis_Mean_rad, dS)
#     plt.plot([], [], ' ', label=str(label[i])+" corr: "+str(round(cor[0][1],3)))
#     plt.plot(np.unique(Mis_Mean_rad), np.poly1d(np.polyfit(Mis_Mean_rad, dS, 1))(np.unique(Mis_Mean_rad)), color=str(color[i]))

# plt.legend(loc='center')
# #plt.title("d(Radial Position) vs. <Misalignment> " + str(dirName) , fontsize=10)
# plt.xlabel("<Misalignment> [um]", fontsize=10)
# plt.ylabel("d(SD(Radial Position)) [mm]", fontsize=10)


# # plt.subplot(224) # SD Beam vs SD Mis  
# # axes=plt.gca()
# # axes.set_xlim(min(Mis_SD_rad)*0.9, max(Mis_SD_rad)*1.1)
# # # axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
# # axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
# # axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

# # for i in range(0, 2):
# #     dS = np.array(RMS_S_rad[i][1:])-RMS_S_rad[i][0] # take away the nominal per statio
# #     plt.scatter(Mis_SD_rad, dS,  marker="+", color=str(color[i])) # plot all cases and reference points 
# #     plt.errorbar(Mis_SD_rad, dS, yerr=RMS_S_rad_error[i][1:],  color=str(color[i]), markersize=14, elinewidth=3, linewidth=0, label=str(label[i])) # add error bar
# #     cor = np.corrcoef(Mis_SD_rad, dS)
# #     plt.plot([], [], ' ', label=str(label[i])+" corr: "+str(round(cor[0][1],3)))
# #     plt.plot(np.unique(Mis_SD_rad), np.poly1d(np.polyfit(Mis_SD_rad, dS, 1))(np.unique(Mis_SD_rad)), color=str(color[i]))
    
# # plt.legend(loc='center')
# # #plt.title("d(Radial Position) vs. <Misalignment> " + str(dirName) , fontsize=10)
# # plt.xlabel("SD(Misalignment) [um]", fontsize=10)
# # plt.ylabel("d(SD(Radial Position)) [mm]", fontsize=10)

# # plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)
# # plt.subplots_adjust(top=0.9)
# # plt.subplots_adjust(bottom=0.1)
# # plt.subplots_adjust(hspace=0.43)
# # plt.subplots_adjust(right=0.9)

# plt.savefig("Corr_"+str(dirName)+".png", dpi=600)

# S12_dm_rad=dM_all_rad[::2]
# S12_ds_rad=dS_all_rad[::2]
# S18_dm_rad=dM_all_rad[1::2]
# S18_ds_rad=dS_all_rad[1::2]

# S12_dm_ver=dM_all_ver[::2]
# S12_ds_ver=dS_all_ver[::2]
# S18_dm_ver=dM_all_ver[1::2]
# S18_ds_ver=dS_all_ver[1::2]


# radVal = [S12_dm_rad, S18_dm_rad, S12_ds_rad,  S18_ds_rad]
# verVal = [S12_dm_ver, S18_dm_ver, S12_ds_ver,  S18_ds_ver]


# radBinsN = 10
# radBinMax = S_rad_max[0]*1e3*1.2
# c_rad = TCanvas("c_rad", "Radial Histograms", 6200, 6200)
# c_rad.Divide(2,2)
# h1r=TH1F("S12: d(Radial Mean)", str(dirName), radBinsN, -radBinMax, radBinMax)
# h2r=TH1F("S18: d(Radial Mean)", str(dirName), radBinsN, -radBinMax, radBinMax)
# h3r=TH1F("S12: d(Radial SD)", str(dirName), radBinsN, -radBinMax, radBinMax)
# h4r=TH1F("S18: d(Radial SD)", str(dirName), radBinsN, -radBinMax, radBinMax)

# histos = [h1r, h2r, h3r, h4r]

# for i, histo in enumerate(histos):
#     for value in radVal[i]:
#         for i_value in value:
#             histo.Fill(float(i_value*1e3))

# for i, histo in enumerate(histos):      
#     c_rad.cd(i+1)
#     histo.Draw()
#     hitso.GetXAxis().SetTitle("[um]")

# # c_rad.Draw()
# c_rad.SaveAs("radCan_"+str(dirName)+".png")

# verBinsN = 20
# verBinMax = 0.5*S_rad_max[0]*1e3*1.2
# c_ver = TCanvas("c_ver", "Vertical Histograms", 6200, 6200)
# c_ver.Divide(2,2)
# h1v=TH1F("S12: d(Ver Mean)", str(dirName), verBinsN, -verBinMax, verBinMax)
# h2v=TH1F("S18: d(Ver Mean)", str(dirName), verBinsN, -verBinMax, verBinMax)
# h3v=TH1F("S12: d(Ver SD)", str(dirName), verBinsN, -verBinMax, verBinMax)
# h4v=TH1F("S18: d(Ver SD)", str(dirName), verBinsN, -verBinMax, verBinMax)

# histos = [h1v, h2v, h3v, h4v]

# for i, histo in enumerate(histos):
#     for value in verVal[i]:
#         for i_value in value:
#             histo.Fill(float(i_value*1e3))

# for i, histo in enumerate(histos):      
#     c_rad.cd(i+1)
#     histo.Draw()
#     hitso.GetXAxis().SetTitle("[um]")

# # c_rad.Draw()
# c_rad.SaveAs("verCan_"+str(dirName)+".png")



# # #Write to the final metric file
# # np.savetxt('metric_'+str(dirName)+'.txt', (Mis_Mean_rad, Mis_SD_rad, S12_dm, S12_ds, S18_dm, S18_ds), delimiter=',' , fmt='%s')   
# # dataPDF = [Mis_Mean_rad, Mis_SD_rad, S12_dm, S12_ds, S18_dm, S18_ds]
# # pd.DataFrame(dataPDF).to_csv('metric_'+str(dirName)+'.csv')

# # f=open("fom_"+dirName+".txt", "w+")
# # f.write(str(misalignment)+"\n")
# # f.write(str(int(S_rad_mean[0]*1e3))+"\n")
# # f.write(str(int(S_rad_mean[1]*1e3))+"\n")
# # f.write(str(int(S_rad_error_mean[0]*1e3))+"\n")
# # f.write(str(int(S_rad_error_mean[1]*1e3))+"\n")
# # f.write(str(int(S_rad_max[0]*1e3))+"\n")
# # f.write(str(int(S_rad_max[1]*1e3))+"\n")
# # f.write(str(int(S_rad_error_max[0]*1e3))+"\n")
# # f.write(str(int(S_rad_error_max[1]*1e3))+"\n")
# # f.write(str(int(RMS_S_rad_mean[0]*1e3))+"\n")
# # f.write(str(int(RMS_S_rad_mean[1]*1e3))+"\n")
# # f.write(str(int(RMS_S_rad_error_mean[0]*1e3))+"\n")
# # f.write(str(int(RMS_S_rad_error_mean[1]*1e3))+"\n")
# # f.write(str(int(RMS_S_rad_max[0]*1e3))+"\n")
# # f.write(str(int(RMS_S_rad_max[1]*1e3))+"\n")
# # f.write(str(int(RMS_S_rad_error_max[0]*1e3))+"\n")
# # f.write(str(int(RMS_S_rad_error_max[1]*1e3))+"\n")
# # f.close()

# print(str(trialN)+" cases analysed!")

# '''