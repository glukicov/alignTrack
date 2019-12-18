################
# Script that produces the final FoM 
#
# Gleb Lukicov 27 December 2018 
############

import argparse, sys, os
import os.path
import subprocess
from ROOT import TFile, TStyle, TCanvas, gStyle, TF1, TPaveStats, gPad, gROOT, TH1, TH1F
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
import math
sys.path.append(os.path.abspath(str(os.environ['MRB_SOURCE'])+"/gm2tracker/align/macros"))
from rootlogon import SetMyStyle, myStyle
SetMyStyle()

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
# trialN = int(args.trialN) 
trialN = 100
# mode = str(args.mode)
# dof=str(args.dof)
# misalignment = int(args.mis)
# shift = int(args.shift) 
# sign=str(args.sign)

# dirName=str(sign)+str(shift)+"_Off_"+str(misalignment)+"_Mis_"+str(mode)+"_"+str(dof)
dirName="P0_Off_100_Mis_U_B"

#Define constant paths and labels 
station12Path = "Extrapolation/vertices/station12/"
station18Path = "Extrapolation/vertices/station18/"
plotPath = [station12Path, station18Path, station12Path, station18Path]
plotName=["h_verticalPos_vs_time", "h_verticalPos_vs_time", "h_radialPos_vs_time", "h_radialPos_vs_time"] # TFile TH2 names  
canvasName = ["s12_vertical", "s18_vertical", "s12_radial", "s18_radial"]
hitTimeHistoPath_s12 = "HitSummary/Station_12/h_hitTime"
hitTimeHistoPath_s18 = "HitSummary/Station_18/h_hitTime"

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
events=[]
hits=[]
tacks=[]
events_s12=[]
events_s18=[]
hits_s12=[]
hits_s18=[]
tracks_s12=[]
tracks_s18=[]
hitsTotal_s12=[] 
hitsTotal_s18=[]

Mis_SD_rad=[]
Mis_Mean_rad=[]
dS_all_rad=[]
dM_all_rad=[]
dS_all_ver=[]
dM_all_ver=[]

#Fil containers in a loop 
for i_trial in range(0, trialN+1):

    print(i_trial)
  
    scr = "/gm2/data/users/glukicov/AUE_plots/"+str(dirName)+"/"+str(i_trial)+"/trackRecoPlots.root"

    #Check file exists
    scrFilePresent=os.path.isfile(str(scr))
       
    if(scrFilePresent):

        scrFile = TFile.Open(scr)
        
        for i_plot in range(0, len(plotPath)):

            #Get the TH2F 
            histo_2D = scrFile.Get(str(plotPath[i_plot])+plotName[i_plot])

            #Apply 30 us time cut 
            first_bin = histo_2D.GetXaxis().FindBin(30.0)

            tmpNameTH1 = "tmpNameTH1_"+str(i_plot)+str(i_trial) # assign a new "name pointer" to the TH1 object for each loop 
            hist_1D = histo_2D.ProjectionY(tmpNameTH1, first_bin, -1)
            hist_1D.GetXaxis().SetRangeUser(-50, 50) # applying a maximum range cut 

            #Do a fit for vertical only 
            if( canvasName[i_plot] == "s12_vertical" or  canvasName[i_plot] == "s18_vertical"):

                # #Define a Gaussian fit function between -30 to 30 mm 
                # tmpNameTF1 = "tmpNameTF1_"+str(i_plot)+str(i_trial)
                # gF = TF1(tmpNameTF1, "gaus", -30.0, 30.0)
                # hist_1D.Fit(gF, "QNR") # quite fit no graphics over the specified range above
                # #Get stats from fit 
                # mean = gF.GetParameter(1)
                # mean_error = gF.GetParError(1)
                # sd = gF.GetParameter(2)
                # sd_error = gF.GetParError(2)

                # #Get stats from hist 
                mean = hist_1D.GetMean()
                mean_error = hist_1D.GetMeanError()
                sd = hist_1D.GetRMS()
                sd_error = hist_1D.GetRMSError()

                if( canvasName[i_plot] == "s12_vertical"):
                    S12_ver.append(mean)
                    S12_ver_error.append(mean_error)
                    RMS_S12_ver.append(sd)
                    RMS_S12_ver_error.append(sd_error)

                    #get normalised # events per station 
                    tracks = hist_1D.GetEntries()
                    tracks_s12.append(tracks)
                    
                    #Total hits 
                    hits_histo = scrFile.Get(hitTimeHistoPath_s12)
                    first_bin = hits_histo.GetXaxis().FindBin(0.03) # get integral for 30us ->
                    last_bin = hits_histo.GetXaxis().FindBin(1.0) # get integral for 30us ->  
                    hitsTotal_s12.append(hits_histo.Integral(first_bin, last_bin) )

                    #Hits in candidates s
                    driftHit_2D = scrFile.Get("HitSummary/Station_12/h_driftTimeVsHitTime")
                    first_bin = driftHit_2D.GetXaxis().FindBin(0.03)
                    tmpNameDrift = "tmpNameDrift_"+str(i_plot)+str(i_trial) # assign a new "name pointer" to the TH1 object for each loop 
                    driftHist = driftHit_2D.ProjectionY(tmpNameDrift, first_bin, -1)
                    hits = driftHist.GetEntries()

                    hits_s12.append(hits)
                    entries = tracks/hits
                    events_s12.append(entries)

                if( canvasName[i_plot] == "s18_vertical"):
                    S18_ver.append(mean)
                    S18_ver_error.append(mean_error)
                    RMS_S18_ver.append(sd)
                    RMS_S18_ver_error.append(sd_error)

                    #get normalised # events per station 
                    tracks = hist_1D.GetEntries()
                    tracks_s18.append(tracks)

                    hits_histo = scrFile.Get(hitTimeHistoPath_s18)
                    first_bin = hits_histo.GetXaxis().FindBin(0.03) # get integral for 30us ->
                    last_bin = hits_histo.GetXaxis().FindBin(1.0) # get integral for 30us ->  
                    hitsTotal_s18.append(hits_histo.Integral(first_bin, last_bin) )

                    driftHit_2D = scrFile.Get("HitSummary/Station_18/h_driftTimeVsHitTime")
                    first_bin = driftHit_2D.GetXaxis().FindBin(0.03)
                    tmpNameDrift = "tmpNameDrift_"+str(i_plot)+str(i_trial) # assign a new "name pointer" to the TH1 object for each loop 
                    driftHist = driftHit_2D.ProjectionY(tmpNameDrift, first_bin, -1)
                    hits = driftHist.GetEntries()
                    hits_s18.append(hits)
                    entries = tracks/hits
                    events_s18.append(entries)

            #Extract mean from a radial histogram only  
            if( canvasName[i_plot] == "s12_radial" or  canvasName[i_plot] == "s18_radial"):

                #Get stats from hist 
                mean = hist_1D.GetMean()
                mean_error = hist_1D.GetMeanError()
                sd = hist_1D.GetRMS()
                sd_error = hist_1D.GetRMSError()

                if( canvasName[i_plot] == "s12_radial"):
                    S12_rad.append(mean)
                    S12_rad_error.append(mean_error)
                    RMS_S12_rad.append(sd)
                    RMS_S12_rad_error.append(sd_error)

                if( canvasName[i_plot] == "s18_radial"):
                    S18_rad.append(mean)
                    S18_rad_error.append(mean_error)
                    RMS_S18_rad.append(sd)
                    RMS_S18_rad_error.append(sd_error)

        if not(i_trial == 0):
            cases.append(int(i_trial)) 

        # if (i_trial != 0):
        #     #Read the misalignments from the FHICL file per case skipping the nominal (0)
        #     file = str(os.environ['MRB_SOURCE'])+"/gm2tracker/align/AUE/"+str(dirName)+"/"+str(i_trial)+"/RunTrackingDAQ_align.fcl"
        #     f=open(file, "r")
        #     radial = (getOffsets(f, "services.Geometry.strawtracker.strawModuleRShift12:"))
        #     radial = radial.replace("services.Geometry.strawtracker.strawModuleRShift12: [", " ") 
        #     radial = radial.replace("]", "") 
        #     radialOff = np.array([float(r) for r in radial.split(',')])
        #     radialOff=radialOff[0:8] 
        #     radialOff=radialOff*1e3 # mm -> um 
        #     Mis_SD_rad.append(int(np.std(radialOff))) # to the nearest um 
        #     Mis_Mean_rad.append(int(np.mean(radialOff)))
    
    #if doesn't exist, move on to the next case 
    else:
        continue

#combine in stations for ease of plotting 
events=[events_s12, events_s18]
tracks=[tracks_s12, tracks_s18]
hits=[hits_s12, hits_s18]
hitsTotal=[hitsTotal_s12, hitsTotal_s18]
S_rad = [S12_rad, S18_rad]
S_rad_error = [S12_rad_error, S18_rad_error]
S_ver = [S12_ver, S18_ver]
S_ver_error = [S12_ver_error, S18_ver_error]
RMS_S_rad = [RMS_S12_rad, RMS_S18_rad]
RMS_S_rad_error = [RMS_S12_rad_error, RMS_S18_rad_error]
RMS_S_ver = [RMS_S12_ver, RMS_S18_ver]
RMS_S_ver_error = [RMS_S12_ver_error, RMS_S18_ver_error]
label=["S12", "S18"]
color=["purple", "orange"]
marker=["+", "*"]

print("len=", len(events[0]))
### "Remove" the nominal



for i_station in range(0, 2):
    events[i_station]=np.array(events[i_station][1:])-events[i_station][0]
    tracks[i_station]=np.array(tracks[i_station][1:])-tracks[i_station][0]
    hits[i_station]=np.array(hits[i_station][1:])-hits[i_station][0]
    hitsTotal[i_station]=np.array(hitsTotal[i_station][1:])-hitsTotal[i_station][0]
    S_rad[i_station]=np.array(S_rad[i_station][1:])-S_rad[i_station][0]
    S_rad_error[i_station]=np.array(S_rad_error[i_station][1:])-S_rad_error[i_station][0]
    S_ver[i_station]=np.array(S_ver[i_station][1:])-S_ver[i_station][0]
    S_ver_error[i_station]=np.array(S_ver_error[i_station][1:])-S_ver_error[i_station][0]
    RMS_S_rad[i_station]=np.array(RMS_S_rad[i_station][1:])-RMS_S_rad[i_station][0]
    RMS_S_rad_error[i_station]=np.array(RMS_S_rad_error[i_station][1:])-RMS_S_rad_error[i_station][0]
    RMS_S_ver[i_station]=np.array(RMS_S_ver[i_station][1:])-RMS_S_ver[i_station][0]
    RMS_S_ver_error[i_station]=np.array(RMS_S_ver_error[i_station][1:])-RMS_S_ver_error[i_station][0]






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
    #line = [[1 ,S_rad[i][0]], [trialN, S_rad[i][0]]]  # put line at mean position 
    #plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = str(color[i]), linestyle="--") # plot the line 
    #plt.text(0, S_rad[i][0]+0.02, str(round(S_rad[i][0], 3))+r" $\pm$ "+str(round(S_rad_error[i][0],3)), fontsize=9) # plot the reference value 
    #plt.text(trialN-1, S_rad[i][trialN-1]+0.02, str(round(S_rad[i][trialN-1], 3))+r" $\pm$ "+str(round(S_rad_error[i][trialN-1],3)), fontsize=9) # plot the reference value 
    plt.scatter(cases, S_rad[i],  marker=str(marker[i]), color=str(color[i]), label=str(label[i])) # plot all cases and reference points 
    plt.errorbar(cases, S_rad[i], yerr=S_rad_error[i],  color=str(color[i]), markersize=14, elinewidth=3, linewidth=0) # add error bar

axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.title("Radial Position case " + str(dirName) , fontsize=10)
plt.xlabel("Sample #", fontsize=10)
plt.ylabel(r"$\Delta$ Mean Radial [mm]", fontsize=10)



plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)

plt.subplot(212) # VERTICAL 
axes=plt.gca()
axes.xaxis.set_major_locator(MaxNLocator(integer=True))
axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')

for i in range(0, 2):
    #line = [[1,S_ver[i][0]], [trialN, S_ver[i][0]]]
    #plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = str(color[i]), linestyle="--")
    #plt.text(0, S_ver[i][0]+0.02, str(round(S_ver[i][0], 3))+r" $\pm$ "+str(round(S_ver_error[i][0],3)), fontsize=9)
    #plt.text(trialN-1, S_ver[i][trialN-1]+0.02, str(round(S_ver[i][trialN-1], 3))+r" $\pm$ "+str(round(S_ver_error[i][trialN-1],3)), fontsize=9)
    plt.scatter(cases, S_ver[i],  marker=str(marker[i]), color=str(color[i]), label=str(label[i]))
    plt.errorbar(cases, S_ver[i], yerr=S_ver_error[i],  color=str(color[i]), markersize=14, elinewidth=3, linewidth=0)

axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.title("Vertical Position case "+str(dirName), fontsize=10)
plt.xlabel("Sample #", fontsize=10)
plt.ylabel(r"$\Delta$ Mean Vertical [mm]", fontsize=10)

plt.subplots_adjust(bottom=0.1)
plt.subplots_adjust(hspace=0.43)
plt.subplots_adjust(right=0.8)
#plt.show()
plt.savefig("xMean_"+str(dirName)+".png", dpi=600)


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
    #line = [[1,RMS_S_rad[i][0]], [trialN, RMS_S_rad[i][0]]]  # put line at mean position 
    #plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = str(color[i]), linestyle="--") 
    #plt.text(0, RMS_S_rad[i][0]+0.02, str(round(RMS_S_rad[i][0], 3))+r" $\pm$ "+str(round(RMS_S_rad_error[i][0],3)), fontsize=9) # plot the reference value 
    #plt.text(trialN-1, RMS_S_rad[i][trialN-1]+0.02, str(round(RMS_S_rad[i][trialN-1], 3))+r" $\pm$ "+str(round(RMS_S_rad_error[i][trialN-1],3)), fontsize=9) # plot the reference value 
    plt.scatter(cases, RMS_S_rad[i],  marker=str(marker[i]), color=str(color[i]), label=str(label[i])) # plot all cases and reference points 
    plt.errorbar(cases, RMS_S_rad[i], yerr=RMS_S_rad_error[i],  color=str(color[i]), markersize=12, elinewidth=2, linewidth=0) # add error bar

axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.title("Radial Position case " + str(dirName) , fontsize=10)
plt.xlabel("Sample #", fontsize=10)
plt.ylabel(r"$\Delta$ Width Radial [mm]", fontsize=10)

plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)

plt.subplot(212) # VERTICAL 
axes=plt.gca()
axes.xaxis.set_major_locator(MaxNLocator(integer=True))
axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')

for i in range(0, 2):
    #line = [[1,RMS_S_ver[i][0]], [trialN, RMS_S_ver[i][0]]]
    #plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = str(color[i]), linestyle="--")
    #plt.text(0, RMS_S_ver[i][0]+0.02, str(round(RMS_S_ver[i][0], 3))+r" $\pm$ "+str(round(RMS_S_ver_error[i][0],3)), fontsize=9)
    #plt.text(trialN-1, RMS_S_ver[i][trialN-1]+0.02, str(round(RMS_S_ver[i][trialN-1], 3))+r" $\pm$ "+str(round(RMS_S_ver_error[i][trialN-1],3)), fontsize=9)
    plt.scatter(cases, RMS_S_ver[i],  marker=str(marker[i]), color=str(color[i]), label=str(label[i]))
    plt.errorbar(cases, RMS_S_ver[i], yerr=RMS_S_ver_error[i],  color=str(color[i]), markersize=12, elinewidth=2, linewidth=0 )

axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.title("Vertical Position case "+str(dirName), fontsize=10)
plt.xlabel("Sample #", fontsize=10)
plt.ylabel(r"$\Delta$ Width Vertical [mm]", fontsize=10)

plt.subplots_adjust(bottom=0.1)
plt.subplots_adjust(hspace=0.43)
plt.subplots_adjust(right=0.8)

# #plt.show()
plt.savefig("xSD_"+str(dirName)+".png", dpi=600)