################
# Script that produces the final FoM 
#
# Gleb Lukicov 27 December 2018 
############

import argparse, sys, os, glob
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

label=("S12", "S18")
color=("purple", "orange")

trialN = len(glob.glob1(".","fom_*.txt"))

list_of_files=[]
for file in glob.glob("*.txt"):
    list_of_files.append(file)

print("Found "+str(trialN)+ " fom files.")
print(list_of_files)


### TODO horrible hack : replace by Pandas Data Frame.... :/ 

misalignment=[]

S12_S_rad_mean=[]
S12_S_rad_error_mean=[]
S12_S_rad_max=[]
S12_S_rad_error_max=[]
S12_RMS_S_rad_mean=[]
S12_RMS_S_rad_error_mean=[]
S12_RMS_S_rad_max=[]
S12_RMS_S_rad_error_max=[]

S18_S_rad_mean=[]
S18_S_rad_error_mean=[]
S18_S_rad_max=[]
S18_S_rad_error_max=[]
S18_RMS_S_rad_mean=[]
S18_RMS_S_rad_error_mean=[]
S18_RMS_S_rad_max=[]
S18_RMS_S_rad_error_max=[]

for i in range(0, trialN):

    f=open(list_of_files[i], "r")
    line = f.readline()
    misalignment.append(int(line))
    
    line = f.readline()
    S12_S_rad_mean.append(int(line))
    line = f.readline()
    S18_S_rad_mean.append(int(line))

    line = f.readline()
    S12_S_rad_error_mean.append(int(line))
    line = f.readline()
    S18_S_rad_error_mean.append(int(line))
    
    line = f.readline()
    S12_S_rad_max.append(int(line))
    line = f.readline()
    S18_S_rad_max.append(int(line))
    
    line = f.readline()
    S12_S_rad_error_max.append(int(line))
    line = f.readline()
    S18_S_rad_error_max.append(int(line))
    
    line = f.readline()
    S12_RMS_S_rad_mean.append(int(line))
    line = f.readline()
    S18_RMS_S_rad_mean.append(int(line))
    
    line = f.readline()
    S12_RMS_S_rad_error_mean.append(int(line))
    line = f.readline()
    S18_RMS_S_rad_error_mean.append(int(line))
    
    line = f.readline()
    S12_RMS_S_rad_max.append(int(line))
    line = f.readline()
    S18_RMS_S_rad_max.append(int(line))
    
    line = f.readline()
    S12_RMS_S_rad_error_max.append(int(line))
    line = f.readline()
    S18_RMS_S_rad_error_max.append(int(line))

    
    f.close()


S_rad_mean=[S12_S_rad_mean, S18_S_rad_mean]
S_rad_error_mean=[S12_S_rad_error_mean, S18_S_rad_error_mean]

S_rad_max=[S12_S_rad_max, S18_S_rad_max]
S_rad_error_max=[S12_S_rad_error_max, S18_S_rad_error_max]

RMS_S_rad_mean=[S12_RMS_S_rad_mean, S18_RMS_S_rad_mean]
RMS_S_rad_error_mean=[S12_RMS_S_rad_error_mean, S18_RMS_S_rad_error_mean]

RMS_S_rad_max=[S12_RMS_S_rad_max, S18_RMS_S_rad_max]
RMS_S_rad_error_max=[S12_RMS_S_rad_error_max, S18_RMS_S_rad_error_max]

S_rad=[S_rad_mean, S_rad_max, RMS_S_rad_mean, RMS_S_rad_max]
S_rad_error=[S_rad_error_mean, S_rad_error_max, RMS_S_rad_error_mean, RMS_S_rad_error_max]

y_title = ["|Vertical Mean| [um]", "|Vertical Max| [um]", "|SD Mean| [um]", "|SD Mean| [um]"]

# Fill Plots 

subplotPos=[221, 222, 223, 224]

plt.figure(1)
for i_case in range(0, 4):

    plt.subplot(subplotPos[i_case]) # RADIAL 
    axes=plt.gca()
    axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
    axes.yaxis.set_major_formatter(FormatStrFormatter("%.0f")) # 1 decimal point of y-axis 
    axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

    for i_station in range(0, 2):
        plt.scatter(misalignment, S_rad[i_case][i_station],  marker="+", color=str(color[i_station])) # plot all cases and reference points 
        plt.errorbar(misalignment, S_rad[i_case][i_station], yerr=S_rad_error[i_case][i_station],  color=str(color[i_station]), markersize=14, elinewidth=3, linewidth=0, label=str(label[i_station])) # add error bar
        plt.plot(np.unique(misalignment), np.poly1d(np.polyfit(misalignment, S_rad[i_case][i_station], 1))(np.unique(misalignment)), color=str(color[i_station]))

    plt.legend(loc='best')
    plt.xlabel("Misalignment [um]", fontsize=8)
    plt.ylabel(y_title[i_case], fontsize=8)

    plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)
    plt.subplots_adjust(bottom=0.1)
    plt.subplots_adjust(hspace=0.43)
    plt.subplots_adjust(right=0.9)

plt.savefig("fom.png", dpi=600)

# def getOffsets(f, name):
#     offsets = [] #tmp storage buffer 
#     for line in f:
#         if re.match(name, line):
#             copy = True
#             offsets=line 

#         else:
#             copy = False

#     return offsets   

# def flatten(input_array):
#     result_array = []
#     for element in input_array:
#         if isinstance(element, int):
#             result_array.append(element)
#         elif isinstance(element, list):
#             result_array += flatten(element)
#     return result_array

# trialN = len(glob.glob1(".","*.csv"))
# print("Found "+str(trialN)+ " metric files.")

# list_of_files=[]
# for file in glob.glob("*.csv"):
#     list_of_files.append(file)


# Mis_Mean_rad=[]
# Mis_SD_rad=[]
# S12_dm=[]
# S12_ds=[]
# S18_dm=[]
# S18_ds=[]

# for i in range(0, trialN):
#     df = pd.read_csv(list_of_files[i])
#     df = df.values
    
#     Mis_Mean_rad.append(df[0][1:])
#     Mis_SD_rad.append(df[1][1:])
#     S12_dm.append(df[2][1:])
#     S12_ds.append(df[3][1:])
#     S18_dm.append(df[4][1:])
#     S18_ds.append(df[5][1:])

# dm=[S12_dm, S18_dm]
# ds=[S12_ds, S18_ds]

# for i_trial in range(0, trialN):
 
#     # Plot for covariance
#     plt.figure(3) 
#     plt.suptitle("Correlation plots" , fontsize=10)
#     plt.subplot(221) # Mean Beam vs Mean Mis  
#     axes=plt.gca()
#    # axes.set_xlim(min(int(Mis_Mean_rad[i_trial]))*0.9, max(int(Mis_Mean_rad[i_trial]))*1.1)
#     # axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
#     axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
#     axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

#     for i in range(0, 2):
#         dM = dm[i][i_trial] # take away the nominal per station 
#         plt.scatter(Mis_Mean_rad[i_trial], dM,  marker="+", color=str(color[i])) # plot all cases and reference points 
#         #plt.errorbar(Mis_Mean_rad[i_trial], dM, yerr=S_rad_error[i][1:],  color=str(color[i]), markersize=14, elinewidth=3, linewidth=0, label=str(label[i])) # add error bar
#         cor = np.corrcoef(Mis_Mean_rad[i_trial], dM)
#         plt.plot([], [], ' ', label=str(label[i])+" corr: "+str(round(cor[0][1],3)))
#         plt.plot(np.unique(Mis_Mean_rad[i_trial]), np.poly1d(np.polyfit(Mis_Mean_rad[i_trial], dM, 1))(np.unique(Mis_Mean_rad[i_trial])), color=str(color[i]))

#     plt.legend(loc='center')
#     #plt.title("d(Radial Position) vs. <Misalignment> " + str(dirName) , fontsize=10)
#     plt.xlabel("<Misalignment> [um]", fontsize=10)
#     plt.ylabel("d(Radial Position) [mm]", fontsize=10)

#     # plt.subplot(222) # Mean Beam vs SD Mis  
#     # axes=plt.gca()
#     # axes.set_xlim(min(Mis_SD_rad)*0.9, max(Mis_SD_rad)*1.1)
#     # # axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
#     # axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
#     # axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

#     # for i in range(0, 2):
#     #     dM = np.array(S_rad[i][1:])-S_rad[i][0] # take away the nominal per station 
#     #     plt.scatter(Mis_SD_rad, dM,  marker="+", color=str(color[i])) # plot all cases and reference points 
#     #     plt.errorbar(Mis_SD_rad, dM, yerr=S_rad_error[i][1:],  color=str(color[i]), markersize=14, elinewidth=3, linewidth=0, label=str(label[i])) # add error bar
#     #     cor = np.corrcoef(Mis_SD_rad, dM)
#     #     plt.plot([], [], ' ', label=str(label[i])+" corr: "+str(round(cor[0][1],3)))
#     #     plt.plot(np.unique(Mis_SD_rad), np.poly1d(np.polyfit(Mis_SD_rad, dM, 1))(np.unique(Mis_SD_rad)), color=str(color[i]))

#     # plt.legend(loc='center')
#     # #plt.title("d(Radial Position) vs. SD(Misalignment) " + str(dirName) , fontsize=10)
#     # plt.xlabel("SD(Misalignment) [um]", fontsize=10)
#     # plt.ylabel("d(Radial Position) [mm]", fontsize=10)


#     # plt.subplot(223) # SD Beam vs Mean Mis  
#     # axes=plt.gca()
#     # axes.set_xlim(min(Mis_Mean_rad)*0.9, max(Mis_Mean_rad)*1.1)
#     # # axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
#     # axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
#     # axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

#     # for i in range(0, 2):
#     #     dS = np.array(RMS_S_rad[i][1:])-RMS_S_rad[i][0] # take away the nominal per station 
#     #     dS_all.append(dS)
#     #     plt.scatter(Mis_Mean_rad, dS,  marker="+", color=str(color[i])) # plot all cases and reference points 
#     #     plt.errorbar(Mis_Mean_rad, dS, yerr=RMS_S_rad_error[i][1:],  color=str(color[i]), markersize=14, elinewidth=3, linewidth=0, label=str(label[i])) # add error bar
#     #     cor = np.corrcoef(Mis_Mean_rad, dS)
#     #     plt.plot([], [], ' ', label=str(label[i])+" corr: "+str(round(cor[0][1],3)))
#     #     plt.plot(np.unique(Mis_Mean_rad), np.poly1d(np.polyfit(Mis_Mean_rad, dS, 1))(np.unique(Mis_Mean_rad)), color=str(color[i]))

#     # plt.legend(loc='center')
#     # #plt.title("d(Radial Position) vs. <Misalignment> " + str(dirName) , fontsize=10)
#     # plt.xlabel("<Misalignment> [um]", fontsize=10)
#     # plt.ylabel("d(SD(Radial Position)) [mm]", fontsize=10)


#     # plt.subplot(224) # SD Beam vs SD Mis  
#     # axes=plt.gca()
#     # axes.set_xlim(min(Mis_SD_rad)*0.9, max(Mis_SD_rad)*1.1)
#     # # axes.xaxis.set_major_locator(MaxNLocator(integer=True)) # int ticks only on x-axis 
#     # axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f")) # 1 decimal point of y-axis 
#     # axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout') #in-out ticks on y axis 

#     # for i in range(0, 2):
#     #     dS = np.array(RMS_S_rad[i][1:])-RMS_S_rad[i][0] # take away the nominal per station
#     #     plt.scatter(Mis_SD_rad, dS,  marker="+", color=str(color[i])) # plot all cases and reference points 
#     #     plt.errorbar(Mis_SD_rad, dS, yerr=RMS_S_rad_error[i][1:],  color=str(color[i]), markersize=14, elinewidth=3, linewidth=0, label=str(label[i])) # add error bar
#     #     cor = np.corrcoef(Mis_SD_rad, dS)
#     #     plt.plot([], [], ' ', label=str(label[i])+" corr: "+str(round(cor[0][1],3)))
#     #     plt.plot(np.unique(Mis_SD_rad), np.poly1d(np.polyfit(Mis_SD_rad, dS, 1))(np.unique(Mis_SD_rad)), color=str(color[i]))
        
#     # plt.legend(loc='center')
#     # #plt.title("d(Radial Position) vs. <Misalignment> " + str(dirName) , fontsize=10)
#     # plt.xlabel("SD(Misalignment) [um]", fontsize=10)
#     # plt.ylabel("d(SD(Radial Position)) [mm]", fontsize=10)

#     # plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)
#     # plt.subplots_adjust(top=0.9)
#     # plt.subplots_adjust(bottom=0.1)
#     # plt.subplots_adjust(hspace=0.43)
#     # # plt.subplots_adjust(right=0.9)

#     plt.savefig("Corr_Final.png", dpi=600)

#     print(str(trialN)+" cases analysed!")

