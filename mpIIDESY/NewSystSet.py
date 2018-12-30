################
# Script that generates a set of random misalignments (X, Y)
# and adds them to the FHICL file for tracking with Data 
#
# Gleb Lukicov 27 December 2018 
############

import argparse, sys
import subprocess
import numpy as np 
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-mis", "--mis") # offset scale (SD) [um]
parser.add_argument("-trialN", "--trialN") # number of iterations
parser.add_argument("-out", "--outFile", help='Output file')
args = parser.parse_args()

# CONSTANTS 
moduleN = 8 #number of movable detectors
stationN = 3 # S0, S12, S18 
globalN = 2 # X, Y
trialN = int(args.trialN) 
misalignment = int(args.mis)
np.random.seed(12345+misalignment)
outFile = str(args.outFile) 

# Global scope containers
MeanMisX=[]
MeanMisY=[] 
SDMisX=[]
SDMisY=[]

#Quickly open the output and check that no previous offsets have been written,
with open(outFile) as f:
     if "services.Geometry.strawtracker.radialOffsetPerModule" in f.read():
         print("The output file already contains offsets - check the correct file is passed/manually backup and delete old offsets!")
         sys.exit()
     if "services.Geometry.strawtracker.verticalOffsetPerModule:" in f.read():
         print("The output file already contains offsets - check the correct file is passed/manually backup and delete old offsets!")
         sys.exit()

#new offsets, dir and FHICL file for each trial 
for i_trial in range(0, trialN):
    
    #Un-smeared for S0 
    misX = np.zeros(moduleN)
    misY = np.zeros(moduleN)

    #Smeared by:
    mu, sigma = 0, 1 # mean and standard deviation * misalignment scale + offset (same for all modules)
    #no overall offsets for first 15 trials 
    offsetX=np.zeros(moduleN)
    offsetY=np.zeros(moduleN)
    if (trialN > 15):
        offsetX=np.ones(moduleN) * np.random.uniform(150, -150)
        offsetY=np.ones(moduleN) * np.random.uniform(150, -150)
    for i in range(0, moduleN):
         misX[i]=offsetX[i]+misalignment*np.random.normal(mu, sigma)
         misY[i]=offsetY[i]+misalignment*np.random.normal(mu, sigma)

    #appends same offsets for S12 and S18 
    a = np.concatenate((np.zeros(moduleN), misX), axis=0)
    b= np.concatenate((np.zeros(moduleN), misY), axis=0)
    radialOffsetPerModule = np.concatenate((a, misX), axis=0)
    verticalOffsetPerModule = np.concatenate((b, misY), axis=0)

    # truncate to the nearest um 
    radialOffsetPerModule = radialOffsetPerModule.astype(int)
    verticalOffsetPerModule = verticalOffsetPerModule.astype(int)

    #store stats for this trial
    MeanMisX.append(np.mean(misX))
    MeanMisY.append(np.mean(misY))
    SDMisX.append(np.std(misX)) # SD of the population 
    SDMisY.append(np.std(misY))

    # print("rad", radialOffsetPerModule)
    # print("ver", verticalOffsetPerModule)

    '''
    #for each trial: (1) create new dir and (2) copy the FHICL file there 
    subprocess.call(["mkdir" , str(i_trial+1)])
    subprocess.call(["cp" , str(outFile), str(i_trial+1)+"/"])

    #write to file: a+ allows to append at the end of the file
    # !! FHICL file expects offsets in mm 
    f = open(str(i_trial+1)+"/"+outFile, 'a+')
    f.write("\n\n//Straw Offsets\n")
    f.write("services.Geometry.strawtracker.radialOffsetPerModule: [ ")
    for item in radialOffsetPerModule[:-1]:
        f.write( str( round(item*1e-3, 3) ) + ", ") # um -> mm 
    f.write( str( round(radialOffsetPerModule[-1]*1e-3, 3) ) )
    f.write(" ]\n")
    f.write("services.Geometry.strawtracker.verticalOffsetPerModule: [ ")
    for item in verticalOffsetPerModule[:-1]:
        f.write( str( round(item*1e-3, 3) ) + ", " )
    f.write( str( round(verticalOffsetPerModule[-1]*1e-3, 3) ) )
    f.write(" ]\n")
    f.close()
    '''

## plot for all trials
plt.subplot(211) 
bins = np.linspace(-200, 200, 50)
plt.hist(MeanMisX, bins, alpha=0.5, label='<X>')
plt.hist(MeanMisY, bins, alpha=0.5, label='<Y>')
plt.legend(loc='upper right', title="<Misalignment>")
plt.text(120, 1, "Entries X: "+str(len(MeanMisX))+"\nMean X: "+str(round(np.mean(MeanMisX),3))+"\nSD X: "+str(round(np.std(MeanMisX),3)), bbox=dict(facecolor='green', alpha=0.5))
plt.text(20, 1, "Entries Y: "+str(len(MeanMisY))+"\nMean Y: "+str(round(np.mean(MeanMisY),3))+"\nSD Y: "+str(round(np.std(MeanMisY),3)), bbox=dict(facecolor='red', alpha=0.5))
plt.xlabel("Mean [um]", fontsize=10)
plt.subplots_adjust(bottom=0.1)

plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=1.0)

plt.subplot(212) 
bins = np.linspace(0, 300, 50)
plt.hist(SDMisX, bins, alpha=0.5, label='SD(X)')
plt.hist(SDMisY, bins, alpha=0.5, label='SD(Y)')
plt.legend(loc='upper right', title="Misalignment SD")
plt.text(240, 2, "Entries X: "+str(len(SDMisX))+"\nMean X: "+str(round(np.mean(SDMisX),3))+"\nSD X: "+str(round(np.std(SDMisX),3)), bbox=dict(facecolor='green', alpha=0.5))
plt.text(165, 2, "Entries Y: "+str(len(SDMisY))+"\nMean Y: "+str(round(np.mean(SDMisY),3))+"\nSD Y: "+str(round(np.std(SDMisY),3)), bbox=dict(facecolor='red', alpha=0.5))
plt.xlabel("Standard Deviation [um]", fontsize=10)

plt.subplots_adjust(bottom=0.1)
plt.subplots_adjust(left=0.1)
plt.subplots_adjust(hspace=0.35)
plt.savefig(str(misalignment)+"_stats.png", dpi=600)
# plt.show()


