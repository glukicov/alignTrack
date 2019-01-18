################
# Script that generates a set of random misalignments (X, Y)
# and adds them to the FHICL file for tracking with Data 
#
# Gleb Lukicov 27 December 2018 
############

import argparse, sys
import subprocess

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-mis", "--mis") # offset scale (SD) [um]
parser.add_argument("-trialN", "--trialN") # number of iterations
args = parser.parse_args()

# CONSTANTS 
trialN = int(args.trialN) 
misalignment = int(args.mis)


for i_trial in range(0, trialN):

    fhiclPath = "/gm2/app/users/glukicov/TrackerAlignment/gm2Dev_v9_14_00/srcs/gm2tracker/align/RunTrackingPlots.fcl"
    scrDir = "/pnfs/GM2/scratch/users/glukicov/Systematics/"+str(misalignment)+"/"+str(i_trial+1)+"/*/data/gm2tracker*.root"
    outDir = "/pnfs/GM2/scratch/users/glukicov/Systematics_Plots/"+str(misalignment)+"/"+str(i_trial+1)+"/"

    subprocess.call(["./gridSetupAndSubmitGM2Data.sh", "--daq", "--ana", "--fhiclFile="+str(fhiclPath), "--localArea", "--output-dir="+str(outDir), "--sam-dataset=Align_"+str(misalignment)+"_"+str(i_trial+1), "--njobs=100", "--offsite" ])    


print(str(trialN)+" jobs submitted to the grid!")