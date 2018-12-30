################
# Script to combine many files into 1 
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

    fhiclPath = "/gm2/app/users/glukicov/TrackerAlignment/gm2Dev_v9_14_00/srcs/gm2tracker/align/macros/mergeTFiles.fcl"
    outDir = "/pnfs/GM2/scratch/users/glukicov/Systematics_Plots_2/"+str(misalignment)+"/"+str(i_trial+1)+"/"

    subprocess.call(["cd", str(outDir)]
    subprocess.call(["gm2", "-c", str(fhiclPath), "-s", str(outDir)]


print(str(trialN)+" jobs submitted to the grid!")