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

    fhiclPath = "/gm2/app/users/glukicov/TrackerAlignment/gm2Dev_v9_14_00/srcs/gm2tracker/align/Systematics/"+str(misalignment)+"/"+str(i_trial+1)+"/RunTrackingDAQ.fcl"
    outDir = "/pnfs/GM2/scratch/users/glukicov/Systematics/"+str(misalignment)+"/"+str(i_trial+1)

    subprocess.call(["./gridSetupAndSubmitGM2Data.sh", "--daq", "--reco", "--fhiclFile="+str(fhiclPath), "--localArea", "--output-dir="+str(outDir), "--sam-dataset=Align_15922_unpacked", "--njobs=100", "--offsite", "--nevents=50" ])    


print(str(trialN)+" jobs submitted to the grid!")