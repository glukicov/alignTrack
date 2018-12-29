################
# Script that produces the final FoM 
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

    scrFile = "/pnfs/GM2/scratch/users/glukicov/Systematics_Plots/"+str(misalignment)+"/"+str(i_trial+1)+"/trackRecoPlots.root"

    subprocess.call(["./gridSetupAndSubmitGM2Data.sh", "--daq", "--ana", "--fhiclFile="+str(fhiclPath), "--localArea", "--output-dir="+str(outDir), "--sam-dataset=Align_"+str(misalignment)+"_"+str(i_trial+1)", "--njobs=1", "--offsite" ])    


print(str(trialN)+" cases analysed!")