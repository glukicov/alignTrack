################
# Script that generates a set of random misalignments (X, Y)
# and adds them to the FHICL file for tracking with Data 
#
# Gleb Lukicov 27 December 2018 
############

import argparse, sys
import subprocess, shlex

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-mis", "--mis") # offset scale (SD) [um]
parser.add_argument("-trialN", "--trialN") # number of iterations
args = parser.parse_args()

# CONSTANTS 
trialN = int(args.trialN) 
misalignment = int(args.mis)


# for i_trial in range(6, trialN):

#     scrDir = "/pnfs/GM2/scratch/users/glukicov/Systematics/"+str(misalignment)+"/"+str(i_trial+1)+"/*/data/"

#     cmd = "sam_add_dataset -d "+str(scrDir)+" -n Align_"+str(misalignment)+"_"+str(i_trial+1)

#     print("cmd=",cmd) 

#     args = shlex.split(cmd)

#     print("args=",args) 
    
#     #subprocess.Popen(str(cmd))

#     subprocess.Popen(args)

# print(str(trialN)+" jobs submitted to the grid!")


cmd = "sam_add_dataset -d /pnfs/GM2/scratch/users/glukicov/Systematics/25/6/2018-12-28-13-21-06/data/ -n Align_25_6"

print("cmd=",cmd) 

args = shlex.split(cmd)

print("args=",args) 

#subprocess.Popen(str(cmd))

subprocess.Popen(args)

