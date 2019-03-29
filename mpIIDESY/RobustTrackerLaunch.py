####################################################################
# Quality control plots for Tracker Alignment: comparison of PEDE results 
#
# Works with data for MC alignment results
# Expects to have offsets written from previous iterations (if used iteratively)
#
# Created: 28 March 2019 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 28 March 2019 by Gleb
#####################################################################
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import argparse 
import re

### HELPER FUNCTIONS ######
## Read offsets from FHICL file 
def getOffsets(f, name):
    offsets = [] #tmp storage buffer 
    for line in f:
        if re.match(name, line):
            copy = True
            offsets=line 
        else:
            copy = False
    return offsets 

#Define and open command line arguments
parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-pF', '--pede_file', help='input pede results (.res) file', default="millepede.res", type=str)  # input PEDE file  
parser.add_argument('-oF', '--offset_file', help='FHICL file used for tracking with offsets', default="RunTrackingDAQ_align.fcl", type=str)  # input PEDE file  
parser.add_argument('-tF', '--truth_file', help='FHICL file used for MDC1 generation', default="runGasGunRing_align.fcl", type=str)  # input PEDE file  
parser.add_argument('-eL', '--extra_label', help='extra plotting label', default="", type=str) # if extra label is passed for plot tittle (e..g iteration #)
args = parser.parse_args()
pede_file_name=args.pede_file
offset_file_name=args.offset_file
truth_file_name=args.truth_file
extra_label_name=args.extra_label

#Define constants (Capital + cammelCase)
ModuleN = 8 # per station 
ModuleArray=np.arange(1, ModuleN+1) #(1, 2,...,8) for plotting  
StrawModuleZPosition = [ 0.0, 134.36, 268.72, 403.08, 537.42, 671.77, 806.10, 940.406] # station CS 
GlobalParNames = ["Radial", "Vertical", "\theta", "\phi", "\psi"] #only ever going to have 5 pars.
GlobalParLabels = [1, 2, 3, 4, 5] # their label
GlobalParDict = dict(zip(GlobalParLabels, GlobalParNames))
FHICLPatchName = ["strawModuleRShift", "strawModuleHShift", None, None, None] # no FHICL patch for angles in art yet  
FHICLServicePath = "services.Geometry.strawtracker." #full path of the tracking FHICL patch 

#Define variables that will be PEDE result-dependent 
globalN = -1 # per module
stationN = -1 # 12 or 18 
useOffsets=False # no offsets unless set
useTruth=False # no truth information unless set 
pars = [] # labels (e.g. 1211 = "S12 M1 Radial")
pede_results = [] #alignment results from PEDE
pede_errors = [] #alignment results errors 
 
# Read PEDE results, skipping the title row 
df = pd.read_csv(pede_file_name, header=None, skiprows=1)

#establish how many alignment parameters were used per module 
globalN = int((len(df[0]))/ModuleN)
parN = globalN * ModuleN  

#Loop through results and store labels, shifts and errors in um 
for i_par in range(0, parN):
    lineString = df[0][i_par] 
    arrayString = [str(i) for i in lineString.split()] # remove spaces 
    pars.append(int(arrayString[0]))
    pede_results.append(round(float(arrayString[1])*1000)) #mm to the nearest um 
    pede_errors.append(round(float(arrayString[4])*1000)) #mm to the nearest um 

# combine into a data structure 
results=[pars, pede_results, pede_errors]
# structure:  data[i_global][i_module][i_results] with result: 0=label, 1=align 2= error 
data = [[[0 for i_result in range(len(results))] for i_module in range(ModuleN)] for i_global in range(globalN)]
for i_global in range(0, globalN):
    for i_module in range(0, ModuleN):
          for i_result in range(0, len(results)):
            data[i_global][i_module][i_result]=results[i_result][i_global::2][i_module] 

#Loop through parameters to establish names ( full label = AA (station) + B (module) + C (par. label) )
stationN=str(pars[0])[0:2]
print("Alignment results for Station:", stationN)
all_labels=[]
for i_par in pars:
    i_par=str(i_par)
    all_labels.append(i_par[3]) # par. label only 
seen = set() 
unique_labels = [x for x in all_labels if x not in seen and not seen.add(x)] # only unseen
sys.stdout.write("Aligning for ")
for global_par in unique_labels:
    sys.stdout.write(GlobalParDict[int(global_par)])
    sys.stdout.write("; ")
sys.stdout.write(" shifts\n")
print("Number of global parameters per module:", globalN)

#Print the PEDE results and errors 
for i_global in range(0, globalN):
    sys.stdout.write(FHICLPatchName[i_global]+stationN+": [ ")
    for i_module in range(0, ModuleN):
        sys.stdout.write(str(data[i_global][i_module][1]))
        if(i_module != 7):
            sys.stdout.write(", ")
    sys.stdout.write(" ] \\\ [um]\n")
for i_global in range(0, globalN):
    sys.stdout.write(FHICLPatchName[i_global]+stationN+" error: [ ")
    for i_module in range(0, ModuleN):
        sys.stdout.write(str(data[i_global][i_module][2]))
        if(i_module != 7):
            sys.stdout.write(", ")
    sys.stdout.write(" ] \\\ [um]\n")

#Check if offsets and/or truth misalignment was loaded from FHICL files, and adjust data accordingly 
offset_file = Path(offset_file_name)
print(offset_file_name)
if offset_file.is_file():
    useOffsets=True
    offset_file=open(offset_file_name, "r")
    for i_global in range(globalN):
        offset_input = (getOffsets(offset_file, FHICLServicePath+FHICLPatchName[i_global]))
        print(offset_input)
        # offset_input = offset_input.replace(FHICLServicePath+FHICLPatchName[i_global]+" [", " ") 
        # offset_input = offset_input.replace("]", "") 
        # offset_input = np.array([float(r) for r in offset_input.split(',')])
        # offset_input=np.array(radialOff[0:8])
        # offset_input=round(float(radialOff[1])*1000)
        # print("offset_input:", offset_input)
print("Offset alignment used for iteration (required if iterating!):", useOffsets)

truth_file = Path(truth_file_name)
if truth_file.is_file():
    useTruth=True
print("Truth alignment used for comparison (simulation only):", useTruth)
