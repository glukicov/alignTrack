####################################################################
# Summary plots for Tracker Alignment: comparison of PEDE results 
#
# Works with data for MC alignment results
#
# Created: 28 March 2019 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 28 March 2019 by Gleb
#####################################################################
import sys # print out to terminal 
from pathlib import Path #check files
import pandas as pd # get data frame from text file 
import numpy as np # arrays 
import argparse # command line inputs sub
import re # to get offsets from file 
import numpy.polynomial.polynomial as poly # linear fit 
import matplotlib.pyplot as plt #for plotting  
import itertools # smart lines in plotting 


### HELPER FUNCTIONS ######
## Read offsets from FHICL file 
def getOffsets(f, name):
    offsets = [] #tmp storage buffer 
    for line in f:
        if re.match(name, line):
            offsets=line 
        else:
            pass
    
    offsets = offsets.replace(name+": [", "") 
    offsets = offsets.replace("]", "") 
    offsets = np.array([float(r) for r in offsets.split(',')])
    offsets = np.array(offsets*1000) # convert to um 
    offsets = offsets.astype(int) 
    return offsets 

#Define and open command line arguments
parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-p', '--path', nargs='+', help="List of directories", type=str)  # iterations  input 
parser.add_argument('-tF', '--truth_file', help='FHICL file used for MDC1 generation', default="runGasGunRing_align.fcl", type=str)  # input PEDE file  
parser.add_argument('-eL', '--extra_label', help='extra plotting label', default="", type=str) # if extra label is passed for plot tittle (e..g iteration #)
args = parser.parse_args()
path=args.path
truth_file_name=args.truth_file
extra_label_name=args.extra_label

print("Summary plots for", len(path), "iterations")
print(path)
colors = ["red", "orange", "purple", "green", "black", "yellow", "pink", "grey"]
labels= ["Truth", "Iter. 1", "Iter. 2", "Iter. 3", "Iter. 4", "Iter. 5", "Iter. 6", "Iter. 7"]

#Get constants from the 1st iteration
ModuleN = 8 # per station 
ModuleArray=np.arange(1, ModuleN+1) #(1, 2,...,8) for plotting  
GlobalParNames = ["Radial", "Vertical", 'Φ', 'ψ', 'θ'] #only ever going to have 5 pars.
units = [r" [$\mathrm{\mu m}$]", r" [$\mathrm{\mu m}$]", " [mrad]", " [mrad]", " [mrad]"]
GlobalParLabels = [1, 2, 3, 4, 5] # their label
GlobalParDict = dict(zip(GlobalParLabels, GlobalParNames))
FHICLPatchName = ["strawModuleRShift", "strawModuleHShift"," strawModuleRotationPhi", "strawModuleRotationPsi", "strawModuleRotationTheta"] # no FHICL patch for angles in art yet  
FHICLServicePath = "services.Geometry.strawtracker." #full path of the tracking FHICL patch 

#Define variables that will be PEDE result-dependent 
globalN = -1 # per module
stationN = -1 # 12 or 18 
pars = [] # labels (e.g. 1211 = "S12 M1 Radial")

# Read PEDE results, skipping the title row 
df = pd.read_csv(path[0]+"/millepede.res", header=None, skiprows=1)
#establish how many alignment parameters were used per module 
globalN = int((len(df[0]))/ModuleN)
parN = globalN * ModuleN  

#Loop through results and store labels, shifts and errors in um 
for i_par in range(0, parN):
    lineString = df[0][i_par] 
    arrayString = [str(i) for i in lineString.split()] # remove spaces 
    pars.append(int(arrayString[0]))

#Loop through parameters to establish names ( full label = AA (station) + B (module) + C (par. label) )
stationN=str(pars[0])[0:2] # XXX check logic 
print("Alignment results for Station:", stationN, " | All results in [um]")
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

 ### Load truth misalignment #####
truth = [[0 for i_module in range(ModuleN)] for i_global in range(globalN)]
for i_global in range(0, 2):
    truth_file=open(path[0]+"/"+truth_file_name, "r")
    truth_input = getOffsets(truth_file, FHICLServicePath+FHICLPatchName[i_global]+stationN)
    truth[i_global]=truth_input
    print("Truth "+FHICLPatchName[i_global]+stationN+" :", truth[i_global])

# For all iterations load the offsets
offset = [[[0 for i_module in range(ModuleN)] for i_global in range(globalN)] for i_iter in range(0, len(path))]
for i_iter in range(0, len(path)):
   
   for i_global in range(0, 2):
    offset_file=open(path[i_iter]+"/"+"OffsetsPerModule"+stationN+".fcl", "r")
    offset_input = getOffsets(offset_file, FHICLServicePath+FHICLPatchName[i_global]+stationN)
    offset[i_iter][i_global]=offset_input
    print("Offset "+FHICLPatchName[i_global]+stationN+" :", offset[i_iter][i_global])


######### Plotting #############
# Plotting "constants"
f = plt.figure(figsize=(7,int(globalN*2+1)))
#TOOD make min max from data
yMax = [120, 120, 25, 25, 25]
yMin = [-120, -120, -25, -25, -25]
bar_width, bar_offset, opacity = 0.15, 0.3, 0.8
#Make subplot for each result 
for i_iter in range(0, len(path)):
    for i_global in range(0, globalN):
        plt.subplot(int( str(globalN)+"1"+str(int(i_global+1)) )) 
        axes = plt.gca()
        axes.set_xlim(ModuleArray[0]-0.5, ModuleArray[-1]+0.5)
        axes.set_ylim(yMin[i_global], yMax[i_global])
        plt.title(GlobalParNames[i_global]+" residuum misalignment in S"+stationN, fontsize=12)
        plt.ylabel(r"Residuum misalignment "+ units[i_global])
        plt.xlabel("Module", fontsize=12)
        plt.xticks(fontsize=10, rotation=0) 
        plt.xticks(ModuleArray + bar_width, ModuleArray)
        plt.yticks(fontsize=10, rotation=0)
        plt.minorticks_on()
        axes.tick_params(axis='x',which='minor',bottom=False)
        axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
        plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)
        # #Plot the 0th line 
        line = [[ModuleArray[0]-0.5,0.0], [ModuleArray[-1]+0.5, 0.0]]
        plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
        # #Plot module lines
        for i_module in range(0, 8):
             line = [[i_module+0.5, yMin[i_global]], [i_module+0.5, yMax[i_global]]]
             plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
        # #Plot data 
        if (i_iter == 0):
            rects1 = plt.bar(ModuleArray-bar_offset, truth[i_global], bar_width, alpha=opacity, color=colors[0],label=labels[0])
        rects2 = plt.bar(ModuleArray + bar_width*(i_iter+1)-bar_offset, offset[i_iter][i_global]-truth[i_global], bar_width, alpha=opacity, color=colors[i_iter+1], label=labels[i_iter+1])
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 


plt.savefig("Summary_S"+stationN+".png", dpi=250)

print("\nSummary results in Summary_S"+stationN+".png\n")