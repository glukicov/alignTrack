####################################################################
# Quality control plots for Tracker Alignment: comparison of PEDE results 
#
# Works with data for MC alignment results
# Expects to have offsets written from previous iterations (if used iteratively)
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
parser.add_argument('-pF', '--pede_file', help='input pede results (.res) file', default="millepede.res", type=str)  # input PEDE file  
parser.add_argument('-oF', '--offset_file', help='FHICL file used for tracking with offsets', default="RunTrackingSim_align.fcl", type=str)  # input PEDE file  
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

#Print the PEDE results and errors 
for i_global in range(0, globalN):
    sys.stdout.write("PEDE result for "+FHICLPatchName[i_global]+stationN+": [ ")
    for i_module in range(0, ModuleN):
        sys.stdout.write(str(data[i_global][i_module][1]))
        if(i_module != 7):
            sys.stdout.write(" ")
    sys.stdout.write(" ] \n")
for i_global in range(0, globalN):
    sys.stdout.write("PEDE result for "+FHICLPatchName[i_global]+stationN+" error: [ ")
    for i_module in range(0, ModuleN):
        sys.stdout.write(str(data[i_global][i_module][2]))
        if(i_module != 7):
            sys.stdout.write(" ")
    sys.stdout.write(" ] \n")

##### Load offsets ######
# Check if offsets were loaded from FHICL file, and adjust data accordingly 
offset_file = Path(offset_file_name)
if offset_file.is_file():
    useOffsets=True
    offsets = [[0 for i_module in range(ModuleN)] for i_global in range(globalN)]
    for i_global in range(0, 2):
        offset_file=open(offset_file_name, "r")
        offset_input = getOffsets(offset_file, FHICLServicePath+FHICLPatchName[i_global]+stationN)
        offsets[i_global]=offset_input
        print("Tracking offset "+FHICLPatchName[i_global]+stationN+" :", offsets[i_global])
        for i_module in range(0, ModuleN):
            data[i_global][i_module][1] += offsets[i_global][i_module]
print("Offset alignment used for iteration (required if iterating!):", useOffsets)   

#Print the PEDE results after offsets to screen and into a FHICL patch
fhicl_out=open("OffsetsPerModule"+str(stationN)+".fcl", "w+") 
for i_global in range(0, globalN):
    sys.stdout.write("PEDE adjusted result for "+FHICLPatchName[i_global]+stationN+": [ ")
    fhicl_out.write(FHICLServicePath+FHICLPatchName[i_global]+stationN+": [ ")
    for i_module in range(0, ModuleN):
        sys.stdout.write(str(data[i_global][i_module][1]))
        fhicl_out.write(str(data[i_global][i_module][1]*1e-3)) # um -> mm 
        if(i_module != 7):
            sys.stdout.write(" ")
            fhicl_out.write(", ")
    sys.stdout.write(" ] \n")
    fhicl_out.write(" ] \n")
fhicl_out.close()

### Load truth misalignment #####
# Check if truth file was loaded 
truth_file = Path(truth_file_name)
if truth_file.is_file():
    useTruth=True
    corrected_truth = []
    truth = [[0 for i_module in range(ModuleN)] for i_global in range(globalN)]
    for i_global in range(0, 2):
        truth_file=open(truth_file_name, "r")
        truth_input = getOffsets(truth_file, FHICLServicePath+FHICLPatchName[i_global]+stationN)
        truth[i_global]=truth_input
        print("Truth "+FHICLPatchName[i_global]+stationN+" :", truth[i_global])
        # Correct the truth offsets to have 0 angle and 0 overall offset via a a linear fit 
        x_new = np.linspace(float(min(ModuleArray)), float(max(ModuleArray)), num=1000) # generate x-points for evaluation 
        coefs = poly.polyfit(ModuleArray, truth[i_global], 1) # straight line 
        corrected_array = []
        for i_module in range(0, ModuleN):
            # corr = Y (in um) + Gradient * X + Intercept 
            correction =  int(truth[i_global][i_module] - coefs[1] * ModuleArray[i_module]  - coefs[0])
            corrected_array.append(correction)
        corrected_truth.append(corrected_array)
        print("Corrected Truth "+FHICLPatchName[i_global]+stationN+" :", corrected_truth[i_global])
print("Truth alignment used for comparison (simulation only):", useTruth)

######### Plotting #############
# Plotting "constants"
yMax = 120
yMin = -120
subplotArray = [211, 212]
#Make subplot for each result 
for i_global in range(0, globalN):
    plt.subplot(subplotArray[i_global]) 
    axes = plt.gca()
    axes.set_xlim(ModuleArray[0]-0.5, ModuleArray[-1]+0.5)
    axes.set_ylim(yMin, yMax)
    plt.title(GlobalParNames[i_global]+" misalignment", fontsize=12)
    plt.ylabel(r"Misalignment [$\mathrm{\mu m}$]")
    plt.xlabel("Module", fontsize=12)
    plt.xticks(fontsize=10, rotation=0) 
    plt.yticks(fontsize=10, rotation=0)
    plt.minorticks_on()
    axes.tick_params(axis='x',which='minor',bottom=False)
    axes.tick_params(axis='y', which='both', left=True, right=True, direction='inout')
    plt.tight_layout(pad=0.4, w_pad=0.1, h_pad=1.0)
    # data containers
    data_points =[] 
    error_points = []
    #Plot the 0th line 
    line = [[ModuleArray[0]-0.5,0.0], [ModuleArray[-1]+0.5, 0.0]]
    plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
    #Plot module lines
    for i_module in range(0, 8):
        line = [[i_module+0.5, yMin], [i_module+0.5, yMax]]
        plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
        data_points.append(data[i_global][i_module][1])
        error_points.append(data[i_global][i_module][2])
    #Plot data 
    plt.errorbar(ModuleArray, data_points, yerr=error_points,  color="purple", markersize=12, elinewidth=1, label="Reco. Mis.\n(this iteration)")
    plt.plot(ModuleArray, data_points, marker="+", color="purple")
    meanAbsReco = np.sum(np.abs(data_points))/len(data_points)
    textstr = '<|Reco|>=%s um'%(int(round(meanAbsReco)))
    plt.text(8.7, yMax*0.8, textstr, fontsize=10, color="purple")
    #Plot previous iteration
    if(useOffsets):
        plt.plot(ModuleArray, offsets[i_global],  marker="+", color="black", label="Reco. Mis.\n(previous iteration)")
    #Plot truth
    if(useTruth):
        plt.plot(ModuleArray, corrected_truth[i_global], marker=".", color="red", label="Truth Mis.")
        meanAbsTruth = np.sum(np.abs(corrected_truth[i_global]))/len(corrected_truth[i_global])
        textstr = '<|Truth|>=%s um'%(int(round(meanAbsTruth)))
        plt.text(8.7, yMax, textstr, fontsize=10, color="red")

    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 
plt.savefig("PEDE_Results_S"+stationN+".png", dpi=250)

print("\nResults plotted in PEDE_Results_S"+stationN+".png\n")
print("\nOffsets for re-tracking are written to OffsetsPerModule"+str(stationN)+".fcl")