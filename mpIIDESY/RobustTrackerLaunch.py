####################################################################
# Quality control plots for Tracker Alignment: comparison of PEDE results 
#
# Works with data for MC alignment results
# Expects to have offsets written from previous iterations (if used iteratively)
#
# Created: 28 March 2019 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 28 March 2019 by Gleb
#####################################################################
import sys, os # print out to terminal 
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
    offsets = offsets.astype(float) 
    return offsets 

#Place all code in main (as we use argparse) to be able to export getOffsets function to other classes
def main():
    #Define and open command line arguments
    parser = argparse.ArgumentParser(description='mode')
    parser.add_argument('-pF', '--pede_file', help='input pede results (.res) file', default="millepede.res", type=str)  # input PEDE file  
    # parser.add_argument('-oF', '--offset_file', help='FHICL file used for tracking with offsets', default="RunGeaneSim_align.fcl", type=str)  # input PEDE file  
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
    tracker=("Tracker 1", "Tracker 2")
    ModuleArray=np.arange(1, ModuleN+1) #(1, 2,...,8) for plotting  
    GlobalParNames = ["Radial", "Vertical", r'$\phi$', r'$\psi$', r'$\theta$'] #only ever going to have 5 pars.
    units = [r" [$\mathrm{\mu m}$]", r" [$\mathrm{\mu m}$]", " [mrad]", " [mrad]", " [mrad]"]
    GlobalParLabels = [1, 2, 3, 4, 5] # their label
    GlobalParDict = dict(zip(GlobalParLabels, GlobalParNames))
    FHICLPatchName = ["strawModuleRShift", "strawModuleHShift"," strawModuleRotationPhi", "strawModuleRotationPsi", "strawModuleRotationTheta"] # no FHICL patch for angles in art yet  
    FHICLServicePath = "services.Geometry.strawtracker." #full path of the tracking FHICL patch 
    round_to = 3 # decimal places 

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
        pede_results.append((float(arrayString[1])*1000)) #mm to the nearest um / rad -> mrad
        pede_errors.append((float(arrayString[4])*1000)) #mm to the nearest um  / rad -> mrad 

    # combine into a data structure 
    results=[pars, pede_results, pede_errors]

    # structure:  data[i_global][i_module][i_results] with result: 0=label, 1=align 2= error 
    data = [[[0 for i_result in range(len(results))] for i_module in range(ModuleN)] for i_global in range(globalN)]
    for i_global in range(0, globalN):
        for i_module in range(0, ModuleN):
              for i_result in range(0, len(results)):
                data[i_global][i_module][i_result]=results[i_result][i_global::globalN][i_module] 

    #Loop through parameters to establish names ( full label = AA (station) + B (module) + C (par. label) )
    stationN=str(pars[0])[0:2]
    if (stationN=="10"):
        stationN="0" #set S0 label  
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
            sys.stdout.write(str(round(data[i_global][i_module][1], round_to)))
            if(i_module != 7):
                sys.stdout.write(" ")
        sys.stdout.write(" ] \n")
    for i_global in range(0, globalN):
        sys.stdout.write("PEDE result for "+FHICLPatchName[i_global]+stationN+" error: [ ")
        for i_module in range(0, ModuleN):
            sys.stdout.write(str(round(data[i_global][i_module][2], round_to)))
            if(i_module != 7):
                sys.stdout.write(" ")
        sys.stdout.write(" ] \n")

    ##### Load offsets ######
    # Check if offsets were loaded from FHICL file, and adjust data accordingly 
    offset_file = os.path.isfile(offset_file_name)
    if offset_file:
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
    fhicl_out=open("OffsetsPerModuleS"+str(stationN)+".fcl", "w+") 
    for i_global in range(0, globalN):
        sys.stdout.write("PEDE adjusted result for "+FHICLPatchName[i_global]+stationN+": [ ")
        fhicl_out.write(FHICLServicePath+FHICLPatchName[i_global]+stationN+": [ ")
        for i_module in range(0, ModuleN):
            sys.stdout.write(str(round(data[i_global][i_module][1], round_to)))
            fhicl_out.write(str(round(data[i_global][i_module][1]*1e-3, (round_to+1) ))) # um -> 1/10 mm for FHICL file 
            if(i_module != 7):
                sys.stdout.write(" ")
                fhicl_out.write(", ")
        sys.stdout.write(" ] \n")
        fhicl_out.write(" ] \n")
    fhicl_out.close()

    ### Load truth misalignment #####
    # if S0 there is no misalignment 
    if (stationN == "0"):
        useTruth=True
        truth = [[0 for i_module in range(ModuleN)] for i_global in range(globalN)]
        print("Using 0th truth offsets for S0")

    # Check if truth file was loaded 
    truth_file = os.path.isfile(truth_file_name)
    if truth_file:
        useTruth=True
        truth = [[0 for i_module in range(ModuleN)] for i_global in range(globalN)]
        for i_global in range(0, 2):
            truth_file=open(truth_file_name, "r")
            truth_input = getOffsets(truth_file, FHICLServicePath+FHICLPatchName[i_global]+stationN)
            truth[i_global]=truth_input
            print("Truth "+FHICLPatchName[i_global]+stationN+" :", truth[i_global])
    print("Truth alignment used for comparison (simulation only):", useTruth)


    ######### Plotting #############
    metric_ficl=open("metricS"+str(stationN)+".txt", "w+") 

    # Plotting "constants"
    f = plt.figure(figsize=(7,int(globalN*2+1)))
    #TOOD make min max from data
    if (stationN=="12"):
        i_station=0
        yMax = [120, 120, 25, 25, 25]
        yMin = [-120, -120, -25, -25, -25]
    if (stationN=="18"):
        i_station=1
        yMax = [120, 220, 25, 25, 25]
        yMin = [-120, -220, -25, -25, -25]
    #Make subplot for each result 
    for i_global in range(0, globalN):
        plt.subplot(int( str(globalN)+"1"+str(int(i_global+1)) )) 
        axes = plt.gca()
        axes.set_xlim(ModuleArray[0]-0.5, ModuleArray[-1]+0.5)
        axes.set_ylim(yMin[i_global], yMax[i_global])
        plt.title(GlobalParNames[i_global]+" misalignment in S"+stationN+" "+extra_label_name, fontsize=12)
        # plt.title(GlobalParNames[i_global]+" misalignment in "+tracker[i_station]+" "+extra_label_name, fontsize=12)
        plt.ylabel(r"Misalignment "+ units[i_global])
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
            line = [[i_module+0.5, yMin[i_global]], [i_module+0.5, yMax[i_global]]]
            plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
            data_points.append(data[i_global][i_module][1])
            error_points.append(data[i_global][i_module][2])
        #Plot data 
        plt.errorbar(ModuleArray, data_points, yerr=error_points,  color="purple", markersize=12, elinewidth=1, label="Reco. Mis.\n(this iteration)")
        plt.plot(ModuleArray, data_points, marker="+", color="purple")
        meanAbsReco = np.sum(np.abs(data_points))/len(data_points)
        textstr = "<|Reco|>="+str(int(meanAbsReco))+str(units[i_global])
        plt.text(8.7, yMax[i_global]*0.8, textstr, fontsize=8, color="purple")
        #Plot previous iteration
        if(useOffsets):
            plt.plot(ModuleArray, offsets[i_global],  marker="+", color="black", label="Reco. Mis.\n(previous iteration)")
        #Plot truth
        if(useTruth and i_global != 2): # ignoring truth for angles 
            plt.plot(ModuleArray, truth[i_global], marker=".", color="red", label="Truth Mis.")
            meanAbsTruth = np.sum(np.abs(truth[i_global]))/len(truth[i_global])
            textstr =  "<|Truth|>="+str(int(meanAbsTruth))+str(units[i_global])
            plt.text(8.7, yMax[i_global], textstr, fontsize=8, color="red")
            metric_ficl.write(str( round(meanAbsReco-meanAbsTruth, round_to) )+" ")
        if(useTruth == False):
            metric_ficl.write(str( round(meanAbsReco, round_to) )+" ")

        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}) # outside (R) of the plot 
    plt.savefig("PEDE_Results_S"+stationN+".png", dpi=250)
    metric_ficl.close()

    print("\nResults plotted in PEDE_Results_S"+stationN+".png\n")
    print("\nOffsets for re-tracking are written to OffsetsPerModuleS"+str(stationN)+".fcl")

#Place all code in main (as we use argparse) to be able to export getOffsets function to other classes
if __name__ == '__main__':
    main()