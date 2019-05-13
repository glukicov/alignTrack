################
# Script that (1) creates multiple dirs with many FHICL files
# that has a changed value for scanning studies in "make" mode 
# (2) runs multiple processes in bkg in "run" mode 
#
# subprocess.call - waits to finish before continues
# subprocess.Popen - does not wait to finish 
# & = bkg process 
# 
# Gleb Lukicov 4 March 2019
############
import glob # glob is awesome! allows regex in process calls 
import subprocess, shlex 

import sys, os # print out to terminal 
import pandas as pd # get data frame from text file 
import numpy as np # arrays 
import argparse # command line inputs sub
import re # to get offsets from file 
import numpy.polynomial.polynomial as poly # linear fit 
import matplotlib.pyplot as plt #for plotting  
import itertools # smart lines in plotting 


parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-mode", "--mode") # create dir or run 
parser.add_argument("-scan", "--scan") # scan study 
parser.add_argument("-dir", "--dir") # scan study 
parser.add_argument('-inputData', "--inputData", type=str, default="/home/gm2/glukicov/defDig_defTrack/*.root")
# /home/gm2/glukicov/defDig_defTrack/*.root
# /pnfs/GM2/scratch/users/glukicov/defDig_defTrack_Iter1/2019-04-19-13-45-25/data/*.root 
# /home/gm2/glukicov/noCutsFullTracks/*.root  - no cuts for dca and min hits 
args = parser.parse_args()

#DS0: /home/gm2/glukicov/defDig_defTrack

mode=str(args.mode)
scan=str(args.scan)
inputData=str(args.inputData)
dirScan=str(args.dir)

print("inputData:", inputData)

#Define tracker constants 
stationName=["S12", "S18"]
# stationName=["S0", "S12", "S18"]
stationN=len(stationName)

#Define input FHICL file(s) 
# Having separation between in/out FHICLs allows to rename these 
# inFile = []
# for i_station in range(0, stationN):
#     inFile_station = str(os.environ['MRB_SOURCE'])+"/gm2tracker/align/fcl/trackerAlignment_"+stationName[i_station]+".fcl" 
#     inFile.append(inFile_station)

# #Define modified FHICL file name(s) 
# outFile=[]
# for i_station in range(0, stationN):
#     outFile_station = "trackerAlignment_"+stationName[i_station]+".fcl"  
#     outFile.append(outFile_station)

# #Plot FHICL file 
# plotFile="RunAlignmentPlots.fcl"
# plotFilePath = str(os.environ['MRB_SOURCE'])+"/gm2tracker/align/fcl/"+str(plotFile)


#Define the scan variables
# cutScans: values to replace the default value 
# defaultValue: default in the FHICL file 
# cutName : must be present in the FHICL file as "cutName : defaultValue" e.g. "cutDCA : 0.5"
if (scan == "dca"):
    cutScans = [0.3, 0.4, 0.5, 0.6, 0.7]
    defaultValue = 0.5
    cutName = "cutDCA"
   
elif (scan == "mom"):
    cutScans = [0, 800, 1000, 1300, 1500, 1700, 2000]
    defaultValue = 0.0
    cutName = "PCut"

elif (scan == "time"):
    cutScans = [0, 10000, 20000, 30000, 35000]
    defaultValue = 0.0
    cutName = "timeCut"

elif (scan == "pval"):
    cutScans = [0.001, 0.002, 0.004]
     # cutScans = [0.0, 0.003, 0.005, 0.007, 0.01, 0.03, 0.05]
    defaultValue = 0.005
    cutName = "pValueCut"

elif (scan == "pzp"):
    cutScans = [0.9, 0.93, 0.95, 0.98, 0.99]
    defaultValue = 0.93
    cutName = "pzCut"

elif (scan == "minHits"):
    # cutScans = [5, 7, 8, 9, 10, 11, 12, 13, 15]
    cutScans = [8, 10]
    defaultValue = 9
    cutName = "minHits"

elif (scan == "t0"):
    # cutScans = [30, 32, 34, 35, 37, 38]
    cutScans = [35, 36, 37]


elif( dirScan != None):
    print("Single Mode ",dirScan)
    cutScans=[dirScan]

else:
    print("Incorrect scan specified!")


if (mode == "make"):
    print("Creating directories and FHICL files for scans...")

    for i_total, i_cut in  enumerate(cutScans):
        # create a top-level cut directory 
        subprocess.call(["mkdir" , str(i_cut)]) 
        # create station dir and copy over the FHICL file 
        for i_station in range(0, stationN):
            print(i_total, "cutName: ", cutName, i_cut, stationName[i_station]) # status printout 

            subprocess.call(["mkdir" , str(i_cut)+"/"+stationName[i_station]])
            subprocess.call(["cp" , str(inFile[i_station]), str(i_cut)+"/"+str(stationName[i_station])+"/"+str(outFile[i_station])]) 
            subprocess.call(["cp" , str(plotFilePath), str(i_cut)+"/"+str(stationName[i_station])+"/"+str(plotFile)])

            #Read the original file and change the scan value by overwriting the file  
            FHICL_nominal = open(str(i_cut)+"/"+str(stationName[i_station])+"/"+str(outFile[i_station])).read()
            FHICL_nominal = FHICL_nominal.replace(cutName+" : "+str(defaultValue), cutName+" : "+str(i_cut))
            FHICL_new = open(str(i_cut)+"/"+str(stationName[i_station])+"/"+str(outFile[i_station]), 'w')
            FHICL_new.write(FHICL_nominal)
            FHICL_new.close()

elif (mode == "run"):
    print("Starting ",str(int(len(cutScans)*stationN))," gm2 ana processes in parallel")

    for i_total, i_cut in  enumerate(cutScans):
        for i_station in range(0, stationN):
        
            os.chdir(str(i_cut)+"/"+str(stationName[i_station])+"/")
            subprocess.call("pwd") # status printout 
            subprocess.Popen(["gm2", "-c", str(outFile[i_station]), "-s"] + glob.glob(inputData))  # use glob to allow for regex in input data definition  
            
            os.chdir("../../") # go back to top dir 

elif (mode == "align"):
    print("Starting ",str(int(len(cutScans)*stationN))," aligning python processes in parallel")

    for i_total, i_cut in  enumerate(cutScans):
        for i_station in range(0, stationN):
        
            os.chdir(str(i_cut)+"/"+str(stationName[i_station])+"/")
            subprocess.call("pwd") # status printout 
            # subprocess.Popen(["python", str(os.environ['MP2'])+"/align.py"])  # use glob to allow for regex in input data definition  
            
            #subprocess.call(["/gm2/app/users/glukicov/PEDE/V04-03-10/pede", "SteeringFile.txt"]) 
            subprocess.call(["python",  str(os.environ['MP2'])+"/RobustTrackerLaunch.py", "--extra_label", str(scan)+": "+str(i_cut)])

            os.chdir("../../") # go back to top dir 

elif (mode == "align-plots"):
    print("Starting ",str(int(len(cutScans)*stationN))," aligning python processes in parallel")

    for i_total, i_cut in  enumerate(cutScans):
        for i_station in range(0, stationN):
        
            os.chdir(str(i_cut)+"/"+str(stationName[i_station])+"/")
            subprocess.call("pwd") # status printout 
            subprocess.call(["python", str(os.environ['MRB_SOURCE'])+"/gm2tracker/align/macros/GetRes.py"]) 

            os.chdir("../../") # go back to top dir 


elif (mode == "run2"):
    print("Starting ",str(int(len(cutScans)*stationN))," gm2 ana processes in parallel")

    for i_station in range(0, stationN):
        for i_total, i_cut in enumerate(cutScans):
        
            os.chdir(str(stationName[i_station])+"/"+str(i_cut)+"/")
            subprocess.call("pwd") # status printout 
            subprocess.Popen(["gm2", "-c", str(outFile[i_station]), "-s"] + glob.glob(inputData))  # use glob to allow for regex in input data definition  
 
            os.chdir("../../") # go back to top dir

elif (mode == "run3"):
    print("Starting ",str(int(stationN))," gm2 ana processes in parallel")

    for i_station in range(0, stationN):

            os.chdir(str(dirScan)+str(stationName[i_station]))
            subprocess.call("pwd") # status printout 
            subprocess.Popen(["gm2", "-c", str(outFile[i_station]), "-s"] + glob.glob(inputData))  # use glob to allow for regex in input data definition  

            os.chdir("../../") # go back to top dir

elif (mode == "make-plot2"):
    print("Copying over ana FHICL files in the scan directories")

    for i_station in range(0, stationN):
        for i_total, i_cut in  enumerate(cutScans):
            print(i_total, "cutName: ", cutName, i_cut, stationName[i_station]) # status printout 
            subprocess.call(["cp" , str(plotFilePath), str(stationName[i_station])+"/"+str(i_cut)+"/"+str(plotFile)]) 


elif (mode == "make-plot"):
    print("Copying over ana FHICL files in the scan directories")

    for i_total, i_cut in  enumerate(cutScans):
        for i_station in range(0, stationN):
            print(i_total, "cutName: ", cutName, i_cut, stationName[i_station]) # status printout 
            subprocess.call(["cp" , str(plotFilePath), str(i_cut)+"/"+str(stationName[i_station])+"/"+str(plotFile)]) 

elif (mode == "run-plot"):
    print("Starting ",str(int(len(cutScans)*stationN))," gm2 ana processes in parallel") 

    for i_total, i_cut in  enumerate(cutScans):
        for i_station in range(0, stationN):
        
            os.chdir(str(i_cut)+"/"+str(stationName[i_station])+"/")
            subprocess.call("pwd") # status printout 
            subprocess.Popen(["gm2", "-c", str(plotFile)])  # assume ana FHICL has a default input TFile   
 
            os.chdir("../../") # go back to top dir 

elif (mode == "run-plot2"):
    print("Starting ",str(int(len(cutScans)*stationN))," gm2 ana processes in parallel") 

    for i_station in range(0, stationN):
        for i_total, i_cut in  enumerate(cutScans):
        
            os.chdir(str(stationName[i_station])+"/"+str(i_cut)+"/")
            subprocess.call("pwd") # status printout 
            subprocess.Popen(["gm2", "-c", str(plotFile)])  # assume ana FHICL has a default input TFile   
 
            os.chdir("../../") # go back to top dir 

else:
    print("Incorrect mode!")