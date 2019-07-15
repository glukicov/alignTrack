################
# Run multiple pede processes in parallel 
# +Prepare data for PEDE 
# 
# Gleb Lukicov 12 July 2019
############
import glob # glob is awesome! allows regex in process calls 
import subprocess, shlex 
import time
import sys, os # print out to terminal 
from os import path
import pandas as pd # get data frame from text file 
import numpy as np # arrays 
import argparse # command line inputs sub
import re # to get offsets from file 
import numpy.polynomial.polynomial as poly # linear fit 
import matplotlib.pyplot as plt #for plotting  
import itertools # smart lines in plotting 

#Input args 
parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-d', "--inputData", type=str)
parser.add_argument('-ignoreRuns', "--ignoreRuns", nargs='+', type=int, default=[])
args = parser.parse_args()
inputData=args.inputData
ignoreRuns=args.ignoreRuns

#Define constants
modulesPerStation = 8 
globalParN = 2 # radial and vertical 
alignmentMethod = "inversion" # PEDE method 
iterationN = 10 # PEDE method 
convergenceRate = 0.01 # PEDE method

#Define global variables
all_runs=set() #no repeating runs 

def main():

    #get station number from the path
    pos = inputData.find("S1")
    stationN=(inputData[pos:pos+3]) 
    print("For station:",stationN)

    prepareData(inputData, stationN)

    runParallel(inputData, all_runs, stationN)


def prepareData(inputData, stationN):
    
    print("Preparing data in",inputData)
    
    #Check all files and runs 
    all_files = next(os.walk(inputData))[2] 
    for file in all_files:
        run=file.strip("gm2tracker_reco_").strip(".bin")[:5]
        all_runs.add(run)
    subRunN=len(all_files)
    runN=len(all_runs)
    print("Total of",subRunN,"subruns in",runN, "runs")

    #create new directory for monitoring
    subprocess.call(["mkdir", inputData+"/Monitoring"])

    #create run subdir and create Steering and Constraints file there
    for run in all_runs:
        fullPath=inputData+"/Monitoring/"+str(run)
        #check if already exists 
        if (os.path.isdir(fullPath)):
            continue
        else:
            subprocess.call(["mkdir", fullPath])
            #grab all subruns for that run 
            run_data=[]
            for file in all_files:
                pos = file.find(run)
                if (pos != -1):
                    run_data.append(inputData+"/"+file)

            #create Steering and Constraint files in the run dir 
            writeConstraintFile(inputData+"/Monitoring/"+str(run)+"/ConstraintFile.txt", stationN)
            writeSteeringFile(inputData+"/Monitoring/"+str(run)+"/SteeringFile.txt", run_data)


def runParallel(inputData, all_runs, stationN):

    print("Starting PEDE running on data in",inputData, ":")

    for run in all_runs:
        fullPath = inputData+"/Monitoring/"+run
        print(fullPath)
        os.chdir(fullPath)

        fileName = "OffsetsPerModule"+stationN+".fcl"
         #check if this run is one the blacklisted runs 
        if ( int(run) in np.array(ignoreRuns) ):
            print(run, "is one of the blacklisted runs, skipping")
            continue
        #check if the PEDE was already run here
        elif (path.exists(fileName)):
           print(fileName,"already exists in", fullPath, "skipping...")
           continue
        else:
            try:
                subprocess.Popen(["/Users/gleb/software/alignTrack/PEDE/pede", "SteeringFile.txt"])
                ########
                # Nice for testing 
                # time.sleep(10) # wait for 2 mins, if still running must be a problem 
                # subprocess.call(["pkill", "-f", "pede" ])
                # print("Killing any remaining PEDE processes..")
                #######
            except Exception:
                continue  
            try:
                subprocess.call(["python3", "/Users/gleb/software/alignTrack/mpIIDESY/RobustTrackerLaunch.py"])
            except Exception:
                continue
            os.chdir("../../")

def writeConstraintFile(constraintFilePath, stationN):
    constraintFile=open(constraintFilePath, "w+")

    stationLabel = str(stationN[1:]) # S12 -> 12 etc.

    #Write a constraint for radial track curvature
    constraintFile.write("!X^2 radial term for track curvature (around the centre of the station):\n")
    constraintFile.write("\nConstraint 0\n")
    for i_module in range(0, modulesPerStation):
        curvedFactor =  str( int( modulesPerStation + 1 - ( 2 * (i_module + 1) ) ) ** 2  )
        labelCN = stationLabel + str(i_module + 1) + str(1) # radial only
        constraintFile.write(labelCN + " " + curvedFactor + "\n")

    constraintFile.write("\n!Radial and Vertical translational DoF (overall movement):\n")
    tanslationFactor = 1 # equal weighting for all
    for i_global in range(0, globalParN):
        constraintFile.write("\nConstraint 0\n")
        for i_module in range(0, modulesPerStation):
            labelCN = stationLabel + str(i_module + 1) + str(i_global + 1)
            constraintFile.write(labelCN + " " + str(tanslationFactor) + "\n")

    constraintFile.write("\n!Radial and Vertical rotational DoF (around the centre of the station):\n")
    for i_global in range(0, globalParN):
        constraintFile.write("\nConstraint 0\n")
        for i_module in range(0, modulesPerStation):
            rotationalFactor =  str(int(modulesPerStation) + 1 - ( 2 * (i_module + 1) ))
            labelCN = stationLabel + str(i_module + 1) + str(i_global + 1)
            constraintFile.write(labelCN + " " + rotationalFactor + "\n")


def writeSteeringFile(steeringFilePath, run_data):
    steeringFile=open(steeringFilePath, "w+")

    pedeMethod = "method " + alignmentMethod +  " " +  str(iterationN) + " " + str(convergenceRate)

    steeringFile.write("* g-2 Tracker Alignment: PEDE Steering File\n")
    steeringFile.write("\n")
    steeringFile.write("ConstraintFile.txt ! constraints text file\n")
    steeringFile.write("Cfiles ! following bin files are Cfiles\n")
    for i_data in run_data:
        steeringFile.write(i_data +" ! binary data file\n")
    steeringFile.write(pedeMethod + "\n")
    steeringFile.write("\n")
    steeringFile.write("end ! optional for end-of-data\n")
   

if __name__ == '__main__':
    main()