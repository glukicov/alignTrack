#!/usr/bin/python


import argparse, sys, glob, os # glob is awesome! allows regex in process calls 
import subprocess, shlex 
import re 

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-scan", "--scan", default="None") # scan study 
parser.add_argument("-dir", "--directory", default="None")  
parser.add_argument('-p', '--path', help='path')
parser.add_argument('-vm', '--virtualMachine', default="Y", type=str)
args = parser.parse_args()

path = str(args.path)
scan=str(args.scan)
directory=str(args.directory)
vm = args.virtualMachine

#Define tracker constants 
# stationName=["S0"]
stationName=["S12", "S18"]
# stationName=["S0", "S12", "S18"]
stationN=len(stationName)

files = ( "*.txt", "T*.root", "gm2tracker_ana.root", "*.fcl", "*.log", "*.bin")
# files = ( "*.txt", "*.fcl", "*.log", "*.png")

machine=""
if (vm == "Y"):
    machine = "gm2gpvm04:"
    print("Coping from VM...")

if (vm == "N"):
    machine = "gm2ucl:"
    print("Coping from gm2ucl")

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
    cutScans = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.007, 0.01, 0.03, 0.05]
    defaultValue = 0.005
    cutName = "pValueCut"

elif (scan == "pzp"):
    cutScans = [0.0]
    # cutScans = [0.9, 0.93, 0.95, 0.98]
    defaultValue = 0.93
    cutName = "pzCut"

elif (scan == "minHits"):
    cutScans = [5, 7, 8, 9, 10, 11]
    defaultValue = 9
    pzCut = "minHits"

elif (scan == "t0"):
    cutScans = [30, 32, 34, 36, 38]

elif (directory != "None"):
    print("Copying over dir: ", directory)
    subprocess.call(["mkdir" , directory])
        #Pet station 
    for i_station in stationName:
        subprocess.call(["mkdir" , str(directory)+"/"+i_station]) 
        fullPath = machine + path+ "/" +i_station+ "/"
        print("Copying from", fullPath)
        #Per file 
        for i in range(0, len(files)):
            command = fullPath+ "/" +str(files[i])
            print("command", command)
            subprocess.call(["scp", str(command), directory+"/"+i_station] )

else:
    print("Incorrect scan specified!")
    sys.exit()


if (scan != "None"):

    # #copy top-level files 
    # print("Copying top-level files")
    # for i in range(0, len(files)):
    #     command = machine  + "/" + str(path)+ "/" + str(files[i])
    #     subprocess.call(["scp", str(command), "." ])
    #     print(command)

    #Per cut 
    for i_total, i_cut in  enumerate(cutScans):
            subprocess.call(["mkdir" , str(i_cut)])
           #Pet station 
            for i_station in stationName:
                subprocess.call(["mkdir" , str(i_cut)+"/"+i_station]) 
                fullPath = machine + path+ "/" +str(i_cut)+ "/" +i_station+ "/"
                print("Copying from", fullPath)
                #Per file 
                for i in range(0, len(files)):
                    command = fullPath+ "/" +str(files[i])
                    #print("command", command)
                    subprocess.call(["scp", str(command), str(i_cut)+"/"+i_station+"/"] )