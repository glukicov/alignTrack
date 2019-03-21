#!/usr/bin/python


import argparse, sys, glob, os # glob is awesome! allows regex in process calls 
import subprocess, shlex 
import re 

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-scan", "--scan") # scan study 
parser.add_argument('-p', '--path', help='path')
args = parser.parse_args()

path = str(args.path)
scan=str(args.scan)

#Define tracker constants 
stationN=2
stationName=["S12", "S18"]

files = ( "*.txt", "T*.root", "gm2tracker_ana.root", "*.fcl", "*.log", "*.bin")

#Define the scan variables
# cutScans: values to replace the default value 
# defaultValue: default in the FHICL file 
# cutName : must be present in the FHICL file as "cutName : defaultValue" e.g. "cutDCA : 0.5"
if (scan == "dca"):
    cutScans = [0.0, 0.250, 0.400, 0.600, 0.700, 0.800, 0.900, 1.000, 1.200, 1.500, 1.700, 2.000, 2.200]
       
elif (scan == "mom"):
    # cutScans = [0, 1000, 1300, 1500, 1700]
    cutScans = [0, 1000, 1300, 1500, 1700, 2000, 2300, 2500, 2800]

elif (scan == "pval"):
    # cutScans = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05]
    cutScans = [0.0, 0.001]

elif (scan == "pzp"):
    cutScans = ["PzP_98", "PzP_95"]

elif (scan == "mixed"):
    cutScans = ["4DLC_0", "4DLC_1000",  "5DLC_0",  "5DLC_1000"]    

else:
    print("Incorrect scan specified!")


#for i_total, i_cut in  enumerate(cutScans):
for i_station in range(0, stationN):
        
        #subprocess.call(["mkdir" , str(i_cut)])
        subprocess.call(["mkdir" , str(stationName[i_station])]) 
        
        #for i_station in range(0, stationN):
        for i_total, i_cut in  enumerate(cutScans):
        
            #subprocess.call(["mkdir" , str(i_cut)+"/"+str(stationName[i_station])]) 
            subprocess.call(["mkdir" , str(stationName[i_station])+"/"+str(i_cut)])
        
            # fullPath = str(path)+"/"+str(i_cut)+"/"+str(stationName[i_station])+"/"
            fullPath = str(path)+"/"+str(stationName[i_station])+"/"+str(i_cut)+"/"
            print("fullPath", fullPath)

            for i in range(0, len(files)):

                command = "gm2gpvm03:"+str(fullPath)+"/"+str(files[i])
                print("command", command)
                # subprocess.call(["scp", str(command), str(i_cut)+"/"+str(stationName[i_station])+"/"] )
                subprocess.call(["scp", str(command), str(stationName[i_station])+"/"+str(i_cut)+"/"] )