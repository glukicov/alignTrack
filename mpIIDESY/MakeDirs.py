############
## 
############
import argparse, sys, glob, os # glob is awesome! allows regex in process calls 
import subprocess, shlex 
import re 

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-name") # create dir or run 
args = parser.parse_args()

name=str(args.name)

#Define tracker constants 
stationName=["S12", "S18"]
# stationName=["S0", "S12", "S18"]
stationN=len(stationName)

subprocess.call(["mkdir" , str(name)])
for i_station in range(0, stationN):
    os.chdir(str(name)+"/")
    subprocess.call(["mkdir" , str(stationName[i_station])])
    os.chdir("../")
