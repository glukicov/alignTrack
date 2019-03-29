import subprocess, shlex 
import argparse, sys, os 

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-s', '--stationN', help='station number') 
parser.add_argument('-c', '--case', help="case study", default=None)
parser.add_argument('-i', '--iteration', help="iteration number", default=0)
args = parser.parse_args()

stationN = str(args.stationN)
case = str(args.case)
iteration = int(args.iteration)

print("Removing old plots")
subprocess.call(["trash", "S"+stationN+".png"])
subprocess.call(["trash", "FoM_Res_S"+str(stationN)+".png"])

print("Calling residual plotter")
subprocess.call(["python3", "/Users/gleb/software/alignTrack/mpIIDESY/GetOVPUV.py", "-s", str(stationN)])

print("Calling pede in the local dir")
subprocess.call(["/Users/gleb/software/alignTrack/PEDE/pede" , "SteeringFile.txt"]) 

print("Calling the alignement plotting script")
subprocess.call(["python3", "/Users/gleb/software/alignTrack/mpIIDESY/FixedTrackerLaunch.py", "-s", str(stationN), "-c", str(case), "-i", str(iteration)])

print("Opening new plots")
subprocess.call(["/Users/gleb/.iterm2/imgcat", "FoM_Res_S"+str(stationN)+".png"])
subprocess.call(["/Users/gleb/.iterm2/imgcat", "S"+stationN+".png"])