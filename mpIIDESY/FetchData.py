#!/usr/bin/python


import argparse, sys
import subprocess

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-p', '--path', type=str)
parser.add_argument('-vm', '--virtualMachine', default="Y", type=str)
args = parser.parse_args()


path = args.path
vm = args.virtualMachine

files = ( "*.txt", "T*.root", "gm2tracker_ana.root", "*.fcl", "*.log", "*.bin")

command=""
for i in range(0, len(files)):
    
    if (vm == "Y"):
        command = "gm2gpvm04:"+str(path)+"/"+str(files[i])
        print("Coping from VM...")
    
    if (vm == "N"):
        command = "gm2ucl:"+str(path)+"/"+str(files[i])
        print("Coping from gm2ucl")
    
    subprocess.call(["scp", str(command), "." ] )