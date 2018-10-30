#!/usr/bin/python


import argparse, sys
import subprocess

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-p', '--path', help='path')
args = parser.parse_args()


path = str(args.path)

files = ( "*.txt", "T*.root", "*.fcl", "*.log", "*.bin")

for i in range(0, len(files)):

	command = "gm2gpvm01:"+str(path)+"/"+str(files[i])
	subprocess.call(["scp", str(command), "." ] )