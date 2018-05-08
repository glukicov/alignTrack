#!/usr/bin/env python

####################################################################
# Sanity plots for Tracker Alignment. 
# Figure-of-Merit for truth misalignment (MC) vs PEDE result.
#
# Created: 16 May 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 8 Jan 2018 by Gleb
#####################################################################

import sys 
import argparse
import glob
import os
import datetime
import time
import subprocess

#Desired number of tracks in the plot
tracksN=(1000, 2000, 3000, 5000, 6000, 7000, 10000)

factor = 66  #Based on track rejection 

tracksN = [x * factor for x in tracksN]

subprocess.call(["clear"])

print "Staring PEDE tests"

print "Parameters from Simulation:"
print "Factor based on track rejection= ",factor

if ( os.path.isfile("PEDE_Mis_art.txt") ):
	os.rename("PEDE_Mis_art.txt", "BK_PEDE_Mis_art.txt")


for i in range(0, len(tracksN)):
	subprocess.call(["gm2", "-c", "fomPzCut.fcl", "-n", str(int(tracksN[i]))])
	
	subprocess.call(["/gm2/app/users/glukicov/PEDE/V04-03-08/pede", "SteeringFile.txt" ])
	
	subprocess.call(["python" ,"ConcatenatePEDE.py"])
	

# subprocess.call(["./TrackerLaunch.py"])
print "PEDE Tests Complete. "

