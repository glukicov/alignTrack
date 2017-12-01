#!/usr/bin/env python
import sys
#import argparse
#import glob
import os
#import datetime
import time
import subprocess
import random
from ROOT import TNtuple, TFile, TTree

# Automate over: Runs, NModules: offsets, histos
NIterations = int(sys.argv[1])
NRuns = int(sys.argv[2])
NTracks = int(sys.argv[3])

# which one are misal.
NModules = 2

dm_m1_run = [0 for number in xrange(NRuns)]
dm_m4_run = [0 for number in xrange(NRuns)]

f = open ("MC_pede_data.txt", "w")
 
print "Staring", NRuns, "PEDE Runs for ", NIterations, "iterations with ", NTracks, "tracks"

# Can add inner loop for NRuns
for i_iter in range(0, NIterations):

	for i_run in range(0, NRuns):

		print "------"
		print "Iteration ", i_iter
		print "------ "
		print "Run ", i_run

		if (i_run == 0): # Reset all containers for the first Run only
			misF = [0.0, 0.0, 0.0, 0.0]
			misC = [0.0, 0.0, 0.0, 0.0]
			Offsets = [0.0, 0.0, 0.0, 0.0]

		randSeed = random.randint(123, 1e6)
		subprocess.call(["./getRandoms.sh", str(NTracks), str(randSeed)], stdout=open(os.devnull, 'wb'))
		subprocess.call(["./AlignTracker", "n", str(NTracks), str(Offsets[0]), str(Offsets[3])], stdout=open(os.devnull, 'wb'))

	    # Rename binary and steering file
		dataFileName = "Tracker_data" + str(int(time.time())) + ".bin"
		strFileName = "Tracker_str" + str(int(time.time())) + ".txt"
		os.rename("Tracker_data.bin", str(dataFileName))

		# Write new steering file for Run 2 specifying new binary file
		f = open(strFileName, 'w')
		f.write("Tracker_con.txt   ! constraints text file (if applicable) \n")
		f.write("Tracker_par.txt   ! parameters (presgima) text file (if applicable)\n")
		f.write("Cfiles ! following bin files are Cfiles\n")
		f.write(str(dataFileName))
		f.write("\n")
		f.write(" method inversion 5 0.001\n")
		f.write("printrecord 2 -1\n")
		f.close()

		subprocess.call(["./pede", str(strFileName)],
	                    stdout=open(os.devnull, 'wb'))

		module_i = 0
		with open("millepede.res") as f:
		    first_line = f.readline()  # skip header
		    for line in f:  # Line is a string
				number_str = line.split()
				misF[module_i] = float(number_str[1])
				module_i += 1

		with open("Tracker_pede_mis.txt") as f:
		    for line in f:  # Line is a string
				number_str = line.split()
				misC[0] = float(number_str[0])
				misC[3] = float(number_str[3])

		dm1 = misF[0] + Offsets[0] - misC[0]
		dm4 = misF[3] + Offsets[3] - misC[3]
		dm_m1_run[i_run] = dm1
		dm_m4_run[i_run] = dm4
		print "RUN", i_run, ":: misF[0]= ", misF[0], "Offsets[0]= ", Offsets[0], "misC[0]= ", misC[0], "dm1=", dm1
		print "RUN", i_run, ":: misF[3]= ", misF[3], "Offsets[1]= ", Offsets[3], "misC[4]= ", misC[0], "dm4=", dm4
		Offsets[0] = misF[0]  # Offset Module 1 (for Run1+)
		Offsets[3] = misF[3]  # Offset Module 4 (for Run1+)

	    #if it is the last run - write arrays into NTuples
	    
	f = open ("MC_pede_data.txt", "a")
	f.write(str(dm_m1_run[0]) + " " + str(dm_m4_run[0]) + " " + str(dm_m1_run[1]) + " " + str(dm_m4_run[1]) + "\n")
	    	

		# end of runs/iterations 

f.close()



