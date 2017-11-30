#!/usr/bin/env python
import sys
import argparse
import glob
import os
import datetime
import time
import subprocess
import random

# Automate over: Runs, NModules: offsets, histos
NIterations = int(sys.argv[1])
NRuns = int(sys.argv[2])
NTracks = int(sys.argv[3])
# NModules = int(sys.argv[4]) #Future TODO: specify how many modules,
# which one are misal.
NModules = 2

dm_m1_run = [[0 for number in xrange(NRuns)]for i_iters in xrange(NIterations)]
dm_m4_run = [[0 for number in xrange(NRuns)]for i_iters in xrange(NIterations)]

print "Staring 2 PEDE Runs for ", NIterations, "iterations with ", NTracks, "tracks"

# Can add inner loop for NRnuns
for i_iter in range(0, NIterations):

	for i_run in range(0, NRnuns):

		print "------"
		print "Iteration ", i_iter
		print "------ "
		print "Run ", i_run

		if (i_run == 0):
            # Reset all containers for the first Run
            misF = [0.0, 0.0, 0.0, 0.0]
            misC = [0.0, 0.0, 0.0, 0.0]
            Offsets = [0.0, 0.0, 0.0, 0.0]

        randSeed = random.randint(123, 1e6)
        subprocess.call(["./getRandoms.sh", str(NTracks),
                            str(randSeed)], stdout=open(os.devnull, 'wb'))
        subprocess.call(["./AlignTracker", "n", str(NTracks),
                            str(Offsets[0]), str(Offsets[3])], stdout=open(os.devnull, 'wb'))

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

	    os.rename("Tracker_data.bin", str(dataFileName))

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

	    dm_m1_run[i_run].append(misF[0] + Offsets[0] - misC[0])
	    dm_m4_run[i_run].append(misF[3] + Offsets[3] - misC[3])
	    print "RUN", i_run, ":: misF[0]= ", misF[0], "Offsets[0]= ", Offsets[0], "misC[0]= ", misC[0]
	    print "RUN", i_run":: misF[3]= ", misF[3], "Offsets[1]= ", Offsets[3], "misC[4]= ", misC[0]
	    Offsets[0] = misF[0]  # Offset Module 1
	    Offsets[3] = misF[3]  # Offset Module 4

	# end of runs/iterations 

Tf = TFile('PEDERuns.root', 'RECREATE')
h_m1_run1 = TH1F("h_m1_run1", "$\detla$M1 in M Run 1", 49, -5, 18)
h_m4_run1 = TH1F("h_m4_run1", "$\detla$M4 in M Run 1", 49, 5, -18)
h_m1_run2 = TH1F("h_m1_run2", "$\detla$M1 in M Run 2", 49, -5, 18)
h_m4_run2 = TH1F("h_m4_run2", "$\detla$M4 in M Run 2", 49, 5, -18)

for i_iter in range(0, NIterations):
	for i_run in range(0, NRuns):
		h_m1_run1.Fill(float(dm_m1_run[1][i_iter]) * 1e4)
		h_m4_run1.Fill(float(dm_m4_run[1][i_iter]) * 1e4)
		h_m1_run2.Fill(float(dm_m1_run[2][i_iter]) * 1e4)
		h_m4_run2.Fill(float(dm_m4_run[2][i_iter]) * 1e4)


f = open('LaunchRuns.txt', 'w')
for i_iter in range(0, NIterations):
	for i_run in range(0, NRuns):
		f.write(str(dm_m1_run[i_run][i_iter])
		f.write(" ")
		f.write(str(dm_m4_run[i_run][i_iter])
		f.write(" ")
	f.write("\n")
	f.close()

	Tf.Write()
	Tf.Close()

	print "ROOT file ready: root PEDERuns.root"
