#!/usr/bin/env python
import sys
import argparse
import glob
import os
import datetime
import time
import subprocess
import random
from ROOT import *

# Automate over: Runs, NModules: offsets, histos
NIterations = int(sys.argv[1])
NTracks = int(sys.argv[2])
# NRuns = int (sys.argv[3])

dm_m1_run1 = []
dm_m4_run1 = []
dm_m1_run2 = []
dm_m4_run2 = []

subprocess.call(["clear"])

print "Staring 2 PEDE Runs"

# Can add inner loop for NRnuns
for i in range(0, NIterations):
    print "------"
    print "Iteration ", i
    print "------ "
    print "Run 1"

    # Reset all containers for new Run
    misF = [0, 0, 0, 0]
    misC = [0, 0, 0, 0]
    Offsets = [0, 0, 0, 0]

    randSeed = random.randint(1, 1e6)
    subprocess.call(["./getRandoms.sh", str(NTracks), str(randSeed)])
    subprocess.call(["./AlignTracker", "n", str(NTracks),
                    str(Offsets[0]), str(Offsets[3])])
    subprocess.call(["./pede", "Tracker_str.txt"])

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

    dm_m1_run1.append(misF[0] + Offsets[0] - misC[0])
    dm_m4_run1.append(misF[3] + Offsets[3] - misC[3])
    print "RUN1:: misF[0]= ", misF[0], "Offsets[0]= ", Offsets[0] , "misC[0]= ", misC[0]
    print "RUN1:: misF[3]= ", misF[3], "Offsets[1]= ", Offsets[3], "misC[4]= ", misC[0]
    Offsets[0] = misF[0]  # Offset Module 1
    Offsets[3] = misF[3]  # Offset Module 4

    # Delete old steering file [just in case!]
    os.remove("Tracker_data.bin")

    # Run 2
    print "Run 2"
    randSeed = random.randint(1, 1e6)
    subprocess.call(["./getRandoms.sh", str(NTracks), str(randSeed)])
    subprocess.call(["./AlignTracker", "n", str(NTracks),
                    str(Offsets[0]), str(Offsets[3])])

    # Rename new binary file
    dataFileName = "Tracker_data" + str(int(time.time())) + ".bin"
    os.rename("Tracker_data.bin", str(dataFileName))

    # Write new steering file for Run 2 specifying new binary file
    f = open('Tracker_str_run2.txt', 'w')
    f.write("Tracker_con.txt   ! constraints text file (if applicable) \n")
    f.write("Tracker_par.txt   ! parameters (presgima) text file (if applicable)\n")
    f.write("Cfiles ! following bin files are Cfiles\n")
    f.write(str(dataFileName))
    f.write("\n")
    f.write(" method inversion 5 0.001\n")
    f.write("printrecord 2 -1\n")
    f.close()

    # how safe it is? TODO name .bin file with timestap as file name
    subprocess.call(["./pede", "Tracker_str_run2.txt"])

    module_i = 0
    with open("millepede.res") as f:
        first_line = f.readline()  # skip header
        for line in f:  # Line is a string
            number_str = line.split()
            misF[module_i] = float(number_str[1])
            module_i += 1

    with open("Tracker_pede_mis.txt") as f:
        for line in f:  #Line is a string
            number_str = line.split()
            misC[0]=float(number_str[0])
            misC[3]=float(number_str[3])

    dm_m1_run2.append(misF[0] + Offsets[0] - misC[0])
    dm_m4_run2.append(misF[3] + Offsets[3] - misC[3])
    print "RUN2:: misF[0]= ", misF[0], "Offsets[0]= ", Offsets[0] , "misC[0]= ", misC[0]
    print "RUN2:: misF[3]= ", misF[3], "Offsets[1]= ", Offsets[3], "misC[4]= ", misC[0]
 
Tf = TFile('PEDERuns.root','RECREATE')
h_m1_run1 = TH1F("h_m1_run1", "$\detla$M1 in M Run 1" , 49, -5, 18)
h_m4_run1 = TH1F("h_m4_run1", "$\detla$M4 in M Run 1" , 49, 5, -18)
h_m1_run2 = TH1F("h_m1_run2", "$\detla$M1 in M Run 2" , 49, -5, 18)
h_m4_run2 = TH1F("h_m4_run2", "$\detla$M4 in M Run 2" , 49, 5, -18)

for n in range(0, len(dm_m1_run1)): 
    h_m1_run1.Fill(float(dm_m1_run1[n])*1e4)
    h_m4_run1.Fill(float(dm_m4_run1[n])*1e4)
    h_m1_run2.Fill(float(dm_m1_run2[n])*1e4)
    h_m4_run2.Fill(float(dm_m4_run2[n])*1e4)

f = open('LaunchRuns.txt', 'w')
for i in range (0, len(dm_m1_run1)):
    f.write(str(dm_m1_run1[i]))
    f.write(" ")
    f.write(str(dm_m4_run1[i]))
    f.write(" ")
    f.write(str(dm_m1_run2[i]))
    f.write(" ")
    f.write(str(dm_m4_run2[i]))
    f.write(" ")
    f.write("\n")
f.close()  

Tf.Write()
Tf.Close()

print "ROOT file ready: root PEDERuns.root"

