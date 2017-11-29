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


#Automate over: Runs, NModules: offsets, histos 

NIterations = int (sys.argv[1]) 
NTracks = int (sys.argv[2]) 
#NRuns = int (sys.argv[3]) 

dm_m1_run1=[]
dm_m4_run1=[]
dm_m1_run2=[]
dm_m4_run2=[]

subprocess.call(["clear"])

print "Staring 2 PEDE Runs"

misF=[] #tmp misalignment holder
misC1=0.0
misC4=0.0

#Can add inner loop for NRnuns
for i in range(0, NIterations):
    print "------"
    print "Iteration " , i
    print "------ "
    print "Run 1"

    Offset1 = 0.0
    Offset2 = 0.0

    randSeed = random.randint(1,1e6)
    subprocess.call( ["./getRandoms.sh", str(NTracks), str(randSeed) ] )
    subprocess.call( ["./AlignTracker", "n", str(NTracks), str(Offset1), str(Offset2) ] ) 
    subprocess.call( ["./pede", "Tracker_str.txt" ] )
    
    with open("millepede.res") as f:
        first_line = f.readline() # skip header
        for line in f:  #Line is a string
            number_str = line.split()
            misF.append(float(number_str[1]))

    with open("Tracker_pede_mis.txt") as f:
        for line in f:  #Line is a string
            number_str = line.split()
            misC1=float(number_str[0])
            misC4=float(number_str[3])

    dm_m1_run1.append(misF[0] + Offset1 - misC1)
    dm_m4_run1.append(misF[3] + Offset2 - misC4)
    #print "RUN1:: misF[0]= ", misF[0], "Offset[0]= ", Offset1, "misC1= ", misC1
    #print "RUN1:: misF[3]= ", misF[3], "Offset[1]= ", Offset2, "misC4= ", misC4
    Offset1 = misF[0]  #Offset Module 1
    Offset2 = misF[3]  #Offset Module 4
    del misF[:] #tmp misalignment holder
    misC1=0.0
    MisC1=0.0

    #Run 2
    print "Run 2"
    randSeed = random.randint(1,1e6)
    subprocess.call( ["./getRandoms.sh", str(NTracks), str(randSeed) ] )
    subprocess.call( ["./AlignTracker", "n", str(NTracks), str(Offset1), str(Offset2) ] )
    subprocess.call( ["./pede", "Tracker_str.txt" ] )
    
    with open("millepede.res") as f:
        first_line = f.readline() # skip header
        for line in f:  #Line is a string
            number_str = line.split()
            misF.append(float(number_str[1]))

    with open("Tracker_pede_mis.txt") as f:
        for line in f:  #Line is a string
            number_str = line.split()
            misC1=float(number_str[0])
            misC4=float(number_str[3])

    dm_m1_run2.append(misF[0] + Offset1 - misC1)
    dm_m4_run2.append(misF[3] + Offset2 - misC4)
    #print "RUN2:: misF[0]= ", misF[0], "Offset[0]= ", Offset1, "misC1= ", misC1
    #print "RUN2:: misF[3]= ", misF[3], "Offset[1]= ", Offset2, "misC4= ", misC4
    del misF[:]

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

f = open('LaunchRuns.txt', 'a')
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

