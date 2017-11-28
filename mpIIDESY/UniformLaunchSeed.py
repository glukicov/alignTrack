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

nEvents = int (sys.argv[1])
UniGen = 2 * nEvents
NIterations = 1000

subprocess.call(["clear"])

print "Staring Uniform tests"

if ( os.path.isfile("Chi2Uni.txt") ):
	os.rename("Chi2Uni.txt", "Chi2Uni.txt.BK")

for i in range(0, 1000):
	randSeed = random.randint(1,1E8)
	subprocess.call( ["python", "randomIntGenerator.py", "-u", "True", "-o", "uniform_ran.txt", "-s", str(randSeed), "-n", str(UniGen) ] )
	subprocess.call( ["./GenerateFitUniform", str(nEvents)] )
	time.sleep(2)  # now we delay for the correct amount of time  

chi2_ndf=[]

with open("Chi2Uni.txt") as f:
    #next(f)
    for lineC in f:  #Line is a string
        #split the string on whitespace, return a list of numbers (as strings)
        number = lineC.split()
        chi2_ndf.append(number[0])

f = TFile('Chi2Uni.root','RECREATE')
h_chi2_ndf  = TH1F("h_chi2_ndf", "h_chi2_ndf", 129, -1, 5)

for n in range(0, len(chi2_ndf)): 
    h_chi2_ndf.Fill(int(chi2_ndf[n]))

f.Write()
f.Close()

print "ROOT file ready: root Chi2Uni.root"

