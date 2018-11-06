#!/usr/bin/env python


####################################################################
# Splits a histogram into 3 ~equal parts
#
# Inputs: ROOT File and Hitoname, # of Slices 
# Output: range values (bin and xvalue )
#
#####################################################################

import os
import string
import time
import decimal
import ROOT
import argparse, sys
ROOT.gROOT.Macro('~/rootlogon.C')
from ROOT import *


parser = argparse.ArgumentParser(description='arguments')
parser.add_argument('-f', '--fileN', help='input ROOT file')
parser.add_argument('-n', '--sliceN', help='slice number')
args = parser.parse_args()

name = "TrackerAlignment/Tracks/Pz"
n = int(args.sliceN)

f = TFile.Open(str(args.fileN))
if f:
    print(str(args.fileN) + "is open")
else:
    print(str(args.fileN) + "not found")

t = f.Get(str(name))

totalN= int(t.GetEntries())
optimalN=totalN/n

print "Slicing", name ,"with",totalN, "entries into slices", n, "of", optimalN

hBinMin = t.FindFirstBinAbove(0, 1)
hBinMax = t.FindLastBinAbove(0, 1)

print "hBinMin= ", hBinMin, " hBinMax ", hBinMax
hBinNumber = hBinMax - hBinMin # number of non-zero bins
print " hBinNumber= ", hBinNumber

valueRange=[]

xaxis = t.GetXaxis()

valueRange.append( int( xaxis.GetBinCenter(hBinMin) ) )  # known starting point 

actionFlag = False

hsum=0
slidingN=0
for i_bin in range(hBinMin, hBinMax):
    hsum += t.GetBinContent(i_bin)
    binCentre = xaxis.GetBinCenter(i_bin)
    
    if (int(hsum/optimalN) > slidingN):
        valueRange.append( int( xaxis.GetBinCenter(i_bin) ) )
        slidingN+=1

valueRange.append( int( xaxis.GetBinCenter(hBinMax) ) )  # known ending point 


print "The equally filled Pz ranges are [MeV]:"
for i in range(0, n):
    print valueRange[i], "< Pz <", valueRange[i+1]



