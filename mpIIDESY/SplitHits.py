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
import argparse, sys
from ROOT import TH1, TFile, TCanvas, TLine, TStyle, gROOT, gStyle, TColor 
sys.path.append("/Users/gleb/")
import rootlogon as rl
rl.SetMyStyle()

parser = argparse.ArgumentParser(description='arguments')
parser.add_argument('-f', '--fileN', help='input ROOT file')
parser.add_argument('-n', '--sliceN', help='slice number')
args = parser.parse_args()

# name = "TrackerAlignment/Tracks/Pz"
# name = "TrackerAlignment/Tracks/pValue"
name = "Extrapolation/vertices/station12/h_mom"
# name = "MomentumSlices/vertices/station12/h_mom"
n = int(args.sliceN)

f = TFile.Open(str(args.fileN))
if f:
    print(str(args.fileN) + "is open")
else:
    print(str(args.fileN) + "not found")

t = f.Get(str(name))

totalN= int(t.GetEntries())
optimalN=totalN/n

print("Slicing", name ,"with",totalN, "entries into slices", n, "of", optimalN)

hBinMin = t.FindFirstBinAbove(0, 1)
hBinMax = t.FindLastBinAbove(0, 1)

print("hBinMin= ", hBinMin, " hBinMax ", hBinMax)
hBinNumber = hBinMax - hBinMin # number of non-zero bins
print(" hBinNumber= ", hBinNumber)

valueRange=[]

xaxis = t.GetXaxis()

valueRange.append( float( xaxis.GetBinCenter(hBinMin) ) )  # known starting point 

actionFlag = False

hsum=0
slidingN=0
for i_bin in range(hBinMin, hBinMax):
    hsum += t.GetBinContent(i_bin)
    binCentre = xaxis.GetBinCenter(i_bin)
    
    if (int(hsum/optimalN) > slidingN):
        valueRange.append( float( xaxis.GetBinCenter(i_bin) ) )
        slidingN+=1

valueRange.append( float( xaxis.GetBinCenter(hBinMax) ) )  # known ending point 

print("The equally filled ranges are:")

canvas = TCanvas("can", "can", 1200, 800)
lines=[]
t.Draw()
t.SetTitle("")
gStyle.SetOptStat(0)
t.GetXaxis().CenterTitle()
nbins = t.GetBinWidth(1)
t.GetYaxis().CenterTitle()
t.GetYaxis().SetTitleOffset(1.0)
t.GetYaxis().SetTitle("Tracks/ "+str(nbins)+ " MeV")
t.GetXaxis().SetTitle("P [MeV]")
for i in range(0, n):
    print(valueRange[i], "< P <", valueRange[i+1])
    l = TLine(valueRange[i], 0, valueRange[i], 1260)
    lines.append(l)
    l.SetLineColor(4)
    l.SetLineWidth(4)
    l.Draw("same")

canvas.Draw()
canvas.Print("Sliced.png")

# print("The equally filled Pz ranges are [MeV]:")
# for i in range(0, n):
#     print(valueRange[i], "< Pz <", valueRange[i+1])



