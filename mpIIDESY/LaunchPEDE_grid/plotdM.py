#!/usr/bin/env python

#Plotter
import glob
from ROOT import *

#These can be inputs from MC 
NModules=2
NRuns=2
dataSize = NModules*NRuns

#How to create a TH1F array in python?
histArray = []
runLabel = ("1", "1", "2", "2")

NBins = 49
lowBinEdge = 0.0
highBinEdge = 20.0

Tf = TFile('PEDERuns.root', 'RECREATE')
c = TCanvas("c", "Results of PEDE runs", 200, 10, 1100, 1100)
c.Divide(2,2)
h_dm1_run1= TH1F("h_dm1_run1", "M1 in Run 1", NBins, -2, 18) 
h_dm4_run1= TH1F("h_dm4_run1", "M4 in Run 1", NBins, -18, 2) 
h_dm1_run2= TH1F("h_dm1_run2", "M1 in Run 2", NBins, -2, 18) 
h_dm4_run2= TH1F("h_dm4_run2", "M4 in Run 2", NBins, -18, 2) 

files = glob.glob("data/*.txt")
for fle in files:
	#Read file and fill histos 
	with open(fle) as f:
	    #next(f)
	    for lineC in f:  #Line is a string
	        #split the string on whitespace, return a list of numbers (as strings)
	        number = lineC.split()
	       	h_dm1_run1.Fill(float(number[0])* 1e4) #cm to um 
	       	h_dm4_run1.Fill(float(number[1])* 1e4) #cm to um 
	       	h_dm1_run2.Fill(float(number[2])* 1e4) #cm to um 
	       	h_dm4_run2.Fill(float(number[3])* 1e4) #cm to um 

myStyle  = TStyle("MyStyle", "My Root Styles")

#Canvas
myStyle.SetCanvasBorderMode(0)  # Transparent
myStyle.SetCanvasColor(0) # Transparent 

#Paper, Pad, Palette, Frame
myStyle.SetPadBorderMode(0) # Transparent 
myStyle.SetPadColor(0) # Transparent 
myStyle.SetPalette(1) # Default 
myStyle.SetFrameBorderMode(1) # Border 

 # Axis 
myStyle.SetLabelSize(0.035, "xyz") # size of axis values
MyStyle.SetTitleSize(0.035, "xyz")
myStyle.SetPadTickX(1)
myStyle.SetPadTickY(1)

# Title 
myStyle.SetTitleColor(1) # Black 
myStyle.SetTitleStyle(0) # Transparent 
myStyle.SetTitleBorderSize(0) # Transparent
myStyle.SetTitleY(0.97) # Set y-position (fraction of pad size)
myStyle.SetTitleX(0.4) # Set x-position (fraction of pad size)

# #Stat box dimensions, position and style 
myStyle.SetStatY(0.89) # Set y-position (fraction of pad size)
myStyle.SetStatX(0.89) # Set x-position (fraction of pad size)
myStyle.SetStatW(0.42) # Set width of stat-box (fraction of pad size)
myStyle.SetStatH(0.09) # Set height of stat-box (fraction of pad size)
myStyle.SetStatStyle(0) # Transparent 
myStyle.SetStatColor(0)  # Transparent
myStyle.SetStatBorderSize(1) # Transparent 


# Histo Filling (visual)
myStyle.SetHistFillColor(kRed)
myStyle.SetHistFillStyle(3012)      

# Stats display options 
myStyle.SetOptStat("ourRmMe") #over/under -flows, Rms and Means with errors, number of entries
myStyle.SetOptFit(1111)  #probability, Chi2, errors, name/values of parameters
myStyle.SetStatFormat("11.5f")  # 4 sig.fig, f=float, no idea what is 6 for? 

gROOT.ForceStyle()

gROOT.SetStyle("MyStyle")

c.cd(1)
h_dm1_run1.GetXaxis().SetTitle("[um]")
h_dm1_run1.Draw()
c.cd(2)
h_dm4_run1.GetXaxis().SetTitle("[um]")
h_dm4_run1.Draw()
c.cd(3)
h_dm1_run2.GetXaxis().SetTitle("[um]")
h_dm1_run2.Draw()
c.cd(4)
h_dm4_run2.GetXaxis().SetTitle("[um]")
h_dm4_run2.Draw()

c.Modified()
c.Update()
c.Print("PEDERuns.png")

Tf.Write()
Tf.Close()

print "ROOT file ready: root PEDERuns.root"