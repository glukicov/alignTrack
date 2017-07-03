#!/usr/bin/env python

####################################################################
# FoM Plots for comparison of actual misalignment vs PEDE results 
#
# 
#
# Created: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 26 June 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
#####################################################################


#import matplotlib.pyplot as plt #for plotting 
import numpy as np  # smart arrays 
import itertools # smart lines
from ROOT import *


# Getting constants from MC
with open("Tracker_p_constants.txt") as f:
     for line in f:  #Line is a string
        number_str = line.split()
        moduleN=int(number_str[0])

# Quickly open the PEDe file and count lines only:
lineN= sum(1 for line in open('PEDE_Mis.txt'))
            

print "Parameters from Simulation and PEDE:"
print "moduleN= ",moduleN
print "PEDE Trials= ",lineN  


mis_C = [0 for i_module in xrange(moduleN)]
# Get 1 set of misalignment from simulation 

with open("Tracker_pede_mis.txt") as f:
	for line in f:  #Line is a string
		number_str = line.split()
        
        for i_module in range(0, moduleN):
        	mis_C[i_module]=float(number_str[i_module])
     

Labels = [[0 for i_module in xrange(moduleN)] for i_lines in xrange(lineN)] 
Misals = [[0 for i_module in xrange(moduleN)] for i_lines in xrange(lineN)] 
Errors = [[0 for i_module in xrange(moduleN)] for i_lines in xrange(lineN)] 
trackN = [] # track count correspond to line number 


with open("PEDE_Mis.txt") as f:
	line_i = 0
	for line in f:  #Line is a string
		number_str = line.split()
		
		for i_module in range(0, moduleN):
			Labels[line_i][i_module]=int(number_str[0+i_module*3])
			Misals[line_i][i_module]=float(number_str[1+i_module*3])
			Errors[line_i][i_module]=float(number_str[2+i_module*3])
	
		if (moduleN==4):
			trackN.append(int(number_str[12]))
		if (moduleN==6):
		    trackN.append(int(number_str[18]))
				#print number_str[18]
		line_i = line_i + 1
		


##################PLOTING##############################

f = TFile('SeedTracker.root','RECREATE')

h_dm0  = TH1F("h_dm0", "dM Module 0", 49, -6, 6)
h_dm1  = TH1F("h_dm1", "dM Module 1", 49, -6, 6)
h_dm2  = TH1F("h_dm2", "dM Module 2", 49, -6, 6)
h_dm3  = TH1F("h_dm3", "dM Module 3", 49, -6, 6)
h_dm4  = TH1F("h_dm4", "dM Module 4", 49, -6, 6)
h_dm5  = TH1F("h_dm5", "dM Module 5", 49, -6, 6)

h_er0  = TH1F("h_er0", "error Module 0", 19, -0.1, 0.1)
h_er1  = TH1F("h_er1", "error Module 1", 19, 1.2, 1.8)
h_er2  = TH1F("h_er2", "error Module 2", 19, 1.2, 1.8)
h_er3  = TH1F("h_er3", "error Module 3", 19, 1.2, 1.8)
h_er4  = TH1F("h_er4", "error Module 4", 1000, -1, 1)
h_er5  = TH1F("h_er5", "error Module 5", 1000, -1, 1)

h_dm0.GetXaxis().SetTitle("[um]");
h_dm1.GetXaxis().SetTitle("[um]");
h_dm2.GetXaxis().SetTitle("[um]");
h_dm3.GetXaxis().SetTitle("[um]");
h_dm4.GetXaxis().SetTitle("[um]");
h_dm5.GetXaxis().SetTitle("[um]");

h_er0.GetXaxis().SetTitle("[um]");
h_er1.GetXaxis().SetTitle("[um]");
h_er2.GetXaxis().SetTitle("[um]");
h_er3.GetXaxis().SetTitle("[um]");
h_er4.GetXaxis().SetTitle("[um]");
h_er5.GetXaxis().SetTitle("[um]");


for i_module in range(0, moduleN):
	for i_lines in range(0, lineN):
		#if (Misals[i_lines][i_module]>=0):
		dM=(Misals[i_lines][i_module]-mis_C[i_module])*1e4
		#if (Misals[i_lines][i_module]<0):
		#dM=(Misals[i_lines][i_module]-mis_C[i_module])*1e4
		errorM=Errors[i_lines][i_module]*1e4
		if (i_module==0):
			h_dm0.Fill(dM)
			h_er0.Fill(errorM)
		if (i_module==1):
			h_dm1.Fill(dM)
			h_er1.Fill(errorM)
		if (i_module==2):
			h_dm2.Fill(dM)
			h_er2.Fill(errorM)
		if (i_module==3):
			h_dm3.Fill(dM)
			h_er3.Fill(errorM)
		if (i_module==4):
			h_dm4.Fill(dM)
			h_er4.Fill(errorM)
		if (i_module==5):
			h_dm5.Fill(dM)
			h_er5.Fill(errorM)


gStyle.SetOptStat(111111)
c1 = TCanvas("c1","c1",700,900)
if (moduleN==4):
	c1.Divide(2,2)
	c1.cd(1) 
	h_er0.Draw()
	c1.cd(2) 
	h_er1.Draw()
	c1.cd(3) 
	h_er2.Draw()
	c1.cd(4) 
	h_er3.Draw()
if (moduleN==6):
	c1.Divide(3,3)
	c1.cd(1) 
	h_er0.Draw()
	c1.cd(2) 
	h_er1.Draw()
	c1.cd(3) 
	h_er2.Draw()
	c1.cd(4) 
	h_er3.Draw()
	c1.cd(5) 
	h_er4.Draw()
	c1.cd(6) 
	h_er5.Draw()
c1.Print("errorsSeed.png")

c2 = TCanvas("c2","c2",700,900)
if (moduleN==4):
	c2.Divide(2,2)
	c2.cd(1) 
	h_dm0.Draw()
	c2.cd(2) 
	h_dm1.Draw()
	c2.cd(3) 
	h_dm2.Draw()
	c2.cd(4) 
	h_dm3.Draw()
if (moduleN==6):
	c2.Divide(3,3)
	c2.cd(1) 
	h_dm0.Draw()
	c2.cd(2) 
	h_dm1.Draw()
	c2.cd(3) 
	h_dm2.Draw()
	c2.cd(4) 
	h_dm3.Draw()
	c2.cd(5) 
	h_dm4.Draw()
	c2.cd(6) 
	h_dm5.Draw()
c2.Print("dMSeed.png")

f.Write()
f.Close()

print "Canvases Printed to files: errorsSeed.png and dMSeed.png"





		
