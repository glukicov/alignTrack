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
     
trackN = [] # track count correspond to line number 
M1=[]
M2=[]


with open("PEDE_Mis.txt") as f:
	line_i = 0
	for line in f:  #Line is a string
		number_str = line.split()
		
		M1.append(float(number_str[4])*1e4)
		M2.append(float(number_str[7])*1e4) 

		if (moduleN==4):
			trackN.append(int(number_str[12]))
		if (moduleN==6):
		    trackN.append(int(number_str[18]))
				#print number_str[18]
		line_i = line_i + 1
		

meanBias=[]
meanBiasError=[]


with open("bias.txt") as f:
	for line in f:  #Line is a string
		number_str = line.split()
		
		meanBias.append(float(number_str[0]))
		meanBiasError.append(float(number_str[1])) 

##################PLOTING##############################

# print M1
# print M2

f = TFile('h_12.root','RECREATE')

h_12  = TH2F("h_12", "dM Module 1 vs 2", 49, -5, 5, 49, -5, 5)
h_12.GetXaxis().SetTitle("dM1 [um]");
h_12.GetYaxis().SetTitle("dM2 [um]");

for i in range(0, len(M1)):
	h_12.Fill(M1[i], M2[i])


h_bias = TH1F("h_bias", "Distribution of Mean 'Line Jitters'", 37, -0.0005, 0.0005)
h_biasError = TH1F("h_biasError", "Distribution of Errors on the Mean 'Line Jitters'", 149, 0, 0.0005)

for i in range(0, len(meanBias)):
	h_bias.Fill(meanBias[i])
	h_biasError.Fill(meanBiasError[i])


print "sigma(dM1)=", h_12.GetStdDev(1)
print "sigma(dM2)=", h_12.GetStdDev(2)
print "Cov(dM1,dM2)=", h_12.GetCovariance(1,2)
print "Correlation Coeff=Cov(dM1,dM2)/(sigma(dM1)*sigma(dM2))=",  h_12.GetCovariance(1,2)/(h_12.GetStdDev(1)*h_12.GetStdDev(2))

gStyle.SetOptStat(1111111111111111111)
c1 = TCanvas("c1","c1",700,900)
h_12.Draw("lego")
c1.Print("h_12.png")

h_12.Draw("lego")
f.Write()
f.Close()



print "Canvases Printed to files: h_12.png"





		
