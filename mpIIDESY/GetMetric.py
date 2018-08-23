#!/usr/bin/env python

#Plotter

from ROOT import *
import matplotlib.pyplot as plt #for plotting 
import numpy as np  # smart arrays 
import itertools # smart lines 
import argparse, sys
from math import log10, floor
from matplotlib.ticker import MaxNLocator
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
import subprocess
def round_sig(x, sig=2):
	return round(x, sig-int(floor(log10(abs(x))))-1)


NModules=8


print "Getting Metric Plots"
f = open("metric.txt", "w") # overwrite old file, if exists 

#Expected dir structure: 0, 1...8 
for i_module in range(0, NModules+1):

	print "Getting into ", i_module 

	#Read the root file for that removed module 
	f = TFile.Open(str(i_module)+'/TrackerAlignment.root')

	#------pValFit Plot---------
	#
	myStyle  =  TStyle("MyStyle", "My Root Styles")
	cUniform = TCanvas("cUnifrom", "cUnifrom", 700, 700)
	cUniform.Divide(1,1)
	name = "TrackerAlignment/Tracks/h_pValue"
	hUniform = f.Get(str(name))
	#Get parameters from the histogram
	hBinMin = hUniform.FindFirstBinAbove(0, 1)
	hBinMax = hUniform.FindLastBinAbove(0, 1)
	#print "hBinMin= ", hBinMin, " hBinMax ", hBinMax
	hBinNumber = hBinMax - hBinMin # number of non-zero bins
	#print " hBinNumber= ", hBinNumber

	#Calculate mean value of all bins
	hsum = 0.0
	for i_bin in range(hBinMin, hBinMax):
		hsum += hUniform.GetBinContent(i_bin)
	hMean = hsum/hBinNumber
	#print "hMean= ", hMean
	xaxis = hUniform.GetXaxis()
	minF = xaxis.GetBinLowEdge(hBinMin)
	maxF = xaxis.GetBinUpEdge(hBinMax)
	#print " minF= ", minF, " maxF= ", maxF

	#Set function to fit
	lineF = TF1("lineF", "pol 0", minF, maxF)
	#lineF->SetParameters(0.0, hMean);
	#Fit function
	cUniform.cd(1)

	#Normalise the histogram
	norm = 1/float(hUniform.GetEntries())
	hUniformNorm = norm * hUniform
	# hUniformNorm.SetMinimum(0.004)
	# hUniformNorm.SetMaximum(0.018)
	hUniformNorm.Draw("E1") #Set errors on all bins
	hUniformNorm.Fit("lineF", "Q") # quite fit 
	pValF = lineF.GetChisquare()/lineF.GetNDF()
	p0=lineF.GetParameter(0) #p0 of the fit 
	p0Error=lineF.GetParError(0) #p0 of the fit
	pMean=hUniformNorm.GetMean() #mean 
	pMeanError=hUniformNorm.GetMeanError() #mean 
	#Largest bin
	binmax = hUniformNorm.GetMaximumBin() 
	largestValue = hUniformNorm.GetXaxis().GetBinCenter(binmax)
	#Bin dist. vs fit
	aboveFitCounter = 0
	for i_bin in range(0, hUniformNorm.GetSize()):
		currentBinValue = hUniformNorm.GetBinContent(i_bin)
		#print "currentBinValue= ", currentBinValue, "p0= ", p0

		if (currentBinValue > p0):
			aboveFitCounter=aboveFitCounter+1
			#print "aboveFitCounter= ", aboveFitCounter

	aboveFitRatio = float(aboveFitCounter)/(hUniformNorm.GetSize()-2)  #divide by the #bin without OF and UF 

	hUniformNorm.SetTitle("Chi^{2}/ndf = " + str(round_sig(pValF,3)) + " p0= " + str(round_sig(p0,3)) + " mean=" + str(round_sig(pMean,3)) + " L.V.= " + str(largestValue) + " >P0= " + str(aboveFitRatio) + " #bins= " +  str(hUniformNorm.GetSize()-2) )
	gStyle.SetOptStat("ourRmMe") #over/under -flows, Rms and Means with errors, number of entries
	gStyle.SetOptFit(1111)  #probability, Chi2, errors, name/values of parameters
	gStyle.SetStatFormat("11.4f")  # 4 sig.fig, f=float

	#Save canvas as .png file
	#cUniform.Modified()
	#cUniform.Update()
	cUniform.Print(str(i_module)+"pValFit.png")
	cUniform.Print(str(i_module)+"pValFit.root")

	#Print data to the text file in the main dir
	f = open("metric.txt", "a")
	f.write(str(pValF) + " ") # hi^{2}/ndf [no error]
	f.write(str(p0) + " ")
	f.write(str(p0Error) + " ")
	f.write(str(pMean) + " ")
	f.write(str(pMeanError) + " ")
	f.write(str(largestValue) + " ")
	f.write(str(aboveFitRatio) + " ")
	f.write("\n")


#Now run over the metric text file
f=open("metric.txt","r")
lines=f.readlines()

data=[]
error=[]

i_type=0 #which data we are on 
#we got 7 columns with 5 types (2 have associated errors)
for i in range(0, 5):

	print "i_type= ", i_type
	for x in lines:
	    data.append(x.split(' ')[i_type])
	    
	    #if there is an error - grab it 
	    if (i_type==1 or i_type==3):
	    	error.append(x.split(' ')[i_type+1])
	    else:
	    	error.append(0.0) 
	   	
	f.close() 	

	if (i_type==1 or i_type==3):
		i_type=i_type+1 #skip over error
	i_type=i_type+1 #go to next type 
	
	
	data = np.array(data)
	error = np.array(error)
	
	#do plots 
	print data, error


	#-----vals vs iteration----
	trialN = len(data)

	yMin = 1.1*float(max(data))
	yMax = 0.9*float(min(data))
	plt.figure(i_type+1)
	axes = plt.gca()

	for i_module in range(0, NModules):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')

	for i in range(0, trialN):
		plt.errorbar(int(i), float(data[i]), yerr=float(error[i]), color="red") 
		plt.plot(int(i), float(data[i]), marker=".", color="red")

	axes.set_xlim(-0.5, trialN-0.5)
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))
	axes.set_ylim(yMax, yMin)
	
	plt.xlabel("Module Removed", fontsize=20)

	if (i_type==1):
		nameStr=r'$\chi^{2}/ndf$ of the fit to p-value dist.'
	if (i_type==3):
		nameStr="p0 of the fit"
		line = [[-0.5,0.01], [8.5, 0.01]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'purple')
		#axes.set_ylim(yMin, 0.012)
	if (i_type==5):
		nameStr="Mean Pvalue"
		line = [[-0.5,0.5], [8.5, 0.5]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'purple')
	if (i_type==6):
		nameStr="Largest Value"
	if (i_type==7):
		nameStr="Fraction above the fit"

	plt.ylabel(str(nameStr), fontsize=20)

	if (i_type==1):
		nameStr="Chi2"
	if (i_type==3):
		nameStr="p0"
	if (i_type==5):
		nameStr="mean"
	if (i_type==6):
		nameStr="LV"
	if (i_type==7):
		nameStr="AF"

	plt.savefig("fig_FoM" + str(nameStr) + ".png")


	#reset for next read 
	data=[]
	error=[]

print "Finished!"
