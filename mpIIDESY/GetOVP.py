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

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-m', '--moduleN', help='mode')
parser.add_argument("-mode", "--mode")
parser.add_argument("-pvals", "--pvals", nargs='+')
parser.add_argument("-mis", "--mis", nargs='+')
args = parser.parse_args()

#These can be inputs from MC 
NModules=int(args.moduleN)
mode = str(args.mode)
NLayers=4
LayerNames = ["U0", "U1", "V0", "V1"]

print "Getting Plots for", NModules, "modules"

if (mode == "plot"):

	f = TFile.Open('TrackerAlignment.root')
	if f:
	    print("is open")
	else:
	    print("Not found")

	#-------layersPz----------
	#
	# totalLayer=0
	# yMin = 0.92
	# yMax = 1.02
	# plt.figure(1)
	# axes = plt.gca()
	# for i_module in range(1, NModules+1):
	# 	line = [[i_module*4+0.5, yMin], [i_module*4+0.5, yMax]]
	# 	plt.plot(
	# 	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	# 	    color = 'green')
	# 	for i_layer in range(0, NLayers):
	# 		name = "TrackerAlignment/UV/h_pzp_M" + str(i_module) + "_" + str(LayerNames[i_layer])
	# 		t = f.Get(str(name))
	# 		#print name
	# 		mean = t.GetMean()
	# 		SD = t.GetRMS()
	# 		#print "mean= ", mean , "SD= ", SD
	# 		totalLayer=totalLayer+1
	# 		plt.errorbar(totalLayer, mean, yerr=SD, color="red") 
	# 		plt.plot(totalLayer, mean, marker="_", color="red")
	# axes.set_xlim(0, totalLayer+1)
	# axes.set_ylim(yMin, yMax)
	# plt.ylabel("<Pz/P> [error = SD]")
	# plt.xlabel("Layer", fontsize=10)
	# plt.savefig("layersPz.png")

	#-------layersResiduals----------
	# #
	# yMin = 0.1
	# yMax = 0.18
	# totalLayer=0
	# plt.figure(2)
	# axes = plt.gca()
	# for i_module in range(1, NModules+1):
	# 	line = [[i_module*4+0.5,yMin], [i_module*4+0.5, yMax]]
	# 	plt.plot(
	# 	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	# 	    color = 'green')

	# 	for i_layer in range(0, NLayers):
	# 		name = "TrackerAlignment/UV/h_Residuals_Module_" + str(i_module) + "_" + str(LayerNames[i_layer])
	# 		t = f.Get(str(name))
	# 		SD = t.GetRMS()
	# 		SDError = t.GetRMSError()
	# 		#print "SD= ", SD , "SDError= ", SDError
	# 		totalLayer=totalLayer+1
	# 		plt.errorbar(totalLayer, SD, yerr=SDError, color="red") 
	# 		plt.plot(totalLayer, SD, marker="_", color="red")
	# axes.set_xlim(0, totalLayer+1)
	# axes.set_ylim(yMin, yMax)
	# plt.ylabel("Residual SD /mm [error = SD error]")
	# plt.xlabel("Layer", fontsize=10)
	# plt.savefig("layersResiduals.png")

	# #-------layersPz(red)----------
	# #
	# totalLayer=0
	# yMin = 0.92
	# yMax = 1.02
	# plt.figure(3)
	# axes = plt.gca()
	# for i_module in range(1, NModules+1):
	# 	line = [[i_module*4+0.5, yMin], [i_module*4+0.5, yMax]]
	# 	plt.plot(
	# 	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	# 	    color = 'green')
	# 	for i_layer in range(0, NLayers):
	# 		name = "TrackerAlignment/UV/h_pzpRed_M" + str(i_module) + "_" + str(LayerNames[i_layer])
	# 		t = f.Get(str(name))
	# 		#print name
	# 		mean = t.GetMean()
	# 		SD = t.GetRMS()
	# 		#print "mean= ", mean , "SD= ", SD
	# 		totalLayer=totalLayer+1
	# 		plt.errorbar(totalLayer, mean, yerr=SD, color="red") 
	# 		plt.plot(totalLayer, mean, marker="_", color="red")
	# axes.set_xlim(0, totalLayer+1)
	# axes.set_ylim(yMin, yMax)
	# plt.ylabel("<Pz/P_Reduced> [error = SD]")
	# plt.xlabel("Layer", fontsize=10)
	# plt.savefig("layersPzP_Reduced.png")


	#-------modulePulls----------
	#
	yMin = -1.1
	yMax = 1.1
	plt.figure(4)
	axes = plt.gca()
	for i_module in range(1, NModules+1):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')
		name = "TrackerAlignment/Modules/h_Pulls_Module_" + str(i_module)
		t = f.Get(str(name))
		mean = t.GetMean()
		SD = t.GetRMS()
		axes.annotate(round_sig(mean, 2), (i_module, mean))
		axes.annotate( "("+str(round_sig(SD, 2))+")", (i_module-0.43, yMin+0.05))
		#print "mean= ", mean , "SD= ", SD
		plt.errorbar(i_module, mean, yerr=SD, color="red") 
		plt.plot(i_module, mean, marker="_", color="red")
	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes.set_xlim(0.5, NModules+1)
	axes.set_ylim(yMin, yMax)
	plt.ylabel("Pulls Mean [error = SD]", fontsize=20)
	plt.xlabel("Module", fontsize=20)
	plt.title("Pulls Means (SD)", fontsize=20)
	plt.savefig("Pulls_M.png")


	#-------modulePulls_Zoom----------
	#
	yMin = -0.1
	yMax = 0.1
	plt.figure(5)
	axes = plt.gca()
	means = []
	MeanErrors=[]
	for i_module in range(1, NModules+1):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')
		name = "TrackerAlignment/Modules/h_Pulls_Module_" + str(i_module)
		t = f.Get(str(name))
		mean = t.GetMean()
		meanError=t.GetMeanError()
		means.append(mean)
		MeanErrors.append(meanError)
		plt.plot(i_module, mean)
		plt.errorbar(i_module, mean, yerr=meanError, color="red") 
		axes.annotate(round_sig(mean), (i_module, mean))
	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	avgMean = sum(means)/float(len(means))
	line = [[0.5,avgMean], [NModules+1.5, avgMean]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'black', linestyle="-")
	plt.text(9.1, avgMean, str(round_sig(avgMean)), fontsize=9)

	for i in range(0, len(means)):
		number = (means[i]-avgMean)/MeanErrors[i]
		if (number != 0):
			#number=int(number*10000)/10000
			number = number
		else:
			number = 0.0
		axes.annotate( "("+str(round_sig(number,2))+")", (i+1-0.4, -0.09))
	axes.set_xlim(0.5, NModules+1)
	axes.set_ylim(yMin, yMax)
	plt.ylabel("Pull Mean [error = Error on the Mean]", fontsize=20)
	plt.xlabel("Module", fontsize=20)
	plt.title("Mean pull ( 'z-value' )", fontsize=18)
	plt.savefig("Pulls_M_Zoom.png")

	#-------moduleResiudals----------
	#
	yMin = -0.2*1e3
	yMax = 0.2*1e3
	plt.figure(6)
	axes = plt.gca()
	for i_module in range(1, NModules+1):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')
		name = "TrackerAlignment/Modules/h_Residuals_Module_" + str(i_module)
		t = f.Get(str(name))
		mean = t.GetMean()
		SD = t.GetRMS()
		#print "mean= ", mean , "SD= ", SD
		# mm to um *1e3 
		plt.errorbar(i_module, mean*1e3, yerr=SD*1e3, color="red") 
		plt.plot(i_module, mean*1e3, marker="_", color="red")
		axes.annotate(round_sig(mean*1e3, 3), (i_module, mean))
		axes.annotate( "("+str(round_sig(SD*1e3, 3))+")", (i_module-0.43, yMin+0.05*1e3))

	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes.set_xlim(0.5, NModules+1)
	axes.set_ylim(yMin, yMax)
	plt.title("Residual Means /um (SD)", fontsize=20)
	plt.ylabel("Residual mean /um [error = SD]", fontsize=20)
	plt.xlabel("Module", fontsize=20)
	plt.savefig("Residuals_M.png")

	#-------moduleResiudals_Zoom----------
	#
	means=[]
	MeanErrors=[]
	yMin = -0.015
	yMax = 0.015
	plt.figure(7)
	axes = plt.gca()
	for i_module in range(1, NModules+1):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')
		name = "TrackerAlignment/Modules/h_Residuals_Module_" + str(i_module)
		t = f.Get(str(name))
		mean = t.GetMean()
		means.append(mean)
		meanError = t.GetMeanError()
		MeanErrors.append(meanError)
		plt.errorbar(i_module, mean, yerr=meanError, color="red") 
		plt.plot(i_module, mean, marker="_", color="red")
		axes.annotate(round_sig(mean), (i_module, mean))

	avgMean = sum(means)/float(len(means))
	line = [[0.5,avgMean], [NModules+1.5, avgMean]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'black', linestyle="-")
	plt.text(9.1, avgMean, str(round_sig(avgMean)), fontsize=9)
	for i in range(0, len(means)):
		number = (means[i]-avgMean)/MeanErrors[i]
		if (number != 0):
			number = number
		else:
			number = 0.0
		axes.annotate( "("+str(round_sig(number,2))+")", (i+1-0.4, -0.014))

	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes.set_xlim(0.5, NModules+1)
	axes.set_ylim(yMin, yMax)
	plt.title("Residual Means ( 'z-value' )", fontsize=20)
	plt.ylabel("Residual Mean /um [error = Error on the Mean]", fontsize=18)
	plt.xlabel("Module", fontsize=20)
	plt.savefig("Residuals_M_Zoom.png")

	#------pValFit---------
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
	hUniform.Draw("E1") #Set errors on all bins
	hUniform.Fit("lineF", "Q") # quite fit 
	pValF = lineF.GetChisquare()/lineF.GetNDF()
	hUniform.SetTitle("p-value fit with Chi^{2}/ndf = " + str(round_sig(pValF,3)))
	gStyle.SetOptStat("ourRmMe") #over/under -flows, Rms and Means with errors, number of entries
	gStyle.SetOptFit(1111)  #probability, Chi2, errors, name/values of parameters
	gStyle.SetStatFormat("11.4f")  # 4 sig.fig, f=float

	#Save canvas as .png file
	#cUniform.Modified()
	#cUniform.Update()
	cUniform.Print("pValFit.png")
	cUniform.Print("pValFit.root")

	subprocess.call(["convert" , "+append", "Residuals_M.png" , "Pulls_M.png", "M.png"])
	subprocess.call(["convert" , "+append", "Residuals_M_Zoom.png" , "Pulls_M_Zoom.png", "M_Zoom.png"])
	subprocess.call(["convert" , "-append", "M.png" , "M_Zoom.png", "Pulls_Res.png"])

	print "ROOT File analysed!"

pvals=[]
#-----pVal vs iteration----
#
if (mode == "pVal"):
	pvals=args.pvals
	pVals=[]
	errors=[]
	pVals=args.pvals[0::2]
	errors=args.pvals[1::2]
	print pVals, errors
	# for i in range (0, trialN):
	# 	pVals.append(args.pvals[i*2])	
	# 	errors.append(args.pvals[i*2+1])
	trialN = len(pVals)

	yMin = 0.28
	yMax = 0.51
	plt.figure(1)
	axes = plt.gca()
	line = [[0.5,0.5], [trialN+0.2, 0.5]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'green')

	for i in range(0, trialN):
		plt.errorbar(int(i), float(pVals[i]), yerr=float(errors[i]), color="red") 
		plt.plot(int(i), float(pVals[i]), marker="_", color="red")

	axes.set_xlim(-0.2, trialN+0.2)
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))
	axes.set_ylim(yMin, yMax)
	plt.ylabel("p-value mean [error = SD]", fontsize=20)
	plt.xlabel("Iteration", fontsize=20)
	plt.savefig("pFoM.png")

	print "pVal FoM produced!"

#-----Truth Misalignment----
#
if (mode == "mis"):
	misXStr=args.mis[0:8]
	misYStr=args.mis[8:16]

	#str in mm to float in um 
	for i in range(0, len(misXStr)):
		misXStr[i] = int(round(float(misXStr[i])*1e3))
		misYStr[i] = int(round(float(misYStr[i])*1e3))

	misX = np.array(misXStr)
	misY = np.array(misYStr)

	yMin = -400
	yMax = 400
	plt.subplot(211)
	axes = plt.gca()
	for i_module in range(1, NModules+1):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')
		plt.plot(i_module, float(misX[i_module-1]), marker="+", color="red")
		mod = read_png('mod.png')
		
		#Add tracker image 
		# imagebox = OffsetImage(mod, zoom=0.15)
		# xy = (i_module-0.5, float(misX[i_module-1])+150) # coordinates to position this image
		# ab = AnnotationBbox(imagebox, xy, xybox=(30., -30.),  xycoords='data', boxcoords="offset points", frameon=False)                                  
		# axes.add_artist(ab)
		axes.annotate(misX[i_module-1], (i_module-0.3, float(misX[i_module-1])+20))
	
	avgMean = sum(misX)/float(len(misX))
	SD = np.sqrt(np.mean(misX**2))
	line = [[0.5,avgMean], [NModules+1.5, avgMean]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'black', linestyle="-")
	plt.text(9.1, avgMean, str(int(round(avgMean))), fontsize=9)

	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes.set_xlim(0.5, NModules+1)
	axes.set_ylim(yMin, yMax)
	plt.title("Misalignment X (<X> = %s um $\sigma$= %s um)" %(int(round(avgMean)), int(round(SD))), fontsize=14)
	plt.ylabel("Misalignment [um]", fontsize=16)

	plt.subplot(212)
	axes2 = plt.gca()
	for i_module in range(1, NModules+1):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')
		plt.plot(i_module, float(misY[i_module-1]), marker="+", color="red")
		mod = read_png('mod.png')
		
		#Add tracker image 
		# imagebox = OffsetImage(mod, zoom=0.15)
		# xy = (i_module-0.5, float(misY[i_module-1])+120) # coordinates to position this image
		# ab = AnnotationBbox(imagebox, xy, xybox=(30., -30.),  xycoords='data', boxcoords="offset points", frameon=False)                                  
		# axes2.add_artist(ab)
		axes2.annotate(misY[i_module-1], (i_module-0.3, float(misY[i_module-1])+20))

	avgMean = sum(misY)/float(len(misY))
	SD = np.sqrt(np.mean(misY**2))
	line = [[0.5,avgMean], [NModules+1.5, avgMean]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'black', linestyle="-")
	plt.text(9.1, avgMean, str(int(round(avgMean))), fontsize=9)
	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes2.set_xlim(0.5, NModules+1)
	axes2.set_ylim(yMin, yMax)
	plt.title("Misalignment Y (<Y> = %s um $\sigma$= %s um)" %(int(round(avgMean)), int(round(SD))), fontsize=14)
	plt.ylabel("Misalignment [um]", fontsize=16)
	plt.xlabel("Module", fontsize=16)

	#plt.show()
	plt.savefig("Misalignment.png")

	print "Misalignment FoM produced!"


