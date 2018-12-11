#!/usr/bin/env python

#Plotter for tracker data 

from ROOT import TFile, TStyle, TCanvas, gStyle, TF1
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

parser = argparse.ArgumentParser(description='arguments')
parser.add_argument('-m', '--moduleN', help='mode')
parser.add_argument('-f', '--fileN', help='input ROOT file')
parser.add_argument("-mode", "--mode") 
args = parser.parse_args()

#Ease-of-life constants
NModules=8
NLayers=4 # per modules
NTotalLayers=NModules*NLayers
mode = str(args.mode)
UVNames = ["U0", "U1", "V0", "V1"]
moduleNamesInitial=np.arange(1,NModules+1) #1-8
layerNamesInitial=np.arange(1, NTotalLayers+1) #1-32

#Ease-of-life with removed modules
# if (int(args.moduleN) != -1):
# 	removedModule=int(args.moduleN)
# 	moduleNames=np.delete(moduleNamesInitial, removedModule-1) # indexing so -1
# 	removedLayers= np.arange(removedModule*4-3,removedModule*4+1)
# 	print "removedPlanes ", removedLayers  # Layers: 0, 1... Planes: 1, 2...
# 	print "layerNamesInitial", layerNamesInitial
# 	layerNames=np.delete(layerNamesInitial, removedLayers-1) # indexing so -1
# else:
# 	removedModule=-1
moduleNames=moduleNamesInitial
layerNames=layerNamesInitial

print("Getting Plots for", len(moduleNames), "modules: ", moduleNames, "and")
print(len(layerNames), "planes: ", layerNames)

if (mode == "plot"):

	f = TFile.Open(str(args.fileN))
	if f:
	    print(str(args.fileN) + "is open")
	else:
	    print(str(args.fileN) + "not found")

		
	####### LAYERS ##############

	#-------LayerPulls----------
	#
	PullsSD=[]
	PullsSDError=[]
	yMin = -4
	yMax = 4
	plt.figure(41)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(layerNames)):
		i_layer=layerNames[i]
		name = "TrackSummary/PerPlane/Plane"+str(i_layer)+"/Measure Pulls/UV Measure Pull Plane "+str(i_layer)
		# print name
		t = f.Get(str(name))
		mean = t.GetMean()
		SD = t.GetRMS()
		SDError = t.GetRMSError()
		PullsSD.append(SD)
		PullsSDError.append(SDError)
		plt.errorbar(i_layer, mean, yerr=SD, color="red") 
		plt.plot(i_layer, mean, marker="_", color="red")
	line = [[0.5,0.0], [NTotalLayers+1, 0.0]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
	axes.set_xlim(0.5, NTotalLayers+1)
	axes.set_ylim(yMin, yMax)
	plt.ylabel("Pulls Mean [error = SD]", fontsize=20)
	plt.xlabel("Layer", fontsize=20)
	plt.title("UV Pulls Means (SD)", fontsize=20)
	plt.savefig("Pulls_L.png")


	#-------LayerPulls_Zoom----------
	#
	yMin = -0.7
	yMax = 0.7
	plt.figure(51)
	axes = plt.gca()
	means = []
	MeanErrors=[]
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(layerNames)):
		i_layer=layerNames[i]
		name = "TrackSummary/PerPlane/Plane"+str(i_layer)+"/Measure Pulls/UV Measure Pull Plane "+str(i_layer)
		t = f.Get(str(name))
		mean = t.GetMean()
		meanError=t.GetMeanError()
		means.append(mean)
		MeanErrors.append(meanError)
		plt.plot(i_layer, mean)
		plt.errorbar(i_layer, mean, yerr=meanError, color="red") 
		#axes.annotate(round_sig(mean), (i_module, mean))
	line = [[0.5,0.0], [NTotalLayers+1, 0.0]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'grey')
	avgMean = sum(means)/float(len(means))
	line = [[0.5,avgMean], [NTotalLayers+1, avgMean]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="-")
	plt.text(32.1, avgMean, str(round_sig(avgMean)), fontsize=9)

	for i in range(0, len(means)):
		number = (means[i]-avgMean)/MeanErrors[i]
		if (number != 0):
			#number=int(number*10000)/10000
			number = number
		else:
			number = 0.0
		#axes.annotate( "("+str(round_sig(number,2))+")", (i+1-0.4, -0.09))
	axes.set_xlim(0.5, NTotalLayers+1)
	axes.set_ylim(yMin, yMax)
	plt.ylabel("Pull Mean [error = Error on the Mean]", fontsize=20)
	plt.xlabel("Layer", fontsize=20)
	plt.title("UV Mean pull ( 'z-value' )", fontsize=18)
	plt.savefig("Pulls_L_Zoom.png")

	#-------LayerResiudals----------
	#
	ResidualRMS=[]
	ResidualRMSError=[]
	yMin =-200
	yMax = 200
	plt.figure(61)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(layerNames)):
			i_layer=layerNames[i]
			name = "TrackSummary/PerPlane/Plane"+str(i_layer)+"/Measure Residuals/UVresidualsMeasPred Plane "+str(i_layer)
			t = f.Get(str(name))
			mean = t.GetMean()
			SD = t.GetRMS()
			SDError = t.GetRMSError()
			ResidualRMS.append(SD*1e3)
			ResidualRMSError.append(SDError*1e3)
			#print "mean= ", mean , "SD= ", SD
			# mm to um *1e3 
			plt.errorbar(i_layer, mean*1e3, yerr=SD*1e3, color="red") 
			plt.plot(i_layer, mean*1e3, marker="_", color="red")
			#axes.annotate(round_sig(mean*1e3, 3), (i_module, mean))
			#axes.annotate( "("+str(round_sig(SD*1e3, 3))+")", (i_module-0.43, yMin+0.05*1e3))

	line = [[0.5,0.0], [NTotalLayers+1, 0.0]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
	axes.set_xlim(0.5, NTotalLayers+1)
	axes.set_ylim(yMin, yMax)
	plt.title("UV Residual Means /um (SD)", fontsize=20)
	plt.ylabel("Residual mean /um [error = SD]", fontsize=20)
	plt.xlabel("Layer", fontsize=20)
	plt.savefig("Residuals_L.png")

	#-------LayerResiudals_Zoom----------
	#
	means=[]
	MeanErrors=[]
	yMin = -80
	yMax = 80
	plt.figure(71)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(layerNames)):
		i_layer=layerNames[i]
		name = "TrackSummary/PerPlane/Plane"+str(i_layer)+"/Measure Residuals/UVresidualsMeasPred Plane "+str(i_layer)
		t = f.Get(str(name))
		mean = t.GetMean()
		means.append(mean*1e3)
		meanError = t.GetMeanError()
		MeanErrors.append(meanError)
		plt.errorbar(i_layer, mean*1e3, yerr=meanError*1e3, color="red") 
		plt.plot(i_layer, mean*1e3, marker="_", color="red")
		#axes.annotate(round_sig(mean*1e3), (i_module, mean*1e3))

	avgMean = sum(means)/float(len(means))
	line = [[0.5,avgMean], [NTotalLayers+1, avgMean]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="-")
	plt.text(32.1, avgMean, str(round_sig(avgMean)), fontsize=9)
	for i in range(0, len(means)):
		number = (means[i]-avgMean)/(MeanErrors[i]*1e3)
		if (number != 0):
			number = number
		else:
			number = 0.0
		#axes.annotate( "("+str(round_sig(number,2))+")", (i+1-0.4, -0.014*1e3))

	line = [[0.5,0.0], [NTotalLayers+1, 0.0]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
	axes.set_xlim(0.5, NTotalLayers+1)
	axes.set_ylim(yMin, yMax)
	plt.title("UV Residual Means ( 'z-value' )", fontsize=20)
	plt.ylabel("Residual Mean /um [error = Error on the Mean]", fontsize=18)
	plt.xlabel("Layer", fontsize=20)
	plt.savefig("Residuals_L_Zoom.png")

	#----Layer Residual SD 
	yMin = 80
	yMax = 180
	means=[]
	plt.figure(81)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(layerNames)):
		i_layer=layerNames[i]
		plt.errorbar(i_layer, ResidualRMS[i], yerr=ResidualRMSError[i], color="red") 
		plt.plot(i_layer, ResidualRMS[i], marker="_", color="red")
		#axes.annotate(int(round_sig( ResidualRMS[i_layer], 3)), (i_module,  ResidualRMS[i_layer]))
		means.append(ResidualRMS[i])

	avgMean = sum(means)/float(len(means))
	line = [[0.5,avgMean], [NTotalLayers+1, avgMean]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="-")
	plt.text(32.1, avgMean, str(int(round_sig(avgMean))), fontsize=9)
	line = [[0.5,0.0], [NTotalLayers+1, 0.0]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
	axes.set_xlim(0.5, NTotalLayers+1)
	axes.set_ylim(yMin, yMax)
	plt.title("UV Residual SD", fontsize=20)
	plt.ylabel("Residual SD /um [error = SD Error]", fontsize=18)
	plt.xlabel("Layer", fontsize=20)
	plt.savefig("ResidualsSD_L.png")

	#----Layer Pull SD 
	yMin = 0.9
	yMax = 2.8
	means=[]
	plt.figure(91)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(layerNames)):
		i_layer=layerNames[i]
		plt.errorbar(i_layer, PullsSD[i], yerr=PullsSDError[i], color="red") 
		plt.plot(i_layer, PullsSD[i], marker="_", color="red")
		means.append(PullsSD[i])

	avgMean = sum(means)/float(len(means))
	line = [[0.5,avgMean], [NTotalLayers+1, avgMean]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="-")
	plt.text(32.1, avgMean, str(int(round_sig(avgMean))), fontsize=9)
	line = [[0.5,0.0], [NTotalLayers+1, 0.0]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
	axes.set_xlim(0.5, NTotalLayers+1)
	axes.set_ylim(yMin, yMax)
	plt.title("UV Pulls SD", fontsize=20)
	plt.ylabel("Pulls SD [error = SD Error]", fontsize=18)
	plt.xlabel("Layer", fontsize=20)
	plt.savefig("PullsSD_L.png")

	subprocess.call(["convert" , "+append", "Residuals_L.png" , "Pulls_L.png", "L.png"])
	subprocess.call(["convert" , "+append", "Residuals_L_Zoom.png" , "Pulls_L_Zoom.png", "L_Zoom.png"])
	subprocess.call(["convert" , "-append", "L.png" , "L_Zoom.png", "Pulls_Res_L.png"])
	subprocess.call(["convert" , "+append", "ResidualsSD_L.png" , "PullsSD_L.png", "SD_L.png"])
	subprocess.call(["convert" , "-append", "SD_L.png", "Pulls_Res_L.png", "L_SD_Pulls_Res_Fom.png"])
	subprocess.call(["trash" , "Residuals_L.png" , "Pulls_L.png", "L.png", "Residuals_L_Zoom.png" , "Pulls_L_Zoom.png", "L_Zoom.png", "L.png" , "L_Zoom.png", "Pulls_Res_L.png", "ResidualsSD_L.png" , "PullsSD_L.png", "SD_L.png"])

	########################################################


	
	print("ROOT File analysed!")