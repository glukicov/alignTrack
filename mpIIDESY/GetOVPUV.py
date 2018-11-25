#!/usr/bin/env python

#Plotter

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

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-m', '--moduleN', help='mode')
parser.add_argument("-mode", "--mode")
parser.add_argument("-pvals", "--pvals", nargs='+')
parser.add_argument("-p0", "--p0", nargs='+')
parser.add_argument("-mis", "--mis", nargs='+')
parser.add_argument("-abs", "--abs")
args = parser.parse_args()


NModules=8
NLayers=4 # per modules
NTotalLayers=32
mode = str(args.mode)
LayerNames = ["U0", "U1", "V0", "V1"]
moduleNamesInitial=np.arange(1,NModules+1) #1-8
layerNamesInitial=np.arange(1, NTotalLayers+1) #1-32


#Metric



if (int(args.moduleN) != -1):
	removedModule=int(args.moduleN)
	moduleNames=np.delete(moduleNamesInitial, removedModule-1) # indexing so -1
	removedLayers= np.arange(removedModule*4-3,removedModule*4+1)
	print("removedPlanes ", removedLayers) # Layers: 0, 1... Planes: 1, 2...
	print("layerNamesInitial", layerNamesInitial)
	layerNames=np.delete(layerNamesInitial, removedLayers-1) # indexing so -1
else:
	removedModule=-1
	moduleNames=moduleNamesInitial
	layerNames=layerNamesInitial

print("Getting Plots for", len(moduleNames), "modules: ", moduleNames, "and")
print(len(layerNames), "planes: ", layerNames)

if (mode == "plot"):

	f = TFile.Open('TrackerAlignment.root')
	if f:
	    print("is open")
	else:
	    print("Not found")

	#-------layersPz----------
	i_totalLayer=0
	yMin = 0.92
	yMax = 1.02
	plt.figure(1)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		for n in range(0, NLayers):
			i_layer=layerNames[i_totalLayer]
			name = "TrackerAlignment/UV/PzoP reduced Module " + str(i_module) + " " + str(LayerNames[n])
			#print name
			t = f.Get(str(name))
			mean = t.GetMean()
			SD = t.GetRMS()
			SDError = t.GetRMSError()
			plt.errorbar(i_layer, mean, yerr=SD, color="red") 
			plt.plot(i_layer, mean, marker="_", color="red")
			i_totalLayer+=1
	
	axes.set_xlim(0, i_totalLayer+1)
	axes.set_ylim(yMin, yMax)
	plt.ylabel("<Pz/P> [error = SD]")
	plt.xlabel("Layer", fontsize=10)
	plt.savefig("layersPz.png")

	#-------layersResiduals----------
	# #
	# yMin = 0.1
	# yMax = 0.18
	# totalLayer=0
	# plt.figure(2)
	# axes = plt.gca()
	# for i in range(0, len(moduleNames)):
	# 	i_module=moduleNames[i]
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
	# for i in range(0, len(moduleNames)):
	# 	i_module=moduleNames[i]
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

	
	####### LAYERS ##############

	#-------LayerPulls----------
	#
	i_totalLayer=0
	PullsSD=[]
	PullsSDError=[]
	yMin = -1.3
	yMax = 1.3
	plt.figure(41)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		for n in range(0, NLayers):
			i_layer=layerNames[i_totalLayer]
			name = "TrackerAlignment/UV/Pulls Module " + str(i_module) + " " + str(LayerNames[n])
			#print name
			t = f.Get(str(name))
			mean = t.GetMean()
			SD = t.GetRMS()
			SDError = t.GetRMSError()
			PullsSD.append(SD)
			PullsSDError.append(SDError)
			plt.errorbar(i_layer, mean, yerr=SD, color="red") 
			plt.plot(i_layer, mean, marker="_", color="red")
			i_totalLayer+=1
	line = [[0.5,0.0], [i_totalLayer+1, 0.0]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
	axes.set_xlim(0.5, NTotalLayers+1)
	axes.set_ylim(yMin, yMax)
	plt.ylabel("Pulls Mean [error = SD]", fontsize=20)
	plt.xlabel("Layer", fontsize=20)
	plt.title("Pulls Means (SD)", fontsize=20)
	plt.savefig("Pulls_L.png")


	#-------LayerPulls_Zoom----------
	#
	i_totalLayer=0
	yMin = -0.5
	yMax = 0.8
	plt.figure(51)
	axes = plt.gca()
	means = []
	MeanErrors=[]
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		for n in range(0, NLayers):
			i_layer=layerNames[i_totalLayer]
			name = "TrackerAlignment/UV/Pulls Module " + str(i_module) + " " + str(LayerNames[n])
			t = f.Get(str(name))
			mean = t.GetMean()
			meanError=t.GetMeanError()
			means.append(mean)
			MeanErrors.append(meanError)
			plt.plot(i_layer, mean)
			plt.errorbar(i_layer, mean, yerr=meanError, color="red") 
			#axes.annotate(round_sig(mean), (i_module, mean))
			i_totalLayer+=1
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
	plt.title("Mean pull ( 'z-value' )", fontsize=18)
	plt.savefig("Pulls_L_Zoom.png")

	#-------LayerResiudals----------
	#
	i_totalLayer=0
	ResidualRMS=[]
	ResidualRMSError=[]
	yMin =-200
	yMax = 200
	plt.figure(61)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		for n in range(0, NLayers):
			i_layer=layerNames[i_totalLayer]
			name = "TrackerAlignment/UV/Residuals UV Module " + str(i_module) + " " + str(LayerNames[n])
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
			i_totalLayer+=1

	line = [[0.5,0.0], [i_totalLayer+1, 0.0]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
	axes.set_xlim(0.5, NTotalLayers+1)
	axes.set_ylim(yMin, yMax)
	plt.title("UV Residual Means (SD)", fontsize=18)
	plt.ylabel("Residual mean [um] (error = SD)", fontsize=18)
	plt.xlabel("Layer", fontsize=18)
	plt.savefig("Residuals_L.png")

	#-------LayerResiudals_Zoom----------
	#
	i_totalLayer=0
	means=[]
	MeanErrors=[]
	yMin = -65
	yMax = 80
	plt.figure(71)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		for n in range(0, NLayers):
			i_layer=layerNames[i_totalLayer]
			name = "TrackerAlignment/UV/Residuals UV Module " + str(i_module) + " " + str(LayerNames[n])
			t = f.Get(str(name))
			mean = t.GetMean()
			means.append(mean*1e3)
			meanError = t.GetMeanError()
			MeanErrors.append(meanError)
			plt.errorbar(i_layer, mean*1e3, yerr=meanError*1e3, color="red") 
			plt.plot(i_layer, mean*1e3, marker="_", color="red")
			#axes.annotate(round_sig(mean*1e3), (i_module, mean*1e3))
			i_totalLayer+=1

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
	plt.title("UV Residual Means ( 'z-value' )", fontsize=18)
	plt.ylabel("Residual Mean [um]", fontsize=18)
	plt.xlabel("Layer", fontsize=18)
	plt.savefig("Residuals_L_Zoom.png")

	#----Layer Residual SD 
	i_totalLayer=0
	yMin = 80
	yMax = 160
	means=[]
	plt.figure(81)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		for n in range(0, NLayers):
			i_layer=layerNames[i_totalLayer]
			plt.errorbar(i_layer, ResidualRMS[i_totalLayer], yerr=ResidualRMSError[i_totalLayer], color="red") 
			plt.plot(i_layer, ResidualRMS[i_totalLayer], marker="_", color="red")
			#axes.annotate(int(round_sig( ResidualRMS[i_totalLayer], 3)), (i_module,  ResidualRMS[i_totalLayer]))
			means.append(ResidualRMS[i_totalLayer])
			i_totalLayer+=1

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
	i_totalLayer=0
	i_totalLayer=0
	yMin = 0.6
	yMax = 1.4
	means=[]
	plt.figure(91)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module*NLayers+0.5,yMin], [i_module*NLayers+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		for n in range(0, NLayers):
			i_layer=layerNames[i_totalLayer]
			plt.errorbar(i_layer, PullsSD[i_totalLayer], yerr=PullsSDError[i_totalLayer], color="red") 
			plt.plot(i_layer, PullsSD[i_totalLayer], marker="_", color="red")
			means.append(PullsSD[i_totalLayer])
			i_totalLayer+=1

	avgMean = sum(means)/float(len(means))
	line = [[0.5,avgMean], [NTotalLayers+1, avgMean]]
	plt.plot( *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'black', linestyle="-")
	plt.text(32.1, avgMean, str(int(round_sig(avgMean))), fontsize=9)
	line = [[0.5,0.0], [NTotalLayers+1, 0.0]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'grey')
	axes.set_xlim(0.5, NTotalLayers+1)
	axes.set_ylim(yMin, yMax)
	plt.title("Pulls SD", fontsize=20)
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



	#-------modulePulls----------
	#
	PullsSD=[]
	PullsSDError=[]
	yMin = -1.1
	yMax = 1.1
	plt.figure(4)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		name = "TrackerAlignment/Modules/Pulls Module " + str(i_module)
		t = f.Get(str(name))
		mean = t.GetMean()
		SD = t.GetRMS()
		SDError = t.GetRMSError()
		PullsSD.append(SD)
		PullsSDError.append(SDError)
		axes.annotate(round_sig(mean, 2), (i_module, mean))
		axes.annotate( "("+str(round_sig(SD, 2))+")", (i_module-0.43, yMin+0.05))
		#print "mean= ", mean , "SD= ", SD
		plt.errorbar(i_module, mean, yerr=SD, color="red") 
		plt.plot(i_module, mean, marker="_", color="red")
	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes.set_xlim(0.5, NModules+0.5)
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
	for i_module in range(0, NModules):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		name = "TrackerAlignment/Modules/Pulls Module " + str(i_module)
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
		i_module=moduleNames[i]
		number = (means[i]-avgMean)/MeanErrors[i]
		if (number != 0):
			#number=int(number*10000)/10000
			number = number
		else:
			number = 0.0
		axes.annotate( "("+str(round_sig(number,2))+")", (i_module-0.4, -0.09))
	axes.set_xlim(0.5, NModules+0.5)
	axes.set_ylim(yMin, yMax)
	plt.ylabel("Pull Mean [error = Error on the Mean]", fontsize=20)
	plt.xlabel("Module", fontsize=20)
	plt.title("Mean pull ( 'z-value' )", fontsize=18)
	plt.savefig("Pulls_M_Zoom.png")

	#-------moduleResiudals----------
	#
	ResidualRMS=[]
	ResidualRMSError=[]
	yMin = -0.2*1e3
	yMax = 0.2*1e3
	plt.figure(6)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		name = "TrackerAlignment/Modules/Residuals UV Module " + str(i_module)
		t = f.Get(str(name))
		mean = t.GetMean()
		SD = t.GetRMS()
		SDError = t.GetRMSError()
		ResidualRMS.append(SD*1e3)
		ResidualRMSError.append(SDError*1e3)
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
	axes.set_xlim(0.5, NModules+0.5)
	axes.set_ylim(yMin, yMax)
	plt.title("UV Residual Means /um (SD)", fontsize=20)
	plt.ylabel("Residual mean /um [error = SD]", fontsize=20)
	plt.xlabel("Module", fontsize=20)
	plt.savefig("Residuals_M.png")

	#-------moduleResiudals_Zoom----------
	#
	means=[]
	MeanErrors=[]
	yMin = -15
	yMax = 15
	plt.figure(7)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		name = "TrackerAlignment/Modules/Residuals UV Module " + str(i_module)
		t = f.Get(str(name))
		mean = t.GetMean()
		means.append(mean*1e3)
		meanError = t.GetMeanError()
		MeanErrors.append(meanError)
		plt.errorbar(i_module, mean*1e3, yerr=meanError*1e3, color="red") 
		plt.plot(i_module, mean*1e3, marker="_", color="red")
		axes.annotate(round_sig(mean*1e3), (i_module, mean*1e3))

	avgMean = sum(means)/float(len(means))
	line = [[0.5,avgMean], [NModules+1.5, avgMean]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'black', linestyle="-")
	plt.text(9.1, avgMean, str(round_sig(avgMean)), fontsize=9)
	for i in range(0, len(means)):
		i_module=moduleNames[i]
		number = (means[i]-avgMean)/(MeanErrors[i]*1e3)
		if (number != 0):
			number = number
		else:
			number = 0.0
		axes.annotate( "("+str(round_sig(number,2))+")", (i_module-0.4, -0.014*1e3))

	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes.set_xlim(0.5, NModules+0.5)
	axes.set_ylim(yMin, yMax)
	plt.title("UV Residual Means ( 'z-value' )", fontsize=20)
	plt.ylabel("Residual Mean /um [error = Error on the Mean]", fontsize=18)
	plt.xlabel("Module", fontsize=20)
	plt.savefig("Residuals_M_Zoom.png")

	#----Residual SD 
	yMin = 100
	yMax = 260
	means=[]
	plt.figure(8)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		plt.errorbar(i_module, ResidualRMS[i], yerr=ResidualRMSError[i], color="red") 
		plt.plot(i_module, ResidualRMS[i], marker="_", color="red")
		axes.annotate(int(round_sig( ResidualRMS[i], 3)), (i_module,  ResidualRMS[i]))
		means.append(ResidualRMS[i])

	avgMean = sum(means)/float(len(means))
	line = [[0.5,avgMean], [NModules+1.5, avgMean]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'black', linestyle="-")
	plt.text(9.1, avgMean, str(int(round_sig(avgMean))), fontsize=9)

	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes.set_xlim(0.5, NModules+0.5)
	axes.set_ylim(yMin, yMax)
	plt.title("UV Residual SD", fontsize=20)
	plt.ylabel("Residual SD /um [error = SD Error]", fontsize=18)
	plt.xlabel("Module", fontsize=20)
	plt.savefig("ResidualsSD_M.png")

	#----Pull SD 
	yMin = 0.6
	yMax = 1.8
	means=[]
	plt.figure(9)
	axes = plt.gca()
	for i_module in range(0, NModules):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
		plt.errorbar(i_module, PullsSD[i], yerr=PullsSDError[i], color="red") 
		plt.plot(i_module, PullsSD[i], marker="_", color="red")
		axes.annotate(round_sig( PullsSD[i], 2), (i_module,  PullsSD[i]))
		means.append(PullsSD[i])

	avgMean = sum(means)/float(len(means))
	line = [[0.5,avgMean], [NModules+1.5, avgMean]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'black', linestyle="-")
	plt.text(9.1, avgMean, str(round_sig(avgMean)), fontsize=9)

	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes.set_xlim(0.5, NModules+0.5)
	axes.set_ylim(yMin, yMax)
	plt.title("Pulls SD", fontsize=20)
	plt.ylabel("Pulls SD [error = SD Error]", fontsize=18)
	plt.xlabel("Module", fontsize=20)
	plt.savefig("PullsSD_M.png")

	subprocess.call(["convert" , "+append", "Residuals_M.png" , "Pulls_M.png", "M.png"])
	subprocess.call(["convert" , "+append", "Residuals_M_Zoom.png" , "Pulls_M_Zoom.png", "M_Zoom.png"])
	subprocess.call(["convert" , "-append", "M.png" , "M_Zoom.png", "Pulls_Res_M.png"])
	subprocess.call(["convert" , "+append", "ResidualsSD_M.png" , "PullsSD_M.png", "SD_M.png"])
	subprocess.call(["convert" , "-append", "SD_M.png", "Pulls_Res_M.png", "M_SD_Pulls_Res_Fom.png"])
	subprocess.call(["trash" , "Residuals_M.png" , "Pulls_M.png", "M.png", "Residuals_M_Zoom.png" , "Pulls_M_Zoom.png", "M_Zoom.png", "M.png" , "M_Zoom.png", "Pulls_Res_M.png", "ResidualsSD_M.png" , "PullsSD_M.png", "SD_M.png"])


	#------PEDE Labels-----
	#
	myStyle  =  TStyle("MyStyle", "My Root Styles")
	cUniform = TCanvas("cLabels", "cLabels", 700, 700)
	cUniform.Divide(1,1)
	name = "TrackerAlignment/Hits/Labels"
	hUniform = f.Get(str(name))
	cUniform.cd(1)
	hUniform.Draw() 
	gStyle.SetOptStat() #over/under -flows, Rms and Means with errors, number of entries
	cUniform.Print("hLabels.png")


	#------pValFit---------
	#
	myStyle  =  TStyle("MyStyle", "My Root Styles")
	cUniform = TCanvas("cUnifrom", "cUnifrom", 700, 700)
	cUniform.Divide(1,1)
	name = "TrackerAlignment/Tracks/pValue"
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

	cutValue= 0.05	
	aboveFractionCounter = 0
	for i_bin in range(11, hBinMax):
		
		currentBinValue = hUniform.GetBinContent(i_bin)
		#print "currentBinValue= ", currentBinValue
		aboveFractionCounter=aboveFractionCounter+currentBinValue

	#print "aboveFractionCounter= ", aboveFractionCounter , "hUniform.GetEntries()= ", hUniform.GetEntries()
	aboveFraction = float(aboveFractionCounter)/(hUniform.GetEntries()) 

	print("Fraction of p-values above ", cutValue,": ", str(aboveFraction), aboveFractionCounter)

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
	
	# hUniformNorm = norm * hUniform
	hUniformNorm=hUniform # TODO normalise 


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

	#Print data to a text file
	f = open("metric.txt", "w")
	f.write(str(pValF) + " ") # hi^{2}/ndf [no error]
	f.write(str(p0) + " ")
	f.write(str(p0Error) + " ")
	f.write(str(pMean) + " ")
	f.write(str(pMeanError) + " ")
	f.write(str(largestValue) + " ")
	f.write(str(aboveFitRatio) + " ")
	

	#Save canvas as .png file
	#cUniform.Modified()
	#cUniform.Update()
	cUniform.Print("pValFit.png")
	cUniform.Print("pValFit.root")

	# for i in range(0, len(moduleNames)):
	# 	i_module=moduleNames[i]
	# 	name = "TrackerAlignment/Modules/h_MCyx_M" + str(i_module)
	# 	hist = f.Get(name)
	# 	c = TCanvas("c", "c", 700, 700)
	# 	c.Divide(1,1)
	# 	c.cd(1)
 #   		hist.Draw("COLZ")
 #   		c.Print("MC"+str(i_module)+".png")
	# subprocess.call(["convert" , "MC*.png", "MCHits_afterCuts.gif"])

	
	print("ROOT File analysed!")

pvals=[]
#-----pVal vs iteration----
#
if (mode == "pVal"):
	pvals=args.pvals
	pVals=[]
	#errors=[]
	pVals=args.pvals
	#errors=args.pvals[1::2]
	print(pVals)
	# for i in range (0, trialN):
	# 	pVals.append(args.pvals[i*2])	
	# 	errors.append(args.pvals[i*2+1])
	trialN = len(pVals)

	yMin = 0.001
	yMax = 0.01
	plt.figure(1)
	axes = plt.gca()

	for i_module in range(0, NModules):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')

	for i in range(0, trialN):
		#plt.errorbar(int(i), float(pVals[i]), yerr=float(errors[i]), color="red") 
		plt.plot(int(i), float(pVals[i]), marker="*", color="red")

	axes.set_xlim(-0.5, trialN-0.5)
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))
	axes.set_ylim(0.2, 0.6)
	line = [[-0.5,0.5], [8+0.5, 0.5]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'purple')
	#plt.ylabel(r'$\chi^{2}/ndf$ of the fit to p-value dist.', fontsize=20)
	plt.ylabel("Mean p-value", fontsize=20)
	#plt.xlabel("Module Removed", fontsize=20)
	plt.xlabel("Iteration", fontsize=20)
	plt.savefig("pFoM.png")

	print("pVal FoM produced!")

p0=[]
#-----p0 vs iteration----
#
if (mode == "p0"):
	p0=args.p0
	print(p0)
	# for i in range (0, trialN):
	# 	pVals.append(args.pvals[i*2])	
	# 	errors.append(args.pvals[i*2+1])
	trialN = len(p0)

	yMin = 0.003
	yMax = 0.01
	plt.figure(1)
	axes = plt.gca()

	for i_module in range(0, NModules):
		line = [[i_module+0.5,yMin], [i_module+0.5, yMax]]
		plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))), color = 'green')

	for i in range(0, trialN):
		#plt.errorbar(int(i), float(pVals[i]), yerr=float(errors[i]), color="red") 
		plt.plot(int(i), float(p0[i]), marker="*", color="red")

	axes.set_xlim(-0.5, trialN-0.5)
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))
	axes.set_ylim(yMin, yMax)
	plt.ylabel("p0 of the fit", fontsize=20)
	plt.xlabel("Module Removed", fontsize=20)
	# plt.xlabel("Iteration", fontsize=20)
	plt.savefig("p0FoM.png")

	print("p0 FoM produced!")

#-----Truth Misalignment----
#
if (mode == "mis"):
	misXStr=args.mis[0:8]
	misYStr=args.mis[8:16]
	print("X Mis:: ", misXStr)
	print("Y Mis:: ",misYStr)

	#absolute misalignment
	if (args.abs == "Y"):
		for i in range(0, len(misXStr)):
			misXStr[i] = abs(int(round(float(misXStr[i])*1e3)))
			misYStr[i] = abs(int(round(float(misYStr[i])*1e3)))
	else:
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
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
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
	SD = np.std(np.array(misX))
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
	for i in range(0, len(moduleNames)):
		i_module=moduleNames[i]
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
	SD = np.std(np.array(misY))
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

	print("Misalignment FoM produced!")


