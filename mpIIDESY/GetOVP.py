#!/usr/bin/env python

#Plotter

import ROOT as r
import matplotlib.pyplot as plt #for plotting 
import numpy as np  # smart arrays 
import itertools # smart lines 
import argparse, sys
from math import log10, floor
from matplotlib.ticker import MaxNLocator

def round_sig(x, sig=2):
	return round(x, sig-int(floor(log10(abs(x))))-1)

parser = argparse.ArgumentParser(description='mode')
parser.add_argument('-m', '--moduleN', help='mode')
parser.add_argument("-mode", "--mode")
parser.add_argument("-pvals", "--pvals", nargs='+',)
args = parser.parse_args()

#These can be inputs from MC 
NModules=int(args.moduleN)
mode = str(args.mode)
NLayers=4
LayerNames = ["U0", "U1", "V0", "V1"]

print "Getting Plots for", NModules, "modules"

if (mode == "plot"):

	f = r.TFile.Open('TrackerAlignment.root')
	if f:
	    print("is open")
	else:
	    print("Not found")

	# totalLayer=0
	# xMin = 0.92
	# xMax = 1.02
	# plt.figure(1)
	# axes = plt.gca()
	# for i_module in range(1, NModules+1):
	# 	line = [[i_module*4+0.5, xMin], [i_module*4+0.5, xMax]]
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
	# axes.set_ylim(xMin, xMax)
	# plt.ylabel("<Pz/P> [error = SD]")
	# plt.xlabel("Layer", fontsize=10)

	# plt.savefig("layersPz.png")

	totalLayer=0
	xMin = 0.92
	xMax = 1.02
	plt.figure(2)
	axes = plt.gca()
	for i_module in range(1, NModules+1):
		line = [[i_module*4+0.5, xMin], [i_module*4+0.5, xMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')

		for i_layer in range(0, NLayers):
			name = "TrackerAlignment/UV/h_pzpRed_M" + str(i_module) + "_" + str(LayerNames[i_layer])
			t = f.Get(str(name))
			#print name
			mean = t.GetMean()
			SD = t.GetRMS()
			#print "mean= ", mean , "SD= ", SD
			totalLayer=totalLayer+1
			plt.errorbar(totalLayer, mean, yerr=SD, color="red") 
			plt.plot(totalLayer, mean, marker="_", color="red")

	axes.set_xlim(0, totalLayer+1)
	axes.set_ylim(xMin, xMax)
	plt.ylabel("<Pz/P_Reduced> [error = SD]")
	plt.xlabel("Layer", fontsize=10)
	plt.savefig("layersPzP_Reduced.png")

	xMin = 0.1
	xMax = 0.18
	totalLayer=0
	plt.figure(3)
	axes = plt.gca()
	for i_module in range(1, NModules+1):
		line = [[i_module*4+0.5,xMin], [i_module*4+0.5, xMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')

		for i_layer in range(0, NLayers):
			name = "TrackerAlignment/UV/h_Residuals_Module_" + str(i_module) + "_" + str(LayerNames[i_layer])
			t = f.Get(str(name))
			SD = t.GetRMS()
			SDError = t.GetRMSError()
			#print "SD= ", SD , "SDError= ", SDError
			totalLayer=totalLayer+1
			plt.errorbar(totalLayer, SD, yerr=SDError, color="red") 
			plt.plot(totalLayer, SD, marker="_", color="red")

	axes.set_xlim(0, totalLayer+1)
	axes.set_ylim(xMin, xMax)
	plt.ylabel("Residual SD [error = SD error]")
	plt.xlabel("Layer", fontsize=10)
	plt.savefig("layersResiduals.png")


	xMin = -1.0
	xMax = 1.0
	plt.figure(4)
	axes = plt.gca()
	for i_module in range(1, NModules+1):
		line = [[i_module+0.5,xMin], [i_module+0.5, xMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')
		name = "TrackerAlignment/Modules/h_Pulls_Module_" + str(i_module)
		t = f.Get(str(name))
		mean = t.GetMean()
		SD = t.GetRMS()
		axes.annotate(round_sig(mean), (i_module, mean))
		#print "mean= ", mean , "SD= ", SD
		plt.errorbar(i_module, mean, yerr=SD, color="red") 
		plt.plot(i_module, mean, marker="_", color="red")

	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes.set_xlim(0.5, NModules+1)
	axes.set_ylim(xMin, xMax)
	plt.ylabel("Pulls Mean [error = SD]")
	plt.xlabel("Module", fontsize=10)
	plt.savefig("layersPulls_M.png")


	xMin = -0.1
	xMax = 0.1
	plt.figure(5)
	axes = plt.gca()
	means = []
	for i_module in range(1, NModules+1):
		line = [[i_module+0.5,xMin], [i_module+0.5, xMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')
		name = "TrackerAlignment/Modules/h_Pulls_Module_" + str(i_module)
		t = f.Get(str(name))
		mean = t.GetMean()
		means.append(mean)
		plt.scatter(i_module, mean,  c='r',s=100)
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
	plt.text(9.1, avgMean, str(round_sig(avgMean)) + " mm", fontsize=9)

	for i in range(0, len(means)):
		number = round_sig(means[i])-round_sig(avgMean)
		if (number != 0):
			#number=int(number*10000)/10000
			number = number
		else:
			number = 0.0
		axes.annotate( "("+str(number)+")", (i+1, -0.01))
		
	axes.set_xlim(0.5, NModules+1)
	axes.set_ylim(xMin, xMax)
	plt.ylabel("Pull Mean")
	plt.xlabel("Module", fontsize=10)
	plt.title("Means pulls per Modules; mean pull (average) ")
	plt.savefig("layersPulls_Zoom_M.png")


	xMin = -0.2
	xMax = 0.2
	plt.figure(6)
	axes = plt.gca()
	for i_module in range(1, NModules+1):
		line = [[i_module+0.5,xMin], [i_module+0.5, xMax]]
		plt.plot(
		    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
		    color = 'green')
		name = "TrackerAlignment/Modules/h_Residuals_Module_" + str(i_module)
		t = f.Get(str(name))
		mean = t.GetMean()
		SD = t.GetRMS()
		#print "mean= ", mean , "SD= ", SD
		plt.errorbar(i_module, mean, yerr=SD, color="red") 
		plt.plot(i_module, mean, marker="_", color="red")

	line = [[0.5,0.0], [NModules+1, 0.0]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'grey')
	axes.set_xlim(0.5, NModules+1)
	axes.set_ylim(xMin, xMax)
	plt.ylabel("Pulls SD [error = SD]")
	plt.xlabel("Module", fontsize=10)
	plt.savefig("layersResiduals_M.png")
	print "ROOT File analysed!"

pvals=[]
#pVal vs iteration
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

	xMin = 0.4
	xMax = 0.51
	plt.figure(6)
	axes = plt.gca()
	line = [[0.5,0.5], [trialN+0.2, 0.5]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'green')

	for i in range(0, trialN):
		plt.errorbar(int(i+1), float(pVals[i]), yerr=float(errors[i]), color="red") 
		plt.plot(int(i+1), float(pVals[i]), marker="_", color="red")

	axes.set_xlim(0.8, trialN+0.2)
	axes.xaxis.set_major_locator(MaxNLocator(integer=True))
	axes.set_ylim(xMin, xMax)
	plt.ylabel("p-value mean [error = SD]")
	plt.xlabel("Iteration", fontsize=10)
	plt.savefig("pFoM.png")


print "pVal FoM produced!"
