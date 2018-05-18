#!/usr/bin/env python

#Plotter

import ROOT as r
import matplotlib.pyplot as plt #for plotting 
import numpy as np  # smart arrays 
import itertools # smart lines 

#These can be inputs from MC 
NModules=8
NLayers=4
LayerNames = ["U0", "U1", "V0", "V1"]


f = r.TFile.Open('TrackerAlignment.root')
if f:
    print("is open")
else:
    print("Not found")

totalLayer=0
xMin = 0.92
xMax = 1.02
plt.figure(1)
axes = plt.gca()
for i_module in range(1, NModules+1):
	line = [[i_module*4+0.5, xMin], [i_module*4+0.5, xMax]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),
	    color = 'green')

	for i_layer in range(0, NLayers):
		name = "TrackerAlignment/UV/h_pzp_M" + str(i_module) + "_" + str(LayerNames[i_layer])
		t = f.Get(str(name))
		mean = t.GetMean()
		SD = t.GetRMS()
		#print "mean= ", mean , "SD= ", SD
		totalLayer=totalLayer+1
		plt.errorbar(totalLayer, mean, yerr=SD, color="red") 
		plt.plot(totalLayer, mean, marker="_", color="red")

axes.set_xlim(0, totalLayer+1)
axes.set_ylim(xMin, xMax)
plt.ylabel("<Pz/P> [error = SD]")
plt.xlabel("Layer", fontsize=10)

plt.savefig("layersPz.png")

xMin = 0.11
xMax = 0.16
totalLayer=0
plt.figure(2)
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

print "layersPz.png and ayersResiduals.png produced"