#!/usr/bin/env python

#Plotter

import ROOT as r
import matplotlib.pyplot as plt #for plotting 
import numpy as np  # smart arrays 
import itertools # smart lines 

#These can be inputs from MC 
NModules=4
NLayers=4
LayerNames = ["U0", "U1", "V0", "V1"]


f = r.TFile.Open('TrackerAlignment.root')
if f:
    print("is open")
else:
    print("Not found")

totalLayer=0
plt.figure(1)
axes = plt.gca()
for i_module in range(1, NModules+1):
	line = [[i_module*4+0.5,0.92], [i_module*4+0.5, 1.02]]
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

axes.set_xlim(0, 17)
axes.set_ylim(0.92, 1.02)
plt.ylabel("<Pz/P> [error = SD]")
plt.xlabel("Layer", fontsize=10)

plt.savefig("layersPz.png")