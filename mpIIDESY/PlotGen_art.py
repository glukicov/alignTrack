#!/usr/bin/python

####################################################################
# Sanity plots for AlginTracker
#
# 
#
# Created: 16 May 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 4 March 2018 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
#####################################################################

#import numpy as np 
import matplotlib.pyplot as plt #for plotting 
import itertools
import csv
import pprint


moduleN=4
viewN=2
layerN=2
strawN=32

# Getting constants from MC
with open("trackN.txt") as f:
	for line in f:  #Line is a string
		number_str = line.split()
		trackN=int(number_str[0])


toalLayerN=layerN*moduleN*viewN

print "Parameters from Simulation:"
print "moduleN= ",moduleN
print "viewN= ",viewN
print "layerN= ",layerN
print "strawN= ",strawN
print "toalLayerN= ",toalLayerN
print "trackN= ",trackN


# X 4D arrays for Mis and Ideal Geom. 
# MisX = [[[[0 for i_straw in xrange(strawN)] for i_layer in xrange(layerN) ] for i_view in xrange(viewN)] for i_module in xrange(moduleN)]
# MisZ = [[[[0 for i_straw in xrange(strawN)] for i_layer in xrange(layerN) ] for i_view in xrange(viewN)] for i_module in xrange(moduleN)]
fit=[[0 for number in xrange(4)] for i_track in xrange(trackN)]


# layerMx=[] #temp storage
# layerMz=[] #temp storage
# with open("Tracker_geom_art.txt") as f:
# 	for line in f:  #Line is a string
# 		layerMx.append(line.split())
# 		nextLine = next(f)
# 		number_str=nextLine.split()
# 		layerMz.append(nextLine.split()) #z

# #Now for straws in X: 
# i_totalLayers=0
# for i_module in range(0, moduleN):
# 	for i_view in range(0, viewN):
# 		for i_layer in range(0, layerN):
# 			for i_straw in range(0, strawN):
# 				dXI= float(layerIx[i_totalLayers][i_straw])
# 				dXM= float(layerMx[i_totalLayers][i_straw])
# 				dZI= float(layerIz[i_totalLayers][i_straw])
# 				dZM= float(layerMz[i_totalLayers][i_straw])
# 				#print "dXI= ", dXI, " dXM= ", dXM, " dZI= ", " dZM= ", dZM
# 				IdealX[i_module][i_view][i_layer][i_straw]=dXI
# 				IdealZ[i_module][i_view][i_layer][i_straw]=dZI
# 				MisX[i_module][i_view][i_layer][i_straw]=dXM
# 				MisZ[i_module][i_view][i_layer][i_straw]=dZM
# 			i_totalLayers+=1 #once we added all straws in that layer -> go to the next absolute layer


i_track = 0
with open("Tracker_tracks_fit_art.txt") as f:  #XXX
    for line in f:  #Line is a string
        number_str = line.split()    
        for i in range(0,4):
        	fit[i_track][i] = float(number_str[i])
        i_track+=1

fit_hitX=[]
fit_hitZ=[]
fit_hitRad=[]

with open("Tracker_hits_fit_art.txt") as f: 
	for line in f:  #Line is a string
		number_str = line.split()
		fit_hitX.append(float(number_str[0]))
		fit_hitZ.append(float(number_str[1]))
		fit_hitRad.append(float(number_str[2]))


##################PLOTING##############################
#Misaligned Geometry and Generated tracks 
plt.figure(1)

#Ideal Geometry and Fitted tracks 
plt.subplot(111)
axes2 = plt.gca()

#First draw tracks and straws with hits
for i_track in range(0, trackN):
	dataI = [[fit[i_track][1],fit[i_track][0]], [fit[i_track][3],fit[i_track][2]]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(dataI, 2))),
	    color = 'red', marker = 'x')
	
#Then draw all other straws 
# i_totalLayers=0
# for i_module in range(0, moduleN):
# 	for i_view in range(0, viewN):
# 		for i_layer in range(0, layerN):
# 			for i_straw in range(0, strawN):
# 				#circle = plt.Circle((IdealZ[i_module][i_view][i_layer][i_straw], IdealX[i_module][i_view][i_layer][i_straw]), 2.5, color='black', fill=False)
# 				#plt.plot(IdealZ[i_module][i_view][i_layer][i_straw], IdealX[i_module][i_view][i_layer][i_straw], color="black", marker = ",")
# 				#axes2.add_artist(circle)	
# 			i_totalLayers+=1 #once we read all straws in that layer -> go to the next absolute layer to get the Z coordinate

for i_hits in range(0, len(fit_hitX)):
	circle = plt.Circle((fit_hitZ[i_hits], fit_hitX[i_hits]), 2.55, color='black', linestyle='-', fill=False)
	# circle = plt.Circle((fit_hitZ[i_hits], fit_hitX[i_hits]), fit_hitRad[i_hits], color='green', linestyle='--', fill=False)
	axes2.add_artist(circle)	
	
axes2.set_ylim([20,50])
axes2.set_xlim([-600,90])
plt.xlabel("z [mm]")
plt.ylabel("x [mm]")
plt.title("Ideal Geometry with Reconstructed Tracks")

plt.show()