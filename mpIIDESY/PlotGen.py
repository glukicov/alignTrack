#!/usr/bin/python

####################################################################
# Sanity plots for AlginTracker
#
# 
#
# Created: 16 May 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
# Modified: 16 May 2017 by Gleb Lukicov (UCL) g.lukicov@ucl.ac.uk
#####################################################################

#import numpy as np 
import matplotlib.pyplot as plt #for plotting 
import itertools
import csv
import pprint


#Ideal geometry Z 
Izs=[]

# Misaligned geometry Z
Mzs=[]

# Getting constants from MC
with open("Tracker_p_constants.txt") as f:
	for line in f:  #Line is a string
		number_str = line.split()
		moduleN=int(number_str[0])
		viewN=int(number_str[1])
		layerN=int(number_str[2])
		strawN=int(number_str[3])
		trackN=int(number_str[4])
		beamX0=float(number_str[5])
		beamZ0=float(number_str[6])
		beamX1=float(number_str[7])
		beamZ1=float(number_str[8])


toalLayerN=layerN*moduleN*viewN

print "Parameters from Simulation:"
print "moduleN= ",moduleN
print "viewN= ",viewN
print "layerN= ",layerN
print "strawN= ",strawN
print "toalLayerN= ",toalLayerN
print "trackN= ",trackN
print "beamX0= ",beamX0
print "beamZ0= ",beamZ0
print "beamX1= ",beamX1
print "beamZ1= ",beamZ1

# X 4D arrays for Mis and Ideal Geom. 
Mis = [[[[0 for i_straw in xrange(strawN)] for i_layer in xrange(layerN) ] for i_view in xrange(viewN)] for i_module in xrange(moduleN)]
Ideal = [[[[0 for i_straw in xrange(strawN)] for i_layer in xrange(layerN) ] for i_view in xrange(viewN)] for i_module in xrange(moduleN)]
# Generated tracks and fitted tracks [x0, x1, z0, z1]
gen=[[0 for number in xrange(4)] for i_track in xrange(trackN)]
fit=[[0 for number in xrange(4)] for i_track in xrange(trackN)]
hitList=[[0 for number in xrange(toalLayerN)] for i_track in xrange(trackN)]


#Read files and store in lists
layerI=[] #temp storage
with open("Tracker_d_geom.txt") as f:
	for line in f:  #Line is a string
		layerI.append(line.split())

for i_layer in range(0, toalLayerN):
	Izs.append(float(layerI[i_layer][strawN])) #Z is the last element in the temp array 

layerM=[] #temp storage
with open("Tracker_d_mis.txt") as f:
	for line in f:  #Line is a string
		layerM.append(line.split())

for i_layer in range(0, toalLayerN):
	Mzs.append(float(layerM[i_layer][strawN]))

#Now for straws in X: 
i_totalLayers=0
for i_module in range(0, moduleN):
	for i_view in range(0, viewN):
		for i_layer in range(0, layerN):
			for i_straw in range(0, strawN):
				dXI= float(layerI[i_totalLayers][i_straw])
				dXM= float(layerM[i_totalLayers][i_straw])
				Ideal[i_module][i_view][i_layer][i_straw]=dXI
				Mis[i_module][i_view][i_layer][i_straw]=dXM
			i_totalLayers+=1 #once we added all straws in that layer -> go to the next absolute layer


#Read file and store in lists for tracks Generated and Fitted:
i_track = 0
with open("Tracker_p_gen.txt") as f:
    for line in f:  #Line is a string
        number_str = line.split()    
        for i in range(0,4):
        	gen[i_track][i] = float(number_str[i])
        #for hit in range(0, toalLayerN):
        #	hitList[i_track][hit] = number_str[4+hit] #hit list starts at 4th position
        i_track+=1

#print hitList

i_track = 0
with open("Tracker_p_fit.txt") as f:  #XXX
#with open("Tracker_p_gen.txt") as f:  #XXX
    for line in f:  #Line is a string
        number_str = line.split()    
        for i in range(0,4):
        	fit[i_track][i] = float(number_str[i])
        i_track+=1

gen_hitX=[]
gen_hitZ=[]
gen_hitRad=[]
fit_hitX=[]
fit_hitZ=[]
fit_hitRad=[]
with open("Tracker_p_hits_gen.txt") as f:
	for line in f:  #Line is a string
		number_str = line.split()
		gen_hitX.append(float(number_str[0]))
		gen_hitZ.append(float(number_str[1]))
		gen_hitRad.append(float(number_str[2]))

with open("Tracker_p_hits_fit.txt") as f: 
	for line in f:  #Line is a string
		number_str = line.split()
		fit_hitX.append(float(number_str[0]))
		fit_hitZ.append(float(number_str[1]))
		fit_hitRad.append(float(number_str[2]))


##################PLOTING##############################
#Misaligned Geometry and Generated tracks 
plt.figure(1)
plt.subplot(211)
axes = plt.gca()

#First draw tracks and straws with hits
for i_track in range(0, trackN):
	dataM = [[gen[i_track][2],gen[i_track][0]], [gen[i_track][3],gen[i_track][1]]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(dataM, 2))),
	    color = 'blue', marker = 'x')
	
#Then draw all other straws 
i_totalLayers=0
for i_module in range(0, moduleN):
	for i_view in range(0, viewN):
		for i_layer in range(0, layerN):
			for i_straw in range(0, strawN):
				circle = plt.Circle((Mzs[i_totalLayers], Mis[i_module][i_view][i_layer][i_straw]), 0.25, color='black', fill=False)
				plt.plot(Mzs[i_totalLayers], Mis[i_module][i_view][i_layer][i_straw], color="black", marker = ",")
				axes.add_artist(circle)
			i_totalLayers+=1 #once we read all straws in that layer -> go to the next absolute layer to get the Z coordinate

for i_hits in range(0, len(gen_hitX)):
	circle3 = plt.Circle((gen_hitZ[i_hits], gen_hitX[i_hits]), gen_hitRad[i_hits], color='red', linestyle='--', fill=False)
	axes.add_artist(circle3)		

axes.set_ylim([beamX0-6,beamX1+3])
axes.set_xlim([beamZ0-1,beamZ1+1])
#plt.xlabel("z [cm]")
plt.ylabel("x [cm]")
plt.title("Misaligned Geometry with True Tracks")

#Ideal Geometry and Fitted tracks 
plt.subplot(212)
axes2 = plt.gca()

#First draw tracks and straws with hits
for i_track in range(0, trackN):
	dataI = [[fit[i_track][2],fit[i_track][0]], [fit[i_track][3],fit[i_track][1]]]
	plt.plot(
	    *zip(*itertools.chain.from_iterable(itertools.combinations(dataI, 2))),
	    color = 'red', marker = 'x')
	
#Then draw all other straws 
i_totalLayers=0
for i_module in range(0, moduleN):
	for i_view in range(0, viewN):
		for i_layer in range(0, layerN):
			for i_straw in range(0, strawN):
				circle = plt.Circle((Izs[i_totalLayers], Ideal[i_module][i_view][i_layer][i_straw]), 0.25, color='black', fill=False)
				plt.plot(Izs[i_totalLayers], Ideal[i_module][i_view][i_layer][i_straw], color="black", marker = ",")
				axes2.add_artist(circle)	
			i_totalLayers+=1 #once we read all straws in that layer -> go to the next absolute layer to get the Z coordinate

for i_hits in range(0, len(fit_hitX)):
	circle = plt.Circle((fit_hitZ[i_hits], fit_hitX[i_hits]), fit_hitRad[i_hits], color='purple', linestyle='--', fill=False)
	axes2.add_artist(circle)	
	
axes2.set_ylim([beamX0-6,beamX1+3])
axes2.set_xlim([beamZ0-1,beamZ1+1])
plt.xlabel("z [cm]")
plt.ylabel("x [cm]")
plt.title("Ideal Geometry with Reconstructed Tracks")

plt.show()
