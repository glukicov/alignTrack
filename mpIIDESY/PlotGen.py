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
import csv
import pprint
import numpy as np  # smart arrays 
import itertools # smart lines 
import argparse, sys
from math import log10, floor
from matplotlib.ticker import MaxNLocator
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
import subprocess
import decimal


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

centresI = [ [0 for i_cord in xrange(2)  ] for i_module in xrange(moduleN)]
centresM = [ [0 for i_cord in xrange(2)  ] for i_module in xrange(moduleN)]

with open("Tracker_p_centre.txt") as f:
	for line in f:
		number_str=line.split()
		for i_module in range(0, moduleN):
			centresI[i_module][0]=float(number_str[2*i_module]) 
			centresI[i_module][1] =float(number_str[2*i_module+1])
			#print 'i_module=', i_module
		nextLine = next(f)
		number_str=nextLine.split()
		for i_module in range(0, moduleN):
			centresM[i_module][0] =float(number_str[2*i_module]) 
			centresM[i_module][1] =float(number_str[2*i_module+1])

#print centresI
#print centresM

####### ROOT Histos ##########

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
MisX = [[[[0 for i_straw in xrange(strawN)] for i_layer in xrange(layerN) ] for i_view in xrange(viewN)] for i_module in xrange(moduleN)]
MisZ = [[[[0 for i_straw in xrange(strawN)] for i_layer in xrange(layerN) ] for i_view in xrange(viewN)] for i_module in xrange(moduleN)]
IdealX = [[[[0 for i_straw in xrange(strawN)] for i_layer in xrange(layerN) ] for i_view in xrange(viewN)] for i_module in xrange(moduleN)]
IdealZ = [[[[0 for i_straw in xrange(strawN)] for i_layer in xrange(layerN) ] for i_view in xrange(viewN)] for i_module in xrange(moduleN)]
# Generated tracks and fitted tracks [x0, x1, z0, z1]
gen=[[0 for number in xrange(4)] for i_track in xrange(trackN)]
fit=[[0 for number in xrange(4)] for i_track in xrange(trackN)]
hitList=[[0 for number in xrange(toalLayerN)] for i_track in xrange(trackN)]


#Read files and store in lists
layerIx=[] #temp storage
layerIz=[] #temp storage
with open("Tracker_d_geom.txt") as f:
	for line in f:  #Line is a string
		layerIx.append(line.split())  # x
		#print "line=", line 
		nextLine = next(f)
		number_str=nextLine.split()
		#print "nextLine=", nextLine 
		layerIz.append(nextLine.split()) #z

layerMx=[] #temp storage
layerMz=[] #temp storage
with open("Tracker_d_mis.txt") as f:
	for line in f:  #Line is a string
		layerMx.append(line.split())
		nextLine = next(f)
		number_str=nextLine.split()
		layerMz.append(nextLine.split()) #z


misXSmeared = []

with open("Tracker_pede_mis.txt") as f:
	for line in f:  #Line is a string
		number_str=line.split()
		for i in range(0, len(number_str)):
			misXSmeared.append(float(number_str[i])) #z

misXSmeared=np.array(misXSmeared)

np.set_printoptions(precision=3) #to nearest um 

print "Misalignments:", misXSmeared*1e2, "[mm]"
misXSmeared = misXSmeared*1e4
print "Misalignments:", misXSmeared, "[um]"

MeanDecayPointZ = -1000 # mm 
print "MeanDecayPointZ= ",MeanDecayPointZ, "mm"


#Now for straws in X: 
i_totalLayers=0
for i_module in range(0, moduleN):
	for i_view in range(0, viewN):
		for i_layer in range(0, layerN):
			for i_straw in range(0, strawN):
				dXI= float(layerIx[i_totalLayers][i_straw])
				dXM= float(layerMx[i_totalLayers][i_straw])
				dZI= float(layerIz[i_totalLayers][i_straw])
				dZM= float(layerMz[i_totalLayers][i_straw])
				#print "dXI= ", dXI, " dXM= ", dXM, " dZI= ", " dZM= ", dZM
				IdealX[i_module][i_view][i_layer][i_straw]=dXI
				IdealZ[i_module][i_view][i_layer][i_straw]=dZI
				MisX[i_module][i_view][i_layer][i_straw]=dXM
				MisZ[i_module][i_view][i_layer][i_straw]=dZM
			i_totalLayers+=1 #once we added all straws in that layer -> go to the next absolute layer

# print "IdealX:: ", IdealX
# print "IdealZ:: ", IdealZ
# print "MisZ:: ", MisZ
# print "MisX:: ", MisX

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
fig=plt.figure(1)

yMin=-5.0
yMax=5.0
# xMin=beamStart-100
# xMax=beamStop-10

xMin=-100
xMax=+150

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
	plt.plot(centresM[i_module][0], centresM[i_module][1], color="red", marker = "*")
	for i_view in range(0, viewN):
		for i_layer in range(0, layerN):
			for i_straw in range(0, strawN):
				circle = plt.Circle((MisZ[i_module][i_view][i_layer][i_straw], MisX[i_module][i_view][i_layer][i_straw]), 0.25, color='black', fill=False)
				plt.plot(MisZ[i_module][i_view][i_layer][i_straw], MisX[i_module][i_view][i_layer][i_straw], color="black", marker = ",")
				axes.add_artist(circle)
			i_totalLayers+=1 #once we read all straws in that layer -> go to the next absolute layer to get the Z coordinate

for i_hits in range(0, len(gen_hitX)):
	circle3 = plt.Circle((gen_hitZ[i_hits], gen_hitX[i_hits]), gen_hitRad[i_hits], color='red', linestyle='--', fill=False)
	axes.add_artist(circle3)


#Do some statistics 
avgMean = ( sum(misXSmeared)/float(len(misXSmeared)) )
SD = np.sqrt(np.mean(misXSmeared**2)) 
#Plot the zeroth line
line = [[xMin,0.0], [xMax, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'grey', linestyle=":")		

axes.set_ylim([yMin, yMax])
axes.set_xlim([xMin, xMax])
#plt.xlabel("z [cm]")
plt.ylabel("x [cm]")
plt.title("Misaligned Geometry with Truth Tracks (Misalignment: <X> = %s um $\sigma$= %s um)" %(int(round(avgMean)), int(round(SD))), fontsize=18)


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
	plt.plot(centresI[i_module][0], centresI[i_module][1], color="red", marker = "*")
	for i_view in range(0, viewN):
		for i_layer in range(0, layerN):
			for i_straw in range(0, strawN):
				circle = plt.Circle((IdealZ[i_module][i_view][i_layer][i_straw], IdealX[i_module][i_view][i_layer][i_straw]), 0.25, color='black', fill=False)
				plt.plot(IdealZ[i_module][i_view][i_layer][i_straw], IdealX[i_module][i_view][i_layer][i_straw], color="black", marker = ",")
				axes2.add_artist(circle)	
			i_totalLayers+=1 #once we read all straws in that layer -> go to the next absolute layer to get the Z coordinate

for i_hits in range(0, len(fit_hitX)):
	circle = plt.Circle((fit_hitZ[i_hits], fit_hitX[i_hits]), fit_hitRad[i_hits], color='purple', linestyle='--', fill=False)
	axes2.add_artist(circle)

#Plot the zeroth line
line = [[xMin,0.0], [xMax, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'grey', linestyle=":")	
	
#And now some plot massaging
axes2.set_ylim([yMin, yMax])
axes2.set_xlim([xMin, xMax])
plt.xlabel("z [cm]", fontsize=18)
plt.ylabel("x [cm]", fontsize=18)
plt.title("Ideal Geometry with Reconstructed Tracks", fontsize=18)

# plt.show()
fig.set_size_inches(26,14)
plt.savefig("AUE.png")