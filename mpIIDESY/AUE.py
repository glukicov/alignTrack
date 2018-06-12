################
# Simple MC to estimate error from Misalignment  
#
# Gleb Lukicov 12th June 2018 
############
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

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-mis", "--mis")
parser.add_argument("-trackN", "--trackN")
args = parser.parse_args()

def getCentre(mod_lyr_strawPositionX, mod_lyr_strawPositionZ):
	CentreXY = [[0 for i in xrange(2)] for i_module in xrange(moduleN)]

	#for a module given by i_global get first and last straws of the first and last layers
	for i_global in range(0, moduleN):
		#Set the centre of a modules as a rotation point
		std::vector<float> U0_x = mod_lyr_strawPositionX[i_global][0][0];  // U0x
		std::vector<float> V1_x = mod_lyr_strawPositionX[i_global][1][1];  // V1x
		float U0_first_x = U0_x[0]; float V1_last_x = V1_x[strawN - 1];
		std::vector<float> U0_z = mod_lyr_strawPositionZ[i_global][0][0];  // z
		std::vector<float> V1_z = mod_lyr_strawPositionZ[i_global][1][1];  // z
		float U0_first_z = U0_z[0]; float V1_last_z = V1_z[strawN - 1];

		float z_c(0), x_c(0);
		x_c = (V1_last_x + U0_first_x) / 2;
		z_c = (V1_last_z + U0_first_z) / 2;

		centre.zCentres.push_back(z_c);
		centre.xCentres.push_back(x_c);
		plot_centres << z_c << " " << x_c << " ";
	}
	plot_centres << endl;
	return centre;


### Set Simulation Parameters #######

misX=float(args.mis)
trackN=int(args.trackN)

print "Generating modules offsets with Misalignment of", misX, "mm in X", misX*1e3, "um"
print "With", trackN, "tracks"

strawModuleZPosition = np.array([0.0, 134.3, 268.6, 402.9, 537.2, 671.5, 805.8, 940.1])
strawModuleXPosition = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

#Convert X position to um and add misalignments
#Smeared by:
mu, sigma = 0, 1 # mean and standard deviation
for i in range(0, len(strawModuleZPosition)):
	strawModuleXPosition[i]=strawModuleXPosition[i]+misX*np.random.normal(mu, sigma)

np.set_printoptions(precision=3) #to nearest um 

print "strawModuleZPosition", strawModuleZPosition, "[mm]"
print "strawModuleXPosition", strawModuleXPosition, "[mm]"

MeanDecayPointZ = 1000 # mm 

print "MeanDecayPointZ= ",MeanDecayPointZ, "mm"


########Construct Geometry [mm]####################
moduleN = 8 #number of movable detectors/module [independent modules]
strawN = 8 #number of measurement elements in x direction  [number of straws per layer]
viewN = 2 #There are two views per module (U and V) 
layerN = 2 #there are 2 layers per view [4 layers per module]
layerTotalN = layerN * viewN * moduleN # total number of layers
startingZDistanceStraw0 = 0.0 # distance of the very first layer [the first straw] relative to the "beam" in z [cm]
startingXDistanceStraw0 = 0.0 # distance of the very first layer [the first straw] in x [cm]
strawSpacing = 6.06 # x distance between straws in a layer
layerSpacing = 5.15 # z distance between layers in a view
viewSpacing = 20.20 # z distance between views in a modules
moduleSpacing = 134.3 # z distance between modules' first layers [first layer of module 1 and first layer of module 2]
layerDisplacement = 3.03 # relative x distance between first straws in adjacent layers in a view [upstream layer is +x shifted]

#Store detector coordinates 
MisX = [[[[0 for i_straw in xrange(strawN)] for i_layer in xrange(layerN) ] for i_view in xrange(viewN)] for i_module in xrange(moduleN)]
MisZ = [[[[0 for i_straw in xrange(strawN)] for i_layer in xrange(layerN) ] for i_view in xrange(viewN)] for i_module in xrange(moduleN)]
IdealX = [[[[0 for i_straw in xrange(strawN)] for i_layer in xrange(layerN) ] for i_view in xrange(viewN)] for i_module in xrange(moduleN)]
IdealZ = [[[[0 for i_straw in xrange(strawN)] for i_layer in xrange(layerN) ] for i_view in xrange(viewN)] for i_module in xrange(moduleN)]

dZ = startingZDistanceStraw0 # the increment in z for consecutive layers
for i_module in range(0, moduleN):
	dX = startingXDistanceStraw0 # starting on the x-axis (z, 0+disp)
	for i_view in range(0, viewN):
		for i_layer in range(0, layerN):
			for i_straw in range(0, strawN):
				# XXX Keeping explicit transformations for addition of rotation later (if needed)
				global_z = dZ
				global_x = dX
				MisX[i_module][i_view][i_layer][i_straw]=(global_x + strawModuleXPosition[i_module])
				#MisX[i_module][i_view][i_layer][i_straw]=(global_x)
				MisZ[i_module][i_view][i_layer][i_straw]=(global_z)
				print "global_x=",global_x,"global_z=",global_z
				dX =  dX - strawSpacing #while we are in the same layer: increment straw spacing in x
			#end of Straws loop
			if (i_view == 0): 
				dX = startingXDistanceStraw0 - layerDisplacement # set displacement in x for the next layer in the view
			if (i_view == 1):
				dX = startingXDistanceStraw0 # set displacement in x for the next layer in the view
			if (i_layer == 0): 
				dZ += layerSpacing # increment spacing between layers in a view once only
		#end of Layers loop
		if (i_view == 0): 
			dZ += viewSpacing # increment spacing between views in a module once only
	# end of View loop
	dZ += moduleSpacing
#Modules


##################PLOTING##############################
#Misaligned Geometry and Generated tracks 
plt.figure(1)
yMin=-50.0
yMax=50.0
xMin=-1100
xMax=1100
axes = plt.gca()

for i_module in range(0, len(strawModuleZPosition)):
	plt.scatter(strawModuleZPosition[i_module], strawModuleXPosition[i_module], color='black', marker="*", s=10)

#Then draw all other straws 
for i_module in range(0, moduleN):
	plt.plot(strawModuleZPosition[i_module], strawModuleXPosition[i_module], color="red", marker = "*")
	for i_view in range(0, viewN):
		for i_layer in range(0, layerN):
			for i_straw in range(0, strawN):
				circle = plt.Circle((MisZ[i_module][i_view][i_layer][i_straw], MisX[i_module][i_view][i_layer][i_straw]), 2.5, color='black', fill=False)
				plt.plot(MisZ[i_module][i_view][i_layer][i_straw], MisX[i_module][i_view][i_layer][i_straw], color="black", marker = ",")
				axes.add_artist(circle)

	
axes.set_ylim([yMin, yMax])
axes.set_xlim([xMin, xMax])
plt.xlabel("z [mm]")
plt.ylabel("x [mm]")
plt.title("Misaligned Modules with Reconstructed Tracks")
plt.savefig("AUE.png")