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
import decimal

parser = argparse.ArgumentParser(description='mode')
parser.add_argument("-mis", "--mis")
parser.add_argument("-trackN", "--trackN")
parser.add_argument('-misCM', '--misCM', nargs='+', help='misalignments')
args = parser.parse_args()

########Geometrical Constants [mm]####################
moduleN = 8 #number of movable detectors/module [independent modules]
strawN = 8 #number of measurement elements in x direction  [number of straws per layer]
viewN = 2 #There are two views per module (U and V) 
layerN = 2 #there are 2 layers per view [4 layers per module]
toalLayerN= layerN * viewN * moduleN # total number of layers
strawRadius = 2.535 # take as the outerRadiusOfTheGas
startingZDistanceStraw0 = 0.0 # distance of the very first layer [the first straw] relative to the "beam" in z [cm]
startingXDistanceStraw0 = strawN/2*6.06-strawRadius # distance of the very first layer [the first straw] in x [cm]
strawSpacing = 6.06 # x distance between straws in a layer
layerSpacing = 5.15 # z distance between layers in a view
viewSpacing = 20.20 # z distance between views in a modules
moduleSpacing = 133.14 # z distance between modules' first layers [first layer of module 1 and first layer of module 2]
layerDisplacement = 3.03 # relative x distance between first straws in adjacent layers in a view [upstream layer is +x shifted]

### Beam and Physics Parameters 
MeanDecayPointZ = -1000 # mm 
print "MeanDecayPointZ= ",MeanDecayPointZ, "mm"

beamPositionLength = 20.0 # x spread 
beamStart = MeanDecayPointZ 
beamStop = 1300 # z 
mX = 0.015  # uniform slope between -0.015 and 0.015 provides nice coverage for 8 straws
trackCut = 0.5 # 500 um
layerCut = 5

resolution = 0.150 #150 um  
useTruthLR = True
pValCut = 0.0

#Define some helper functions 

CentreXY = [[0 for i in xrange(2)] for i_module in xrange(moduleN)] # global scope
def getCentre(mod_lyr_strawPositionX, mod_lyr_strawPositionZ):
	

	#for a module given by i_global get first and last straws of the first and last layers
	for i_module in range(0, moduleN):
		#Set the centre of a modules as a rotation point
		U0_x = mod_lyr_strawPositionX[i_module][0][0];  # U0x
		V1_x = mod_lyr_strawPositionX[i_module][1][1];  # V1x
		U0_first_x = U0_x[0] 
		V1_last_x = V1_x[strawN - 1]
		U0_z = mod_lyr_strawPositionZ[i_module][0][0]  # z
		V1_z = mod_lyr_strawPositionZ[i_module][1][1]  # z
		U0_first_z = U0_z[0]
		V1_last_z = V1_z[strawN - 1]

		x_c = (V1_last_x + U0_first_x) / 2.0
		z_c = (V1_last_z + U0_first_z) / 2.0

		CentreXY[i_module][0]=(z_c)
		CentreXY[i_module][1]=(x_c)

def pointToLineDCA(zStraw, xStraw, xSlope, Intercept):

	#http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
	#converting from x=mz+c -> mz+x+c=0
	a = -xSlope;
	b = 1.0;
	c = -xIntercept;
	x0 = zStraw;
	y0 = xStraw;

	dca = np.fabs( a * x0 + b * y0 + c ) / np.sqrt( a * a + b * b  ) ;
	return dca

### Set Simulation Parameters #######

misXInput=float(args.mis)
trackN=int(args.trackN)

print "Generating modules offsets with smeared Misalignment of", misXInput, "mm in X", misXInput*1e3, "um"
print "With", trackN, "tracks"

misXSmeared = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
#Smeared by:
mu, sigma = 0, 1 # mean and standard deviation

if (misXInput != -1):
	
	#Convert X position to um and add misalignments
	
	for i in range(0, len(misXSmeared)):
		misXSmeared[i]=misXSmeared[i]+misXInput*np.random.normal(mu, sigma)

	np.set_printoptions(precision=3) #to nearest um 


	print "Misalignments:", repr(misXSmeared*0.1), "[cm]"
	print "Misalignments:", misXSmeared, "[mm]"
	print "Misalignments:", misXSmeared*1e3, "[um]"
	avgMean = ( sum(misXSmeared)/float(len(misXSmeared)) )
	SD = np.sqrt(np.mean(misXSmeared**2)) 
	print  "<Mis>", np.around(avgMean*1e3, decimals=3) ," [um] Mis SD:", np.around(SD*1e3, decimals=3) , " [um]"
else:
	for i in range(0, 8):
		misXSmeared[i]=args.misCM[i]
	
	np.set_printoptions(precision=3) #to nearest um 
	print "Misalignments:", misXSmeared, "[cm]"
	print "Misalignments:", misXSmeared*10.0, "[mm]"
	print "Misalignments:", misXSmeared*1e4, "[um]"
	avgMean = ( sum(misXSmeared)/float(len(misXSmeared)) )
	SD = np.sqrt(np.mean(misXSmeared**2)) 
	print  "<Mis>", np.around(avgMean*1e4, decimals=1) ," [um] Mis SD:", np.around(SD*1e4, decimals=1) , " [um]"


####### ROOT Histos ##########
Tf = TFile('Extrapolation.root', 'RECREATE')
h_residualRecon = TH1F("h_residualRecon", "Recon Residual (Track - Hit DCA); Residual [mm]",  99, -1.0, 1.0) 
h_residualTruth = TH1F("h_residualTruth", "Truth Track DCA - Hit DCA; Residual [mm]",  99, -1.0, 1.0)
h_hitDCA = TH1F("h_hitDCA", "Hit DCA; DCA [mm]",  90,  -0.5, 4.0) 
h_trackDCA = TH1F("h_trackDCA", "Track DCA; DCA [mm]",  90,  -0.5, 4.0)
h_trackTruthDCA = TH1F("h_trackTruthDCA", "Truth Track DCA; DCA [mm]",  90,  -0.5, 4.0) 

h_pValue = TH1F("h_pValue", "p-value", 100, 0.0, 1.0)
h_vertexTruth = TH1F("h_vertexTruth", "Vertex: Truth; [mm]", 100, -20.0, 20.0)
h_vertexReco = TH1F("h_vertexReco", "Vertex: Reco; [mm]", 100, -20.0, 20.0)
h_vertexResidual = TH1F("h_vertexResidual", "Vertex: #Delta (Recon - Truth)", 100, -20.0, 20.0)


########Construct Geometry [mm]####################

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
				#Create Misaligned and Ideal geometries
				MisX[i_module][i_view][i_layer][i_straw]=(global_x + misXSmeared[i_module])
				MisZ[i_module][i_view][i_layer][i_straw]=(global_z)
				IdealX[i_module][i_view][i_layer][i_straw]=(global_x)
				IdealZ[i_module][i_view][i_layer][i_straw]=(global_z)
				#print "global_x=",global_x,"global_z=",global_z
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


############ GENERATE HITS ########################

hitCount=0
trackCount=0
trackCutCount=0

gen=[[0 for number in xrange(4)] for i_track in xrange(trackN)]
fit=[[0 for number in xrange(4)] for i_track in xrange(trackN)]
hitList=[[0 for number in xrange(toalLayerN)] for i_track in xrange(trackN)]

for i_track in range(0, trackN):
	#print ""
	#print "i_track=", i_track
	#Track parameters for rand-generated line MC [start and end positions outside of detectors]
	x0 = beamPositionLength *np.random.uniform(0,1) - beamPositionLength/2.0;
	xIntercept = x0 
	h_vertexTruth.Fill(xIntercept)
	x1 = x0 # for parallel lines only
	signXSlope = 0.0 # for parallel lines only
	generalLines = True  # XXX quick hack
	if (generalLines == True):
		if (np.random.uniform(0,1) >= 0.5):
			signXSlope = 1.0
		else:
			signXSlope = -1.0

	xSlope = (np.random.uniform(0,1) * signXSlope) * mX
	x1 = xSlope * beamStop + xIntercept # "xExit"

	xRecon = []
	zRecon = []
	radRecon = []
	DCAs = []
	LRTruth = []

	hitsPerTrack=0
	i_hit=0
	cutTriggered = False
	for i_module in range(0, moduleN):
		for i_view in range(0, viewN):
			for i_layer in range(0, layerN):
				# true track position in-line with the layer [from line x=ym+c]
				xTrack = xSlope * MisX[i_module][i_view][i_layer][int(strawN / 2)] + xIntercept
				# position on the detector [from dca]
				#print "i_hit=", i_hit
				#Find the two closest x straw coordinates given x on the line - with same z [the vector of straw x is naturally in an ascending order]
				# xTrack is the point of the line "in-line with the layers"
				lower=0 
				upper=0 
				hitDistance=0
				lastID = strawN - 1 # the ID of the very last straw in the vector
				LR=0 
				index=0 #L= +1 or R=-1 hit looking downstream
				xLayer = MisX[i_module][i_view][i_layer]
				zLayer = MisZ[i_module][i_view][i_layer]
				# // If the very first straw [straw closest to the beam path], this straw has the highest x coordinate value in the distance-ordered descending straw-vector
				# // compares less than the hit, the hit must be from the L (+1)
				if (xTrack > xLayer[0]):
					#hit distance is then the dca from the L of the straw with highest x coordinate
					hitDistance = pointToLineDCA(zLayer[0], xLayer[0], xSlope, xIntercept)
					# print "hitDistance in above", hitDistance
					LR = +1 # assign truth L
					index = 0 #
					
				#Another scenario, is that the hit is smaller than the x coordinate of the last straw
				elif (xTrack < xLayer[lastID]):
					# hit distance is the dca from the R of the last straw
					hitDistance = pointToLineDCA(zLayer[0], xLayer[lastID], xSlope, xIntercept)
					# print "hitDistance in below", hitDistance
					LR = -1
					index = lastID

				# // All other cases will have a hit between two straws
				# // to check which DCA is actually shorter we need the calculate and compare
				# // the two [the straw-position vector iterator only checks the vertical distance] and assign LRma
				else:
					#we have already taken care of the end-straws, so now just need to look in between
					z_counter = 0;
					for i_counter in range(0, len(xLayer) ):
						if (xLayer[i_counter] < xTrack):
							lower = xLayer[i_counter];
							upper = xLayer[i_counter - 1];
							#print "xLayer[i_counter], xTrack, lower, upper, z_counter, i_counter", xLayer[i_counter], xTrack, lower, upper, z_counter, i_counter
							z_counter+=1
							break
							

					hit_distance_low = pointToLineDCA(zLayer[z_counter], lower, xSlope, xIntercept);
					hit_distance_up = pointToLineDCA(zLayer[z_counter - 1], upper, xSlope, xIntercept);
					#print "hit_distance_low ", hit_distance_low
					#print "hit_distance_up ", hit_distance_up

					if (hit_distance_low < strawRadius and hit_distance_up < strawRadius):
						sys.stderr.write("Multiple Layers Hit")

					#if DCA in higher straw (lower ID) is bigger, select straw ID with smaller DCA
					if (hit_distance_up > hit_distance_low):
						hitDistance = hit_distance_low
						LR = +1
						#unique and ordered straw positions in vector guarantee correct id
						index = [i for i,x in enumerate(xLayer) if x == lower]
				
					if (hit_distance_up < hit_distance_low):
						hitDistance = hit_distance_up;
						LR = -1
						index = [i for i,x in enumerate(xLayer) if x == upper]
					
					#if DCAs are equal, drop the dice
					if (hit_distance_up == hit_distance_low):
						sys.stderr.write("Ambiguity which straw registered hit")

				 # Finish looking for the correct hit 

				#Now smear the DCA data
				dcaUnsmeared = hitDistance;
				#hitDistanceSmeared = hitDistance ### TODO 
				hitDistanceSmeared = hitDistance + resolution * np.random.normal(mu, sigma);
				
				dca = hitDistanceSmeared
				DCAs.append(dca)
				LRTruth.append(LR)
				residualTruth = dca - hitDistance;
	
				#proceed only if hit is valid
				if (dca < strawRadius):
					#Find the truth ID, and LR hit for a straw
					xRecon.append(IdealX[i_module][i_view][i_layer][index[0]]); 
					zRecon.append(IdealZ[i_module][i_view][i_layer][index[0]]); 
					radRecon.append(dca);
					h_hitDCA.Fill(dca)
					#print "dca= ", dca

					hitsPerTrack+=1
					z_counter+=1;  # incrementing distance of planes
					# missed hit check

				# DCA cut if
				i_hit+=1
			#end of Layers loop
	    # end of View loop
	# end of looping over modules

	#//Apply a 500 um cut on the WHOLE track
	for i in range(0, len(DCAs)):
		if (DCAs[i] < trackCut ):
			cutTriggered = True
			trackCutCount +=1
			#print "Cut triggered"
			break
	#print DCAs
	#print hitsPerTrack
	#Now for the whole track, if not cut
	if (cutTriggered == False and hitsPerTrack >= layerCut):
	
		dataSize = len(radRecon)

		#LSR fit from: http://codesam.blogspot.com/2011/06/least-square-linear-regression-of-data.html
		SUMx = 0.0;     #//sum of x values
		SUMy = 0.0;     #//sum of y values
		SUMxy = 0.0;    #//sum of x * y
		SUMxx = 0.0;    #//sum of x^2
		SUMres = 0.0;   #//sum of squared residue
		res = 0.0;      #//residue squared
		slope = 0.0;    #//slope of regression line
		y_intercept = 0.0; #//y intercept of regression line
		SUM_Yres = 0.0; #//sum of squared of the discrepancies
		AVGy = 0.0;     #//mean of y
		AVGx = 0.0;     #//mean of x
		Yres = 0.0;     #//squared of the discrepancies
		Rsqr = 0.0;     #//coefficient of determina		
		#//calculate various sums 
		for i in range(0, dataSize):
		    #//sum of x
		    SUMx = SUMx + zRecon[i];
		    #//sum of y
		    SUMy = SUMy + xRecon[i];
		    #//sum of squared x*y
		    SUMxy = SUMxy +zRecon[i] * xRecon[i];
		   # //sum of squared x
		    SUMxx = SUMxx + zRecon[i] * zRecon[i];
		
		#//calculate the means of x and y
		AVGy = SUMy / dataSize;
		VGx = SUMx / dataSize;		
		#//slope or a1
		slope = (dataSize * SUMxy - SUMx * SUMy) / (dataSize * SUMxx - SUMx*SUMx);	
		#//y itercept or a0
		y_intercept = AVGy - slope * AVGx;		
		if (1==1):
		    print("x mean(AVGx) = %0.5E\n", AVGx);
		    print("y mean(AVGy) = %0.5E\n", AVGy);
		    print ("\n");
		    print ("The linear equation that best fits the given data:\n");
		    print ("       y = %2.8lfx + %2.8f\n", slope, y_intercept);
		    print ("------------------------------------------------------------\n");
		    print ("   Original (x,y)   (y_i - y_avg)^2     (y_i - a_o - a_1*x_i)^2\n");
		    print ("------------------------------------------------------------\n")
		#//calculate squared residues, their sum etc.
		for i in range(0, dataSize):
		    #//current (y_i - a0 - a1 * x_i)^2
		    Yres = np.power((xRecon[i] - y_intercept - (slope * zRecon[i])), 2);
		    
		    xFit = slope*zRecon[i]+y_intercept; 
		    #resData.x_fitted.append(xFit);
		    xRes = xFit - xRecon[i];
		    #resData.residuals.push_back(xRes); 
		    
		    #sum of (y_i - a0 - a1 * x_i)^2
		    SUM_Yres += Yres;	
		    #current residue squared (y_i - AVGy)^2
		    rres = np.power(xRecon[i] - AVGy, 2);	
		    #sum of squared residues
		    SUMres += res;
			#calculate r^2 coefficient of determination
    		Rsqr = (SUMres - SUM_Yres) / SUMres;	
		if (1==1):
		    print("--------------------------------------------------\n");
		    print("Sum of (y_i - y_avg)^2 = %0.5E\t\n", SUMres);
		    print("Sum of (y_i - a_o - a_1*x_i)^2 = %0.5E\t\n", SUM_Yres);
		    print("Standard deviation(St) = %0.5E\n", np.sqrt(SUMres / (dataSize - 1)));
		    print("Standard error of the estimate(Sr) = %0.5E\t\n", np.sqrt(SUM_Yres / (dataSize-2)));
		    print("Coefficient of determination(r^2) = %0.5E\t\n", (SUMres - SUM_Yres)/SUMres);
		    print("Correlation coefficient(r) = %0.5E\t\n", np.sqrt(Rsqr));		
    
    

    	# if (1==1):
    	# 	print slope*beamStart+y_intercept , slope*beamStop+y_intercept  , beamStart  , beamStop 
 
		gen[trackCount][0]=beamStart
		gen[trackCount][1]= x0
		gen[trackCount][2]=beamStop
		gen[trackCount][3]=x1

		fit[trackCount][0]=beamStart
		x0Recon = slope*beamStart + y_intercept
		fit[trackCount][1]= x0Recon
		fit[trackCount][2]=beamStop
		fit[trackCount][3]= slope*beamStop + y_intercept
			
		h_vertexReco.Fill(x0Recon)
		h_vertexResidual.Fill(xIntercept-x0Recon)

		hitCount+=len(radRecon);
		trackCount+=1

		# MC.residualsRecon = resData.residualsRecon;
		# MC.slopeRecon = resData.slopeRecon;
		# MC.interceptRecon = resData.interceptRecon;
		# MC.xStraw = xRecon;
		# MC.zStraw = zRecon;
		# MC.pValue = resData.pValue;
		# MC.chi2Circle = resData.chi2Circle;
		# MC.driftRad = radRecon; // vector

print "Total number of passed tracks ", trackCount, "with", hitCount, "hits"
print "Tracks cut due to DCA", trackCutCount

##################PLOTING##############################
fig = plt.figure(1)

#Misaligned Geometry and Generated tracks 
getCentre(MisX, MisZ)
plt.subplot(211)
yMin=-50.0
yMax=50.0
xMin=beamStart-100
xMax=beamStop-10
axes = plt.gca()

#Centre of Module (Z,X)
for i_module in range(0, moduleN):
	plt.scatter(CentreXY[i_module][0], CentreXY[i_module][1], color='green', marker="*", s=10)

#Draw all other straws 
for i_module in range(0, moduleN):
	for i_view in range(0, viewN):
		for i_layer in range(0, layerN):
			for i_straw in range(0, strawN):
				circle = plt.Circle((MisZ[i_module][i_view][i_layer][i_straw], MisX[i_module][i_view][i_layer][i_straw]), strawRadius, color='black', fill=False)
				plt.plot(MisZ[i_module][i_view][i_layer][i_straw], MisX[i_module][i_view][i_layer][i_straw], color="black", marker = ",")
				axes.add_artist(circle)

#Draw Truth tracks 
for i_track in range(0, trackCount):
	dataM = [[gen[i_track][0],gen[i_track][1]], [gen[i_track][2],gen[i_track][3]]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(dataM, 2))), color = 'blue', marker = 'x')

#Do some statistics 
#Plot the zeroth line
line = [[xMin,0.0], [xMax, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'grey', linestyle=":")

#And now some plot massaging
axes.set_ylim([yMin, yMax])
axes.set_xlim([xMin, xMax])
plt.xlabel("z [mm]", fontsize=18)
plt.ylabel("x [mm]", fontsize=18)
plt.title("Misaligned Geometry with Truth Tracks (Misalignment: <X> = %s um $\sigma$= %s um)" %(int(round(avgMean*1e3)), int(round(SD*1e3))), fontsize=18)


#Ideal Geometry and Reco tracks 
getCentre(IdealX, IdealZ)
plt.subplot(212)
axes = plt.gca()

#Centre of Module (Z,X)
for i_module in range(0, moduleN):
	plt.scatter(CentreXY[i_module][0], CentreXY[i_module][1], color='green', marker="*", s=10)

#Draw all other straws 
for i_module in range(0, moduleN):
	for i_view in range(0, viewN):
		for i_layer in range(0, layerN):
			for i_straw in range(0, strawN):
				circle = plt.Circle((IdealZ[i_module][i_view][i_layer][i_straw], IdealX[i_module][i_view][i_layer][i_straw]), strawRadius, color='black', fill=False)
				plt.plot(IdealZ[i_module][i_view][i_layer][i_straw], IdealZ[i_module][i_view][i_layer][i_straw], color="black", marker = ",")
				axes.add_artist(circle)

#Draw Recon tracks 
for i_track in range(0, trackCount):
	dataM = [[fit[i_track][0],fit[i_track][1]], [fit[i_track][2],fit[i_track][3]]]
	plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(dataM, 2))), color = 'red', marker = 'x')

#Do some statistics 
avgMean = ( sum(misXSmeared)/float(len(misXSmeared)) )
SD = np.sqrt(np.mean(misXSmeared**2)) 
#Plot the zeroth line
line = [[xMin,0.0], [xMax, 0.0]]
plt.plot(*zip(*itertools.chain.from_iterable(itertools.combinations(line, 2))),color = 'grey', linestyle=":")

#And now some plot massaging
axes.set_ylim([yMin, yMax])
axes.set_xlim([xMin, xMax])
plt.xlabel("z [mm]", fontsize=18)
plt.ylabel("x [mm]", fontsize=18)
plt.title("Ideal Geometry with Reconstructed Tracks (Extrapolation: <X> = %s um $\sigma$= %s um)" %(int(round(avgMean*1e3)), int(round(SD*1e3))), fontsize=18)



# plt.show()
fig.set_size_inches(26,14)
plt.savefig("AUE.png")


Tf.Write()
Tf.Close()