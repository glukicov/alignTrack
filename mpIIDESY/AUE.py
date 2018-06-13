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

########Geometrical Constants [mm]####################
moduleN = 8 #number of movable detectors/module [independent modules]
strawN = 8 #number of measurement elements in x direction  [number of straws per layer]
viewN = 2 #There are two views per module (U and V) 
layerN = 2 #there are 2 layers per view [4 layers per module]
layerTotalN = layerN * viewN * moduleN # total number of layers
strawRadius = 2.535 # take as the outerRadiusOfTheGas
startingZDistanceStraw0 = 0.0 # distance of the very first layer [the first straw] relative to the "beam" in z [cm]
startingXDistanceStraw0 = strawN/2*6.06-strawRadius # distance of the very first layer [the first straw] in x [cm]
strawSpacing = 6.06 # x distance between straws in a layer
layerSpacing = 5.15 # z distance between layers in a view
viewSpacing = 20.20 # z distance between views in a modules
moduleSpacing = 133.14 # z distance between modules' first layers [first layer of module 1 and first layer of module 2]
layerDisplacement = 3.03 # relative x distance between first straws in adjacent layers in a view [upstream layer is +x shifted]

### Beam and Physics Parameters 
beamPositionLength = 20.0 # x spread 
beamStop = (moduleSpacing + viewSpacing + layerSpacing*float(layerN)) * float(moduleN) # z 
mX = 0.015  # uniform slope between -0.015 and 0.015 provides nice coverage for 8 straws
trackCut = 0.5 # 500 um

resolution = 0.150 #150 um  
useTruthLR = True

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

#Convert X position to um and add misalignments
#Smeared by:
mu, sigma = 0, 1 # mean and standard deviation
for i in range(0, len(misXSmeared)):
	misXSmeared[i]=misXSmeared[i]+misXInput*np.random.normal(mu, sigma)

np.set_printoptions(precision=3) #to nearest um 

print "Misalignments:", misXSmeared, "[mm]"
print "Misalignments:", misXSmeared*1e3, "[um]"

MeanDecayPointZ = 1000 # mm 

print "MeanDecayPointZ= ",MeanDecayPointZ, "mm"


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

for i_track in range(0, trackN):
	hitsPerTrack=0
	#Track parameters for rand-generated line MC [start and end positions outside of detectors]
	x0 = beamPositionLength *np.random.uniform(0,1) - beamPositionLength/2.0; 
	xIntercept = x0 
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

	for i_module in range(0, moduleN):
		for i_view in range(0, viewN):
			for i_layer in range(0, layerN):
				# true track position in-line with the layer [from line x=ym+c]
				xTrack = xSlope * MisX[i_module][i_view][i_layer][int(strawN / 2)] + xIntercept
				# position on the detector [from dca]

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
					LR = +1 # assign truth L
					index = 0 #
					
				#Another scenario, is that the hit is smaller than the x coordinate of the last straw
				elif (xTrack < xLayer[lastID]):
					# hit distance is the dca from the R of the last straw
					hitDistance = pointToLineDCA(zLayer[0], xLayer[lastID], xSlope, xIntercept)
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
							z_counter+=1
							break
							#goto jmp;
					#jmp:

					hit_distance_low = pointToLineDCA(zLayer[z_counter], lower, xSlope, xIntercept);
					hit_distance_up = pointToLineDCA(zLayer[z_counter - 1], upper, xSlope, xIntercept);

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
				hitDistanceSmeared = hitDistance ### TODO 
				#hitDistanceSmeared = hitDistance + resolution * randomFacility->Gaus(0.0, 1.0);
				#//Apply a 500 um cut on the WHOLE track
				cutTriggered = False
				if (abs(hitDistanceSmeared) < trackCut ):
					cutTriggered = True # // will be checked in the MC_Launch on DCA return
					print "Track Cut"
				dca = hitDistanceSmeared
				residualTruth = dca - hitDistance;

				missedHit=False
				if (cutTriggered == False):
					if (dca > strawRadius):
						print "Hit Rejected: outside of straw layer with dca =", dca 
						z_counter+=1
						missedHit==True
					
					#proceed only if hit is valid
					if (missedHit == False):
						#Find the truth ID, and LR hit for a straw
						xRecon.append(IdealX[i_module][i_view][i_layer][index[0]]); 
						zRecon.append(IdealZ[i_module][i_view][i_layer][index[0]]); 
						radRecon.append(dca); 

						hitCount+=1;
						hitsPerTrack+=1
						z_counter+=1;  # incrementing distance of planes
					# missed hit check

				# DCA cut if
			#end of Layers loop
	    # end of View loop
	# end of looping over modules

	#Now for the whole track, if not cut
	if (cutTriggered == False):
		#from James's code: https://cdcvs.fnal.gov/redmine/projects/gm2tracker/repository/entry/teststand/StraightLineTracker_module.cc?utf8=%E2%9C%93&rev=feature%2FtrackDevelop
		nHits = hitsPerTrack # same for no hit rejection

		#These sums are parameters for the analytic results that don't change between LR combos (use U here but equally applicable to V coordinate)
		S, Sz, Su, Szz, Suu, Suz =0 # // also good declaration style fur custom types
		for i_hit in range(0, nHits): 
			z = zRecon[i_hit];
			u = xRecon[i_hit];
			err2 = np.power(resolution, 2); # the error is determined by the resolution [constant]
			S   += 1. / err2;
			Sz  += z / err2;
			Su  += u / err2;
			Szz += z * z / err2;
			Suu += u * u / err2;
			Suz += u * z / err2;
		# // hits

		# Number of LR combinations (2^N or 1 if using truth)
		nLRCombos = np.pow(2, nHits);
		if (useTruthLR):
			nLRCombos = 1;

		#Loop over all LR combinations and produce line fit for each one
		for LRCombo  in range (0, nLRCombos):
			#These sums are the other parameters for the analytic results
			Sr, Sru, Srz = 0
			for i_hit in range (0, nHits):
				z = zRecon[i_hit];
				u = xRecon[i_hit];
				err2 = pow(resolution, 2); # the error is determined by the resolution
				r = radRecon[i_hit];

				# Set r based on whether it's left (+ve r) or right (-ve r)
				if (LRTruth[i_hit] == -1): 
					r = -r;
				Sr  += r / err2;
				Sru += r * u / err2;
				Srz += r * z / err2;
			#hits

			#Make function of derivate of Chi-Squared w.r.t. gradient - quite an algebraically intensive calculation
			#Range assumes that we've hit more than one layer.
			dX2_dm = TF1("dX2_dm", "-x*x + [0]*x*sqrt(x*x+1) + [1]*x + [2]*sqrt(x*x+1) + 1", -40, 40)
			dX2_dm.SetParameter(0, (Sr * Su / S - Sru) / (Su * Sz / S - Suz));
			dX2_dm.SetParameter(1, ( (Su * Su / S - Suu) - (Sz * Sz / S - Szz) ) / (Su * Sz / S - Suz) );
			dX2_dm.SetParameter(2, (Sr * Sz / S - Srz) / (Su * Sz / S - Suz));

			# // Roots of this function are minima or maxima of Chi2
			# // Finding one that has positive derivative doesn't work, so fill in all roots and take one with best Chi2
			# // TF1::GetX(0) isn't too clever at finding roots - so we'll loop over the function range and give tighter range for root finding
			# // Holders for intercepts & gradients that satisfy Chi2 minimisation/maximisation
			gradients=[]
			intercepts=[]

			#Set some step size for loop - needs to be small enough that we don't miss roots where function crosses and re-crosses zero within this range
			stepSize = 0.1;

			#Loop over range and push back all roots to vectors
			prevValue = dX2_dm.Eval(dX2_dm.GetXmin());
			for mVal in range(dX2_dm.GetXmin() + stepSize, dX2_dm.GetXmax(), stepSize):
				newValue = dX2_dm.Eval(mVal);
				if (np.signbit(prevValue) != np.signbit(newValue)):
					m_tmp = dX2_dm.GetX(0, mVal - stepSize, mVal);
					gradients.append(m_tmp);
					intercepts.append( (Su - m_tmp * Sz + sqrt(m_tmp * m_tmp + 1)*Sr) / S );
				prevValue = newValue;
			# for mVal loop

			#Throw if we didn't find a root - something went wrong
			if (len(gradients) == 0):
				sys.stderr.write("No roots found from chi-squared minimisation function. Check step size and function range.")

			# Holders for final gradient/intercept result - initialised to 0 here to keep compiler happy, but should never make it through logic with these
			gradient = 0;
			intercept = 0;

			#Loop over possible gradient/intercepts and calculate chi2 value - then take lowest value as best gradient/intercept
			chi2ValMin = float("inf")
			for grad in range(0, len(gradients)):

				chi2Val = 0;
				for i_hit in range(0, nHits):

					z = zRecon[i_hit];
					u = xRecon[i_hit];
					r = radRecon[i_hit];
					err2 = pow(resolution, 2); # the error is determined by the resolution

					# Set r based on whether it's left (+ve r) or right (-ve r)
					if (LRTruth[i_hit] == -1):
					 	r = -r;

					#Calculate distance of track from wire and use it for Chi2 calculation
					d = (gradients[grad] * z + intercepts[grad] - u) / sqrt(gradients.at(grad) * gradients.at(grad) + 1);
					#double d = Tracker::pointToLineDCA(z, u, gradient, intercept);
					chi2Val += pow(d - r, 2) / err2;

					if (debugBool && StrongDebugBool) {
						cout << "grad= " << grad << " gradients.size()= " << gradients.size() << " chi2Val " << chi2Val
						     << " z= " << z << " u= "  << u <<  " r= " << r << " LRTruth= " << LRTruth[i_hit] << " err2= " << err2
						     << " d= " << d << " intercepts.at(grad)= " << intercepts.at(grad)
						     << " gradients.at(grad)= " << gradients.at(grad) << endl;
					}

				} // hits

				// Store gradient/intercept for lowest chi2 val;
				if (chi2Val < chi2ValMin) {
					gradient  = gradients.at(grad);
					intercept = intercepts.at(grad);
					chi2ValMin = chi2Val;
				}
			} // end of grad/inter.

			// Convert to p-value and add track to vector if it passes p-value cut
			double pVal = TMath::Prob(chi2ValMin, nHits - 2); //Two fit parameters
			resData.pValue = pVal;
			resData.chi2Circle = chi2ValMin;

			if (debugBool && StrongDebugBool) {cout << "pVal=" << pVal << " chi2ValMin= " << chi2ValMin << endl << endl;}
			if (pVal > pValCut) {
				// We'll want to store left/right hits so set these
				for (int i_hit = 0; i_hit < nHits; i_hit++) {
					//XXX Now that we know the slope and the gradient of the best fit line through drift circles,
					// we can calculate the residual for each drift circle, which is
					// "DCA from that line to the straw centre" - "Radius of the drift circle"
					float residualRecon = Tracker::pointToLineDCA(zRecon[i_hit],  xRecon[i_hit], gradient, intercept) - radRecon[i_hit];
					resData.residualsRecon.push_back(residualRecon);   // residual between the (centre of the straw and the fitted line [pointToLineDCA]) and radius of the fit circle;
				} // hits

				// Passing recon track parameters to MC
				resData.slopeRecon = gradient;     // slope of the best fit line
				resData.interceptRecon = intercept; // intercept
				// Python plotting
				if (debugBool && cutTriggered == false) plot_fit <<  gradient*beamStart + intercept << " "  << gradient*beamStop + intercept  <<   " " <<  beamStart  << " " << beamStop << endl;
			} // p-value cut
		} // LRCombinations [once for useTruthLR]

		return resData;

		trackCount+=1
		# MC.residualsRecon = resData.residualsRecon;
		# MC.slopeRecon = resData.slopeRecon;
		# MC.interceptRecon = resData.interceptRecon;
		# MC.xStraw = xRecon;
		# MC.zStraw = zRecon;
		# MC.pValue = resData.pValue;
		# MC.chi2Circle = resData.chi2Circle;
		# MC.driftRad = radRecon; // vector
	if (cutTriggered == True):
		break  
			


##################PLOTING##############################
fig = plt.figure(1)

#Misaligned Geometry and Generated tracks 
getCentre(MisX, MisZ)
plt.subplot(211)
yMin=-50.0
yMax=50.0
xMin=-1200
xMax=1300
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
