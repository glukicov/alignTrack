/*
This source file contains definitions of various functions in the method class.
*/

#include "AlignTracker_methods.h"

using namespace std;

// Set up empty pointer for instance of the class.
Tracker* Tracker::s_instance = NULL;

/**
Constructor for the tracker class.
 */
Tracker::Tracker() {
	// Set mapping for U0...V1
	string tempMapping[4] = {"U0", "U1", "V0", "V1"};
	for (int i_view = 0; i_view < viewN; i_view++) {
		UVmapping.push_back(vector<string> ()); //initialize the first index with a 2D vector
		for (int i_layer = 0; i_layer < layerN; i_layer++) {
			if (i_view == 0) {UVmapping[i_view].push_back(tempMapping[i_layer]);}
			if (i_view == 1) {UVmapping[i_view].push_back(tempMapping[i_layer + 2]);}
		} //layers
	} // views

} // constructor

/**
    Empty destructor for tracker class.
*/
Tracker::~Tracker() {
}

/**
   Get pointer to only instance of tracker class, creating this instance if it doesn't already exist.
   instance is a static member function of Tracker (singleton)
   @return Pointer to tracker instance.
 */
Tracker* Tracker::instance() {
	// Create pointer to class instance if one doesn't exist already, then return that pointer.
	if (s_instance == NULL) s_instance = new Tracker();
	return s_instance;
}

/**
   Write a steering file to the supplied file-stream.
    @param steering_file reference to ofstream to write steering file to.
 */
void Tracker::writeSteeringFile(ofstream& steering_file, ofstream& metric) {
	if (steering_file.is_open()) {

		std::stringstream pede_method; pede_method.str(""); pede_method << "method inversion 5 0.001";
		metric << "| " << pede_method.str().c_str();
		std::stringstream msg_method;
		msg_method << Logger::yellow() << pede_method.str().c_str();
		Logger::Instance()->write(Logger::NOTE, msg_method.str());

		steering_file <<  "* g-2 Tracker Alignment: PEDE Steering File" << endl
		              << " "  << endl
		              << "Tracker_con.txt   ! constraints text file (if applicable) " << endl
		              << "Tracker_par.txt   ! parameters (presgima) text file (if applicable)" << endl
		              << "Cfiles ! following bin files are Cfiles" << endl
		              << "Tracker_data.bin   ! binary data file" << endl
		              << "method inversion 5 0.001" << endl
		              << "printrecord 2 -1 ! produces mpdebug.txt for record 2 with the largest value of χ2/Ndf" << endl
		              << " "  << endl
		              << "end ! optional for end-of-data" << endl;
	} // steering file open
} // steering function

/**
   Write a constraint file to the supplied file-stream.
    @param constraint_file reference to ofstream to write constraint file to.
 */
void Tracker::writeConstraintFile(ofstream& constraint_file, ofstream& debug_con, bool debugBool, ofstream& metric) {
	// Check constraints file is open, then write.
	if (constraint_file.is_open()) {
		//Evaluation of constraints
		float one = 1.0;
		metric << " | C: ";
		stringstream labelt;
		// Given number of constraints
		int i_module = 0; //select first module
		for (int i_NC = 0; i_NC < matNC ; i_NC++) {

			constraint_file << "Constraint 0.0" << endl;
			//labelt << "-; "; // == no constraintss // XXX
			int labelt = i_module + 1; // Millepede accepts +ive labels only
			constraint_file << labelt << " " << fixed << setprecision(5) << one << endl;
			i_module = moduleN - 1 ; // select last module

		} // end of NC loop
		metric << labelt.str().c_str();
	} // constrain file open
	cout << "Memory space requirement (inversion method, i.e. upper bound) = " << ( matN * matN + matN ) / 2 + matN * matNC + ( matNC * matNC + matNC ) / 2 << endl;
} // end of writing cons file

/**
   Write a presigma parameter file to the supplied file-stream.
    @param presigma_file reference to ofstream to write presigma file to.
 */
void Tracker::writePresigmaFile(ofstream& presigma_file, ofstream& metric) {
	// Check constraints file is open, then write.
	if (presigma_file.is_open()) {
		presigma_file << "PARAMETERS" << endl;
		metric << " | P: ";
		//Fixing module 0 and the last module
		for (int i_module = 0; i_module < moduleN; i_module++) {
			for (int i_par = 0; i_par < ngl; i_par++) {
				if (i_module == 0 || i_module == moduleN - 1) {
					//if (i_module == 1 || i_module == 2) {
					float initialValue = 0.0; //modules at x0
					float preSigma = -1.0;
					int labelM = i_module + 1; // Millepede accepts +ive labels only
					int labelP = i_par + 1; // Millepede accepts +ive labels only
					presigma_file << labelM << labelP << " " << fixed << setprecision(5) << initialValue << fixed << setprecision(5) << " " << preSigma << endl;
					metric << labelM << labelP << " " << initialValue << " " << preSigma << "; ";
				} // end of fixed modules
			} //end of global parameter loop
		} // end of detectors loop
	} // presigma file open
} // end of writing presigma file

/**
   Define centre of rotation
   @return rotation centres for each module
 */
RotationCentres Tracker::getCentre(std::vector< std::vector< std::vector< std::vector< float > > > > mod_lyr_strawPositionX, std::vector< std::vector< std::vector< std::vector< float > > > > mod_lyr_strawPositionZ, ofstream& plot_centres) {
	RotationCentres centre;
	// for a module given by i_global get first and last straws of the first and last layers
	for (int i_global = 0; i_global < moduleN; i_global++) {

		//Set the centre of a modules as a rotation point
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
}

/**
   DCA for a point to a line in 2D space
   @return dca
*/
float Tracker::pointToLineDCA(float zStraw, float xStraw, float xSlope, float xIntercept) {

	//http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
	//converting from x=mz+c -> mz+x+c=0
	float a = -xSlope;
	float b = 1.0;
	float c = -xIntercept;
	float x0 = zStraw;
	float y0 = xStraw;

	float dca = abs( a * x0 + b * y0 + c ) / sqrt( a * a + b * b  ) ;
	return dca;
}

/** Uses DCA function to find the shortest dca between straws (i.e. which straw was hit in that layer)
    @Inputs (see DCA method) + vector of straws' x coordinates in a layer, and a return type: "dca_hit"  or "x_line"
    @return hit_distance - dca of the straw that was hit
*/
DCAData Tracker::DCAHit(vector<float> xLayer, vector<float> zLayer, float xTrack, float xSlpoe, float xIntercept, bool debugBool) {

	DCAData dcaData;

	bool StrongDebugBool = false; //quick hack XXX for even more debug output

	//Find the two closest x straw coordinates given x on the line - with same z [the vector of straw x is naturally in an ascending order]
	// xTrack is the point of the line "in-line with the layers"
	float lower, upper, hitDistance;
	int lastID = strawN - 1; // the ID of the very last straw in the vector
	int LR, index; //L= +1 or R=-1 hit looking downstream
	vector<float>::iterator it;

	// If the very first straw [straw closest to the beam path], this straw has the highest x coordinate value in the distance-ordered descending straw-vector
	// compares less than the hit, the hit must be from the L (+1)
	if (xTrack > xLayer[0]) {
		// hit distance is then the dca from the L of the straw with highest x coordinate
		hitDistance = pointToLineDCA(zLayer[0], xLayer[0], xSlpoe, xIntercept);
		LR = +1; // assign truth L
		index = 0;
		if (debugBool && StrongDebugBool) {
			cout << "Track in line at " << xTrack << " The first straw is closest at " << upper <<  "; with DCA " << hitDistance << endl;
		}
	}
	// Another scenario, is that the hit is smaller than the x coordinate of the last straw
	else if (xTrack < xLayer[lastID]) {
		// hit distance is the dca from the R of the last straw
		hitDistance = pointToLineDCA(zLayer[0], xLayer[lastID], xSlpoe, xIntercept);
		LR = -1;
		index = lastID;
		if (debugBool && StrongDebugBool) {
			cout << "Track in line at " << xTrack << " The last straw is closest at " << lower <<  "; with DCA " << hitDistance << endl;
		}
	}
	// All other cases will have a hit between two straws
	// to check which DCA is actually shorter we need the calculate and compare
	// the two [the straw-position vector iterator only checks the vertical distance] and assign LRma
	else {
		//we have already taken care of the end-straws, so now just need to look in between
		int z_counter = 0;
		for (int i_counter = 0; i_counter < xLayer.size(); i_counter++) {
			if (xLayer[i_counter] < xTrack) {
				lower = xLayer[i_counter];
				upper = xLayer[i_counter - 1];
				z_counter++;
				if (debugBool && StrongDebugBool) cout << "lower= " << lower << " upper " << upper << endl;
				//as soon as we find a single straw that is at lower x than the hit,
				// our search is over
				goto jmp;
			}
		}
jmp:

		float hit_distance_low = pointToLineDCA(zLayer[z_counter], lower, xSlpoe, xIntercept);
		float hit_distance_up = pointToLineDCA(zLayer[z_counter - 1], upper, xSlpoe, xIntercept);

		if (hit_distance_low < strawRadius && hit_distance_up < strawRadius) {
			if (debugBool) {cout << "Multiple straws in layer were hit!" << endl;}
			incMultipleHitsLayer();
		}

		// if DCA in higher straw (lower ID) is bigger, select straw ID with smaller DCA
		if (hit_distance_up > hit_distance_low) {
			hitDistance = hit_distance_low;
			LR = +1;
			// unique and ordered straw positions in vector guarantee correct id
			it = std::find(xLayer.begin(), xLayer.end(), lower);
			index = std::distance(xLayer.begin(), it);
		}
		if (hit_distance_up < hit_distance_low) {
			hitDistance = hit_distance_up;
			LR = -1;
			it = std::find(xLayer.begin(), xLayer.end(), upper);
			index = std::distance(xLayer.begin(), it);
		}
		//if DCAs are equal, drop the dice... XXX won't be relevant with hit rejection
		if (hit_distance_up == hit_distance_low) {
			float random = randomFacility->Uniform(0.0, 1.0);
			if (random < 0.5) {
				hitDistance = hit_distance_low;
			}
			if (random > 0.5) {
				hitDistance = hit_distance_up;
			}
			cout << "Ambiguity which straw registered hit" << endl;
			incAmbiguityHit();
		}
		if (debugBool && StrongDebugBool) {
			cout <<  "Track in Line " << xTrack << "Two straws closest to that point are " << lower << ", and " << upper <<  "; with DCAs " << hit_distance_low <<  " and " << hit_distance_up << ", respectively." << endl;
		}
	} // end of iterator to find straws between hits

	//Now smear the DCA data
	dcaData.dcaUnsmeared = hitDistance;
	float hitDistanceSmeared = hitDistance + Tracker::instance()->getResolution() * randomFacility->Gaus(0.0, 1.0);
	//Apply a 500 um cut on the WHOLE track
	if (abs(hitDistanceSmeared) < trackCut && trackCutBool==true) {
		cutTriggered = true; // will be checked in the MC_Launch on DCA return
	}

	dcaData.dca = hitDistanceSmeared;
	float residualTruth = hitDistanceSmeared - hitDistance;
	dcaData.residualTruth = residualTruth;
	dcaData.LRSign = LR;
	dcaData.strawID = index;

	if (debugBool && StrongDebugBool) {
		cout << "Selected DCA as the correct hit distance is " << hitDistance << ". Straw ID: " << dcaData.strawID;
		if (LR  < 0) {cout << ". The straw was hit from the right" << endl;}
		if (LR > 0) {cout << ". The straw was hit from the left" << endl;}
	}//debug
	return dcaData; // as dca and id of the closest straw
}


// Adds DCA to the ideal geometry
ReconData Tracker::HitRecon(int detID, float detDca, vector<float> xLayer, vector<float> zLayer) {

	ReconData reconData;
	bool StrongDebugBool = false; //quick hack XXX for even more debug output
	float xIdealStraw = xLayer[detID];
	float reconDca = detDca; //adds dca to the assumed straw x coordinate
	reconData.z = zLayer[detID];
	reconData.x = xIdealStraw;
	reconData.dcaRecon = reconDca;
	return reconData;
}

// Function to return residuals to the fitted line (due to dca point scatter + resolution of the detector)
// @ Inputs: z,x position of straws, radius of drift circle, # circles, plotting file for input, verbosity flag, truth flag
// Truth can be used as an optional input [e.g. LR information]
// @return vector of residuals for each straw
ResidualData Tracker::GetResiduals(vector<float> zRecon, vector<float> xRecon, vector<float> radRecon, int dataSize, ofstream& plot_fit, bool debugBool, bool useTruthLR, vector<int> LRTruth) {

	ResidualData resData; // return slope, intercept of the fitted track

	bool StrongDebugBool = false; //quick hack XXX for even more debug output

	// from James's code: https://cdcvs.fnal.gov/redmine/projects/gm2tracker/repository/entry/teststand/StraightLineTracker_module.cc?utf8=%E2%9C%93&rev=feature%2FtrackDevelop
	// line 392 onwards, inputs to the original function: vector<DriftCircle>& circles, double pValCut, long long truthLRCombo
	int nHits = dataSize; // same for no hit rejection

	// These sums are parameters for the analytic results that don't change between LR combos (use U here but equally applicable to V coordinate)
	double S(0), Sz(0), Su(0), Szz(0), Suu(0), Suz(0); // also good declaration style fur custom types
	for (int i_hit = 0; i_hit < nHits; i_hit++) {
		double z = zRecon[i_hit];
		double u = xRecon[i_hit];
		double err2 = pow(Tracker::instance()->getResolution(), 2); // the error is determined by the resolution [constant]
		S   += 1. / err2;
		Sz  += z / err2;
		Su  += u / err2;
		Szz += z * z / err2;
		Suu += u * u / err2;
		Suz += u * z / err2;
	} // hits

	// Number of LR combinations (2^N or 1 if using truth)
	int nLRCombos = pow(2, nHits);
	if (useTruthLR) nLRCombos = 1;

	// Loop over all LR combinations and produce line fit for each one
	for (int LRCombo = 0; LRCombo < nLRCombos; LRCombo++) {
		// These sums are the other parameters for the analytic results
		double Sr(0), Sru(0), Srz(0);
		for (int i_hit = 0; i_hit < nHits; i_hit++) {
			double z = zRecon[i_hit];
			double u = xRecon[i_hit];
			double err2 = pow(Tracker::instance()->getResolution(), 2); // the error is determined by the resolution
			double r = radRecon[i_hit];

			// Set r based on whether it's left (+ve r) or right (-ve r)
			if (LRTruth[i_hit] == -1) r = -r;
			Sr  += r / err2;
			Sru += r * u / err2;
			Srz += r * z / err2;
		} // hits

		// Make function of derivate of Chi-Squared w.r.t. gradient - quite an algebraically intensive calculation
		// Range assumes that we've hit more than one layer.
		TF1* dX2_dm = new TF1("dX2_dm", "-x*x + [0]*x*sqrt(x*x+1) + [1]*x + [2]*sqrt(x*x+1) + 1", -40, 40);
		dX2_dm->SetParameter(0, (Sr * Su / S - Sru) / (Su * Sz / S - Suz));
		dX2_dm->SetParameter(1, ( (Su * Su / S - Suu) - (Sz * Sz / S - Szz) ) / (Su * Sz / S - Suz) );
		dX2_dm->SetParameter(2, (Sr * Sz / S - Srz) / (Su * Sz / S - Suz));

		// Roots of this function are minima or maxima of Chi2
		// Finding one that has positive derivative doesn't work, so fill in all roots and take one with best Chi2
		// TF1::GetX(0) isn't too clever at finding roots - so we'll loop over the function range and give tighter range for root finding
		// Holders for intercepts & gradients that satisfy Chi2 minimisation/maximisation
		vector<double> gradients, intercepts;

		// Set some step size for loop - needs to be small enough that we don't miss roots where function crosses and re-crosses zero within this range
		double stepSize = 0.1;

		// Loop over range and push back all roots to vectors
		double prevValue = dX2_dm->Eval(dX2_dm->GetXmin());
		for (double mVal = dX2_dm->GetXmin() + stepSize; mVal <= dX2_dm->GetXmax(); mVal += stepSize) {
			double newValue = dX2_dm->Eval(mVal);
			if (signbit(prevValue) != signbit(newValue)) {  // if sign doesn't match
				double m_tmp = dX2_dm->GetX(0, mVal - stepSize, mVal);
				gradients.push_back(m_tmp);
				intercepts.push_back( (Su - m_tmp * Sz + sqrt(m_tmp * m_tmp + 1)*Sr) / S );
				if (debugBool && StrongDebugBool) {cout << "m_tmp= " << m_tmp << " mVal= " << mVal << " prevValue= " << prevValue << " newValue= " << newValue << endl;}
			}
			prevValue = newValue;
			if (debugBool && StrongDebugBool) {cout << "newValue= " << newValue << endl;}
		} // for mVal loop
		delete dX2_dm;

		// Throw if we didn't find a root - something went wrong
		if (gradients.size() == 0) {
			stringstream exception2;
			exception2 << "StraightLineTracker::calculateUVLineFits" << "No roots found from chi-squared minimisation function. Check step size and function range.\n";
			Logger::Instance()->write(Logger::WARNING, exception2.str());
		}

		// Holders for final gradient/intercept result - initialised to 0 here to keep compiler happy, but should never make it through logic with these
		double gradient = 0;
		double intercept = 0;

		// Loop over possible gradient/intercepts and calculate chi2 value - then take lowest value as best gradient/intercept
		double chi2ValMin = std::numeric_limits<double>::infinity();
		for (unsigned int grad = 0; grad < gradients.size(); grad++) {

			double chi2Val = 0;
			for (int i_hit = 0; i_hit < nHits; i_hit++) {

				double z = zRecon[i_hit];
				double u = xRecon[i_hit];
				double r = radRecon[i_hit];
				double err2 = pow(Tracker::instance()->getResolution(), 2); // the error is determined by the resolution

				// Set r based on whether it's left (+ve r) or right (-ve r)
				if (LRTruth[i_hit] == -1) r = -r;

				// Calculate distance of track from wire and use it for Chi2 calculation
				double d = (gradients.at(grad) * z + intercepts.at(grad) - u) / sqrt(gradients.at(grad) * gradients.at(grad) + 1);
				//double d = Tracker::pointToLineDCA(z, u, gradient, intercept);
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
}

/**
   Simulate the passage of a randomly generated linear track through the detector, calculating properties of hits by the track on detector planes.
   @return MCData struct containing data about detector hits.
*/
MCData Tracker::MCLaunch(float scatterError, ofstream& debug_calc, ofstream& debug_off, ofstream& debug_mc, ofstream& plot_fit, ofstream& plot_gen, ofstream& plot_hits_gen, ofstream& plot_hits_fit, bool debugBool) {

	bool StrongDebugBool = false; // HACK XXX set by hand for now
	// Set up new container for track data, with hit count set to zero
	MCData MC;
	int layerCount = 0; //count only counts through detector
	cutTriggered = false; //set trigger to OFF until DCA is less than 500 um, then kill this track

	//clearing containers for next track
	xRecon.clear();  // ideal straw x
	zRecon.clear();  // ideal straw z
	radRecon.clear();

	//Track parameters for rand-generated line MC [start and end positions outside of detectors]
	// redefining the track as x=ym+c
	float x0 = beamPositionLength * randomFacility->Uniform(0.0, 1.0) - 1.0; //uniform vertex
	float xIntercept = x0; // by definition

	float x1 = x0; // for parallel lines only
	float xSlope = 0.0; // for parallel lines only

	bool generalLines = true;    // XXX quick hack
	if (generalLines == true) {

		float signXSlope;
		if (randomFacility->Uniform(0.0, 1.0) >= 0.5) {
			signXSlope = 1.0;
		}
		else {
			signXSlope = -1.0;
		}
		xSlope = (randomFacility->Uniform(0.0, 1.0) * signXSlope) * 0.015; // unifrom slope between -0.015 and 0.015: provides nice coverage for 8 straws
		x1 = xSlope * beamStop + xIntercept; // "xExit"

	} // end of generalLines == true HACK

	float xTrack; //true track position x=zm+c

	MC.x0 = x0;
	MC.x1 = x1;
	MC.slopeTruth = xSlope;
	MC.interceptTruth = xIntercept;
	//std::cout << " " << "\n";
	//Now get the centres of the ideal modules (not truth) and pass to main
	RotationCentres idealCentres = Tracker::getCentre(mod_lyr_strawIdealPositionX, mod_lyr_strawIdealPositionZ, debug_calc);

	//The main loop is for modules [they produce label of Global Parameters]:
	// Then looping over layers and views
	// within layers we will have loops over straws in that layer [to find out which straw was hit and dca]
	int z_counter = 0; // planes in z are incremented from 0 to TotalLayerN
	if (debugBool && StrongDebugBool) {cout << "Calculating hits/dca:" << endl;}
	//loops over modules, views, layers
	for (int i_module = 0; i_module < moduleN; i_module++) {
		for (int i_view = 0; i_view < viewN; i_view++) {
			for (int i_layer = 0; i_layer < layerN; i_layer++) { //per layer
				//missedHit = false; // set flag to false until triggered
				//std::cout << "layerCount= " << layerCount << "\n";
				// TRUTH PARAMETERS
				//The registered hit position on the misaligned detector is smeared by its resolution
				// zTrack is just the distance[z_counter] - in-line with the layer [chose the middle straw]
				xTrack = xSlope * mod_lyr_strawIdealPositionZ[i_module][i_view][i_layer][int(strawN / 2)] + xIntercept; // true track position in-line with the layer [from line x=ym+c]
				// DCA will return the pointToLineDCA (smeared by resolution) and strawID [+truth parameters: LR sign etc. xHit and zHit]
				DCAData MisDetector = Tracker::DCAHit(mod_lyr_strawMisPositionX[i_module][i_view][i_layer], mod_lyr_strawMisPositionZ[i_module][i_view][i_layer], xTrack, xSlope, xIntercept, debugBool); // position on the detector [from dca]
				//Check for trigger on DCA < 500 um, and only continue calculations not triggered
				if (cutTriggered) {MC.cut = true;} // set cut to true for the main, and return MC.cut as true only.
				if (cutTriggered == false) {
					float dca =  MisDetector.dca;  //dca (smeared) of the hit straw
					// proceed only if hit is valid
					if (dca <= strawRadius) {
						//std::cout << "DCA= " << dca << " strawRadius= " << strawRadius << "\n";
						//Find the truth ID, and LR hit for a straw
						float residualTruth = MisDetector.residualTruth;
						int misID =  MisDetector.strawID; // ID of the hit straw [to identify the correct straw in the Fit function]
						int misLRSign = MisDetector.LRSign;
						MC.dcaUnsmeared.push_back(MisDetector.dcaUnsmeared);
						MC.strawID.push_back(misID);
						MC.LR.push_back(misLRSign);
						//if (1 == 1) {cout << "DCA is= " << dca << " for straw ID= " << misID << " was hit from " << misLRSign << endl;}
						// RECONSTRUCTED PARAMETERS
						//Reconstructing the hit as seen from the ideal detector [shift circle x coordinate by misalignment]

						ReconData ReconDetector = Tracker::HitRecon(misID, dca, mod_lyr_strawIdealPositionX[i_module][i_view][i_layer], mod_lyr_strawIdealPositionZ[i_module][i_view][i_layer]);
						float zCircle = ReconDetector.z;
						float xCircle = ReconDetector.x;
						float radCircle = ReconDetector.dcaRecon;
						xRecon.push_back(xCircle); // vector to store x coordinates of circles as seen from the ideal detector
						zRecon.push_back(zCircle); // vector to store z coordinates of circles as seen from the ideal detector
						radRecon.push_back(radCircle); // vector to store radius of circles as seen from the ideal detector

						//Module number [for labelling] - after (if) passing the rejection.
						// Millepede accepts only positive non-zero integers as labels
						// M0=10, M1=20... // x=1, z=2, 𝛉=3
						int label1 = (i_module + 1) * 10 + 1; // x
						int label2 = (i_module + 1) * 10 + 2; // z
						int label3 = (i_module + 1) * 10 + 3; // 𝛉
						MC.label1.push_back(label1);
						MC.label2.push_back(label2);
						MC.label3.push_back(label3);

						//Sanity Plots: Hits
						MC.dca.push_back(dca); // DCA [sanity]
						MC.residualsTruth.push_back(residualTruth); // Truth residual
						MC.Module_i.push_back(i_module);
						MC.View_i.push_back(i_view);
						MC.Layer_i.push_back(i_layer);
						MC.Straw_i.push_back(misID);

						//push back the tracker centre for the hit
						MC.zCentreModule.push_back(idealCentres.zCentres[i_module]);
						MC.xCentreModule.push_back(idealCentres.xCentres[i_module]);

						layerCount++;
						//std::cout << "layerCount2= " << layerCount << "\n";
						z_counter++;  // incrementing distance of planes
					} // missed hit check
					//Rejection of hits due to geometry (i.e. missed hits)
					//No signal in a straw = no signal in the layer
					else { //[between (0,0.25]
						incRejectedHitsDCA();
						//cout << "Hit Rejected: outside of straw layer with dca =" << dca << endl;
						z_counter++;  // incrementing distance of planes
						//missedHit = true;
					}
					//std::cout << "layerCount3= " << layerCount << "\n";
				}// DCA cut if
			}//end of Layers loop
		}// end of View loop
	}// end of looping over modules
	//std::cout << "layerCount4= " << layerCount << "\n";
	// Now for the whole track, if not cut
	if (cutTriggered == false) {
		MC.totalLayerHits = layerCount;
		//std::cout << "layerCount5= " << layerCount << "\n";
		if (debugBool && StrongDebugBool) {cout << "Calculating residuals..." << endl;}
		//This happens once per MC function call [as we now accumulated x coordinates of "ideal" points for all hits
		// and need to do a simultaneous fit once - to return #hits worth of residuals]
		ResidualData resData = GetResiduals(zRecon, xRecon, radRecon, layerCount, plot_fit, debugBool, useTruthLR, MC.LR);

		MC.residualsRecon = resData.residualsRecon;
		MC.slopeRecon = resData.slopeRecon;
		MC.interceptRecon = resData.interceptRecon;
		MC.xStraw = xRecon;
		MC.zStraw = zRecon;
		MC.pValue = resData.pValue;
		MC.chi2Circle = resData.chi2Circle;
		MC.driftRad = radRecon; // vector
		//cout << "MC.residualsRecon.size()= " << MC.residualsRecon.size() << "\n";

		//plotting files
		if (debugBool) {
			plot_gen << x0 << " " << x1 << " " << beamStart << " " << beamStop << " " << endl;
			int i_counter = 0;
			for (int i_module = 0; i_module < moduleN; i_module++) {
				for (int i_view = 0; i_view < viewN; i_view++) {
					for (int i_layer = 0; i_layer < layerN; i_layer++) { //per layer
						plot_hits_gen << mod_lyr_strawMisPositionX[i_module][i_view][i_layer][MC.strawID[i_counter]] << " " << mod_lyr_strawMisPositionZ[i_module][i_view][i_layer][MC.strawID[i_counter]] << " "  << MC.dca[i_counter] << endl;
						plot_hits_fit << xRecon[i_counter] << " " << zRecon[i_counter] << " "  << radRecon[i_counter] << endl;
						i_counter++;
					} // layers
				} // views
			} // modules
		} // plotting

	} // dca cut

	return MC; // Return data from simulated track
} // end of MC

// MC misalignment of detectors and setting assumed geometry
void Tracker::setGeometry(ofstream& debug_geom, ofstream& debug_mis, ofstream& pede_mis, ofstream& plot_centres, bool debugBool, ofstream& metric) {

	//Geometry of detector arrangement (Ideal Geometry)
	float dZ = startingZDistanceStraw0; // the increment in z for consecutive layers
	int layer_n = 0;  // layer label
	float dX = startingXDistanceStraw0; // the increment in x for straws in a layers
	for (int i_module = 0; i_module < moduleN; i_module++) {
		mod_lyr_strawIdealPositionX.push_back(vector<vector<vector<float> > >()); //initialize the first index with a 2D vector
		mod_lyr_strawIdealPositionZ.push_back(vector<vector<vector<float> > >()); //initialize the first index with a 2D vector
		for (int i_view = 0; i_view < viewN; i_view++) {
			mod_lyr_strawIdealPositionX[i_module].push_back(vector<vector<float> >());
			mod_lyr_strawIdealPositionZ[i_module].push_back(vector<vector<float> >());
			for (int i_layer = 0; i_layer < layerN; i_layer++) {
				mod_lyr_strawIdealPositionX[i_module][i_view].push_back(vector<float> ());
				mod_lyr_strawIdealPositionZ[i_module][i_view].push_back(vector<float> ());
				layer.push_back(layer_n);  // layer label array [starting from 0th layer]
				for (int i_straw = 0; i_straw < strawN; i_straw++) {
					mod_lyr_strawIdealPositionZ[i_module][i_view][i_layer].push_back(dZ + offsetZ[i_module]); // vector will contain all z coordinates of layers
					mod_lyr_strawIdealPositionX[i_module][i_view][i_layer].push_back(dX + offsetX[i_module]);
					dX = dX - strawSpacing; //while we are in the same layer: increment straw spacing in x
				} //end of Straws loop
				if (i_view == 0) { dX = startingXDistanceStraw0 - layerDisplacement; } //set displacement in x for the next layer in the view
				if (i_view == 1) { dX = startingXDistanceStraw0; } //set displacement in x for the next layer in the view
				layer_n++;
				if (i_layer == 0) { dZ += layerSpacing; } // increment spacing between layers in a view once only
			} // layers
			if (i_view == 0) { dZ += viewSpacing; } // increment spacing between views in a module once only
		} // views
		dZ += moduleSpacing;
	} // modules

	//Now get rotational centres for the assumed geometry, add offsets if necessary
	RotationCentres idealCentres = Tracker::getCentre(mod_lyr_strawIdealPositionX, mod_lyr_strawIdealPositionZ, plot_centres);
	dZ = startingZDistanceStraw0; // the increment in z for consecutive layers
	dX = startingXDistanceStraw0; // the increment in x for straws in a layers
	for (int i_module = 0; i_module < moduleN; i_module++) {
		for (int i_view = 0; i_view < viewN; i_view++) {
			for (int i_layer = 0; i_layer < layerN; i_layer++) {
				for (int i_straw = 0; i_straw < strawN; i_straw++) {
					float local_z = dZ - idealCentres.zCentres[i_module]; // translating to local module coordinates
					float local_x = dX - idealCentres.xCentres[i_module];
					//2D rotation matrix with CC convention from (z=0 horizontal/beam axis)
					float rotated_z = local_z * cos(offsetTheta[i_module]) - local_x * sin (offsetTheta[i_module]);
					float rotated_x = local_z * sin(offsetTheta[i_module]) + local_x * cos (offsetTheta[i_module]);
					float global_z = rotated_z + idealCentres.zCentres[i_module];
					float global_x = rotated_x + idealCentres.xCentres[i_module];
					mod_lyr_strawIdealPositionZ[i_module][i_view][i_layer][i_straw] = global_z; // vector will contain all z coordinates of layers
					mod_lyr_strawIdealPositionX[i_module][i_view][i_layer][i_straw] = global_x;
					dX = dX - strawSpacing; //while we are in the same layer: increment straw spacing in x
				} //end of Straws loop
				if (i_view == 0) { dX = startingXDistanceStraw0 - layerDisplacement; } //set displacement in x for the next layer in the view
				if (i_view == 1) { dX = startingXDistanceStraw0; } //set displacement in x for the next layer in the view
				layer_n++;
				if (i_layer == 0) { dZ += layerSpacing; } // increment spacing between layers in a view once only
			} // layers
			if (i_view == 0) { dZ += viewSpacing; } // increment spacing between views in a module once only
		} // views
		dZ += moduleSpacing;
	} // modules

	// Print out:
	if (debugBool) {
		stringstream gm1;
		gm1 << endl << "Ideal Detector Position:";
		Logger::Instance()->write(Logger::WARNING, gm1.str());
		int Zcounter = 0;
		for (int i_module = 0; i_module < moduleN; i_module++) {
			for (int i_view = 0; i_view < viewN; i_view++) {
				for (int i_layer = 0; i_layer < layerN; i_layer++) { //per module
					cout << "IDEAL M" << i_module + 1 << noshowpos << UVmapping[i_view][i_layer] << " X : ";
					for (int i_straw = 0; i_straw < strawN; i_straw++) {
						cout << showpos << mod_lyr_strawIdealPositionX[i_module][i_view][i_layer][i_straw] << " ";
						debug_geom << mod_lyr_strawIdealPositionX[i_module][i_view][i_layer][i_straw] << " ";
					} // straws
					cout << endl;
					debug_geom << endl;
					cout << "IDEAL M" << i_module + 1 << noshowpos << UVmapping[i_view][i_layer] << " Z : ";
					for (int i_straw = 0; i_straw < strawN; i_straw++) {
						cout << showpos << mod_lyr_strawIdealPositionZ[i_module][i_view][i_layer][i_straw] << " ";
						debug_geom << mod_lyr_strawIdealPositionZ[i_module][i_view][i_layer][i_straw] << " ";
					} // straws
					cout << endl;
					debug_geom << endl;
				} // end of Layers
			}// end of View loop
			cout << endl;
		}//Modules
	} // debug


	//Now misaligning detectors
	dZ = startingZDistanceStraw0; // the increment in z for consecutive layers
	layer_n = 0;  // layer label
	//metric << "MTheta: ";
	for (int i_module = 0; i_module < moduleN; i_module++) {
		float dX = startingXDistanceStraw0; // starting on the x-axis (z, 0+disp)
		mod_lyr_strawMisPositionX.push_back(vector<vector<vector<float> > >()); //initialize the first index with a 2D vector
		mod_lyr_strawMisPositionZ.push_back(vector<vector<vector<float> > >()); //initialize the first index with a 2D vector
		for (int i_view = 0; i_view < viewN; i_view++) {
			mod_lyr_strawMisPositionX[i_module].push_back(vector<vector<float> >()); //initialize the first index with a 2D vector
			mod_lyr_strawMisPositionZ[i_module].push_back(vector<vector<float> >()); //initialize the first index with a 2D vector
			for (int i_layer = 0; i_layer < layerN; i_layer++) {
				mod_lyr_strawMisPositionX[i_module][i_view].push_back(vector<float> ()); //initialize the first index with a 2D vector
				mod_lyr_strawMisPositionZ[i_module][i_view].push_back(vector<float> ()); //initialize the first index with a 2D vector
				layer.push_back(layer_n);  // layer label array [starting from 0th layer]
				for (int i_straw = 0; i_straw < strawN; i_straw++) {
					float local_z = dZ - idealCentres.zCentres[i_module]; // translating to local module coordinates
					float local_x = dX - idealCentres.xCentres[i_module];
					//2D rotation matrix with CC convention from (z=0 horizontal/beam axis)
					float rotated_z = local_z * cos(dispTheta[i_module]) - local_x * sin (dispTheta[i_module]);
					float rotated_x = local_z * sin(dispTheta[i_module]) + local_x * cos (dispTheta[i_module]);
					float global_z = rotated_z + idealCentres.zCentres[i_module];
					float global_x = rotated_x + idealCentres.xCentres[i_module];
					mod_lyr_strawMisPositionZ[i_module][i_view][i_layer].push_back(global_z + dispZ[i_module]); // vector will contain all z coordinates of layers
					mod_lyr_strawMisPositionX[i_module][i_view][i_layer].push_back(global_x + dispX[i_module]);
					dX =  dX - strawSpacing; //while we are in the same layer: increment straw spacing in x
				} //end of Straws loop
				if (i_view == 0) { dX = startingXDistanceStraw0 - layerDisplacement; } //set displacement in x for the next layer in the view
				if (i_view == 1) { dX = startingXDistanceStraw0; } //set displacement in x for the next layer in the view
				layer_n++;
				if (i_layer == 0) { dZ += layerSpacing; } // increment spacing between layers in a view once only
			}//end of Layers loop
			if (i_view == 0) { dZ += viewSpacing; } // increment spacing between views in a module once only
		}// end of View loop
		dZ += moduleSpacing;
	}//Modules

	//Now get rotational centres for the misaligned geometry (sanity check)
	RotationCentres misCentres = Tracker::getCentre(mod_lyr_strawMisPositionX, mod_lyr_strawMisPositionZ, plot_centres);

	// Print out:
	if (debugBool) {
		stringstream gm2;
		gm2 << "Misaligned Detector Position:";
		Logger::Instance()->write(Logger::WARNING, gm2.str());
		int Zcounter = 0;
		for (int i_module = 0; i_module < moduleN; i_module++) {
			for (int i_view = 0; i_view < viewN; i_view++) {
				for (int i_layer = 0; i_layer < layerN; i_layer++) { //per module
					cout << "MIS M" << noshowpos << i_module + 1 << UVmapping[i_view][i_layer] << " X : ";
					for (int i_straw = 0; i_straw < strawN; i_straw++) {
						cout << showpos  << mod_lyr_strawMisPositionX[i_module][i_view][i_layer][i_straw] << " ";
						debug_mis << mod_lyr_strawMisPositionX[i_module][i_view][i_layer][i_straw] << " ";
					} //end of Straws loop
					cout << endl;
					debug_mis << endl;
					cout << "MIS M" << noshowpos << i_module + 1 << UVmapping[i_view][i_layer] << " Z : ";
					for (int i_straw = 0; i_straw < strawN; i_straw++) {
						cout << showpos  << mod_lyr_strawMisPositionZ[i_module][i_view][i_layer][i_straw] << " ";
						debug_mis << mod_lyr_strawMisPositionZ[i_module][i_view][i_layer][i_straw] << " ";
					} //end of Straws loop
					cout << endl;
					debug_mis << endl;
				} // end of Layers
			}// end of View loop
			cout << endl;
		}//Modules
	}

	// Print out misalignment per module
	cout << "Misalignment(M)" << endl;
	for (int i_module = 0; i_module < moduleN; i_module++) {
		cout << "M" << noshowpos << i_module + 1  << " x :: "  << showpos << dispX[i_module] << " cm. " << "z :: "  << showpos << dispZ[i_module] << " cm. " << "𝛉 :: " << showpos << dispTheta[i_module] << " rad. " << endl;
		//if (i_module == 1 || i_module == 2) { // TODO move to pre-sigma function
		pede_mis << dispX[i_module] << " ";
		//}
	} // modules
	cout << noshowpos;

}//end of misalign and set geometry