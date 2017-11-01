/*
This source file contains definitions of various functions in the method class.
*/

#include "AlignTracker_methods.h"

using namespace std;

// Set up empty pointer for instance of the class.
Tracker* Tracker::s_instance = NULL;

/**
Constructor for tracker class.
 */
Tracker::Tracker() {
	// non-static variables definition here

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

void Tracker::write_steering_file(ofstream& steering_file) {
	if (steering_file.is_open()) {
		steering_file <<  "* g-2 Tracker Alignment: PEDE Steering File" << endl
		              << " "  << endl
		              << "Tracker_con.txt   ! constraints text file " << endl
		              << "Cfiles ! following bin files are Cfiles" << endl
		              << "Tracker_data.bin   ! binary data file" << endl
		              << "method inversion 1 0.01" << endl
		              << "printrecord  -1 -1      ! debug printout for bad data records" << endl
		              << "printrecord 1 -1 ! produces mpdebug.txt" << endl    //
		              << " "  << endl
		              << "end ! optional for end-of-data" << endl;
	} // steering file open
} // steering function
/**
   Write a constraint file to the supplied file-stream.
    @param constraint_file Reference to ofstream to write constraint file to.
 */
// XXX constraints are ignored with HIP method
void Tracker::write_constraint_file(ofstream& constraint_file, ofstream& debug_con, bool debugBool) {
	// Check constraints file is open, then write.
	if (constraint_file.is_open()) {
		//Evaluation of constraints
		float one = 1.0;
		//Fixing module 0 and the last module
		for (int i_module = 0; i_module < moduleN; i_module++) {
			if (i_module == 0 || i_module == moduleN - 1) {
				constraint_file << "Constraint 0.0" << endl;
				int labelt = i_module + 1; // Millepede accepts +ive labels only
				constraint_file << labelt << " " << fixed << setprecision(5) << one << endl;
			} // end of fixed modules
		} // end of detectors loop
	} // constrain file open
} // end of writing cons file

/**
   DCA for a point to a line in 2D space
*/
float Tracker::pointToLineDCA(float z_straw, float x_straw, float x_slope, float x_intercept) {

	//http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html

	//converting from x=mz+c -> mz+x+c=0
	float a = -x_slope;
	float b = 1.0;
	float c = -x_intercept;
	float x0 = z_straw;
	float y0 = x_straw;

	float dca = abs( a * x0 + b * y0 + c ) / sqrt( a * a + b * b  ) ;

	//cout << "m= " << a << " c= " << c << " z= " << x0 << " -x= " << y0 << " dca= " << dca <<  endl;

	return dca;
}


/** Uses DCA function to find the shortest dca between straws (i.e. which straw was hit in that layer)
    @Inputs (see DCA method) + vector of straws' x coordinates in a layer, and a return type: "dca_hit"  or "x_line"
    @return hit_distance - dca of the straw that was hit, or x coordinate of the line on the same z as a straw
*/
DCAData Tracker::DCAHit(vector<float> xLayer, float zStraw, float xTrack, float xSlpoe, float xIntercept, bool debugBool) {

	DCAData dca_data;

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
		hitDistance = pointToLineDCA(zStraw, xLayer[0], xSlpoe, xIntercept);
		LR = +1; // assign truth L
		index = 0;
		if (debugBool && StrongDebugBool) {
			cout << "Track in line at " << xTrack << " The first straw is closest at " << upper <<  "; with DCA " << hitDistance << endl;
		}
	}
	// Another scenario, is that the hit is smaller than the x coordinate of the last straw
	else if (xTrack < xLayer[lastID]) {
		// hit distance is the dca from the R of the last straw
		hitDistance = pointToLineDCA(zStraw, xLayer[lastID], xSlpoe, xIntercept);
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
		for (int i_counter = 0; i_counter < xLayer.size(); i_counter++) {
			if (xLayer[i_counter] < xTrack) {
				lower = xLayer[i_counter];
				upper = xLayer[i_counter - 1];
				if (debugBool && StrongDebugBool) cout << "lower= " << lower << " upper " << upper << endl;
				//as soon as we find a single straw that is at lower x than the hit,
				// our search is over
				goto jmp;
			}
		}
jmp:

		float hit_distance_low = pointToLineDCA(zStraw, lower, xSlpoe, xIntercept);
		float hit_distance_up = pointToLineDCA(zStraw, upper, xSlpoe, xIntercept);

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
			float random = Tracker::generate_uniform();
			if (random < 0.5) {
				hitDistance = hit_distance_low;
			}
			if (random > 0.5) {
				hitDistance = hit_distance_up;
			}
			if (1 == 1) {cout << "Ambiguity which straw registered hit" << endl;}
			cout << "Ambiguity which straw registered hit" << endl;
			incAmbiguityHit();
		}
		if (debugBool && StrongDebugBool) {
			cout <<  "Track in Line " << xTrack << "Two straws closest to that point are " << lower << ", and " << upper <<  "; with DCAs " << hit_distance_low <<  " and " << hit_distance_up << ", respectively." << endl;
		}
	} // end of iterator to find straws between hits

	//Now smear the DCA data
	dca_data.dcaUnsmeared = hitDistance;
	float hitDistanceSmeared = hitDistance + Tracker::instance()->getResolution() * Tracker::generate_gaus();
	//Apply a 500 um cut on the WHOLE track

	if (abs(hitDistanceSmeared) < trackCut && trackCutBool) {
		cutTriggered = true; // will be checked in the MC_Launch on DCA return
	}

	dca_data.dca = hitDistanceSmeared;
	float residualTruth = hitDistanceSmeared - hitDistance;
	dca_data.residualTruth = residualTruth;
	dca_data.LRSign = LR;
	dca_data.strawID = index;

	if (debugBool && StrongDebugBool) {
		cout << "Selected DCA as the correct hit distance is " << hitDistance << ". Straw ID: " << dca_data.strawID;
		if (LR  < 0) {cout << ". The straw was hit from the right" << endl;}
		if (LR > 0) {cout << ". The straw was hit from the left" << endl;}
	}//debug
	return dca_data; // as dca and id of the closest straw
}


// Adds dca to the ideal geometry
ReconData Tracker::HitRecon(int det_ID, float det_dca, vector<float> xLayer, float z_distance) {

	ReconData recon_data;
	bool StrongDebugBool = false; //quick hack XXX for even more debug output
	float x_IdealStraw = xLayer[det_ID];
	float recon_dca = det_dca; //adds dca to the assumed straw x coordinate
	if (StrongDebugBool) { cout << "recon_dca= " << recon_dca << "; x_IdealStraw= " << x_IdealStraw <<  " det_dca= " << det_dca << endl;}
	recon_data.z = z_distance;
	recon_data.x = x_IdealStraw;
	recon_data.dcaRecon = recon_dca;
	return recon_data;
}

// Function to return residuals to the fitted line (due to dca point scatter + resolution of the detector)
// @ Inputs: ideal points (x,z)

//Truth can be used as an optional input

// z,x position of straws, radius of drift circle, # circles, plotting file for input, verbosity flag, truth flag
ResidualData Tracker::GetResiduals(vector<float> zRecon, vector<float> xRecon, vector<float> radRecon, int dataSize, ofstream& plot_fit, bool debugBool, bool useTruthLR, vector<int> LR_truth) {

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
			if (LR_truth[i_hit] == -1) r = -r;
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
		double stepSize = 0.1; // XXX

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
				if (LR_truth[i_hit] == -1) r = -r;

				// Calculate distance of track from wire and use it for Chi2 calculation
				double d = (gradients.at(grad) * z + intercepts.at(grad) - u) / sqrt(gradients.at(grad) * gradients.at(grad) + 1);
				//double d = Tracker::pointToLineDCA(z, u, gradient, intercept);
				chi2Val += pow(d - r, 2) / err2;

				if (debugBool && StrongDebugBool) {
					cout << "grad= " << grad << " gradients.size()= " << gradients.size() << " chi2Val " << chi2Val
					     << " z= " << z << " u= "  << u <<  " r= " << r << " LR_truth= " << LR_truth[i_hit] << " err2= " << err2
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
		resData.p_value = pVal;
		resData.chi2_circle = chi2ValMin;

		if (debugBool && StrongDebugBool) {cout << "pVal=" << pVal << " chi2ValMin= " << chi2ValMin << endl << endl;}
		if (pVal > pValCut) {
			// We'll want to store left/right hits so set these
			for (int i_hit = 0; i_hit < nHits; i_hit++) {
				//XXX Now that we know the slope and the gradient of the best fit line through drift circles,
				// we can calculate the residual for each drift circle, which is
				// "DCA from that line to the straw centre" - "Radius of the drift circle"
				float Res = Tracker::pointToLineDCA(zRecon[i_hit],  xRecon[i_hit], gradient, intercept) - radRecon[i_hit];
				resData.residuals.push_back(Res);   // residual between the (centre of the straw and the fitted line [pointToLineDCA]) and radius of the fit circle;

				float xTrack_recon = gradient * zRecon[i_hit] + intercept; // Recon track position in-line with the layer [from line x=ym+c]
				resData.x_track_recon.push_back(xTrack_recon);  // Recon in-line (gen.) track position

			} // hits

			// Passing recon track parameters to MC
			resData.slope_recon = gradient;     // slope of the best fit line
			resData.intercept_recon = intercept; // intercept
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
MCData Tracker::MC_launch(float scatterError, ofstream& debug_calc, ofstream& debug_off, ofstream& debug_mc, ofstream& plot_fit, ofstream& plot_gen, ofstream& plot_hits_gen, ofstream& plot_hits_fit, bool debugBool) {

	bool StrongDebugBool = false; // HACK XXX set by hand for now
	// Set up new container for track data, with hit count set to zero
	MCData MC;
	MC.hit_count = 0; //count only counts through detector
	cutTriggered = false; //set trigger to OFF until DCA is less than 500 um, then kill this track

	//clearing containers for next track
	xRecon.clear();  // ideal straw x
	zRecon.clear();  // ideal straw z
	radRecon.clear();

	//Track parameters for rand-generated line MC [start and end positions outside of detectors]
	// redefining the track as x=ym+c
	float x0 = beamPositionLength * Tracker::generate_uniform() - 1.0; //uniform vertex
	float signXSlope;
	if (Tracker::generate_uniform() >= 0.5) {
		signXSlope = 1.0;
	}
	else {
		signXSlope = -1.0;
	}
	float xSlope = (Tracker::generate_uniform() * signXSlope) * (0.5 * beamPositionLength / beamStop);
	float xIntercept = x0; // by definition
	float x1 = xSlope * beamStop + xIntercept;

	float xTrack; //true track position x=zm+c

	MC.x0 = x0;
	MC.x1 = x1;
	MC.slope_truth = xSlope;
	MC.intercept_truth = xIntercept;

	//The main loop is for modules [they produce label of Global Parameters]:
	// Then looping over layers and views
	// within layers we will have loops over straws in that layer [to find out which straw was hit and dca]
	int z_counter = 0; // distance in z is incremented from 0 to TotalLayerN
	if (debugBool && StrongDebugBool) {cout << "Calculating hits/dca:" << endl;}
	//loops over modules, views, layers
	for (int i_module = 0; i_module < moduleN; i_module++) {
		for (int i_view = 0; i_view < viewN; i_view++) {
			for (int i_layer = 0; i_layer < layerN; i_layer++) { //per layer
				missedHit = false; // set flag to false until triggered
				// TRUTH PARAMETERS
				//The registered hit position on the misaligned detector is smeared by its resolution
				// zTrack is just the distance[z_counter] - in-line with the layer
				xTrack = xSlope * distance[z_counter] + xIntercept; // true track position in-line with the layer [from line x=ym+c]
				MC.x_track_true.push_back(xTrack);  // True in-line (gen.) track position

				// DCA will return the pointToLineDCA (smeared by resolution) and strawID [+truth parameters: LR sign etc. xHit and zHit]
				DCAData MisDetector = Tracker::DCAHit(mod_lyr_strawMisPosition[i_module][i_view][i_layer], distance[z_counter], xTrack, xSlope, xIntercept, debugBool); // position on the detector [from dca]
				//Check for trigger on DCA < 500 um, and only continue calculations not triggered
				if (cutTriggered) {MC.cut = true;} // set cut to true for the main, and return MC.cut as true only.
				if (cutTriggered == false) {
					float dca =  MisDetector.dca;  //dca (smeared) of the hit straw
					//Rejection of hits due to geometry (i.e. missed hits)
					//No signal in a straw = no signal in the layer
					if (dca > strawRadius && hitCut == true) { //[between (0,0.25]
						MC.hit_bool.push_back(0);
						incRejectedHitsDCA();
						if (debugBool) {cout << "Hit Rejected: outside of straw layer with dca =" << dca << endl;}
						z_counter++;  // incrementing distance of planes
						missedHit = true;
					}
					// proceed only if hit is valid
					if (missedHit == false) {
						//Find the truth ID, and LR hit for a straw
						float residualTruth = MisDetector.residualTruth;
						int mis_ID =  MisDetector.strawID; // ID of the hit straw [to identify the correct straw in the Fit function]
						int mis_LRSign = MisDetector.LRSign;
						MC.dca_unsmeared.push_back(MisDetector.dcaUnsmeared);
						MC.strawID.push_back(mis_ID);
						MC.LR.push_back(mis_LRSign);
						//Recording hit information
						MC.hit_list.push_back(z_counter); //push back the absolute layer # into the hit list
						MC.hit_bool.push_back(1);  // 1 = hit
						stringstream absolute_hit;
						//Making absolute strawID for hit for plotting:
						absolute_hit << mis_ID;
						MC.absolute_straw_hit.push_back(absolute_hit.str().c_str());
						if (debugBool && StrongDebugBool) {cout << "DCA is= " << dca << " for straw ID= " << mis_ID << " was hit from " << mis_LRSign << endl;}
						// RECONSTRUCTED PARAMETERS
						//Reconstructing the hit as seen from the ideal detector [shift circle x coordinate by misalignment]

						ReconData ReconDetector = Tracker::HitRecon(mis_ID, dca, mod_lyr_strawIdealPosition[i_module][i_view][i_layer], distance[z_counter]);
						float zCircle = ReconDetector.z;
						float xCircle = ReconDetector.x;
						float radCircle = ReconDetector.dcaRecon;
						xRecon.push_back(xCircle); // vector to store x coordinates of circles as seen from the ideal detector
						zRecon.push_back(zCircle); // vector to store z coordinates of circles as seen from the ideal detector
						radRecon.push_back(radCircle); // vector to store radius of circles as seen from the ideal detector

						//Module number [for labelling] - after (if) passing the rejection.
						// Millepede accepts only positive non-zero integers as labels
						ostringstream oss; oss << i_module + 1;
						istringstream iss(oss.str()); int label_int; iss >> label_int;
						MC.i_hits.push_back(label_int); // vector of modules that were actually hit [after passing rejection test: for MP2 labelling]

						//Z-coordinate of hits [it was always known from geometry - no z-misalignment for now...]
						MC.z_hits.push_back(distance[z_counter]);
						MC.hit_sigmas.push_back(Tracker::instance()->getResolution());

						//Sanity Plots: Hits
						MC.dca.push_back(dca); // DCA [sanity]
						MC.residuals_gen.push_back(residualTruth); // Truth residual
						MC.Module_i.push_back(i_module);
						MC.View_i.push_back(i_view);
						MC.Layer_i.push_back(i_layer);
						MC.Straw_i.push_back(mis_ID);

						MC.hit_count++;
						z_counter++;  // incrementing distance of planes
					} // missed hit check

				}// DCA cut if
			}//end of Layers loop
		}// end of View loop
	}// end of looping over modules

	if (cutTriggered == false) {

		if (debugBool && StrongDebugBool) {cout << "Calculating residuals..." << endl;}
		//This happens once per MC function call [as we now accumulated x coordinates of "ideal" points for all hits
		// and need to do a simultaneous fit once - to return #hits worth of residuals]
		ResidualData res_Data = GetResiduals(zRecon, xRecon, radRecon, MC.hit_count, plot_fit, debugBool, useTruthLR, MC.LR);

		MC.residuals = res_Data.residuals;
		MC.slope_recon = res_Data.slope_recon;
		MC.intercept_recon = res_Data.intercept_recon;
		MC.x_straw = xRecon;
		MC.z_straw = zRecon;
		MC.x_track_recon = res_Data.x_track_recon;
		MC.p_value = res_Data.p_value;
		MC.chi2_circle = res_Data.chi2_circle;
		MC.driftRad = radRecon; // vector

		if (debugBool) {
			plot_gen << x0 << " " << x1 << " " << beamStart << " " << beamStop << " " << endl;
			int i_counter = 0;
			for (int i_module = 0; i_module < moduleN; i_module++) {
				for (int i_view = 0; i_view < viewN; i_view++) {
					for (int i_layer = 0; i_layer < layerN; i_layer++) { //per layer
						plot_hits_gen << mod_lyr_strawMisPosition[i_module][i_view][i_layer][MC.strawID[i_counter]] << " " << distance[i_counter] << " "  << MC.dca[i_counter] << endl;
						plot_hits_fit << xRecon[i_counter] << " " << zRecon[i_counter] << " "  << radRecon[i_counter] << endl;
						i_counter++;
					}
				}
			}
		}

	} // dca cut
	return MC; // Return data from simulated track

} // end of MC

//Geometry of detector arrangement (Ideal Geometry)
void Tracker::setGeometry(ofstream& debug_geom,  bool debugBool) {
	// float dZ=startingZDistanceStraw0; // to add distance in z

	float dZ = startingZDistanceStraw0; // the increment in z for consecutive layers
	int layer_n = 0;  // layer label
	for (int i_module = 0; i_module < moduleN; i_module++) {
		for (int i_view = 0; i_view < viewN; i_view++) {
			for (int i_layer = 0; i_layer < layerN; i_layer++) {
				distance.push_back(dZ); // vector will contain all z coordinates of layers
				layer.push_back(layer_n);  // layer label array [starting from 0th layer]
				//resolutionLayer.push_back(Tracker::instance()->getResolution()); //resolution in each layer
				projectionX.push_back(float(1.0));  // x projection of hits in each layer
				layer_n++;
				if (i_layer == 0) { dZ += layerSpacing; } // increment spacing between layers in a view once only
			} // layers
			if (i_view == 0) { dZ += viewSpacing; } // increment spacing between views in a module once only
		} // views
		dZ += moduleSpacing;
	} // modules

	// Geometry of detector arrangement in X
	float dX = startingXDistanceStraw0; // starting on the x-axis (z, 0)
	for (int i_module = 0; i_module < moduleN; i_module++) {
		mod_lyr_strawIdealPosition.push_back(vector<vector<vector<float> > >()); //initialize the first index with a 2D vector
		for (int i_view = 0; i_view < viewN; i_view++) {
			mod_lyr_strawIdealPosition[i_module].push_back(vector<vector<float> >()); //initialize the first index with a 2D vector
			for (int i_layer = 0; i_layer < layerN; i_layer++) {
				mod_lyr_strawIdealPosition[i_module][i_view].push_back(vector<float> ()); //initialize the first index with a 2D vector
				for (int i_straw = 0; i_straw < strawN; i_straw++) {
					mod_lyr_strawIdealPosition[i_module][i_view][i_layer].push_back(dX);
					dX = dX - strawSpacing; //while we are in the same layer: increment straw spacing in x
				} //end of Straws loop
				if (i_view == 0) { dX = startingXDistanceStraw0 - layerDisplacement; } //set displacement in x for the next layer in the view
				if (i_view == 1) { dX = startingXDistanceStraw0; } //set displacement in x for the next layer in the view
			}//end of Layers loop
		}// end of View loop
	}//Modules

	// Print out:
	if (debugBool) {
		stringstream gm1;
		gm1 << endl << "Ideal Detector Position:";
		Logger::Instance()->write(Logger::WARNING, gm1.str());
		int Zcounter = 0;
		for (int i_module = 0; i_module < moduleN; i_module++) {
			for (int i_view = 0; i_view < viewN; i_view++) {
				for (int i_layer = 0; i_layer < layerN; i_layer++) { //per module
					cout << "IDEAL Mod " << i_module  << " " << UVmapping[i_view][i_layer] << " X : ";
					for (int i_straw = 0; i_straw < strawN; i_straw++) {
						cout << mod_lyr_strawIdealPosition[i_module][i_view][i_layer][i_straw] << " ";
						debug_geom << mod_lyr_strawIdealPosition[i_module][i_view][i_layer][i_straw] << " ";
					} // straws
					cout << "  | Z= " << distance[Zcounter] << " [cm]" << endl;  // TODO align the cout better
					debug_geom << distance[Zcounter] << endl;
					Zcounter++;
				} // end of Layers

			}// end of View loop
			cout << endl;
		}//Modules
	}

} // end of geom

// MC misalignment of detectors
void Tracker::misalign(ofstream& debug_mis, ofstream& pede_mis, bool debugBool) {

	//Now misaligning detectors in x
	float misDispX; // effective misalignment
	for (int i_module = 0; i_module < moduleN; i_module++) {
		misDispX = dispX[i_module];

		float dX = startingXDistanceStraw0 + (misDispX); // starting on the x-axis (z, 0+disp)
		sdevX.push_back(misDispX);  // vector to store the actual of misalignment
		mod_lyr_strawMisPosition.push_back(vector<vector<vector<float> > >()); //initialize the first index with a 2D vector
		for (int i_view = 0; i_view < viewN; i_view++) {
			mod_lyr_strawMisPosition[i_module].push_back(vector<vector<float> >()); //initialize the first index with a 2D vector
			for (int i_layer = 0; i_layer < layerN; i_layer++) {
				mod_lyr_strawMisPosition[i_module][i_view].push_back(vector<float> ()); //initialize the first index with a 2D vector
				for (int i_straw = 0; i_straw < strawN; i_straw++) {
					mod_lyr_strawMisPosition[i_module][i_view][i_layer].push_back(dX);
					dX =  dX - strawSpacing; //while we are in the same layer: increment straw spacing in x
				} //end of Straws loop
				if (i_view == 0) { dX = startingXDistanceStraw0 - layerDisplacement + (misDispX); } //set displacement in x for the next layer in the view
				if (i_view == 1) { dX = startingXDistanceStraw0 + (misDispX); } //set displacement in x for the next layer in the view
			}//end of Layers loop
		}// end of View loop

	}//Modules

	// Print out:
	if (debugBool) {
		stringstream gm2;
		gm2 << "Misaligned Detector Position:";
		Logger::Instance()->write(Logger::WARNING, gm2.str());
		int Zcounter = 0;
		for (int i_module = 0; i_module < moduleN; i_module++) {
			for (int i_view = 0; i_view < viewN; i_view++) {
				for (int i_layer = 0; i_layer < layerN; i_layer++) { //per module
					cout << "MIS Mod " << i_module  << " " << UVmapping[i_view][i_layer] << " X : ";
					for (int i_straw = 0; i_straw < strawN; i_straw++) {
						cout << mod_lyr_strawMisPosition[i_module][i_view][i_layer][i_straw] << " ";
						debug_mis << mod_lyr_strawMisPosition[i_module][i_view][i_layer][i_straw] << " ";

					} //end of Straws loop
					cout << " | Z= " << distance[Zcounter] << " [cm]" << endl;
					debug_mis << distance[Zcounter] << endl;
					Zcounter++;
				} // end of Layers

			}// end of View loop
			cout << endl;
		}//Modules
	}

// Estimating misalignment parameters from geometry and assumed constants:

	//--- For Modules only [MC misalignment set per module -> transfer to each layer] -- //

	//Overall Misalignment [w.r.t to other modules]
	float sum_of_elems = 0.0;
	for (vector<float>::iterator it = sdevX.begin(); it != sdevX.end(); ++it) {
		sum_of_elems += *it;
	}
	overallMis = sum_of_elems / moduleN;
	cout << "Manual Misalignment: " << endl;
	for (int i_module = 0; i_module < moduleN; i_module++) {
		cout << showpos << "Module " << i_module << " :: Characteristic:  " << sdevX[i_module] << " cm. "; // absolute misalignment [as set by MC]
		if (debugBool) {pede_mis << sdevX[i_module] << " "; }
		float relMisTmp = sdevX[i_module] - overallMis;
		// now push these misalignment parameters for all layers in the module [for use later]
		for (int i_view = 0; i_view < viewN; i_view++) {
			for (int i_layer = 0; i_layer < layerN; i_layer++) {
				relMis.push_back(relMisTmp);           //relative misalignment per layer [w.r.t to other modules]
				charMis.push_back(sdevX[i_module]);  // absolute misalignment per layer
			} // view
		} // layer
		cout << showpos << "Relative: " << relMisTmp << " cm." << endl;
	} // modules
	cout << noshowpos;
	cout << "The overall misalignment was " << overallMis << " cm" <<  endl << endl;

	//--- Calculations for Each Layer -- //
	// XXX ONLY VALID FOR 2-parameter line fit (i.e. not circle fit)
	// XXX Verbose calculation via separate loops on purpose - to demonstrate the alignment methods
	// For straight-line case only.

	// constants
	float N = float(layerTotalN); // total N of layers
	float D = Tracker::instance()->getResolution(); // Original detector resolution
	int i_totalLayers; // dummy layer counter

	// Pivot Point
	i_totalLayers = 0;
	for (int i_module = 0; i_module < moduleN; i_module++) {
		for (int i_view = 0; i_view < viewN; i_view++) {
			for (int i_layer = 0; i_layer < layerN; i_layer++) {
				pivotPoint_estimated += distance[i_totalLayers];
				i_totalLayers++;
			}// layer
		} // view
	} // modules
	pivotPoint_estimated = pivotPoint_estimated / N;

	// Centred Distances and squared distances
	i_totalLayers = 0;
	for (int i_module = 0; i_module < moduleN; i_module++) {
		for (int i_view = 0; i_view < viewN; i_view++) {
			for (int i_layer = 0; i_layer < layerN; i_layer++) {
				float tmpZDistance = distance[i_totalLayers] - pivotPoint_estimated;
				zDistance_centered.push_back(tmpZDistance);
				squaredZSum += pow(tmpZDistance, 2);
				i_totalLayers++;
			}// layer
		} // view
	} // modules

	// Estimated SD of residuals, Sum of mis., and squared sum of mis.
	i_totalLayers = 0;
	for (int i_module = 0; i_module < moduleN; i_module++) {
		sigma_recon_estimated.push_back(vector< vector< float > > ());
		for (int i_view = 0; i_view < viewN; i_view++) {
			sigma_recon_estimated[i_module].push_back( vector<float>  () );
			for (int i_layer = 0; i_layer < layerN; i_layer++) {
				float simga_est = D * sqrt( (N - 1) / N - (pow(zDistance_centered[i_totalLayers], 2) / squaredZSum) );
				sigma_recon_estimated[i_module][i_view].push_back(simga_est);
				MisZdistanceSum += charMis[i_totalLayers] * zDistance_centered[i_totalLayers];
				MisSum += charMis[i_totalLayers];
				i_totalLayers++;
			}// layer
		} // view
	} // modules

	// Predicted means of residuals [a.k.a "shear misalignment"]
	i_totalLayers = 0;
	for (int i_module = 0; i_module < moduleN; i_module++) {
		for (int i_view = 0; i_view < viewN; i_view++) {
			for (int i_layer = 0; i_layer < layerN; i_layer++) {
				float charM = charMis[i_totalLayers];
				float z = zDistance_centered[i_totalLayers];
				float shearMistmp = charM - MisSum / N - (z * MisZdistanceSum) / squaredZSum;
				shearMis.push_back(shearMistmp);
				i_totalLayers++;
			}// layer
		} // view
	} // modules

	// The Chi2 estimation requires misalignment parameters per layer
	Chi2_recon_estimated = N - 2.0; // 2 parameter fit
	i_totalLayers = 0;
	for (int i_module = 0; i_module < moduleN; i_module++) {
		for (int i_view = 0; i_view < viewN; i_view++) {
			for (int i_layer = 0; i_layer < layerN; i_layer++) {
				float M = shearMis[i_totalLayers];
				float z = zDistance_centered[i_totalLayers];
				Chi2_recon_estimated += ( (M * M) - (M) / (N) - (2.0 * M * z) / (squaredZSum) ) / (D * D) ;
				i_totalLayers++;
			}// layer
		} // view
	} // modules

}//end of misalign

/**
   Set filename to read uniform random numbers from.
   @param uniform_filename String for filename of file of uniform random numbers.
 */
void Tracker::set_uniform_file(string uniform_filename) {
	RandomBuffer::instance()->open_uniform_file(uniform_filename);
}

/**
   Set filename to read Gaussian random numbers from.
   @param uniform_filename String for filename of file of Gaussian random numbers.
 */
void Tracker::set_gaussian_file(string gaussian_filename) {
	RandomBuffer::instance()->open_gaussian_file(gaussian_filename);
}

float Tracker::generate_gaus() {
	float gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
	return gaus;
}

float Tracker::generate_uniform() {
	float uniform = (( RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()));
	return uniform;
}


/**
 * Returns the peak (maximum so far) resident set size [RSS] (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS. By Dr. David R. Nadeau:
 * http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use
 */
size_t Tracker::getPeakRSS( ) {
#if defined(_WIN32)
	/* Windows -------------------------------------------------- */
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
	/* AIX and Solaris ------------------------------------------ */
	struct psinfo psinfo;
	int fd = -1;
	if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
		return (size_t)0L;      /* Can't open? */
	if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
	{
		close( fd );
		return (size_t)0L;      /* Can't read? */
	}
	close( fd );
	return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	/* BSD, Linux, and OSX -------------------------------------- */
	struct rusage rusage;
	getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
	return (size_t)rusage.ru_maxrss;
#else
	return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
	/* Unknown OS ----------------------------------------------- */
	return (size_t)0L;          /* Unsupported. */
#endif
}