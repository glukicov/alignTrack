#ifndef ALIGNTRACKER_M
#define ALIGNTRACKER_M

/* This header file contains definitions of constant variables used in
* the method class, as well as function declarations, and definitions of functions.
*/
#include "random_buffer.h" // courtesy of John Smeaton (UCL)
///XXX some includes may become redundant
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <bitset>
#include <TF1.h>
#include <TMath.h>

/**
   Structure to contain data of a generated track, with the number of hits, their positions, the uncertainty in the positions, and the plane number hit.
*/
struct MCData {
	int hit_count; /** Number of hits in detector */
	//Hits
	std::vector<float> z_true; /** Z-positions of generated hits in detector */ //distances
	std::vector<float> z_recon; /** Z-positions of recon hits in detector */ //distances
	std::vector<float> residuals; /** X-positions of residuals between ideal and real geometry */
	std::vector<float> hit_sigmas; /** Resolution for hits in detector */
	std::vector<int> i_hits; /* vector of modules that were actually hit [after passing rejection test] */
	std::vector<int> label_1;
	std::vector<int> label_2;
	std::vector<int> hit_list;  // same of layers (absolute)
	std::vector<int> hit_bool;  // same of layers (absolute) +1 = hit 0=no hit
	std::vector<string> absolute_straw_hit; // for plotting
	std::vector<int> LR; // +1 = L -1=R
	std::vector<float> residuals_gen;
	std::vector<int> residuals_fit;
	std::vector<float> strawID;
	std::vector<float> dca; /** X-positions of recorded hits in a real detector */
	std::vector<float> dca_unsmeared; /** X-positions of recorded hits in a real detector */
	//Detector coordinates
	std::vector<int> Module_i;
	std::vector<int> View_i;
	std::vector<int> Layer_i;
	std::vector<int> Straw_i;
	//track parameters
	float x0;
	float x1;
	float slope_truth;
	float intercept_truth;
	float slope_recon;
	float intercept_recon;
	float meanXReconTrack;
	float meanZReconTrack;
	std::vector<float> driftRad;
	std::vector<float> z_straw;
	std::vector<float> x_straw;
	std::vector<float> x_track_true; /** X-positions of true hits (generated line x coordinate in line with a layer) in detector */
	std::vector<float> x_track_recon; // reconstructed x position of the line
	float p_value;
	float chi2_circle;
	bool cut = false; // cut trigger to kill the track
};

// DCA structure - calculated for each hit
struct DCAData {
	int strawID;
	float dcaUnsmeared;
	float dca; //smeared dca
	float residualTruth;
	int LRSign; // L=+1, R=-1
};

struct ReconData {
	float z;
	float x;
	float dcaRecon;
};


// Returned - calculated once per track
struct ResidualData {
	std::vector<float> residuals; // residual between the (centre of the straw and the fitted line [pointToLineDCA]) and radius of the fit circle;
	float slope_recon;
	float intercept_recon;
	float meanXReconTrack;
	float meanZReconTrack;
	std::vector<float> x_track_recon;
	float p_value;
	float chi2_circle;
};

/**
   Singleton class to represent the detector
 */
class Tracker {

private:

	static Tracker* s_instance; // Pointer to instance of class

	int trackNumber; // Number of tracks (i.e. records) to be simulated passing through detector - passed as command line argument

	static constexpr float twoR = 2.0; //For normalisation of uniform random numbers [0,1] : (MAX+RND)/(twoR*MAX)

	// **** COUNTERS ****  //
	int rejectedHitsDCA = 0; // rejected hits due to DCA > straw radius
	int multipleHitsLayer = 0; // passed over from DCAData [if >1 hit per layer]
	int ambiguityHit = 0; //Exactly in the middle of 2 straws
	int hitLayerCounter; // absolute layer ID for the hit
	bool cutTriggered; // FLAG: set as false at each track generation, triggered if smeared DCA < trackCut [below]
	bool missedHit; // FLAG: set as false at each track generation, triggered if smeared DCA > strawRadius

	//initialising physics variables
	// MF + inhomogeneity, E_loss, MS

	//float dispX[8] = {0.0, 0.03, -0.03, 0.00, 0.0, 0.0, 0.0, 0.0}; // manual misalignment [relative misalignment per module]
	//float dispZ[8] = {0.0, 0.005, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0}; // manual misalignment [relative misalignment per module]
	//float offsetX[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // To the ideal detector, for second pede iteration
	//float offsetZ[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // To the ideal detector, for second pede iteration
	//float dispTheta[8] = {0.0, 0.611, -0.610, 0.0, 0.0, 0.0, 0.0}; // radians
	float dispTheta[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // radians
	float offsetTheta[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //radians

	static constexpr float resolution = 0.015; // 150um = 0.015 cm for hit smearing
	static constexpr float trackCut = 0.05; //500 um = 0.5 mm for dca cut on tracks

	static const int nlc = 2; // dR/dc dR/dm
	// static const int ngl = 2; // dR/dx dR/dz 
	static const int ngl = 1; // dR/d𝛉 

	float pValCut = 0.00; // from 0->1
	bool trackCutBool = true; // if true, tracks will be rejected if DCA > trackCut
	bool hitCut = false; // if true, hits will be rejected if DCA > strawRadius
	bool useTruthLR = true;


	// **** GEOMETRIC CONSTANTS ****  //
	// define detector geometry [all distances are in cm]
	static const int moduleN = 4; //number of movable detectors/module [independent modules]
	static const int strawN = 8; //number of measurement elements in x direction  [number of straws per layer]
	static const int viewN = 2; //There are two views per module (U and V) XXX
	static const int layerN = 2; //there are 2 layers per view [4 layers per module]
	static const int layerTotalN = layerN * viewN * moduleN; // total number of layers
	static constexpr float startingZDistanceStraw0 = 5.0; // distance of the very first layer [the first straw] relative to the "beam" in z // [cm]
	static constexpr float startingXDistanceStraw0 = 2.0; // distance of the very first layer [the first straw] in x // [cm]
	static constexpr float strawSpacing = 0.606;  // x distance between straws in a layer
	static constexpr float layerSpacing = 0.515; // z distance between layers in a view
	static constexpr float viewSpacing = 2.020; // z distance between views in a modules
	static constexpr float moduleSpacing = 13.735; // z distance between modules' first layers [first layer of module 1 and first layer of module 2]
	static constexpr float layerDisplacement = 0.303; // relative x distance between first straws in adjacent layers in a view [upstream layer is +x shifted]
	//std::vector<float> staircaseXDisplacment = [0.0, 0.723, 0.979, 1.238, 1.503, 1.755, 2.008, 2.259];  // TODO staircase for MF in future [cm]
	//[ 6880.52, 6873.29, 6863.50, 6851.12, 6836.09, 6818.54, 6798.46, 6775.87]; // from gm2geom [mm]
	//Area/volume/width required for MS (later on), and for rejection of "missed" hits [dca > strawRadius]
	static constexpr float strawRadius = 0.2535; // takes as the outerRadiusOfTheGas from gm2geom/strawtracker/strawtracker.fcl // [cm]
	static constexpr float stereoTheta = 0.1309;  // stereo angle [rad]  // [rad] (7.5000 deg = 0.1309...rad)   // XXX for later 3D versions

	// **** BEAM PARAMETERS ****  // [all distances are in cm]
	//static constexpr float beamPositionLength = strawN*strawSpacing+strawSpacing; // max x coordinate = beamPositionLength - beamOffset; mix x = -dispX
	static constexpr float beamPositionLength = 2.0;
	static constexpr float beamOffset = 1.0; // offset from 0 in x
	static constexpr float beamStart = startingZDistanceStraw0 - 5.0; // z
	static constexpr float beamStop = (moduleSpacing + viewSpacing + layerSpacing*float(layerN)) * float(moduleN); // z

	// **** MC CALCULATION CONTAINERS ****  //
	// Hits-based
	std::vector<int> layer; // record of layers that were hit
	std::vector<float> projectionX; //projection of measurement direction in (X)
	std::vector<float> distanceIdealZ;  // Z distance between planes [this is set in geometry]
	std::vector<float> distanceMisZ;  // Z distance between misaligned planes [this is set in misalignment]
	std::vector<float> resolutionLayer;   //resolution [to record vector of resolution per layer if not constant] //XXX [this is not used at the moment]
	std::vector<float> sdevX;// shift in x due to the imposed misalignment (alignment parameter)
	std::vector<float> sdevZ;// shift in x due to the imposed misalignment (alignment parameter)

	// Vectors to hold Ideal and Misaligned (true) X positions [moduleN 0-7][viewN 0-1][layerN 0-1][strawN 0-31]
	std::vector< std::vector< std::vector< std::vector< float > > > > mod_lyr_strawIdealPosition;
	std::vector< std::vector< std::vector< std::vector< float > > > > mod_lyr_strawMisPosition;
	std::vector< std::vector< std::vector< float > > > sigma_recon_estimated;

	// Vector to store the mapping of Views and Layers in a module U0, U1, V0, V1
	std::vector< std::vector< string > > UVmapping; // set in the constructor

	//vector to store x coordinates of the track as seen from the ideal detector
	vector<float> xRecon;  // ideal straw x
	vector<float> zRecon;  // ideal straw z
	vector<float> radRecon; // reconstructed fit circle radius (DCA)

	//Misalignment
	vector<float> charMis;  // The alignment parameter: absolute misalignment of a plane
	vector<float> relMis;  // Relative misalignment (w.r.t to overall mis. - per layer)
	vector<float> shearMis; // vector to hold the shear misalignment for each plane [mean of the residuals per plane]
	vector<float> zDistance_centered; // SD calculations assumes mean z of 0 [vector to hold z distances w.r.t to pivot point]
	float overallMis; // the overall misalignment - calculated in the misalignment method  [same for all modules]
	float Chi2_recon_estimated = 0.0; // Estimated (recon.) Chi2 of the fit
	float pivotPoint_estimated = 0.0; // Mean z point (pivot point)
	float squaredZSum = 0.0; // Squared sum of centred Z distances
	float MisZdistanceSum = 0.0; // Sum of Mis * z distance
	float MisSum = 0.0; // Sum of all misalignments [per layer]

	//Matrix memory space
	int mat_n = moduleN; // # number of global parameters
	int mat_nc = 0; // # number of constraints

	// Class constructor and destructor
	Tracker();
	~Tracker();

public:

	static Tracker* instance(); // Function to return pointer to the class instance

	// Function to simulate a track through the detector, then return data for plane hits.
	// Uses MC method to reject hits going outside of the detector

	float generate_gaus(); // using the RandomBuffer class

	float generate_uniform(); // using the RandomBuffer class

	float pointToLineDCA(float z_straw, float x_straw, float x_slope, float x_intercept); //simple 2D DCA

	DCAData DCAHit(std::vector<float> xLayer, float zStraw, float xTrack, float xSlpoe, float xIntercept, bool debugBool); // calls DCA to chose the right straw

	ReconData HitRecon(int det_ID, float det_dca, std::vector<float> xLayer, float z_distance); //return the dca to ideal geometry, and centre of the circle

	ResidualData GetResiduals(std::vector<float> zRecon, std::vector<float> xRecon, std::vector<float> radRecon, int dataSize, std::ofstream& plot_fit, bool debugBool, bool useTruth, std::vector<int> LR_truth);

	MCData MC_launch(float, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, bool);

	void setGeometry(std::ofstream&, bool); //Geometry of detector arrangement

	void misalign(std::ofstream&, std::ofstream&, bool, std::ofstream& metric); // MC misalignment of detectors

	void write_constraint_file(std::ofstream&, std::ofstream&, bool, std::ofstream& metric);  // Writes a constraint file for use with PEDE.

	void write_presigma_file(std::ofstream&, std::ofstream& metric);  // Writes a pre-sigma parameter file for use with PEDE.

	void write_steering_file(std::ofstream&, std::ofstream& metric); // Steering file for PEDE.

	void set_uniform_file(std::string); // Set filename for uniform random numbers [randomIntGenerator.py - see https://github.com/glukicov/alignTrack]

	void set_gaussian_file(std::string); // Set filename for Gaussian random numbers

	//size_t getPeakRSS( ); // Peak Dynamic Memory used

	//size_t getCurrentRSS( );

	//
	// Setter methods
	//
	void setTrackNumber(int tracks) {
		trackNumber = tracks;
	}

	// void setXOffset1(float off1) {
	// 	offsetX[1] = off1;
	// }

	// void setXOffset2(float off2) {
	// 	offsetX[2] = off2;
	// }

	// void setZOffset1(float off1) {
	// 	offsetZ[1] = off1;
	// }

	// void setZOffset2(float off2) {
	// 	offsetZ[2] = off2;
	// }

	void setThetaOffset1(float off1) {
		offsetTheta[1] = off1;
	}

	void setThetaOffset2(float off2) {
		offsetTheta[2] = off2;
	}

	void incRejectedHitsDCA() {
		rejectedHitsDCA = rejectedHitsDCA + 1;
	}

	void incMultipleHitsLayer() {
		multipleHitsLayer = multipleHitsLayer + 1;
	}

	void incAmbiguityHit() {
		ambiguityHit = ambiguityHit + 1;
	}

	//
	// Getter methods
	//

	float get_Chi2_recon_estimated() {
		return Chi2_recon_estimated;
	}

	float get_pivotPoint_estimated() {
		return pivotPoint_estimated;
	}

	vector<float> get_sigma_recon_estimatedVector() {
		vector<float> result;
		for (int i_module = 0; i_module < moduleN; i_module++) {
			for (int i_view = 0; i_view < viewN; i_view++) {
				for (int i_layer = 0; i_layer < layerN; i_layer++) {
					result.push_back(sigma_recon_estimated[i_module][i_view][i_layer]);
				}
			}
		}
		return result;
	}

	float get_shearMis(int i) {
		return shearMis[i];
	}


	float get_sigma_recon_estimated(int i, int j, int k) {
		return sigma_recon_estimated[i][j][k];
	}

	string getUVmapping(int i, int j) {
		return UVmapping[i][j];
	}

	std::vector<float> getIdealZDistanceVector() {
		return distanceIdealZ;
	}

	float getIdealZDistance(int i) {
		return distanceIdealZ[i];
	}

	int getLayer(int i) {
		return layer[i];
	}

	float getProjectionX(int i) {
		return projectionX[i];
	}

	float getSdevX(int i) {
		return sdevX[i];
	}

	float getSdevZ(int i) {
		return sdevZ[i];
	}

	float getOverallMis() {
		return overallMis;
	}

	int getRejectedHitsDCA() {
		return rejectedHitsDCA;
	}

	int getMultipleHitsLayer() {
		return multipleHitsLayer;
	}

	int getAmbiguityHit() {
		return ambiguityHit;
	}


	int getStrawN() {
		return strawN;
	}

	int getViewN() {
		return viewN;
	}

	int getTrackNumber() {
		return trackNumber;
	}

	float getStrawRadius() {
		return strawRadius;
	}

	float getStrawSpacing() {
		return strawSpacing;
	}

	float getResolution() {
		return resolution;
	}

	int getLayerN() {
		return layerN;
	}

	int getLayerTotalN() {
		return layerTotalN;
	}

	int getModuleN() {
		return moduleN;
	}


	float getBeamStart() {
		return beamStart;
	}

	float getBeamStop() {
		return beamStop;
	}

	float getBeamOffset() {
		return beamOffset;
	}

	float getBeamPositionLength() {
		return beamPositionLength;
	}

	bool getLRStatus() {
		return useTruthLR;
	}

	bool getHitCutStatus() {
		return hitCut;
	}

	bool getTrackCutBool() {
		return trackCutBool;
	}

	float getPValCut() {
		return pValCut;
	}

	float getTrackCut() {
		return trackCut;
	}

};


#endif