#ifndef ALIGNTRACKER_M
#define ALIGNTRACKER_M

/* This header file contains definitions of constant variables used in
* the method class, as well as function declarations, and definitions of functions.
*/

///XXX some includes may be redundant
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
#include <TRandom3.h> // Random number facility XXX use art random? 
#include "Logger.hh"

/**
   Structure to contain data of a generated track, with the number of hits, their positions, the uncertainty in the positions, and the plane number hit.
*/
struct MCData {
	int totalLayerHits; /** Number of hits in detector */
	// *** Hit Parameters *** // 
	std::vector<float> residualsRecon; // Recon residuals between the fitted track and a drift circle 
	std::vector<float> residualsTruth; // Truth residuals
	std::vector<int> label1; // dR/dx
	std::vector<int> label2; // dR/dz 
	std::vector<int> label3; // dR/d𝛉
	std::vector<int> LR; // +1 = L -1=R - truth LR information 
	std::vector<float> strawID; // for reconstruction
	std::vector<float> dca; // truth dca 
	std::vector<float> driftRad; // recon dca of the hit [should be the same as truth!]
	std::vector<float> dcaUnsmeared;  
	std::vector<float> zStraw; // the z position of straw centre
	std::vector<float> xStraw; // the x position of straw centre
	std::vector<float> zCentreModule; // rotational centre of a module in z 
	std::vector<float> xCentreModule; // rotational centre of a module in x 
	//Detector coordinates [for a hit]
	std::vector<int> Module_i;
	std::vector<int> View_i;
	std::vector<int> Layer_i;
	std::vector<int> Straw_i;
	// ** Track parameters ** //
	float x0; // entrance of beam in x 
	float x1; // exit 
	float slopeTruth; 
	float interceptTruth;
	float slopeRecon;
	float interceptRecon;
	float pValue;
	float chi2Circle;
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
	// Straw positions:
	float z;
	float x;
	float dcaRecon; // should be the same as truth
};

// Returned once per track
struct ResidualData {
	std::vector<float> residualsRecon; // residual between the fitted track and radius of the fit circle;
	float slopeRecon;
	float interceptRecon;
	float pValue;
	float chi2Circle;
};

// Define and return the centre of rotational symmetry
struct RotationCentres {
	std::vector<float> zCentres;
	std::vector<float> xCentres;
};

/**
   Singleton class to represent the detector
 */
class Tracker {

private:

	static Tracker* s_instance; // Pointer to instance of class
	int trackNumber; // Number of tracks (i.e. records) to be simulated passing through detector - passed as command line argument

	// XXX use as input for grid jobs
	int randomSeed = 1234;
	TRandom3* randomFacility = new TRandom3(randomSeed);

	// **** COUNTERS ****  //
	int rejectedHitsDCA = 0; // rejected hits due to DCA > straw radius
	int multipleHitsLayer = 0; // passed over from DCAData [if >1 hit per layer]
	int ambiguityHit = 0; //Exactly in the middle of 2 straws
	bool cutTriggered; // FLAG: set as false at each track generation, triggered if smeared DCA < trackCut [below]
	bool missedHit; // FLAG: set as false at each track generation, triggered if smeared DCA > strawRadius

	//initialising physics variables
	// MF + inhomogeneity, E_loss, MS
	static constexpr float resolution =0.015; // 150um = 0.015 cm for hit smearing
	static constexpr float trackCut = 0.05; //500 um = 0.5 mm for dca cut on tracks
	static constexpr int layerCut = 11;
	static constexpr int nlc = 2; // dR/dc dR/dm
	static constexpr int ngl = 1; //  dR/dx dR/dz  dR/d𝛉
	//Matrix memory space
	int matN = moduleN; // # number of global parameters
	int matNC = 0; // # number of constraints

	float pValCut = 0.0; // from 0.0->1.0
	bool trackCutBool = true; // if true, tracks will be rejected if DCA < trackCut
	bool useTruthLR = true; // use LR information from generated tracks [requires DCA cut!]
	bool hitCut = true; // if true, hits will be rejected if DCA > strawRadius
	
	//Set truth misalignment of modules 
	// float dispX[4] =     { 0.00, 0.01, -0.01, 0.00}; // manual misalignment [relative misalignment per module]
	// float dispX[4] =     { 0.00, 0.00, 	0.020, 0.00}; // manual misalignment [relative misalignment per module]
	float dispX[4] =     { 0.00, 0.05, 	-0.05, 0.00}; // manual misalignment [relative misalignment per module]
	float dispZ[4] =     {0.00,  0.00,  0.0,  0.00}; // manual misalignment [relative misalignment per module]
	float dispTheta[4] = {0.00,  0.00,  0.0,  0.00}; // radians

	// **** GEOMETRIC CONSTANTS ****  // XXX will be taken from gm2geom in the future
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
	// The offset is a "software" fix to get the ideal (assumed) geometry closer to the truth geometry
	float offsetX[8] =     {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // To the ideal detector, for second pede iteration
	float offsetZ[8] =     {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // To the ideal detector, for second pede iteration
	float offsetTheta[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //radians
	
	// **** BEAM PARAMETERS ****  // [all distances are in cm]
	//static constexpr float beamPositionLength = strawN*strawSpacing+strawSpacing; // max x coordinate = beamPositionLength - beamOffset; mix x = -dispX
	static constexpr float beamPositionLength = 2.0;
	static constexpr float beamOffset = 1.0; // offset from 0 in x
	static constexpr float beamStart =-100.0; // z
	static constexpr float beamStop = 160.0; // z

	// **** MC CALCULATION CONTAINERS ****  //
	// Hits-based
	std::vector<int> layer; // record of layers that were hit

	// Vectors to hold Ideal and Misaligned (true) X positions [moduleN 0-7][viewN 0-1][layerN 0-1][strawN 0-31]
	std::vector< std::vector< std::vector< std::vector< float > > > > mod_lyr_strawIdealPositionX;
	std::vector< std::vector< std::vector< std::vector< float > > > > mod_lyr_strawMisPositionX;
	std::vector< std::vector< std::vector< std::vector< float > > > > mod_lyr_strawIdealPositionZ;  // Z distance between planes [this is set in geometry]
	std::vector< std::vector< std::vector< std::vector< float > > > > mod_lyr_strawMisPositionZ;  // Z distance between misaligned planes [this is set in misalignment]

	// Vector to store the mapping of Views and Layers in a module U0, U1, V0, V1
	std::vector< std::vector< std::string > > UVmapping; // set in the constructor

	//vector to store x coordinates of the track as seen from the ideal detector
	std::vector<float> xRecon;  // ideal straw x
	std::vector<float> zRecon;  // ideal straw z
	std::vector<float> radRecon; // reconstructed fit circle radius (DCA)

	// Class constructor and destructor
	Tracker();
	~Tracker();

public:

	static Tracker* instance(); // Function to return pointer to the class instance

	// Function to simulate a track through the detector, then return data for plane hits.
	// Uses MC method to reject hits going outside of the detector

	float pointToLineDCA(float zStraw, float xStraw, float xSlopeTruth, float xInterceptTruth); //simple 2D DCA

	DCAData DCAHit(std::vector<float> xLayer, std::vector<float> zLayer, float xTrack, float xSlpoeTruth, float xInterceptTruth, bool debugBool); // calls DCA to chose the right straw

	ReconData HitRecon(int detID, float detDca, std::vector<float> xLayer, std::vector<float> zLayer); //return the dca to ideal geometry, and centre of the circle

	ResidualData GetResiduals(std::vector<float> zRecon, std::vector<float> xRecon, std::vector<float> radRecon, int dataSize, std::ofstream& plot_fit, bool debugBool, bool useTruth, std::vector<int> LR_truth);

	MCData MCLaunch(float scatterError, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, bool);

	void setGeometry(std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, bool, std::ofstream&); // MC misalignment of detectors

	//This function will return a vector for the centre of rotation for all modules
	RotationCentres getCentre(std::vector< std::vector< std::vector< std::vector< float > > > > mod_lyr_strawPositionX, std::vector< std::vector< std::vector< std::vector< float > > > > mod_lyr_strawPositionZ, std::ofstream& plot_centres);

	void writeConstraintFile(std::ofstream&, std::ofstream&, bool, std::ofstream& metric);  // Writes a constraint file for use with PEDE.

	void writePresigmaFile(std::ofstream&, std::ofstream& metric);  // Writes a pre-sigma parameter file for use with PEDE.

	void writeSteeringFile(std::ofstream&, std::ofstream& metric); // Steering file for PEDE.

	//
	// Setter methods
	//
	void setTrackNumber(int tracks) {
		trackNumber = tracks;
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
	std::string getUVmapping(int i, int j) {
		return UVmapping[i][j];
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

	const int getNLC() {
		return nlc;
	}
	
	const int getNGL() {
		return ngl;
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