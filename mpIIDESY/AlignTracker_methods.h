#ifndef ALIGNTRACKER_M
#define ALIGNTRACKER_M

/* This header file contains definitions of constant variables used in 
* the method class, as well as function declarations, and definitions of functions.
*/  

#include "random_buffer.h" // courtesy of John Smeaton (UCL)
///XXX some includes may become redundant
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <TMath.h>


/**
   Structure to contain data of a generated track, with the number of hits, their positions, the uncertainty in the positions, and the plane number hit.
*/
struct MCData {
	int hit_count; /** Number of hits in detector */
	//Hits
	std::vector<float> z_hits; /** Z-positions of generated hits in detector */ //distances 
	std::vector<float> x_residuals; /** X-positions of residuals between ideal and real geometry */
	std::vector<float> hit_sigmas; /** Resolution for hits in detector */
	std::vector<int> i_hits; /* vector of modules that were actually hit [after passing rejection test] */
	std::vector<int> hit_list;  // same of layers (absolute)
	std::vector<int> hit_bool;  // same of layers (absolute) +1 = hit 0=no hit
 	std::vector<string> absolute_straw_hit; // for plotting
 	std::vector<float> LR;
	std::vector<float> residuals_gen;
	std::vector<int> residuals_fit;
 	std::vector<float> strawID;
  	std::vector<float> x_mis_dca; /** X-positions of recorded hits in a real detector */
	std::vector<float> x_hit_recon; // ideal straw hit position + dca (from mis.)
	std::vector<float> x_hit_true; // true xHit position
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
	std::vector<float> x_track_true; /** X-positions of true hits (generated line x coordinate in line with a layer) in detector */
	std::vector<float> x_track_recon; // reconstructed x position of the line 
	float p_value; 
	float chi2;
};

// DCA structure - calculated for each hit
struct DCAData{
	int strawID; 
	float dca;
	float LRSign; // L=-ive, R=+ive 
};

// Returned - calculated once per track
struct ResidualData{
	std::vector<float> residuals; // residual between  
	std::vector<float> x_fitted;   // the x-coordinate of the fitted line in a layer [z-coordinate corresponds to distance vector]
	float slope_recon;
	float intercept_recon;
	float meanXReconTrack;
	float meanZReconTrack;
	float p_value; 
	float chi2; 
};

/**
   Singleton class to represent the detector
 */
class Tracker {

 private:

	static Tracker* s_instance; // Pointer to instance of class

	int trackNumber; // Number of tracks (i.e. records) to be simulated passing through detector - passed as command line argument
	
	static constexpr float twoR=2.0; //For normalisation of uniform random numbers [0,1] : (MAX+RND)/(twoR*MAX)

 	// **** COUNTERS ****  //
    int rejectedHitsDCA=0;  // rejected hits due to DCA > straw radius 
    int multipleHitsLayer=0; // passed over from DCAData [if >1 hit per layer]
    int ambiguityHit=0; //Exactly in the middle of 2 straws
    int hitLayerCounter; // absolute layer ID for the hit
	
	//initialising physics variables
 	// MF + inhomogeneity, E_loss, MS

    float dispX[8] = {0.0, 0.0, 0.0, 0.00, 0.0, 0.0, 0.0}; // manual misalignment [relative misalignment per module]
    
    static constexpr float resolution=0.015;  // 150um = 0.015 cm for hit smearing
 	  
	// **** GEOMETRIC CONSTANTS ****  //
	// define detector geometry [all distances are in cm]
	static const int moduleN = 4; //number of movable detectors/module [independent modules]
	static const int strawN = 8; //number of measurement elements in x direction  [number of straws per layer]
	static const int viewN = 2; //There are two views per module (U and V) XXX
	static const int layerN = 2; //there are 2 layers per view [4 layers per module]
	static const int layerTotalN = layerN*viewN*moduleN; // total number of layers 
	static constexpr float startingZDistanceStraw0 = 5.0; // distance of the very first layer [the first straw] relative to the "beam" in z // [cm]
	static constexpr float startingXDistanceStraw0 = -2.0; // distance of the very first layer [the first straw] in x // [cm]
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
	static constexpr float beamOffset=1.0; // offset from 0 in x
	static constexpr float beamStart = startingZDistanceStraw0-5.0; // z 
	static constexpr float beamStop = (moduleSpacing+viewSpacing+layerSpacing*float(layerN))*float(moduleN);  // z  

	// **** MC CALCULATION CONTAINERS ****  // 
	// Hits-based
	std::vector<int> layer; // record of layers that were hit
    std::vector<float> projectionX; //projection of measurement direction in (X)
    std::vector<float> distance;  // Z distance between planes [this is set in geometry]
    std::vector<float> resolutionLayer;   //resolution [to record vector of resolution per layer if not constant] //XXX [this is not used at the moment]
    std::vector<float> sdevX;// shift in x due to the imposed misalignment (alignment parameter)

    // Vectors to hold Ideal and Misaligned (true) positions [moduleN 0-7][viewN 0-1][layerN 0-1][strawN 0-31]
    std::vector< std::vector< std::vector< std::vector< float > > > > mod_lyr_strawIdealPosition; 
    std::vector< std::vector< std::vector< std::vector< float > > > > mod_lyr_strawMisPosition;  
    std::vector< std::vector< std::vector< float > > > sigma_recon_estimated; 

    // Vector to store the mapping of Views and Layers in a module U0, U1, V0, V1
    std::vector< std::vector< string > > UVmapping; // set in the constructor

    //Misalignment 
    vector<float> charMis;  // The alignment parameter: absolute misalignment of a plane 
    vector<float> relMis;  // Relative misalignment (w.r.t to overall mis. - per layer) 
    vector<float> shearMis; // vector to hold the shear misalignment for each plane [mean of the residuals per plane]
    vector<float> zDistance_centered; // SD calculations assumes mean z of 0 [vector to hold z distances w.r.t to pivot point]
    float overallMis; // the overall misalignment - calculated in the misalignment method  [same for all modules]
    float Chi2_recon_estimated=0.0;  // Estimated (recon.) Chi2 of the fit 
    float pivotPoint_estimated=0.0;  // Mean z point (pivot point)
    float squaredZSum = 0.0; // Squared sum of centred Z distances 
    float MisZdistanceSum=0.0; // Sum of Mis * z distance
    float MisSum = 0.0; // Sum of all misalignments [per layer]
    	
	// Class constructor and destructor
	Tracker();
	~Tracker();

 public:

	static Tracker* instance(); // Function to return pointer to the class instance

	// Function to simulate a track through the detector, then return data for plane hits.
	// Uses MC method to reject hits going outside of the detector
	
	float generate_gaus(); // using the RandomBuffer class

	float generate_uniform(); // using the RandomBuffer class

	float DCA(float, float);

	DCAData DCAHit(std::vector<float>, float, float, bool);

	float HitRecon(int, float, float, std::vector<float>);

	//ResidualData GetResiduals(std::vector<float>,  std::vector<float>, std::ofstream&, bool);

	ResidualData GetResiduals(std::vector<float>,  std::vector<float>, int, std::ofstream&, bool);

	MCData MC_launch(float, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, bool); 

	void setGeometry(std::ofstream&, bool); //Geometry of detector arrangement 

	void misalign(std::ofstream&, std::ofstream&, bool); // MC misalignment of detectors 
		
    void write_constraint_file(std::ofstream&, std::ofstream&, bool);  // Writes a constraint file for use with PEDE. 

    void write_steering_file(std::ofstream&); // Steering file for PEDE.

	void set_uniform_file(std::string); // Set filename for uniform random numbers [randomIntGenerator.py - see https://github.com/glukicov/alignTrack]

	void set_gaussian_file(std::string); // Set filename for Gaussian random numbers

	size_t getPeakRSS( ); // Peak Dynamic Memory used 

	size_t getCurrentRSS( );

	//
	// Setter methods
	//
	void setTrackNumber(int tracks){
		trackNumber=tracks; 
	}

	void incRejectedHitsDCA(){
		rejectedHitsDCA=rejectedHitsDCA+1;
	}

	void incMultipleHitsLayer(){
		multipleHitsLayer=multipleHitsLayer+1;
	}

	void incAmbiguityHit(){
		ambiguityHit=ambiguityHit+1;
	}
    
	//
	// Getter methods
	//

	float get_Chi2_recon_estimated(){
		return Chi2_recon_estimated;
	}

	float get_pivotPoint_estimated(){
		return pivotPoint_estimated;
	}

	vector<float> get_sigma_recon_estimatedVector(){
		vector<float> result;
		for (int i_module=0; i_module<moduleN; i_module++){
        	for (int i_view=0; i_view<viewN; i_view++){
            	for (int i_layer=0; i_layer<layerN; i_layer++){		
					result.push_back(sigma_recon_estimated[i_module][i_view][i_layer]);
      			}
  	 		}
  	 	}	
		return result;
	}

	float get_shearMis(int i){
		return shearMis[i];
	}


	float get_sigma_recon_estimated(int i, int j, int k){
		return sigma_recon_estimated[i][j][k];
	}

	string getUVmapping(int i, int j){
		return UVmapping[i][j];
	}

	std::vector<float> getZDistanceVector() {
		return distance;
	}
	
	float getZDistance(int i) {
		return distance[i];
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

	float getOverallMis() {
		return overallMis;
	}

	int getRejectedHitsDCA(){
		return rejectedHitsDCA;
	}

	int getMultipleHitsLayer(){
		return multipleHitsLayer;
	}

	int getAmbiguityHit(){
		return ambiguityHit;
	}


	int getStrawN(){
		return strawN;
	}

	int getViewN(){
		return viewN;
	}

	int getTrackNumber() {
		return trackNumber;
	}

	float getStrawRadius() {
		return strawRadius;
	}

	float getStrawSpacing(){
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

    	
	float getBeamStart(){
		return beamStart;
	}
	
	float getBeamStop(){
		return beamStop;
	}

	float getBeamOffset(){
		return beamOffset;
	}

	float getBeamPositionLength(){
		return beamPositionLength;
	}

};


#endif