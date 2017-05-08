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


/**
   Structure to contain data of a generated track, with the number of hits, their positions, the uncertainty in the positions, and the plane number hit.
*/
struct LineData {
	int hit_count; /** Number of hits in detector */
	std::vector<float> z_hits; /** Z-positions of generated hits in detector */ //distances 
	std::vector<float> x_hits; /** X-positions of generated hits in detector */
	std::vector<float> x_true; /** X-positions of true hits in detector */
	std::vector<float> x_det; /** X-positions of recorded hits in detector */
	std::vector<float> x_mis; /** X-positions of misalignment in hits in detector */
	std::vector<float> hit_sigmas; /** Resolution for hits in detector */
	std::vector<int> i_hits; /** Number for plane struck in detector hits, with the plane numbers starting at 1, and increasing by one for each adjacent plane */
	std::vector<float> x0_gen; // generated points of the track
	std::vector<float> z0_gen; // generated points of the track
	std::vector<float> x1_gen; // generated points of the track
	std::vector<float> z1_gen; // generated points of the track
	std::vector<float> x_m; 
	std::vector<float> x_c;
};

/**
   Singleton class to represent a detector made up of an array of planar drift chambers, simulating the passage of linear tracks through the detector.
 */
class Tracker {

 private:

	static Tracker* s_instance; // Pointer to instance of class

	static const int trackCount=1000; /** Number of tracks (i.e. records) to be simulated passing through detector */
	//[all distances are in cm]
	static const int beamPositionLength = 10.0;  // max x position of beam origin [0, 10]
	static const float beamStart = 0.0; // z 
	static const float beamStop = 25.0;  // z 

	float resolution; //TODO decide if this is a constant  [value is assigned in the constructor]

	float dispX;
	//--------------------//
 
	///initialising physics variables
	static const int moduleN = 2; //number of movable detectors/module [independent modules]
	static const int viewN = 2; //There are two views per module (U and V)
	static const int layerN = 2; //there are 2 layers per view
	static const int strawN = 16; //number of measurement elements in x direction  [number of straws per layer]

	static const float twoR=2.0; //For normalisation of uniform random numbers [0,1] : (MAX+RND)/(twoR*MAX)

	//  define detector geometry [all distances are in cm]
	static const float startingDistanceModule0=5.0; // distance of first layer relative to the "beam" // [cm]
	static const float strawSpacing= 0.6;  // x distance between straws in a layer
	static const float layerSpacing = 0.5; // z distance between layers in a view
	static const float viewSpaing = 2.0; // z distance between views in a modules
	static const float moduleSpacing= 11; // z distance between modules
	static const float layerDisplacement = 0.3; // relative x distance between first straws in adjacent layers in a view [upstream layer is +x shifted] 

	//Area/volume/width required for MS (later on), and for rejection of "missed" hits [dca > strawRadius]
	static const float strawRadius =0.25; //thickness/width of a plane (X0) // [cm]
	//static const float stereoTheta=0.1309;  // stereo angle [rad]  // [rad] (7.5000 deg = 0.1309...rad)   // for later 3D versions
	
	std::vector<int> layer; // record of layers that were hit
    std::vector<float> projectionX; //projection of measurement direction in (X)
    
    std::vector<float> distance;  // distance between planes [this is set in geometry]
    std::vector<float> resolutionLayer;   //resolution [to record vector of resolution per layer if not constant] //XXX [this is not used at the moment]
    
    float sdevX[moduleN];// shift in x due to the imposed misalignment (alignment parameter)
    	
	// Class constructor and destructor
	Tracker();
	~Tracker();

 public:

	static Tracker* instance(); // Function to return pointer to the class instance

	// Function to simulate a track through the detector, then return data for plane hits.
	// Uses MC method to reject hits going outside of the detector
	
	float DCA(float, float, float, float, float, float);

	LineData MC(float, std::ofstream&, std::ofstream&, std::ofstream&, bool); 

	void setGeometry(std::ofstream&, bool); //Geometry of detector arrangement 

	void misalign(std::ofstream&, bool); // MC misalignment of detectors 
		
    void write_constraint_file(std::ofstream&, std::ofstream&, bool);  // Writes a constraint file for use with PEDE. 

	void set_uniform_file(std::string); // Set filename for uniform random numbers [randomIntGenerator.py - see https://github.com/glukicov/alignTrack]

	void set_gaussian_file(std::string); // Set filename for Gaussian random numbers


	//
	// Getter methods
	//

	std::vector<int> getLayer() {
		return layer;
	}
	
	std::vector<float> getProjectionX() {
		return projectionX;
	}

	int getStrawN(){
		return strawN;
	}


	int getTrackCount() {
		return trackCount;
	}

	float getStrawRadius() {
		return strawRadius;
	}

	float getResolution() {
		return resolution;
	}

	int getLayerN() {
		return layerN;
	}

	int getModuleN() {
		return moduleN;
	}

    float getDispX(){
		return dispX;
	}
	
	float getBeamStart(){
		return beamStart;
	}
	
	float getBeamStop(){
		return beamStop;
	}

};


#endif