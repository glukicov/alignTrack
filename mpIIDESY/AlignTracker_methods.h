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
	std::vector<float> x_track; /** X-positions of true hits in detector */
	std::vector<float> x_mis_dca; /** X-positions of recorded hits in detector */
	std::vector<float> x_ideal; // ideal straw potion + dca (from mis.)
	std::vector<float> hit_sigmas; /** Resolution for hits in detector */
	std::vector<int> i_hits; /** Number for plane struck in detector hits, with the plane numbers starting at 1, and increasing by one for each adjacent plane */
	std::vector<float> x_m; // generated slopes
	std::vector<float> x_c; //generated intercepts

	std::vector<float> x0_gen; 
	std::vector<float> x1_gen; 
	std::vector<float> z0_gen; 
	std::vector<float> z1_gen; 
};


struct DCAData{
	int strawID; 
	float dca;
	float LRSign; // L=-ive, R=+ive 
};

/**
   Singleton class to represent a detector made up of an array of planar drift chambers, simulating the passage of linear tracks through the detector.
 */
class Tracker {

 private:

	static Tracker* s_instance; // Pointer to instance of class

	//[all distances are in cm]
	static constexpr float beamPositionLength = 1.3; // max x coordinate = beamPositionLength + beamOffset
	static constexpr float beamOffset=0.6; // offset from 0 in x
	static constexpr float beamStart = 0.0; // z 
	static constexpr float beamStop = 25.0;  // z 

	float resolution; //TODO decide if this is a constant  [value is assigned in the constructor]
	float dispX; // [value is assigned in the constructor]
 	
 	static constexpr float twoR=2.0; //For normalisation of uniform random numbers [0,1] : (MAX+RND)/(twoR*MAX)
	
	//initialising physics variables
 	//XXX
	
	//  define detector geometry [all distances are in cm]
	static const int moduleN = 2; //number of movable detectors/module [independent modules]
	//static const int viewN = 2; //There are two views per module (U and V) XXX
	static const int layerN = 2; //there are 2 layers per view [4 layers per module]
	static const int layerTotalN = layerN*moduleN; // total number of layers 
	static const int strawN = 4; //number of measurement elements in x direction  [number of straws per layer]
	static constexpr float startingZDistanceStraw0 = 5.0; // distance of the very first layer [the first straw] relative to the "beam" in z // [cm]
	static constexpr float startingXDistanceStraw0 = 0.0; // distance of the very first layer [the first straw] in x // [cm]
	static constexpr float strawSpacing = 0.606;  // x distance between straws in a layer
	static constexpr float layerSpacing = 0.515; // z distance between layers in a view
	static constexpr float viewSpaing = 2.020; // z distance between views in a modules
	static constexpr float moduleSpacing = 13.735; // z distance between modules' first layers [first layer of module 1 and first layer of module 2]
	static constexpr float layerDisplacement = 0.303; // relative x distance between first straws in adjacent layers in a view [upstream layer is +x shifted] 
	
	//Area/volume/width required for MS (later on), and for rejection of "missed" hits [dca > strawRadius]
	static constexpr float strawRadius = 0.25; //thickness/width of a plane (X0) // [cm]
	static constexpr float stereoTheta = 0.1309;  // stereo angle [rad]  // [rad] (7.5000 deg = 0.1309...rad)   // for later 3D versions
	
	std::vector<int> layer; // record of layers that were hit
    std::vector<float> projectionX; //projection of measurement direction in (X)
    
    std::vector<float> distance;  // Z distance between planes [this is set in geometry]
    std::vector<float> resolutionLayer;   //resolution [to record vector of resolution per layer if not constant] //XXX [this is not used at the moment]
    std::vector<float> sdevX;// shift in x due to the imposed misalignment (alignment parameter)
    int trackNumber; /** Number of tracks (i.e. records) to be simulated passing through detector */


    // Vectors to hold Ideal and Misaligned (true) positions [moduleN 0-7][viewN 0-1][layerN 0-1][strawN 0-31]
    std::vector< std::vector< std::vector< float > > > mod_lyr_strawIdealPosition; 
    std::vector< std::vector< std::vector< float > > > mod_lyr_strawMisPosition;  
    
    //vector to store x coordinates of the track as seen from the ideal detector 
     std::vector<float> x_idealPoints;  
 
    	
	// Class constructor and destructor
	Tracker();
	~Tracker();

 public:

	static Tracker* instance(); // Function to return pointer to the class instance

	// Function to simulate a track through the detector, then return data for plane hits.
	// Uses MC method to reject hits going outside of the detector
	
	

	float DCA(float ,float ,float ,float ,float ,float);
	
	//TODO rewrite like DCA_simple struc. 
	float DCAHit(std::vector<float>, float, float, float, float, float, float, float, bool);

	float DCA_simple(float, float);

	DCAData DCAHit_simple(std::vector<float>, float, float, bool);

	float GetIdealPoint(int, float, float, std::vector<float>);

	vector<float> GetResiduals(std::vector<float>, std::vector<float>);

	vector<float> GetResiduals_simple(std::vector<float>, std::ofstream&);

	LineData MC(float, std::ofstream&, std::ofstream&, std::ofstream&, std::ofstream&, bool); 

	void setGeometry(std::ofstream&, bool); //Geometry of detector arrangement 

	void misalign(std::ofstream&, bool); // MC misalignment of detectors 
		
    void write_constraint_file(std::ofstream&, std::ofstream&, bool);  // Writes a constraint file for use with PEDE. 

	void set_uniform_file(std::string); // Set filename for uniform random numbers [randomIntGenerator.py - see https://github.com/glukicov/alignTrack]

	void set_gaussian_file(std::string); // Set filename for Gaussian random numbers



	//
	// Setter methods
	//

	void setTrackNumber(int tracks){
		trackNumber=tracks; 
	}


	//
	// Getter methods
	//


	int getLayer(int i) {
		return layer[i];
	}
	
	float getProjectionX(int i) {
		return projectionX[i];
	}

	float getSdevX(int i) {
		return sdevX[i];
	}

	// float getX_idealPoints(int i){
	// 	return x_idealPoints[i];
	// }


	int getStrawN(){
		return strawN;
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