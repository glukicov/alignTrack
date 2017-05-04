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
	std::vector<float> x_hits; /** X-positions of generated hits in detector */
	std::vector<float> y_hits; /** Y-positions of generated hits in detector */
	std::vector<float> x_true; /** X-positions of true hits in detector */
	std::vector<float> y_true; /** Y-positions of true hits in detector */
	std::vector<float> x_det; /** X-positions of recorded hits in detector */
	std::vector<float> y_det; /** Y-positions of recorded hits in detector */
	std::vector<float> x_mis; /** X-positions of misalignment in hits in detector */
	std::vector<float> y_mis; /** Y-positions of misalignment in hits in detector */
	std::vector<float> hit_sigmas; /** Resolution for hits in detector */
	std::vector<int> i_hits; /** Number for plane struck in detector hits, with the plane numbers starting at 1, and increasing by one for each adjacent plane */
};

/**
   Singleton class to represent a detector made up of an array of planar drift chambers, simulating the passage of linear tracks through the detector.
 */
class Tracker {

 private:

	static Tracker* s_instance; // Pointer to instance of class

	static const int trackCount=100; /** Number of tracks (i.e. records) to be simulated passing through detector */

	static const int beamPositionLength = 10; 

	float resolution; //TODO decide if this is a constant  [value is assigned in the constructor]

	float dispX;
	//--------------------//
 
	///initialising physics variables
	static const int detectorN = 6; //number of detector layers [independent layers]
	static const int layerN = 8; //number of measurement layers  [total number of layers] This will be layers 
	static const int pixelXN = 16 ; //number of measurement elements in x direction this will be straws 16

	static const float twoR=2.0; //For normalisation of uniform random numbers [0,1] : (MAX+RND)/(twoR*MAX)

	static const int pixelTotalN = detectorN*pixelXN; //total number of measurement elements 
	//  define detector geometry
	static const float startingDistancePlane1=50.0; // distance of first layer relative to the "beam" // [cm]
	static const float planeDistance=10.0; // distance between planes //[cm] 
	static const float width =0.02; //thickness/width of a plane (X0) // [cm]
	static const float offset=0.5;  // offset of stereo pixels [cm] 
	static const float stereoTheta=0.08727;  // stereo angle  // [rad] (5 deg = 0.087.. rad)
	//static const float stereoTheta=0.1309;  // stereo angle [rad]  // [rad] (7.5000 deg = 0.1309...rad)   // XXX
	static const float layerSize=20.0; //length of layers // [cm] 


	std::vector<int> layer; // record of layers that were hit
    std::vector<float> projectionX; //projection of measurement direction in (X)
    std::vector<float> projectionY; //projection of measurement direction in (Y)

    std::vector<float> distance;  // distance between planes [this is set in geometry]
    std::vector<float> resolutionLayer;   //resolution [to record vector of resolution per layer if not constant] //XXX [this is not used at the moment]
    
    float sdevX[detectorN][pixelXN];// shift in x due to the imposed misalignment (alignment parameter)
    float sdevY[detectorN][pixelXN] ; //shift in y due to the imposed misalignment (alignment GLOBAL parameter)
	
	// Class constructor and destructor
	Tracker();
	~Tracker();

 public:

	static Tracker* instance(); // Function to return pointer to the class instance

	// Function to simulate a track through the detector, then return data for plane hits.
	// Uses MC method to reject hits going outside of the detector
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

	int getPixelXN(){
		return pixelXN;
	}


	int getTrackCount() {
		return trackCount;
	}

	float getWidth() {
		return width;
	}

	float getResolution() {
		return resolution;
	}

	int getLayerN() {
		return layerN;
	}


	int getPixelTotalN() {
		return pixelTotalN;
	}

    float getDispX(){
		return dispX;
	}

};


#endif