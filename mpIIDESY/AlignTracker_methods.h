#ifndef ALIGNTRACKER_M
#define ALIGNTRACKER_M

/* This header file contains definitions of constant variables used in 
* method class, as well as function declarations, and definitions of functions.
*/  

#include "random_buffer.h" // courtesy of John Smeaton (UCL)

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

	static const int trackCount=2000; /** Number of tracks (i.e. records) to be simulated passing through detector */
 
	///initialising physics variables
	static const int detectorN = 10; //number of detector layers [independent layers]
	static const int layerN = 14; //number of measurement layers  [total number of layers]
	static const int pixelXN = 10 ; //number of measurement elements in x direction
	static const int pixelYN = 5; //number of measurement elements in y direction  
	static const int pixelXYN = pixelXN*pixelYN;

	static const float twoR=2.0; //For normalisation of uniform random numbers [0,1] : (MAX+RND)/(twoR*MAX)

	static const int pixelTotalN = detectorN*pixelYN*pixelXN; //total number of measurement elements 
	//  define detector geometry
	static const float startingDistancePlane1=10.0; // distance of first layer relative to the "beam"
	static const float planeDistance=10.0; // distance between planes [cm] 
	static const float width =0.02; //thickness/width of a plane (X0)
	static const float offset=0.5;  // offset of stereo pixels
	static const float stereoTheta=0.1309;  // stereo angle [rad]  // radians (7.5000 deg = 0.1309...rad)  
	static const float layerSize=20.0; //length of layers  [cm] 

	float resolution; //TODO decide if this is a constant  [value is assigned in the constructor]
	
		
	std::vector<int> layer; // record of layers that were hit
    std::vector<float> projectionX; //projection of measurement direction in (X)
    std::vector<float> projectionY; //projection of measurement direction in (Y)

    std::vector<float> distance;  // distance between planes [this is set in geometry]
    std::vector<float> resolutionLayer;   //resolution [to record vector of resolution per layer if not constant] //XXX [this is not used at the moment]
    
    float sdevX[detectorN][pixelYN][pixelXN];// shift in x due to the imposed misalignment (alignment parameter)
    float sdevY[detectorN][pixelYN][pixelXN] ; //shift in y due to the imposed misalignment (alignment GLOBAL parameter)
	
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
		
    void write_constraint_file(std::ofstream&); // Writes a constraint file for use with PEDE. 

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

	std::vector<float> getProjectionY() {
		return projectionY;
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

	int getPixelXYN() {
		return pixelXYN;
	}

	int getPixelTotalN() {
		return pixelTotalN;
	}

};


#endif