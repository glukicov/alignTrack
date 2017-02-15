#ifndef MPTEST2DETECTOR_H
#define MPTEST2DETECTOR_H

// This header file contains definitions of constant variables used in 
// Detector class, as well as function declarations, and 
//definitions of some inline functions. 

#include "random_buffer.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath> //math class`

/**
   Structure to contain data of a generated track, with the number of hits, their positions, the uncertainty in the positions, and the plane number hit.
*/
struct LineData {
	int hit_count; /** Number of hits in detector */
	std::vector<float> x_hits; /** X-positions of hits in detector */
	std::vector<float> y_hits; /** Y-positions of hits in detector */
	std::vector<float> hit_sigmas; /** Resolution for hits in detector */
	std::vector<int> i_hits; /** Number for plane struck in detector hits, with the plane numbers starting at zero, and increasing by one for each adjacent plane */
};

/**
   Singleton class to represent a detector made up of an array of planar drift chambers, simulating the passage of linear tracks through the detector.
 */
class Detector {

 private:

	static Detector* s_instance; // Pointer to instance of class

	const int track_count = 10000; /** Number of tracks to be simulated passing through detector */

	///initialsing physics varibles
	static const int detectorN = 10; //number of detector layers
	static const int layerN = 14; //number of measurement layers   //XXX why 14 is to do with stereo-angles for 1,4,7,10
	static const int moduleXN = 10; //number of modules in x direction
	static const int moduleYN = 5; //number of modules in y direction   //total 50 modules moduleXN * moduleYN = 50
	static const int moduleXYN = moduleXN*moduleYN;

	static const int modulesTotalN=detectorN*moduleYN*moduleXN; //total number of modules
	//  define detector geometry
	float arcLength_Plane1= 10.0; // arclength of first plane
	float planeDistance= 10.0; // distance between planes //cm / Pede works in cm
	float width= 0.02; //thickness/width of plane (X0)
	float offset=  0.5;  // offset of stereo modules
	float stereoTheta=0.08727;  // stereo angle  // radians (5 deg = 0.087.. rad)  
	float layerSize= 20.0; //size of layers  //cm 
	float resolution =0.002;  // <resolution  // 20um = 0.002 cm 
	float scatterError = 0; // multiple scattering error
	
	
    std::vector<int> layer; 
    std::vector<float> projectionX; //projection of measurent direction in (X)
    std::vector<float> projectionY; //projection of measurent direction in (Y)
    /// TODO rewrite those as vectors for detector class 
    float sdevX[moduleXYN*detectorN];// shift in x (alignment parameter)
    float sdevY[moduleXYN*detectorN] ; //shift in y (alignment GLOBAL parameter)
    float arcLength[layerN];  // arc length
    float resolutionLayer[layerN];   //resolution
	
	// Class constructor and destructor
	Detector();
	~Detector();

 public:
 	
    

/// XXX need more methods? 

	static Detector* instance(); // Function to return pointer to class instance

	LineData genlin2(); // Function to simulate a track through the detector, then return data for plane hits. 

	void setGeometry(); //Geometry of detecor arrangement 

	void misalign(); // MC misalignment of detecors 
		
    void write_constraint_file(std::ofstream&); // Writes a constraint file to the provided file stream, for use with pede. 

	void set_uniform_file(std::string); // Set filename for uniform random numbers

	void set_gaussian_file(std::string); // Set filename for gaussian random numbers


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

	/**
	   Get the number of tracks to be generated in the detector.

	   @return Number of tracks.
	 */
	int get_track_count() {
		return track_count;
	}

	int getWidth() {
		return width;
	}

	int getLayerN() {
		return layerN;
	}

	int getModuleXYN() {
		return moduleXYN;
	}

	int getModulesTotalN() {
		return modulesTotalN;
	}

	
	// XXX more getter methods!!!
	


};


#endif