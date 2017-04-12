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
	std::vector<double> x_hits; /** X-positions of hits in detector */
	std::vector<double> y_hits; /** Y-positions of hits in detector */
	std::vector<double> x_true; /** X-positions of hits in detector */
	std::vector<double> y_true; /** Y-positions of hits in detector */
	std::vector<double> x_det; /** X-positions of hits in detector */
	std::vector<double> y_det; /** Y-positions of hits in detector */
	std::vector<double> x_mis; /** X-positions of hits in detector */
	std::vector<double> y_mis; /** Y-positions of hits in detector */
	std::vector<double> hit_sigmas; /** Resolution for hits in detector */
	std::vector<int> i_hits; /** Number for plane struck in detector hits, with the plane numbers starting at 1, and increasing by one for each adjacent plane */
};

/**
   Singleton class to represent a detector made up of an array of planar drift chambers, simulating the passage of linear tracks through the detector.
 */
class Detector {

 private:

	static Detector* s_instance; // Pointer to instance of class

	const int trackCount = 10000; /** Number of tracks to be simulated passing through detector */

	///initialsing physics varibles
	static const int detectorN = 10; //number of detector layers
	static const int layerN = 14; //number of measurement layers  [extra modules at 1, 4 ,7 10]
	static const int pixelXN = 10; //number of pixels in x direction
	static const int pixelYN = 5; //number of pixels in y direction   //total 50 pixels pixelXN * pixelYN = 50
	static const int pixelXYN = pixelXN*pixelYN;

	static const int pixelTotalN=detectorN*pixelYN*pixelXN; //total number of pixels extra modules at 1, 4 ,7 10 have no pixels]
	//  define detector geometry
	double startingDistancePlane1= 10.0; // arclength of first plane
	double planeDistance= 10.0; // distance between planes //cm / Pede works in cm
	double width= 0.02; //thickness/width of plane (X0)
	double offset=  0.5;  // offset of stereo pixels
	double stereoTheta=0.08727;  // stereo angle  // radians (5 deg = 0.087.. rad)  
	double layerSize= 20.0; //size of layers  //cm 
	double resolution =0.002;  // <resolution  // 20um = 0.002 cm 
	double scatterError = 0; // multiple scattering error
	
    std::vector<int> layer; 
    std::vector<double> projectionX; //projection of measurent direction in (X)
    std::vector<double> projectionY; //projection of measurent direction in (Y)

    std::vector<double> distance;  // arc length/distance between planes 
    std::vector<double> resolutionLayer;   //resolution
   
    double sdevX[detectorN][pixelYN][pixelXN];// shift in x (alignment parameter)
    double sdevY[detectorN][pixelYN][pixelXN] ; //shift in y (alignment GLOBAL parameter)
    
	
	// Class constructor and destructor
	Detector();
	~Detector();

 public:

	static Detector* instance(); // Function to return pointer to class instance

	LineData genlin2(std::ofstream&, std::ofstream&, bool); // Function to simulate a track through the detector, then return data for plane hits. 

	void setGeometry(std::ofstream&, bool); //Geometry of detecor arrangement 

	void misalign(std::ofstream&, bool); // MC misalignment of detecors 
		
    void write_constraint_file(std::ofstream&); // Writes a constraint file to the provided file stream, for use with pede. 

	void set_uniform_file(std::string); // Set filename for uniform random numbers

	void set_gaussian_file(std::string); // Set filename for gaussian random numbers


	//
	// Getter methods
	//

	std::vector<int> getLayer() {
		return layer;
	}
	
	std::vector<double> getProjectionX() {
		return projectionX;
	}

	std::vector<double> getProjectionY() {
		return projectionY;
	}

	/**
	   Get the number of tracks to be generated in the detector.

	   @return Number of tracks.
	 */
	int getTrackCount() {
		return trackCount;
	}

	int getWidth() {
		return width;
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

	int getScatterError() {
		return scatterError;
	}


};


#endif