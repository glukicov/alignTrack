/** 
	mptest1_port_detector.h

	Purpose: Simulate linear tracks passing through a plane drift chamber detector, with misaligned plane positions uncalibrated drift velocities, in order to generate the necessary data for the correct plane positions and drift velocities to be calculated using pede. This header file contains definitions of constant variables used in Detector class, as well as function declarations, and definitions of some inline functions. 

	@author John Smeaton
	@version 03/02/2017

 */

#ifndef MPTEST1DETECTOR_H
#define MPTEST1DETECTOR_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "random_buffer.h"

/**
   Structure to contain data of a generated track, with the number of hits, their positions, the uncertainty in the positions, and the plane number hit.
*/
struct LineData {
	int hit_count; /** Number of hits in detector */
	std::vector<float> x_hits; /** X-positions of hits in detector */
	std::vector<float> y_hits; /** Y-positions of hits in detector */
	std::vector<float> hit_sigmas; /** Resolution for hits in detector */
	std::vector<float> y_drifts; /** Drift distance from hit position to closest wire */
	std::vector<int> i_hits; /** Number for plane struck in detector hits, with the plane numbers starting at zero, and increasing by one for each adjacent plane */
};

/**
   Singleton class to represent a detector made up of an array of planar drift chambers, simulating the passage of linear tracks through the detector.
 */
class Detector {

 private:

	static Detector* s_instance; // Pointer to instance of class

	// Numbers of planes, tracks
	const int PLANE_COUNT = 100; /** Number of detector planes */
	const int TRACK_COUNT = 10000; /** Number of tracks to be simulated passing through detector */
	
	// Parameters for detector plane properties 
	const float PLANE_X_BEGIN = 10.0; /** Beginning x-position of detector */
	const float PLANE_X_SEP = 10.0; /** Separation distance between planes */
	const float PLANE_THICKNESS = 2.0; /** Thickness of planes (note: This is not used)*/
	const float PLANE_HEIGHT = 100.0; /** Height of planes (distance covered by planes in y-direction) */
	const float PLANE_EFF = 0.90; /** Default efficiency of planes (fraction of hits recorded to planes a track passes through) */
	const float MEAS_SIGMA = 0.0150; /** Default resolution of planes (smearing distance for y-position of hits) */

	const float DISPL_SIGMA = 0.1; /** Standard deviation of plane hit displacement distribution */
	const float DRIFT_SIGMA = 0.02; /** Standard deviation of plane drift velocity fractional deviation distribution */

	std::vector<float> plane_pos_devs; /** Vector of plane position deviations */
	std::vector<float> drift_vel_devs; /** Vector of plane drift velocity fractional deviations */

	std::vector<float> true_plane_effs; /** Vector of plane efficiencies */
	std::vector<float> true_meas_sigmas; /** Vector of plane resolutions */

	// Class constructor and destructor
	Detector();
	~Detector();

 public:

	static Detector* instance(); // Function to return pointer to class instance

	LineData gen_lin(); // Function to simulate a track through the detector, then return data for plane hits. 
		
	void set_plane_properties(); // Sets up plane position, velocity deviations, using random number generator. 

	void write_constraint_file(std::ofstream&); // Writes a constraint file to the provided file stream, for use with pede. 

	void set_uniform_file(std::string); // Set filename for uniform random numbers

	void set_gaussian_file(std::string); // Set filename for gaussian random numbers


	//
	// Getter methods
	//

	/**
	   Get the number of tracks to be generated in the detector.

	   @return Number of tracks.
	 */
	int get_track_count() {return TRACK_COUNT;}


	/**
	   Get the number of planes in the detector.

	   @return Detector plane count.
	 */
	int get_plane_count() {return PLANE_COUNT;}


	/**
	   Get the drift velocity deviations for the detector planes.

	   @return Vector of plane drift velocity deviations.
	*/
	std::vector<float> get_drift_vel_devs() {return drift_vel_devs;}


	/**
	   Get the plane position deviations from zero for the detector planes.

	   @return Vector of plane position deviations.
	 */
	std::vector<float> get_plane_pos_devs() {return plane_pos_devs;}

};


#endif
