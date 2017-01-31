#ifndef MPTEST1DETECTOR_H
#define MPTEST1DETECTOR_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "TRandom3.h"

// Structure to contain data of a generated line, with the number of hits, their positions, the uncertainty in the positions, and the plane number hit.
struct LineData {
	int hit_count;
	std::vector<float> x_hits;
	std::vector<float> y_hits;
	std::vector<float> hit_sigmas;
	std::vector<float> y_drifts;
	std::vector<int> i_hits;
};


class Detector {

 private:

	static Detector* s_instance;

	TRandom3* rand_gen;
	int seed;

	// Numbers of planes, tracks
	int PLANE_COUNT = 100;
	int TRACK_COUNT = 10000;
	
	// Parameters for detector plane properties 
	float PLANE_X_BEGIN = 10.0; // Beginning of detector
	float PLANE_X_SEP = 10.0; // Separation between planes
	float PLANE_THICKNESS = 2.0; // Thickness of planes
	float PLANE_HEIGHT = 100.0; // Height of planes
	float PLANE_EFF = 0.90; // Default efficiency of planes
	float MEAS_SIGMA = 0.0150; // Default resolution of planes

	// Standard deviations of distributions of plane displacement and drift velocity.
	float DISPL_SIGMA = 0.1;
	float DRIFT_SIGMA = 0.02;

	// Arrays of the plane displacement and velocity deviations.
	std::vector<float> plane_pos_devs;
	std::vector<float> drift_vel_devs;

	// Arrays of the plane efficiencies and resolutions
	std::vector<float> true_plane_effs;
	std::vector<float> true_meas_sigmas;

	Detector();

 public:

	static Detector* instance();
	LineData gen_lin();
	void set_seed(int new_seed);
	
	void set_plane_properties();
	void write_constraint_file(std::ofstream&);

	int get_track_count() {return TRACK_COUNT;}
	int get_plane_count() {return PLANE_COUNT;}
	int get_seed() {return seed;}
	std::vector<float> get_drift_vel_devs() {return drift_vel_devs;}
	std::vector<float> get_plane_pos_devs() {return plane_pos_devs;}
};



#endif
