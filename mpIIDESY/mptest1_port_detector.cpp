/** 
	mptest1_port_detector.cpp

	Purpose: Simulate linear tracks passing through a plane drift chamber detector, with misaligned plane positions uncalibrated drift velocities, in order to generate the necessary data for the correct plane positions and drift velocities to be calculated using pede. This source file contains definitions of various functions in Detector class. 

	@author John Smeaton
	@version 03/02/2017

 */

#include "mptest1_port_detector.h"

using namespace std;

// Set up empty pointer for instance of class.
Detector* Detector::s_instance = NULL; 

/** 
	Empty constructor for detector class.
 */
Detector::Detector() {
}

/** 
	Empty destructor for detector class.
*/
Detector::~Detector() {
}


/**
   Get pointer to only instance of detector class, creating this instance if it doesn't already exist.

   @return Pointer to detector instance.
 */
Detector* Detector::instance() {

	// Create pointer to class instance if one doesn't exist already, then return that pointer.
	if (s_instance == NULL) s_instance = new Detector();
	return s_instance;
}


/**
   Simulate the passage of a randomly generated linear track through the detector, calculating properties of hits by the track on detector planes.
  
   @return LineData struct containing data about detector hits.
*/
LineData Detector::gen_lin() {

	// Set up new container for track data, with hit count set to zero`
	LineData line;
	line.hit_count = 0;

	// Generate random values of track intercept and gradient.
	float y_intercept = (0.5 * PLANE_HEIGHT) + (0.1 * PLANE_HEIGHT * (RandomBuffer::instance()->get_uniform_number() - 0.5)); 
	float gradient = ((RandomBuffer::instance()->get_uniform_number() - 0.5) * PLANE_HEIGHT) / ((PLANE_COUNT - 1) * PLANE_X_SEP);

	// Iterate across planes
	for (int i=0; i<PLANE_COUNT; i++) {

		// Get position of plane
		float x = PLANE_X_BEGIN + (i * PLANE_X_SEP);

		// Check if hit is registered, due to limited plane efficiency.
		if (RandomBuffer::instance()->get_uniform_number() < true_plane_effs[i]) {

			// Calculate true value of y where line intercects plane, and biased value where hit is recorded, due to plane displacement
			float y_true = y_intercept + (gradient * x);
			float y_biased = y_true - plane_pos_devs[i];

			// Calculate number of struck wire. Do not continue simulating this track if it passes outside range of wire values.
			int wire_num = int(1 + (y_biased / 4));
			if (wire_num <=0 || wire_num > 25) break;

			// Record x-position, and plane number, of hit.
			line.x_hits.push_back(x);
			line.i_hits.push_back(i);

			// Calculate smear value from detector resolution.			
			float y_smear = true_meas_sigmas[i] * RandomBuffer::instance()->get_gaussian_number();

			// Calculate y-position of hit wire, then calculate drift distance.
			float y_wire = (float) (wire_num * 4.0) - 2.0;
			line.y_drifts.push_back(y_biased - y_wire);

			// Calculate deviation in recorded y position due to drift velocity deviation
			float y_dvd = line.y_drifts[line.hit_count] * drift_vel_devs[i];

			// Calculate recorded hit y-position, with deviations due to smearing and drift velocity deviation.
			line.y_hits.push_back(y_biased + y_smear - y_dvd);

			// Record uncertainty in hit y-position, and increment number of hits.
			line.hit_sigmas.push_back(true_meas_sigmas[i]);
			line.hit_count++;

		}
	}
	return line;
}


/**
   Set up plane position deviations, and fractional drift velocity deviations, using random gaussian distributions.
 */
void Detector::set_plane_properties() {

	// Loop across all planes in detector
	for (int i=0; i<PLANE_COUNT; i++) {

		// Set up vectors of plane efficiencies and resolutions.
		true_plane_effs.push_back(PLANE_EFF);
		true_meas_sigmas.push_back(MEAS_SIGMA);

		// Set up random plane position deviations, and velocity deviations
		plane_pos_devs.push_back(DISPL_SIGMA * RandomBuffer::instance()->get_gaussian_number());
		drift_vel_devs.push_back(DRIFT_SIGMA * RandomBuffer::instance()->get_gaussian_number());
	}

	// Set up bad plane with low efficiency, and bad resolution.
	true_plane_effs[6] = 0.1;
	true_meas_sigmas[6] = 0.0400;

	// Set two planes to have no deviation in position, to constrain fit.
	plane_pos_devs[9] = 0.0;
	plane_pos_devs[89] = 0.0;

}

/**
   Write a constraint file to the supplied file-stream.

   TODO: Use abstractions instead of specifying ofstream?

   @param constraint_file Reference to ofstream to write constraint file to. 
 */
void Detector::write_constraint_file(ofstream& constraint_file) {

	// Check constraints file is open, then write. [Note - Don't yet understand these]
	if (constraint_file.is_open()) {
		
		cout << "Writing constraint file..." << endl;

		constraint_file << "Constraint 0.0" << endl;
		for (int i=0; i<PLANE_COUNT; i++) {
			int labelt = 10 + (i + 1) * 2;
			constraint_file << labelt << " " << fixed << setprecision(7) << 1.0 << endl;
		}


		float d_bar = 0.5 * (PLANE_COUNT - 1) * PLANE_X_SEP; 
		float x_bar = PLANE_X_BEGIN + (0.5 * (PLANE_COUNT - 1) * PLANE_X_SEP);
		constraint_file << "Constraint 0.0" << endl;
		for (int i=0; i<PLANE_COUNT; i++) {
			
			int labelt = 10 + (i + 1) * 2;
			
			float x = PLANE_X_BEGIN + (i * PLANE_X_SEP);
			float ww = (x - x_bar) / d_bar;

			constraint_file << labelt << " " << fixed << setprecision(7) << ww << endl;
		}	
	}
}


/**
   Set filename to read uniform random numbers from.

   @param uniform_filename String for filename of file of uniform random numbers.
 */
void Detector::set_uniform_file(string uniform_filename) {
	RandomBuffer::instance()->open_uniform_file(uniform_filename);
}


/**
   Set filename to read gaussian random numbers from.

   @param uniform_filename String for filename of file of gaussian random numbers.
 */
void Detector::set_gaussian_file(string gaussian_filename) {
	RandomBuffer::instance()->open_gaussian_file(gaussian_filename);
}

