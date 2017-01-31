#include "mptest1_port_detector.h"

using namespace std;

Detector* Detector::s_instance = NULL; 

Detector::Detector() {

	seed = 453032763; 
	rand_gen = new TRandom3(seed);

	s_instance = 0;

}

Detector* Detector::instance() {

	if (s_instance == NULL) s_instance = new Detector();
	return s_instance;
}


LineData Detector::gen_lin() {

	// Set up new container for track data, with hit count set to zero`
	LineData line;
	line.hit_count = 0;

	// Generate random values of track intercept and gradient.
	float y_intercept = (0.5 * PLANE_HEIGHT) + (0.1 * PLANE_HEIGHT * (rand_gen->Rndm() - 0.5)); 
	float gradient = ((rand_gen->Rndm() - 0.5) * PLANE_HEIGHT) / ((PLANE_COUNT - 1) * PLANE_X_SEP);

	// Iterate across planes
	for (int i=0; i<PLANE_COUNT; i++) {

		// Get position of plane
		float x = PLANE_X_BEGIN + (i * PLANE_X_SEP);

		// Check if hit is registered, due to limited plane efficiency.
		if (rand_gen->Rndm() < true_plane_effs[i]) {

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
			float y_smear = true_meas_sigmas[i] * rand_gen->Gaus(1,0);

			// Calculate y-position of hit wire, then calculate drift distance.
			float y_wire = (float) (wire_num * 4.0) - 2.0;
			line.y_drifts.push_back(y_biased - y_wire);

			// Calculate deviation in recorded y position due to drift velocity deviation
			float y_dvd = line.y_drifts[line.hit_count] * drift_vel_devs[i];

			// Calculate recorded hit y-position, with deviations due to smearing and drift velocity deviation.
			line.y_hits.push_back(y_biased + y_smear + y_dvd);

			// Record uncertainty in hit y-position, and increment number of hits.
			line.hit_sigmas.push_back(true_meas_sigmas[i]);
			line.hit_count++;

		}
	}
	return line;
}


void Detector::set_seed(int new_seed) {
	seed = new_seed;
	rand_gen = new TRandom3(seed);
	cout << "New Seed Set: " << seed << endl;
}


void Detector::set_plane_properties() {

	for (int i=0; i<PLANE_COUNT; i++) {

		true_plane_effs.push_back(PLANE_EFF);
		true_meas_sigmas.push_back(MEAS_SIGMA);

		plane_pos_devs.push_back(DISPL_SIGMA * rand_gen->Gaus(0,1));
		drift_vel_devs.push_back(DRIFT_SIGMA * rand_gen->Gaus(0,1));
	}

	true_plane_effs[6] = 0.1;
	true_meas_sigmas[6] = 0.0400;

	plane_pos_devs[9] = 0.0;
	plane_pos_devs[89] = 0.0;
}


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

