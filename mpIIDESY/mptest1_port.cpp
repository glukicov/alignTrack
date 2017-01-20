/**
 * Port of mptest1.f90 Mille test program. Simulates a plane drift chamber, with variable
 * plane offset and drift velocity. Writes global and local derivatives to a binary file,
 * and writes appropriate steering and constraint files.
 *
 * John Smeaton 14/01/2017
 *
 **/

#include "mptest1_port.h"

using namespace std;

// Number of planes and tracks.
const int plane_count = 100;
int track_count = 10000;

// Parameters for detector plane properties 
float plane_x_begin = 10.0; // Beginning of detector
float plane_x_sep = 10.0; // Separation between planes
float plane_thickness = 2.0; // Thickness of planes
float plane_height = 100.0; // Height of planes
float plane_eff = 0.90; // Default efficiency of planes
float meas_sigma = 0.0150; // Default resolution of planes
	
// Standard deviations of distributions of plane displacement and drift velocity.
float displ_sigma = 0.1;
float drift_sigma = 0.02;

// Arrays of the plane displacement and velocity deviations.
float plane_pos_devs[plane_count];
float drift_vel_devs[plane_count];

// Arrays of the plane efficiencies and resolutions
float true_plane_effs[plane_count];
float true_meas_sigmas[plane_count];

// Structure to contain data of a generated line, with the number of hits, their positions, the uncertainty in the positions, and the plane number hit.
struct Line_data {
	int hit_count;
	vector<float> x_hits;
	vector<float> y_hits;
	vector<float> hit_sigmas;
	vector<float> y_drifts;
	vector<int> i_hits;
};

// Random number generators and distributions, for uniform and gaussian distribution

random_device uniform_device;
random_device gaus_device;

uniform_real_distribution<float> uniform_dist(0.0, 1.0);
normal_distribution<float> gaus_dist(0.0, 1.0);

// Function to simulate a linear track through the detector, returning data about detector hits.
Line_data genlin() {

	seed_seq uniform_seeds{uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device()}; 
	seed_seq gaus_seeds{gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device()}; 

	mt19937 uniform_generator(uniform_seeds);
	mt19937 gaus_generator(gaus_seeds);


	// Set up new container for track data, with hit count set to zero
	Line_data line;
	line.hit_count = 0;

	// Generate random values of track intercept and gradient.
	float y_intercept = (0.5 * plane_height) + (0.1 * plane_height * (uniform_dist(uniform_generator) - 0.5));
	float gradient = ((uniform_dist(uniform_generator) - 0.5) * plane_height) / ((plane_count - 1) * plane_x_sep);
	   
	// Iterate across planes
	for (int i=0; i < plane_count; i++) {

		// Get position of plane
		float x = plane_x_begin + (i * plane_x_sep);

		// Check if hit is registered, due to limited plane efficiency.
		if (uniform_dist(uniform_generator) < true_plane_effs[i]) {
			
			// Calculate true value of y where line intercects plane, and biased value where hit is recorded, due to plane displacement
			float true_y = y_intercept + (gradient * x);
			float biased_y = true_y - plane_pos_devs[i];
			
			// Calculate number of struck wire. Do not continue simulating this track if it passes outside range of wire values.
			int wire_num = int(1 + (biased_y / 4));
			if(wire_num <= 0 || wire_num > 25) break;		

			// Record x-position, and plane number, of hit.
			line.x_hits.push_back(x);
			line.i_hits.push_back(i);

			// Calculate smear value from detector resolution.
			float smear_y = true_meas_sigmas[i] * gaus_dist(gaus_generator); 
				
			// Calculate y-position of hit wire, then calculate drift distance.
			float y_wire = (float) (wire_num * 4.0) - 2.0;
			line.y_drifts.push_back(biased_y - y_wire);			

			// Calculate deviation in recorded y position due to drift velocity deviation
			float y_dvds = line.y_drifts[line.hit_count] * drift_vel_devs[i];

			// Calculate recorded hit y-position, with deviations due to smearing and drift velocity deviation.
			line.y_hits.push_back(biased_y + smear_y - y_dvds);

			// Record uncertainty in hit y-position, and increment number of hits.
			line.hit_sigmas.push_back(true_meas_sigmas[i]);
			line.hit_count++;

		}
 
	} 

	return line; // Return data from simulated track

}



int main() {

	seed_seq uniform_seeds{uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device()}; 
	seed_seq gaus_seeds{gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device()}; 

	mt19937 uniform_generator(uniform_seeds);
	mt19937 gaus_generator(gaus_seeds);

	// Name and properties of binary output file
	string binary_file_name = "mp2tst.bin_c";
	bool as_binary = true;
	bool write_zero = false;

	// Create instance of class for writing binary data.
	Mille m (binary_file_name.c_str(), as_binary, write_zero);

	// Names of constraint and steering files
	string constraint_file_name = "mp2con.txt";
	string steering_file_name = "mp2str.txt";
	string true_params_file_name = "mp2test1_true_params.txt";

	cout << "" << endl;
	cout << "Generating test data for mp II..." << endl;
	cout << "" << endl;

	// Open file streams for constraint and steering files, overwriting any original files
	ofstream constraint_file(constraint_file_name);
	ofstream steering_file(steering_file_name);
	ofstream true_params_file(true_params_file_name);

	// Iterate across planes, setting plane efficiencies, resolutions, deviations in position and drift velocity.
	for (int i=0; i<plane_count; i++) {
		true_plane_effs[i] = plane_eff;
		true_meas_sigmas[i] = meas_sigma;

		plane_pos_devs[i] = displ_sigma * gaus_dist(gaus_generator);
		drift_vel_devs[i] = drift_sigma * gaus_dist(gaus_generator);

	}

	// Print plane labels, and plane displacements to file
	for (int i=0; i<plane_count; i++) 
		true_params_file << 10 + (2 * (i + 1)) << " " << -plane_pos_devs[i] << endl;

	true_params_file << endl; // Insert blank line

	// Print drift velocity labels, and drift velocity displacements to file.
	for (int i=0; i<plane_count; i++) 
		true_params_file << 500 + i + 1 << " " << -drift_vel_devs[i] << endl; 


	// To constrain measurement
	plane_pos_devs[9] = 0.0;
	plane_pos_devs[89] = 0.0;
				
	// Set bad plane, with poor efficiency and resolution.
	true_plane_effs[6] = 0.1;
	true_meas_sigmas[6] = 0.0400;
	

	// Check steering file is open, then write
	if (steering_file.is_open()) {

		cout << "" << endl;
		cout << "Writing Steering File" << endl;
		cout << "" << endl;


		steering_file << "*            Default test steering file" << endl
					  << "fortranfiles ! following bin files are fortran" << endl
					  << "mp2con.txt   ! constraints text file " << endl
					  << "Cfiles       ! following bin files are Cfiles" << endl
					  << "mp2tst.bin   ! binary data file" << endl
					  << "*hugecut 50.0     !cut factor in iteration 0" << endl
					  << "*chisqcut 1.0 1.0 ! cut factor in iterations 1 and 2" << endl
					  << "*entries  10 ! lower limit on number of entries/parameter" << endl
					  << "*pairentries 10 ! lower limit on number of parameter pairs" << endl
					  << "                ! (not yet!)" << endl
					  << "*printrecord   1  2      ! debug printout for records" << endl
					  << "*printrecord  -1 -1      ! debug printout for bad data records" << endl
					  << "*outlierdownweighting  2 ! number of internal iterations (> 1)" << endl
					  << "*dwfractioncut      0.2  ! 0 < value < 0.5" << endl
					  << "*presigma           0.01 ! default value for presigma" << endl
					  << "*regularisation 1.0      ! regularisation factor" << endl
					  << "*regularisation 1.0 0.01 ! regularisation factor, pre-sigma" << endl
					  << " " << endl
					  << "*bandwidth 0         ! width of precond. band matrix" << endl
					  << "method diagonalization 3 0.001 ! diagonalization      " << endl
					  << "method fullMINRES       3 0.01 ! minimal residual     " << endl
					  << "method sparseMINRES     3 0.01 ! minimal residual     " << endl
					  << "*mrestol      1.0D-8          ! epsilon for MINRES" << endl
					  << "method inversion       3 0.001 ! Gauss matrix inversion" << endl
					  << "* last method is applied" << endl
					  << "*matiter      3  ! recalculate matrix in iterations" << endl
					  << " " << endl
					  << "end ! optional for end-of-data" << endl;
 
	}	

	// Check constraints file is open, then write. [Note - Don't yet understand these]
	if (constraint_file.is_open()) {
		
		cout << "" << endl;
		cout << "Writing Constraint File" << endl;
		cout << "" << endl;

		constraint_file << "Constraint 0.0" << endl;
		for (int i=0; i<plane_count; i++) {
			int labelt = 10 + (i + 1) * 2;
			constraint_file << labelt << " " << fixed << setprecision(7) << 1.0 << endl;
		}


		float d_bar = 0.5 * (plane_count - 1) * plane_x_sep; 
		float x_bar = plane_x_begin + (0.5 * (plane_count - 1) * plane_x_sep);
		constraint_file << "Constraint 0.0" << endl;
		for (int i=0; i<plane_count; i++) {
			
			int labelt = 10 + (i + 1) * 2;
			
			float x = plane_x_begin + (i * plane_x_sep);
			float ww = (x - x_bar) / d_bar;

			constraint_file << labelt << " " << fixed << setprecision(7) << ww << endl;
		}
		
	}

	
	// Set up counters for number of hits, and tracks
	int all_hit_count = 0;
	int all_record_count = 0;

	// Iterate over number of tracks
	for (int i=0; i<track_count; i++) {

		// Simulate track, and get data
		Line_data generated_line = genlin();
		
		// Iterate over hits in detector
		for (int j=0; j<generated_line.hit_count; j++) {
			
			// Create arrays of local and global derivatives.
			float local_derivs[2] {1.0, generated_line.x_hits[j]};
			float global_derivs[2] {1.0, generated_line.y_drifts[j]};

			// Labels for plane displacement, and velcity deviation. 
			int labels[2] {10 + (2 * (generated_line.i_hits[j] + 1)), 500 + generated_line.i_hits[j] + 1};
			
			int local_deriv_count = 2;
			int global_deriv_count = 2;

			// Write to binary file.
			m.mille(2, local_derivs, 2, global_derivs, labels, generated_line.y_hits[j], generated_line.hit_sigmas[j]);

			all_hit_count++; // Increment total number of recorded hits

		}

		m.end(); // End of this record

		all_record_count++; // Increment total number of records
	}

	cout << " " << endl;
	cout << " " << endl;
	cout << track_count << " tracks generated with " << all_hit_count << " hits." << endl;
	cout << all_record_count << " records written." << endl;
	cout << " " << endl; 

	constraint_file.close();
	steering_file.close();
	true_params_file.close();

	// Terminate program.
	return 0;

}

									   
