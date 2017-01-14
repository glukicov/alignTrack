#include "mptest1.h"

using namespace std;

// Parameters for detector plane properties 
const int plane_count = 100;
int track_count = 10000;

float plane_x_begin = 10.0;
float plane_x_sep = 10.0;
float plane_thickness = 2.0;
float plane_height = 100.0;
float plane_eff = 0.90;
float meas_sigma = 0.0150;
	
// Standard deviations of distributions of plane displacement and drift velocit.
float displ_sigma = 0.1;
float drift_sigma = 0.02;

float plane_pos_devs[plane_count];
float drift_vel_devs[plane_count];

float true_plane_effs[plane_count];
float true_meas_sigmas[plane_count];
float hit_sigmas[plane_count];

float x_hits[plane_count];
float y_hits[plane_count];
float y_drifts[plane_count];

int i_hits[plane_count];

// Random number generators and distributions, for uniform and gaussian distribution
default_random_engine uniform_generator;
uniform_real_distribution<float> uniform_dist(0.0,1.0);

default_random_engine gaus_generator;
normal_distribution<float> gaus_dist(0.0, 1.0);


vector<float> genlin(int& hit_count, vector<float>& x_hits, vector<float>& hit_sigmas, vector<float>& y_drifts, vector<int>& i_hits) {

	vector<float> y_hits;

	float y_intercept = (0.5 * plane_height) + (0.1 * plane_height * (uniform_dist(uniform_generator) - 0.5));
	float gradient = ((uniform_dist(uniform_generator) - 0.5) * plane_height) / ((plane_count - 1) * plane_x_sep);
	   
	for (int i=0; i < plane_count; i++) {
		float x = plane_x_begin + (i * plane_x_sep);

		if (uniform_dist(uniform_generator) < true_plane_effs[i]) {
			
			float true_y = y_intercept + (gradient * x);
			float biased_y = true_y - plane_pos_devs[i];
			int wire_num = int(1 + (biased_y / 4));
			if(wire_num <= 0 || wire_num > 25) break;		

			x_hits.push_back(x);
			i_hits.push_back(i);

			float smear_y = true_meas_sigmas[i] * gaus_dist(gaus_generator); 
			float y_dvds = 0.0;
				
			y_hits.push_back(smear_y + biased_y + y_dvds);
			
			float y_wire = (float) (wire_num * 4.0) - 2.0;
			y_drifts.push_back(biased_y - y_wire);
			
			y_dvds = y_drifts[hit_count] * drift_vel_devs[i];

			y_hits[hit_count] = biased_y + smear_y - y_dvds;
			hit_sigmas.push_back(true_meas_sigmas[i]);

			hit_count++;

		}
 
	} 

	return y_hits;

}



void mptest1() {

	// Name and properties of binary output file
	string binary_file_name = "mp2tst.bin";
	bool as_binary = true;
	bool write_zero = false;

	Mille m (binary_file_name.c_str(), as_binary, write_zero);

	// Names of constraint and steering files
	string constraint_file_name = "mp2con.txt";
	string steering_file_name = "mp2str.txt";

	cout << "" << endl;
	cout << "Generating test data for mp II..." << endl;
	cout << "" << endl;

	// Open file streams for constraint and steering files, overwriting any original files
	ofstream constraint_file(constraint_file_name);
	ofstream steering_file(steering_file_name);

	// Iterate across planes, setting plane efficiencies, resolutions, deviations in position and drift velocity.
	for (int i=0; i<plane_count; i++) {
		true_plane_effs[i] = plane_eff;
		true_meas_sigmas[i] = meas_sigma;

		plane_pos_devs[i] = displ_sigma * gaus_dist(gaus_generator);
		drift_vel_devs[i] = drift_sigma * gaus_dist(gaus_generator);

	}

	// To constrain measurement
	plane_pos_devs[9] = 0.0;
	plane_pos_devs[89] = 0.0;
				
	// Set bad plane, with poor efficiency and resolution.
	true_plane_effs[6] = 0.1;
	true_meas_sigmas[6] = 0.0400;
	

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

	
	int all_hit_count = 0;
	int all_record_count = 0;

	for (int i=0; i<track_count; i++) {

		int hit_count = 0;
		vector<float> x_hits;
		vector<float> y_drifts;
		vector<float> hit_sigmas;
		vector<int> i_hits;

		vector<float> y_hits = genlin(hit_count, x_hits, y_drifts, hit_sigmas, i_hits);
		
		for (int j=0; j<hit_count; j++) {
			float local_derivs[2] {1.0, x_hits[j]};
			float global_derivs[2] {1.0, y_drifts[j]};
			int labels[2] {10 + (2 * (i_hits[j] + 1)), 500 + i_hits[j] + 1};
			
			// cout << local_derivs[0] << ", " << local_derivs[1] << endl;
			// cout << global_derivs[0] << ", " << global_derivs[1] << endl;
			// cout << labels[0] << ", " << labels[1] << endl;

			m.mille(2, local_derivs, 2, global_derivs, labels, y_hits[j], hit_sigmas[j]);

			all_hit_count++;
		}

		m.end();

		all_record_count++;
	}

	cout << " " << endl;
	cout << " " << endl;
	cout << track_count << " tracks generated with " << all_hit_count << " hits." << endl;
	cout << all_record_count << " records written." << endl;
	cout << " " << endl; 

}

int main() {

	mptest1();

	return 0;

}
									   
