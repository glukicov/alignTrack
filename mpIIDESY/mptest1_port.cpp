/**
   mptest1_port.cpp
   
   Purpose: Port of mptest1.f90 Mille test program. Simulates a plane drift chamber, with variable plane offset and drift velocity. Writes global and local derivatives to a binary file, and writes appropriate steering and constraint files. This source file defines the main function, and a structure used to define a TTree structure, allowing parameter values to be output.
  
   @author John Smeaton
   @version 03/02/2017
 
 */

#include "mptest1_port.h"

using namespace std;

/**
   Structure used to write data to TTree, with variable for each branch in tree.
 */
struct Parameter_data {
	int label; /** Parameter label, ranging from 12-210 for plane displacements, and 501-600 for velocity deviations.*/
	int fitType; /** Index for type of fit, with this defined as 0 for true parameter values */ 
	float paramValue; /** Value of parameter */
	float paramError; /** Error on parameter value (zero for true parameters) */
};


/**
   Main function for test data generation. Simulates tracks passing through plane detector, then writing data to a binary file, and generating appropriate constraint and steering files so plane displacements and drift velocity deviations may be fitted using the pede executable
 */
int main(int argc, char* argv[]) {

	// Minimum level of log entries to show - set to preference
	Logger::Instance()->setLogLevel(Logger::NOTE);

	// Allow logger to throw critical errors
	Logger::Instance()->enableCriticalErrorThrow();

	// Write program header
	Logger::Instance()->write(Logger::NOTE, "");
	Logger::Instance()->write(Logger::NOTE, "********************************************");
	Logger::Instance()->write(Logger::NOTE, "*                 MPTEST 1                 *");
	Logger::Instance()->write(Logger::NOTE, "********************************************");
	Logger::Instance()->write(Logger::NOTE, "");
	Logger::Instance()->write(Logger::NOTE, "    _____________________________  \\  /");
	Logger::Instance()->write(Logger::NOTE, "   {_|_|_|_|_|_|_|_|_|_|_|_|_|_|_( ͡° ͜ʖ ͡°) ");
    Logger::Instance()->write(Logger::NOTE, "    /\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/");
	Logger::Instance()->write(Logger::NOTE, ""); 


	try {

		// Check if correct number of arguments specified, logging error if not
		if (argc > 3) {
			Logger::Instance()->write(Logger::ERROR, "Too many arguments - please specify filename for file of uniform random numbers, followed by filename for file of gaussian random numbers. ");
		} else if (argc < 3) {
			Logger::Instance()->write(Logger::ERROR, "Too few arguments - please specify filename for file of uniform random numbers, followed by filename for file of gaussian random numbers. ");
		} else {

			// Set filenames to read random numbers from, using arguments. Log error if files don't exist.
			try {
				Detector::instance()->set_uniform_file(argv[1]);
				Detector::instance()->set_gaussian_file(argv[2]);
			} catch (ios_base::failure& e) {
				Logger::Instance()->write(Logger::ERROR, "Filestream exception caught: " + string(e.what()) + "\nPlease ensure valid filenames are specified!");
			}	
		}

		// Try to set plane properties, logging error if exception is caught.
		try {
			Detector::instance()->set_plane_properties();
		} catch (invalid_argument& e) {
			Logger::Instance()->write(Logger::ERROR, "Invalid argument error! " + string(e.what()));
		}
	}
	
	// If critical error occurs, then terminate programme.
	catch (CriticalError& e) {
		cerr << "Terminating..." << endl;
		return 1;
	}


	// Name and properties of binary output file
	string binary_file_name = "mp2tst1_c.bin";
	bool as_binary = true;
	bool write_zero = false;

	// Create instance of class for writing binary data.
	Mille m (binary_file_name.c_str(), as_binary, write_zero);

	// Names of constraint and steering files
	string constraint_file_name = "mp2test1con_c.txt";
	string steering_file_name = "mp2test1str_c.txt";
	string true_params_file_name = "mp2test1_true_params_c.txt";

	// Write to log
	Logger::Instance()->write(Logger::INFO, "");
	Logger::Instance()->write(Logger::INFO, "Generating test parameters...");

	// Open file streams for constraint and steering files, overwriting any original files
	ofstream constraint_file(constraint_file_name);
	ofstream steering_file(steering_file_name);
	ofstream true_params_file(true_params_file_name);

	// Initialise structure to hold TTree data
	Parameter_data true_params;

	// Create file, tree for parameters
	TFile f ("mptest1_parameters_c.root", "recreate");
	TTree t ("paramTree", "Tree to contain true, fitted parameter values");

	// Set up tree branches
	t.Branch("fitType", &true_params.fitType, "fitType/F");
	t.Branch("label", &true_params.label, "label/I");
	t.Branch("paramValue", &true_params.paramValue, "paramValue/F");
	t.Branch("paramError", &true_params.paramError, "paramError/F");

	true_params_file << endl; // Insert blank line

	// Print plane labels, and plane displacements to file
	for (int i=0; i<Detector::instance()->get_plane_count(); i++) { 
		true_params_file << 10 + (2 * (i + 1)) << " " << -Detector::instance()->get_plane_pos_y_devs()[i] << endl;
		
		// Add parameters, with labels, to TTree. Should fitType = 0 denoting true parameter values.
		true_params.fitType = 0;
		true_params.paramError = 0;
		true_params.label = 10 + (2 * (i + 1));
		true_params.paramValue = -Detector::instance()->get_plane_pos_y_devs()[i];
		t.Fill();
	}

		
	// Print drift velocity labels, and drift velocity displacements to file.
	for (int i=0; i<Detector::instance()->get_plane_count(); i++) { 
		true_params_file << 500 + i + 1 << " " << -Detector::instance()->get_drift_vel_devs()[i] << endl; 

		// Add parameters, with labels, to TTree. Should fitType = 0 denoting true parameter values.
		true_params.fitType = 0;
		true_params.label = 500 + i + 1;
		true_params.paramValue = -Detector::instance()->get_drift_vel_devs()[i];
		true_params.paramError = 0;
		t.Fill();
	}

	// Write values to TTree.
	t.Write();

	// Write constraint file, and initial parameter file, for use with pede
	Detector::instance()->write_constraint_file(constraint_file);

	// Check steering file is open, then write
	if (steering_file.is_open()) {

		Logger::Instance()->write(Logger::INFO, "Writing steering file...");

		steering_file << "*            Default test steering file" << endl
					  << "fortranfiles ! following bin files are fortran" << endl
					  << "mp2test1con_c.txt   ! constraints text file " << endl
					  << "Cfiles       ! following bin files are Cfiles" << endl
					  << "mp2tst1_c.bin   ! binary data file" << endl
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


	Logger::Instance()->write(Logger::INFO, "Generating tracks, writing hit data to binary...");
	
	// Set up counters for number of hits, and tracks
	int all_hit_count = 0;
	int all_record_count = 0;

	// Iterate over number of tracks
	for (int i=0; i<Detector::instance()->get_track_count(); i++) {

		// Simulate track, and get data
		LineData generated_line = Detector::instance()->gen_lin();

		// Get local gradient of this track, from recorded hit distances
		float local_gradient = generated_line.gradient; // (generated_line.y_hits.back() - generated_line.y_hits[0]) / (generated_line.x_hits.back() - generated_line.x_hits[0]);

		
		// Iterate over hits in detector
		for (int j=0; j<generated_line.hit_count; j++) {
						
			// Create arrays of local and global derivatives.
			float local_derivs[2] {1.0, generated_line.x_hits[j]};
			float global_derivs[2] {1.0, generated_line.y_drifts[j]};

			// Labels for plane displacements, and velcity deviation. 
			int labels[2] {10 + (2 * (generated_line.i_hits[j] + 1)), 500 + generated_line.i_hits[j] + 1};
			
			// Number of local, global derivatives for use in mille.
			int local_deriv_count = 2;
			int global_deriv_count = 2;

			// Write to binary file.
			m.mille(local_deriv_count, local_derivs, global_deriv_count, global_derivs, labels, generated_line.y_hits[j], generated_line.hit_sigmas[j]);

			// For debugging
			if (i==0) {
				true_params_file << "Hit " << j << endl 
								 << "Local: " << local_derivs[0] << " " << local_derivs[1] << endl
								 << "Global: " << global_derivs[0] << " " << global_derivs[1] << endl
								 << "Label: " << labels[0] << " " << labels[1] << endl
								 << "Hit: " << generated_line.y_hits[j] << endl
								 << "Sigma: " << generated_line.hit_sigmas[j] << endl << endl;

			}

			all_hit_count++; // Increment total number of recorded hits
		}

		m.end(); // End of this record

		all_record_count++; // Increment total number of records
	}

	Logger::Instance()->write(Logger::INFO, "");
	Logger::Instance()->write(Logger::INFO, to_string(Detector::instance()->get_track_count()) + " tracks generated with " + to_string(all_hit_count) + " hits.");
	Logger::Instance()->write(Logger::INFO, to_string(all_record_count) + " records written.");
	Logger::Instance()->write(Logger::INFO, "");

	// Close text files
	constraint_file.close();
	steering_file.close();
	true_params_file.close();

	// Terminate program.
	return 0;
}

									   
