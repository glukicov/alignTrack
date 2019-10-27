/*
*   Gleb Lukicov (g.lukicov@ucl.ac.uk) @ Fermilab
*   Created: 3 January 2018
*   Modified: 5 January 2018 
----------------------------------------------------------------
This programme uses MC methods to produce a .bin data file for the
PEDE routine, to align the tracking detector for the g-2
experiment.
Methods and functions are contained in AlignTracker_methods.cpp (Tracker Singleton class)

***Python and Bash Scripts:
*PlotGen.py - plot generated tracks and modules [after e.g. ./AlignTrack p 1]
*randomIntGenerator.py - generates random numbers
*plotdM.py - plot difference between PEDE and MC
*timeLaunch.py - measures time of MC execution for different number of tracks
*PlotFunc.py - plots relevant derivative functions analytically 

Python Iterative Test Scripts for establishing PEDE performance:
*Seed:
	PEDE_tests_seed.bash - starts tests for same number of tracks with varying seed.

*FoM:
	PEDE_tests.py - starts iterative comparison of PEDE vs MC for different number of tracks
	[using TrackerLaunch.py for plotting and ConcatenatePEDE.py for combining results]
*MisProfile.py - FoM for Chi2 

Grid Submission Scripts:
in LaunchPEDE_grid/ : 
*MC_pede/ - clean source files to tarball on a grid node for compilation
*send.sh - specifies the grid conditions and input files
*grid_submit.sh - specifies the job execution

Function Derivations in alignment.tex

ROOT Plotting macro rootlogon.C [loaded automatically e.g.: rootbrowse Tracker.root]

============================= Test Model v1.0 =======================================
*Simple 2D case (2D alignment as 2D translation and rotation along detector centre) with Circle Fit.
* No B-field, straight tracks, no MS, 100% efficiency.
*(!) NATURAL UNITS (PEDE): cm, rad, GeV [natural units will be omitted]

(!) Module, Straw, View, Layer labelling starts from 0
[+1 is put by hand as MP2 doesn't accept 0 as a label]

Terminology:
-detector=module (==the alignment object)
-Local Parameters: track parameters (e.g slope and intercept)
-Global Parameters: x and z position of a detector (and 𝛉 as a rotation through centre)

Gaussian (gaus) Smearing: Mean=0, RMS=1  (on hits via detector resolution)
Beam originates (in z direction) from z=0, straight tracks, no scattering,
    tracks are lines between z=beamStart and z=beamStop, with uniformly generated x1 and slope.
Resolution (i.e. tracker systematic) is implemented.
Misalignment is implemented "manually" for each module in +/- x/z direction, +/- 𝛉
[all layers in module are dependent (i.e. move together)]

Features [can be on or off - see AlignTracker.h]:
Hit rejection (>0.5 straw radius)
DCA rejection (>500 um) on the whole track
p-value cut
Truth LR information
Offsets [i.e. alignment "software" fix] as inputs

Constrains:
Pre-sigmas: -1 for "fixed" modules
Steering options: e.g. method inversion 5 0.001 [5 iterations, 0.01 dF convergence;
inversion method provides uncertainties on all alignment parameters]

                       			--Tracker Geometry [cm]--:
Spacing (in x) between straws in layer is 0.606 [strawSpacing]
Spacing between layers in same view is 0.515 [layerSpacing]
Spacing between U and V views (last and first layers in views) is 2.020 [viewSpaing]
Spacing between modules (last and first layers in modules) 13.735 [moduleSpacing]
Layers in a view are (e.g. U0 vs. U1) have an extra relative x displacement of 0.303 [layerDisplacement]
The staircaseXDisplacment of modules is not implemented yet XXX

For more information on the Tracker conventions see:
http://gm2-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=4375&filename=Orientation.pdf&version=1

=======================================================================================

#TODO Model type [physics] will be implemented via imodel = 0, 1 .....

#Verbosity level is implemented via command line command: n = normal, d = debug, p=plotting (./PlotGen.py) [e.g. './AlignTracker d']

Private (for now) Git repository: https://github.com/glukicov/alignTrack [further instructions are there]
Millepede II Manual and SC: http://www.desy.de/~kleinwrt/MP2/doc/html/index.html
*
*
**/
#include "AlignTracker.h"

using namespace std;
// for compiler version determination:
string ver_string(int a, int b, int c) {
	ostringstream ss;
	ss << a << '.' << b << '.' << c;
	return ss.str();
}

//**********************************************MAIN*************************************//
int main(int argc, char* argv[]) {
	//Starting a clock
	clock_t t_cpu; // CPU ticks for programme execution
	t_cpu = clock();
	auto t_start = chrono::high_resolution_clock::now(); // Wall clock ticks

	//Determining compiler used:
	string true_cxx =
#ifdef __clang__
	    "clang++";
#else
	    "g++";
#endif
	string true_cxx_ver =
#ifdef __clang__
	    ver_string(__clang_major__, __clang_minor__, __clang_patchlevel__);
#else
	    ver_string(__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

//------------------------------------------Global Variables Initialisation---------------------------------------------------------//

	int imodel = 0;  //Model type (see above) TODO implement this as an argument to main [for broken lines, MF, etc.]
	int setPrecision = 7; // precision (# decimal points) of printout for debug text files and cout
	string compareStr; //for debug vs. normal output as specified by the user
	int tracksInput; // number of tracks to generate as specified by the user
	bool debugBool = false; // './AlignTracker n' - for normal, of ./AlignTracker d' - for verbose debug output
	bool plotBool = false; // './AlignTracker p' - for plotting with PlotGen.py
	//Set up counters for hits and records (tracks)
	int hitsN = 0; // actually recorded (i.e. non-rejected hits)
	int recordN = 0; //records = tracks
	float scatterError; // multiple scattering error [calculated here and passed back to the Tracker class]
	// Those variable are set inside the hits loop, and used outside - hence global scope
	float residualsTrueSum_2 = 0.0; // squared sum of the true residuals (from measurements)
	float residualsReconSum_2 = 0.0; // squared sum of the recon. residuals (from measurements)
	float Chi2ReconActual = 0.0; // Reconstructed Chi2 (from measurements)
	int negDCA = 0; // counting negatively smeared DCAs
	const Color_t colourVector[] = {kMagenta, kOrange, kBlue, kGreen, kYellow, kRed, kGray, kBlack}; //8 colours for up to 8 modules
	gErrorIgnoreLevel = kWarning; // Display ROOT Warning and above messages [i.e. suppress info]
	gROOT->Macro("rootlogon.C");
	// Simple LR mapping for ROOT plots
	char nameLR[] = {'L', 'R'};
	int valueLR[] = { -1, 1};
	const char* boolYN[2] = {"No", "Yes"}; // True/False -> Yes/No mapping

	// use overloading to duplicate cout into the log file
	ofstream LOG; LOG.open("Tracker_log.txt"); MyStreamingHelper helper(LOG, cout);
	cout << fixed << setprecision(setPrecision); // set precision of standard screen output
	helper << fixed << setprecision(setPrecision); // set precision of log file output

	Logger::Instance()->setLogLevel(Logger::NOTE); //Tell the logger to only show message at INFO level or above
	Logger::Instance()->enableCriticalErrorThrow(); //Tell the logger to throw exceptions when ERROR messages are received

	// Check if correct number of arguments specified, exiting if not
	if (argc > 3) { Logger::Instance()->write(Logger::ERROR, "Too many arguments -  please specify verbosity flag. #tracks, offset1, offset2.  e.g. ./AlignTracker n 100000 0.0 0.0");}
	else if (argc < 3) {Logger::Instance()->write(Logger::ERROR, "Too few arguments - please specify verbosity flag, #tracks, offset1, offset2. e.g. ./AlignTracker n 100000 0.0 0.0");}
	else { // Set filenames to read random numbers from, using arguments. Catch exception if these files do not exist.
		try {
			compareStr = argv[1];
			tracksInput = stoi(argv[2]);
		}
		catch (ios_base::failure& e) {
			Logger::Instance()->write(Logger::ERROR, "Exception caught: " + string(e.what()) + "\nPlease ensure valid verbosity level specified!");
		}
	} // end of 2nd else [correct # arguments]
	Tracker::instance()->setTrackNumber(tracksInput*10); // request more tracks in case of rejection 
	//this is also passed to Tracker functions, with debug file names
	if (compareStr == "d") {
		debugBool = true; // print out to debug files [and verbose cout output]
		Logger::Instance()->write(Logger::WARNING,  "******DEBUG MODE*****");
	}
	else if (compareStr == "p") {
		plotBool = true;
		debugBool = true; // print out to debug files and plotting files - use with low track #
		Logger::Instance()->write(Logger::WARNING,  "******PLOTTING MODE*****");
	}
	else if (compareStr == "n" || compareStr == "a") {
		debugBool = false; // print out to debug files
		plotBool = false;  // print out to plotting files
	}
	else {
		Logger::Instance()->write(Logger::ERROR, "Please specify verbosity flag. (e.g. debug [d], plot[p] or align/normal [a/n])");
	}

	Logger::Instance()->setUseColor(false); // will be re-enabled below [to use custom colour output to terminal]
	stringstream msg0, msg01, msg02, msg1;
	Logger::Instance()->write(Logger::NOTE, "");
	msg0 << Logger::blue() <<  "*********************************************************************" << Logger::def();
	Logger::Instance()->write(Logger::NOTE, msg0.str());
	msg01 << Logger::yellow() << "  g-2 Tracker Alignment (v1.0) - Gleb Lukicov (UCL) - January 2018         " << Logger::def();
	Logger::Instance()->write(Logger::NOTE, msg01.str());
	msg1 << Logger::blue() <<  "*********************************************************************" << Logger::def();
	Logger::Instance()->write(Logger::NOTE, msg1.str());
	Logger::Instance()->setUseColor(true); // back to default colours
	helper << endl;
	helper << "Alignment Model with " << Tracker::instance()->getModuleN() << " tracker modules, having " << Tracker::instance()->getStrawN()
	       << " straws per layer." << endl;
	helper << "[" << Tracker::instance()->getLayerN() << " layers per module; " << Tracker::instance()->getViewN() << " views per module]." << endl;
	helper << "Total of " << Tracker::instance()->getLayerTotalN() << " measurement layers." << endl;
	helper << "No B-field, Straight Tracks (general lines), 100% efficiency." << endl;
	helper << "Hit rejection for (DCA > StrawRadius) is used: " << boolYN[Tracker::instance()->getHitCutStatus()] << endl;
	helper << "Tracks are rejected if one of hits have a (DCA < " << Tracker::instance()->getTrackCut() << ")" << " : " << boolYN[Tracker::instance()->getTrackCutBool()] << endl;
	helper << "Truth LR information is used in the reconstruction: " << boolYN[Tracker::instance()->getLRStatus()] << endl;
	helper << "p-value cut (<) for reconstructed tracks: " << Tracker::instance()->getPValCut() << endl;
	helper << "Straight Tracks with Circle Fit: single hit per layer allowed" << endl;
	helper << "Resolution is " << Tracker::instance()->getResolution() << " cm  [hit smearing]." << endl;
	helper << endl;

//------------------------------------------Setting Output Files---------------------------------------------------------//

	//arguments for Mille constructor:
	const char* outFileName = "Tracker_data.bin";
	bool asBinary = true; // set false for debugging
	bool writeZero = false; // to write 0 LC/DLC labels - not accepted by pede

	//Constraints, Steering and Parameter files are the inputs to Pede (together with binary file).
	ofstream constraint_file("Tracker_con.txt");
	ofstream steering_file("Tracker_str.txt");
	ofstream presigma_file("Tracker_par.txt");

	//Debug files [only filled with "d" option]: file streams for debug files
	// Setting fixed precision for floating point values
	ofstream debug_calc("Tracker_d_calc.txt");  debug_calc << fixed << setprecision(setPrecision);     // intermediate calculation
	ofstream debug_mis("Tracker_d_mis.txt");    debug_mis << fixed << setprecision(setPrecision);    // Misalignment
	ofstream debug_geom("Tracker_d_geom.txt");  debug_geom << fixed << setprecision(setPrecision);    // Geometry
	ofstream debug_off("Tracker_d_off.txt");    debug_off << fixed << setprecision(setPrecision); // Offsets/Missed hits
	ofstream debug_mc("Tracker_d_MC.txt");      debug_mc << fixed << setprecision(setPrecision);  // Final results from MC
	ofstream debug_con("Tracker_d_con.txt");    debug_con << fixed << setprecision(setPrecision);   //Constraints
	ofstream plot_gen("Tracker_p_gen.txt");     plot_gen << fixed << setprecision(setPrecision);  //Truth Track points
	ofstream plot_fit("Tracker_p_fit.txt");     plot_fit << fixed << setprecision(setPrecision);   //Reconstructed Track points
	ofstream contsants_plot("Tracker_p_constants.txt");  contsants_plot << fixed << setprecision(setPrecision);   // passing constants (e.g. strawN to python script)
	ofstream plot_centres("Tracker_p_centre.txt");  plot_centres << fixed << setprecision(setPrecision);   // passing constants (e.g. strawN to python script)
	ofstream plot_hits_gen("Tracker_p_hits_gen.txt");  plot_hits_gen << fixed << setprecision(setPrecision); // Truth Hits points
	ofstream plot_hits_fit("Tracker_p_hits_fit.txt");   plot_hits_fit << fixed << setprecision(setPrecision); // Recon Hits points
	ofstream pede_mis("Tracker_pede_mis.txt");  pede_mis << fixed << setprecision(setPrecision);  // Misalignments
	ofstream metric("Tracker_metric.txt"); stringstream metricStr; metricStr.str(""); //metric << fixed << setprecision(setPrecision);  // Misalignments

	//------------------------------------------ROOT: Booking etc.---------------------------------------------------------//

	//output ROOT file
	TFile* file = new TFile("Tracker.root", "recreate");
	//create a subdirectories
	TDirectory* cd_All_Hits = file->mkdir("All_Hits");
	TDirectory* cd_UV = file->mkdir("UV");
	TDirectory* cd_Modules = file->mkdir("Modules");
	TDirectory* cd_Tracks = file->mkdir("Tracks");
	TDirectory* cd_Straws = file->mkdir("Straws");
	TDirectory* cd_PEDE = file->mkdir("PEDE");

	// Book histograms [once only] - Key quantities
	TH1F* h_sigma_MP2 = new TH1F("h_sigma_MP2", "Resolution (#sigma)",  49,  Tracker::instance()->getResolution() - 0.001,
	                             Tracker::instance()->getResolution() + 0.001); // F=float bins, name, title, nBins, Min, Max
	TH1F* h_res_MP2 = new TH1F("h_res_MP2", "Residuals: Recon",  199, -0.08, 0.08);
	TH1F* h_dca = new TH1F("h_dca", "DCA",  149,  -0.1, Tracker::instance()->getStrawRadius() + 0.25);
	TH1F* h_driftRad = new TH1F("h_driftRad", "Drift Rad: circle fit",  149,  -0.1, Tracker::instance()->getStrawRadius() + 0.25);
	TH1F* h_dca_unsmeared = new TH1F("h_dca_unsmeared", "Unsmeared DCA",  98,  -0.1, Tracker::instance()->getStrawRadius() + 0.25);
	TH1I* h_id = new TH1I("h_id", "Straw IDs", Tracker::instance()->getStrawN(), 0, Tracker::instance()->getStrawN());
	// Track-generation-based
	TH1F* h_truth_slope = new TH1F("h_truth_slope", "Slope: Truth",  170,  -0.017, 0.017);

	TH1F* h_recon_slope = new TH1F("h_recon_slope", "Slope: Recon", 170,  -0.017, 0.017);
	
	TH1F* h_recon_intercept = new TH1F("h_recon_intercept", "Decay Vertex: Recon; [mm]",  104,  -13, 13);
	TH1F* h_truth_intercept = new TH1F("h_truth_intercept", "Decay Vertex: Truth; [mm]",  104,  -13, 13);
	TH1F* h_reconMinusTrue_track_intercept = new TH1F("h_reconMinusTrue_track_intercept", " #Delta (Recon - True) Vertex; [um]",  519,  -2000.0, 2000.0);
	
	TH1F* h_x0 = new TH1F("h_x0", "Truth #x_{0}",  99,  -2, 2);
	TH1F* h_x1 = new TH1F("h_x1", "Truth #x_{1}",  99,  -3, 3);
	//Track/Hits-based
	TH1I* h_labels = new TH1I("h_labels", "Labels in PEDE",  72, 10, 46);
	TH1F* h_residual_true = new TH1F("h_residual_true", "Residuals: Truth", 500, -0.06, 0.06);
	TH1F* h_chi2_true = new TH1F("h_chi2_true", "#Chi^{2}: Truth", 40, -1, 50);
	TH1F* h_residual_recon = new TH1F("h_residual_recon", "Residuals: Recon", 199, -0.2, 0.2);
	TH1F* h_chi2_recon = new TH1F("h_chi2_recon", "#Chi^{2}: Recon", 250, 0, 250);
	TH1I* h_hitCount = new TH1I("h_hitCount", "Hit count", 33 , 0, 33);
	TH1F* h_reconMinusTrue_track_slope = new TH1F("h_reconMinusTrue_track_slope", "#Delta (Recon - True) Slope",  119,  -0.002, 0.002);

	TH1F* h_pval = new TH1F("p_value", "p-value", 49, -0.01, 1.1);
	TH1F* h_chi2_circle = new TH1F("h_chi2_circle", "#Chi^{2}: circle-fit", 149, -0.1, 90);
	TH1F* h_chi2_circle_ndf = new TH1F("h_chi2_circle_ndf", "#Chi^{2}/ndf: circle-fit", 149, -0.1, 10.0);
	//PEDE: Derivatives
	TH1F* h_DLC1 = new TH1F("h_DLC1", "DLC1: All Modules",  149,  -1.1, 1.1); h_DLC1->SetDirectory(cd_PEDE);
	TH1F* h_DLC2 = new TH1F("h_DLC2", "DLC2: All Modules",  879,  -65.0, 65.0); h_DLC2->SetDirectory(cd_PEDE);
	TH1F* h_DGL1 = new TH1F("h_DGL1", "DGL1: All Modules",  349,  -1.1, 1.1); h_DGL1->SetDirectory(cd_PEDE);
	TH1F* h_DGL2 = new TH1F("h_DGL2", "DGL2: All Modules",  149,  -0.017, 0.017); h_DGL2->SetDirectory(cd_PEDE);
	TH1F* h_DGL3 = new TH1F("h_DGL3", "DGL3: All Modules",  149,  -1.1, 1.1); h_DGL3->SetDirectory(cd_PEDE);

	//Use array of pointer of type TH1x to set axis titles and directories
	TH1F* cmTitle[] = {h_sigma_MP2, h_res_MP2, h_dca,
	                   h_truth_intercept, h_x0, h_x1, h_residual_true, h_residual_recon,
	                   h_driftRad, h_dca_unsmeared
	                  };
	for (int i = 0; i < (int) sizeof( cmTitle ) / sizeof( cmTitle[0] ); i++) {
		TH1F* temp = cmTitle[i];
		cmTitle[i]->SetXTitle("[cm]");
	}
	TH1F* cdAllHits_F[] = {h_sigma_MP2, h_res_MP2, h_dca, h_chi2_true, h_residual_recon, h_residual_true,
	                       h_chi2_recon, h_driftRad, h_dca_unsmeared };
	TH1F* cdTracks_F[] = {h_truth_slope, h_x0, h_x1, h_reconMinusTrue_track_slope,
	                      h_recon_slope, h_pval, h_chi2_circle, h_chi2_circle_ndf };
	TH1I* cdAllHits_I[] = {h_labels, h_hitCount, h_id};
	for (int i = 0; i < (int) sizeof( cdAllHits_F ) / sizeof( cdAllHits_F[0] ); i++) {
		cdAllHits_F[i]->SetDirectory(cd_All_Hits);
	}
	for (int i = 0; i < (int) sizeof( cdTracks_F ) / sizeof( cdTracks_F[0] ); i++) {
		cdTracks_F[i]->SetDirectory(cd_Tracks);
	}
	for (int i = 0; i < (int) sizeof( cdAllHits_I ) / sizeof( cdAllHits_I[0] ); i++) {
		cdAllHits_I[i]->SetDirectory(cd_All_Hits);
	}

	stringstream h_name; //to book hist. in a loop
	stringstream h_title;
	// Modules, views, layers
	for (int i_module = 0; i_module < Tracker::instance()->getModuleN(); i_module++) {
		for (int i_view = 0; i_view < Tracker::instance()->getViewN(); i_view++) {
			for (int i_layer = 0; i_layer < Tracker::instance()->getLayerN(); i_layer++) {
				string UV = Tracker::instance()->getUVmapping(i_view, i_layer); // converting view/layer ID into conventional labels
				h_name.str(""); h_name << "h_residual_recon_M_" << i_module << "_" << UV;
				h_title.str(""); h_title << "Residuals Recon M" << i_module << " " << UV ;
				auto hl4 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  149, -0.2, 0.2);
				hl4->GetXaxis()->SetTitle("[cm]"); hl4->SetDirectory(cd_UV);
			} // layers
		} // views
	} // modules

	//Modules
	for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++) {
		h_name.str(""); h_name << "h_Residuals_Module_" << i_module;
		h_title.str(""); h_title << "Residuals Recon M" << i_module;
		auto hm4 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  199, -0.2, 0.2);
		hm4->GetXaxis()->SetTitle("[cm]"); hm4->SetDirectory(cd_Modules);
		// LR splitting
		for (int i_LR = 0; i_LR < 2; i_LR++) {
			h_name.str(""); h_name << "h_Residuals_Module_" << i_module << "_" <<  valueLR[i_LR];
			h_title.str(""); h_title << "Residuals Recon M" << i_module  << " " << nameLR[i_LR];
			auto hm6 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  199, -0.2, 0.2);
			hm6->GetXaxis()->SetTitle("[cm]"); hm6->SetDirectory(cd_Modules);
		}
	} // end of modules loop

//------------------------------------------Mille Routines---------------------------------------------------------//

	// Creating .bin, steering, presigma and constrain files
	Mille M(outFileName, asBinary, writeZero);  // call to Mille.cc to create a .bin file
	
	helper << "Generating test data for g-2 Tracker Alignment in PEDE:" << endl;

	helper << fixed << setprecision(4);
	// XXX: definition of broken lines here in the future

	// SETTING GEOMETRY and MISALIGNMENT
	Tracker::instance()->setGeometry(debug_geom, debug_mis, pede_mis, plot_centres, debugBool, metric);
	helper << "Misalignment is complete!" << endl;
	helper << fixed << setprecision(setPrecision);

	// Write a constraint file, for use with pede
	Tracker::instance()->writeConstraintFile(constraint_file, debug_con, debugBool, metric);
	helper << "Constraints are written! [see Tracker_con.txt]" << endl;

	//Now writing the steering file
	Tracker::instance()->writeSteeringFile(steering_file, metric);
	helper << "Steering file was generated! [see Tracker_con.txt]" << endl;

	// Now set pre-sigma for know global parameters
	Tracker::instance()->writePresigmaFile(presigma_file, metric);
	helper << "Presigma Parameter file was generated! [see Tracker_par.txt]" << endl;

	helper << "Calculating residuals..." << endl;

//------------------------------------------Main Mille Track Loop---------------------------------------------------------//

	bool StrongDebugBool = true;
	//Generating tracks
	for (int trackCount = 0; trackCount < Tracker::instance()->getTrackNumber(); trackCount++) {
		//float p=pow(10.0, 1+Tracker::instance()->generate_uniform()); // XXX
		//scatterError=sqrt(Tracker::instance()->getWidth())*0.014/p; // XXX
		scatterError = 0; // set no scatterError for now
		//if (debugBool && StrongDebugBool) { helper << "Track: " << trackCount << endl; }
		
		//Generating tracks
		MCData generatedMC = Tracker::instance()->MCLaunch(scatterError, debug_calc, debug_off, debug_mc, plot_fit, plot_gen, plot_hits_gen,
		                     plot_hits_fit, debugBool);

		int totalHitLayers = generatedMC.totalLayerHits;
		//std::cout << "totalHitLayers= " <<totalHitLayers << "\n"; 

		//First of all, check if track has not failed the cut
		if (generatedMC.cut == false) {

			for (int hitCount = 0; hitCount < totalHitLayers; hitCount++) { //counting only hits going though detector

				//******************PEDE INPUTS******************************
				float resiudalRecon =  generatedMC.residualsRecon[hitCount]; //Reconstructed Residual
				float resolution =  Tracker::instance()->getResolution(); //Resolution
				//Number of local and global parameters
				const int nalc = Tracker::instance()->getNLC();
				const int nagl = Tracker::instance()->getNGL();

				//Variables for derivative calculations:
				float z = generatedMC.zStraw[hitCount];
				float x = generatedMC.xStraw[hitCount];
				float m = generatedMC.slopeRecon;
				float c = generatedMC.interceptRecon;
				float zc = generatedMC.zCentreModule[hitCount]; // rotational centre
				float xc = generatedMC.xCentreModule[hitCount]; // rotational centre

				//Local derivatives: see alignment.tex for derivations
				float dlc1 = ( c + m * z - x ) / ( sqrt(m * m + 1) * abs(c + m * z - x) ) ; // "DCA magnitude" dR/dc
				float dlc2 = ( (m * m + 1) * z * (c + m * z - x) - m * pow(abs(c + m * z - x), 2) ) / ( pow(m * m + 1, 1.5) * abs(c + m * z - x)  ) ; //dR/dm
				float derlc[] = {dlc1, dlc2};
				//Global derivatives
				float dgl1 = ( c + m * z - x ) / ( sqrt(m * m + 1) * abs(c + m * z - x) );  //dR/dx
				float dgl2 = ( m * ( c + m * z - x ) ) / ( sqrt(m * m + 1) * abs(c + m * z - x) ); //dR/dz
				float dgl3 = ( ( m * ( c + m * z - x ) ) / ( sqrt(m * m + 1) * abs(c + m * z - x) ) * (-x + xc )  )  +  ( ( c + m * z - x ) / ( sqrt(m * m + 1) * abs(c + m * z - x) ) * (z - zc) ); //dR/d𝛉
				float dergl[] = {dgl1, dgl2, dgl3};
				//Labels
				int l1 = generatedMC.label1[hitCount]; //Mx
				int l2 = generatedMC.label2[hitCount]; //Mz
				int l3 = generatedMC.label3[hitCount]; //M𝛉
				int label[] = {l1, l2, l3};

				//TODO multiple scattering errors (no correlations) (for imodel == 1)
				//add break points multiple scattering later XXX (for imodel == 2)
				//! add 'broken lines' offsets for multiple scattering XXX (for imodel == 3)

				M.mille(nalc, derlc, nagl, dergl, label, resiudalRecon, resolution);

				//***********************************Sanity Plots: Hits******************************************//
				//Getting the "detector coordinates"
				int moduleN = generatedMC.Module_i[hitCount];
				int viewN = generatedMC.View_i[hitCount];
				int layerN = generatedMC.Layer_i[hitCount];
				int strawID = generatedMC.Straw_i[hitCount];
				string UV = Tracker::instance()->getUVmapping(viewN, layerN); // converting view/layer ID into conventional labels

				//Fill for all hits
				h_res_MP2 -> Fill (resiudalRecon); // residuals
				h_sigma_MP2 -> Fill(resolution); // errors
				h_dca->Fill(generatedMC.dca[hitCount]); // truth DCA
				h_driftRad->Fill(generatedMC.driftRad[hitCount]); // recon DCA
				h_dca_unsmeared->Fill(generatedMC.dcaUnsmeared[hitCount]);
				h_labels->Fill(l1); h_labels->Fill(l2); h_labels->Fill(l3);
				h_id ->Fill(strawID);
				if (generatedMC.driftRad[hitCount] < 0) negDCA++;
				h_DLC1->Fill(dlc1); h_DLC2->Fill(dlc2); h_DGL1->Fill(dgl1); h_DGL2->Fill(dgl2); h_DGL3->Fill(dgl3);

				//Calculating Chi2:
				float residualTruth = generatedMC.residualsTruth[hitCount];
				h_residual_true->Fill(residualTruth);
				residualsTrueSum_2 += pow(residualTruth / resolution, 2);
				h_residual_recon->Fill(resiudalRecon); //already used as input to mille
				residualsReconSum_2 += pow(resiudalRecon / resolution, 2);

				//Fill for hits in modules/layers/straws
				h_name.str(""); h_name << "Modules/h_Residuals_Module_" << moduleN;
				TH1F* h8 = (TH1F*)file->Get( h_name.str().c_str() );
				h8->Fill(resiudalRecon);

				h_name.str(""); h_name << "UV/h_residual_recon_M_" << moduleN << "_" << UV;
				TH1F* h9 = (TH1F*)file->Get( h_name.str().c_str() );
				h9->Fill(resiudalRecon);

				h_name.str(""); h_name << "Modules/h_Residuals_Module_" << moduleN << "_" << generatedMC.LR[hitCount];
				TH1F* h13 = (TH1F*)file->Get( h_name.str().c_str() );
				h13->Fill(resiudalRecon);

				hitsN++; //count hits
				//std::cout << "hitsN=" << hitsN << "\n";
			} // end of hits loop

			//***********************************Sanity Plots: Tracks******************************************/
			//For true tracks
			float chi2True = residualsTrueSum_2;
			h_chi2_true->Fill(chi2True);
			//For recon tracks
			float chi2Recon = residualsReconSum_2;
			h_chi2_recon->Fill(chi2Recon);
			//Resetting counters for next track
			residualsTrueSum_2 = 0;
			residualsReconSum_2 = 0;
			h_hitCount->Fill(totalHitLayers);

			//Filling Track-based plots
			h_truth_slope->Fill(generatedMC.slopeTruth);
			h_truth_intercept->Fill(generatedMC.interceptTruth*10);
			h_x0->Fill(generatedMC.x0);
			h_x1->Fill(generatedMC.x1);
			h_reconMinusTrue_track_intercept->Fill((generatedMC.interceptTruth - generatedMC.interceptRecon)*10000.0); // um
			h_reconMinusTrue_track_slope->Fill(generatedMC.slopeTruth - generatedMC.slopeRecon);
			h_recon_slope->Fill(generatedMC.slopeRecon);
			h_recon_intercept->Fill(generatedMC.interceptRecon*10.0);
			h_pval->Fill(generatedMC.pValue);
			h_chi2_circle->Fill(generatedMC.chi2Circle);
			h_chi2_circle_ndf->Fill(generatedMC.chi2Circle / (totalHitLayers - 2));

			// XXX additional measurements from MS IF (imodel == 2) THEN
			//IF (imodel >= 3) THEN

			M.end(); // Write buffer (set of derivatives with same local parameters) to file.
			recordN++; // count records (i.e. written tracks);

            //Stop when requested number of tracks is reached 
            if (recordN >= Tracker::instance()->getTrackNumber()) goto stop;

		} // cut on DCA check
	} // end of track count // End of Mille // End of collecting residual records
	//Passing constants to plotting scripts
    stop:
    std::cout<<"stopped\n";
	contsants_plot << Tracker::instance()->getModuleN() << " " << Tracker::instance()->getViewN() << " "
	               << Tracker::instance()->getLayerN() << " " << Tracker::instance()->getStrawN() << " " << recordN << " "
	               << Tracker::instance()->getBeamOffset()   << " " << Tracker::instance()->getBeamStart() << " " <<  Tracker::instance()->getBeamPositionLength()
	               << "  " << Tracker::instance()->getBeamStop() << " " 

	               << endl;

	helper << "Mille residual-accumulation routine completed! [see Tracker_data.bin]" << endl;

	//------------------------------------------ROOT: Fitting Functions---------------------------------------------------------//

	helper << endl;
	helper << "-------------------------------------------------------------------------" << endl;
	helper << "ROOT fitting parameters and output:" << endl;

	Chi2ReconActual = h_chi2_recon->GetMean(); // Get reconstructed Chi2

	helper << "-------------------------------------------------------------------------" << endl;
	helper << " " << endl;
	helper << Tracker::instance()->getTrackNumber() << " tracks requested; " << recordN << " generated with " << hitsN << " hits." << endl;
	float rejectedTracks = (Tracker::instance()->getTrackNumber() - recordN) / float(Tracker::instance()->getTrackNumber()) * 100.0;
	helper << Tracker::instance()->getTrackNumber() - recordN << " records rejected (" << rejectedTracks << " %)." << endl;
	helper << "Number of DCA (==drift radii) smeared below 0 is " << negDCA << endl;
	stringstream out1, out2, out3, out4, out5;
	out2 << "Measured Mean Chi2 (circle fit) " << Chi2ReconActual;
	Logger::Instance()->write(Logger::WARNING, out2.str());
	LOG << out1.rdbuf() << endl; LOG << out2.rdbuf() << endl;
	float rejectsFrac = Tracker::instance()->getRejectedHitsDCA();
	rejectsFrac = rejectsFrac / (Tracker::instance()->getLayerTotalN() * recordN);
	helper << fixed << setprecision(1);
	out3 << "Hits that missed a straw (DCA rejection for all layers): " << Tracker::instance()->getRejectedHitsDCA() << " (" << rejectsFrac * 100 << "%).";
	out4 << "Multiple hits (for all layers): " << Tracker::instance()->getMultipleHitsLayer() << ".";
	out5 << "Ambiguity Hits Resolved (rand.): " << Tracker::instance()->getAmbiguityHit() << ".";
	LOG << out3.rdbuf() << endl; LOG << out4.rdbuf() << endl; LOG << out5.rdbuf() << endl;
	Logger::Instance()->write(Logger::WARNING, out3.str());
	Logger::Instance()->write(Logger::WARNING, out4.str());
	Logger::Instance()->write(Logger::WARNING, out5.str());
	helper << " " << endl;
	Logger::Instance()->setUseColor(false); // will be re-enabled below
	stringstream msg2, msg3, msg4, msgA, msgB;
	msgA <<  Logger::blue() << "Ready for PEDE algorithm: ./pede Tracker_str.txt" << Logger::def();
	msgB << Logger::blue() << "Sanity Plots: root Tracker.root" << Logger::def();
	Logger::Instance()->write(Logger::NOTE, msgA.str()); Logger::Instance()->write(Logger::NOTE, msgB.str());
	// Millepede courtesy of John
	msg2 << Logger::green() << "    _____________________________  \\  /" << Logger::def();
	Logger::Instance()->write(Logger::NOTE, msg2.str());
	msg3 << Logger::yellow() << "   {_|_|_|_|_|_|_|_|_|_|_|_|_|_|_( ͡° ͜ʖ ͡°) " << Logger::def();
	Logger::Instance()->write(Logger::NOTE, msg3.str());
	msg4 << Logger::red() << "    /\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/" << Logger::def();
	Logger::Instance()->write(Logger::NOTE, msg4.str());
	Logger::Instance()->setUseColor(true); // back to default colours
	helper << "<d(vertex)>= " << h_reconMinusTrue_track_intercept->GetMean() << " +- " << h_reconMinusTrue_track_intercept->GetMeanError() << "\n";
	helper << "d(vertex) SD = " << h_reconMinusTrue_track_intercept->GetRMS() << " +- " << h_reconMinusTrue_track_intercept->GetRMSError() << "\n";
	helper << fixed << setprecision(4);
	helper << "p-value= " << h_pval->GetMean() << " +- " << h_pval->GetMeanError() << "\n";
	if (debugBool) {
		Logger::Instance()->write(Logger::WARNING, "Text debug files were produced: ls Tracker_d_*.txt");
	}
	if (plotBool) {
		Logger::Instance()->write(Logger::WARNING, "./PlotMP2.py to draw geometry with hits");
	}

	// Close text files
	constraint_file.close();
	steering_file.close();
	debug_mis.close();
	debug_geom.close();
	debug_calc.close();
	debug_off.close();
	debug_con.close();
	// Close ROOT file
	file->Write();
	file->Close();

	helper << endl;
	helper << "Programme log written to: Tracker_log.txt" << endl;
	helper << fixed << setprecision(4);
	t_cpu = clock() - t_cpu;
	auto t_end = chrono::high_resolution_clock::now();
	helper << "Programme execution took " <<  t_cpu << " CPU clicks (" << ((float)t_cpu) / CLOCKS_PER_SEC << " s)." << " Wall clock time passed: "
	       << chrono::duration<double>(t_end - t_start).count() << " s." << endl;
	time_t now = time(0);
	char* dt = ctime(&now);
	helper << "The C++ compiler used: " << true_cxx << " " << true_cxx_ver
	       << " Job finished on: " << dt << endl;
	return 0;
} //end of main