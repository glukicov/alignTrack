/*
*   Gleb Lukicov (g.lukicov@ucl.ac.uk) @ Fermilab
*   Created: 17 April 2017
*   Modified: 9 October 2017
----------------------------------------------------------------
This programme uses MC methods to produce a .bin data file for the
PEDE routine, to align the tracking detector for the g-2
experiment.
Methods and functions are contained in AlignTracker_methods.cpp (Tracker class)

============================= Test Model v0.4 =======================================
*Simple 2D case (1D alignment along x).
* No B-field, straight tracks, no MS, 100% efficiency.
*(!) NATURAL UNITS (same as PEDE): cm, rad, GeV [natural units will be omitted]

(!) Module, Straw, View, Layer labelling starts from 0 [+1 is put by hand as MP2 doesn't accept 0 as a label]

Terminology: detector=module (==the alignment object)

Gaussian (gaus) Smearing: Mean=0, RMS=1  (on hits, resolution, etc.)
Beam originates (in z direction) from z=0, straight tracks, no scattering, only resolution (gaus) smearing,
    tracks are lines between z=beamStart and z=beamStop, with uniformly generated x1 and slope.
Resolution (i.e. tracker systematic) is implemented.
Misalignment is implemented "manually" for each module in +/- x-direction [all layers in module are dependent (i.e. move together)]

Hit rejection: (>0.5 straw radius)  -not used yet

Constrains: TODO
Steering options: method inversion 1 0.01 [1 iteration, 0.01 dF convergence; inversion method provides uncertainties on all alignment parameters]

                       			--Tracker Geometry [cm]--:
Spacing (in x) between straws in layer is 0.606 [strawSpacing]
Spacing between layers in same view is 0.515 [layerSpacing]
Spacing between U and V views (last and first layers in views) is 2.020 [viewSpaing]
Spacing between modules (last and first layers in modules) 13.735 [moduleSpacing]
Layers in a view are (e.g. U0 vs. U1) have an extra relative x displacement of 0.303 [layerDisplacement]
The staircaseXDisplacment of modules is not implemented yet

For more information see:
http://gm2-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=4375&filename=Orientation.pdf&version=1

=======================================================================================

#TODO Model type [physics] will be implemented via imodel = 0, 1 .....

#Verbosity level implemented via command line command: n = normal, d = debug, p=plotting (./PlotGen.py) [e.g. './AlignTracker d']

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

//------------------------------------------Variable Initialisation---------------------------------------------------------//

	int imodel = 0;  //Model type (see above) TODO implement this as an argument to main [for broken lines, MF, etc.]
	int setPrecision = 7; // precision (# decimal points) of printout for debug text files and cout
	string compareStr; //for debug vs. normal output as specified by the user
	int tracksInput; // number of tracks to generate as specified by the user
	float offset1;
	float offset2;
	bool debugBool = false; // './AlignTracker n' - for normal, of ./AlignTracker d' - for verbose debug output
	bool plotBool = false; // './AlignTracker p' - for plotting with PlotGen.py
	//Set up counters for hits and records (tracks)
	int hitsN = 0; // actually recorded (i.e. non-rejected hits)
	int recordN = 0; //records = tracks
	float scatterError; // multiple scattering error [calculated here and passed back to the Tracker class]
	// Those variable are set inside the hits loop, and used outside - hence global scope
	float residuals_true_sum_2 = 0.0; // squared sum of the true residuals (from measurements)
	float residuals_recon_sum_2 = 0.0; // squared sum of the recon. residuals (from measurements)
	float pivotPoint_actual = 0.0; // central z point from measurements
	float Chi2_recon_actual = 0.0; // Reconstructed Chi2 (from measurements)
	int negDCA = 0; // counting negatively smeared DCAs
	const Color_t colourVector[] = {kMagenta, kOrange, kBlue, kGreen, kYellow, kRed, kGray, kBlack}; //8 colours for up to 8 modules
	gErrorIgnoreLevel = kWarning; // Display ROOT Warning and above messages [i.e. suppress info]
	// Simple LR mapping for ROOT plots
	char nameLR[] = {'L', 'R'};
	char nameResSign[] = {'P', 'N'};
	const char* nameLC[] = {"DLC1", "DLC2"};
	int valueLR[] = { -1, 1};
	// True/False -> Yes/No mapping
	const char* boolYN[2] = {"No", "Yes"};

	// use overloading to duplicate cout into the log file
	ofstream LOG; LOG.open("Tracker_log.txt"); MyStreamingHelper helper(LOG, cout);
	cout << fixed << setprecision(setPrecision); // set precision of standard screen output
	helper << fixed << setprecision(setPrecision); // set precision of log file output

	//Tell the logger to only show message at INFO level or above
	Logger::Instance()->setLogLevel(Logger::NOTE);
	//Tell the logger to throw exceptions when ERROR messages are received
	Logger::Instance()->enableCriticalErrorThrow();

	// Check if correct number of arguments specified, exiting if not
	if (argc > 5) { Logger::Instance()->write(Logger::ERROR, "Too many arguments -  please specify verbosity flag. (e.g. debug [d], plot[p] or align/normal [n])");}
	else if (argc < 5) {Logger::Instance()->write(Logger::ERROR, "Too few arguments - please specify verbosity flag. (e.g. e.g. debug [d], plot[p] or align/normal [n])");}
	else { // Set filenames to read random numbers from, using arguments. Catch exception if these files do not exist.
		try {
			compareStr = argv[1];
			tracksInput = stoi(argv[2]);
			offset1 = stof(argv[3]);
			offset2 = stof(argv[4]);
		}
		catch (ios_base::failure& e) {
			Logger::Instance()->write(Logger::ERROR, "Exception caught: " + string(e.what()) + "\nPlease ensure valid verbosity level specified!");
		}
	} // end of 2nd else [correct # arguments]

	//this is also passed to Tracker functions, with debug file names
	if (compareStr == "d") {
		debugBool = true; // print out to debug files [and verbose cout output]
		Tracker::instance()->setTrackNumber(tracksInput);
		Logger::Instance()->write(Logger::WARNING,  "******DEBUG MODE*****");
		Tracker::instance()->setOffset1(offset1);
		Tracker::instance()->setOffset2(offset2);
	}
	else if (compareStr == "p") {
		plotBool = true;
		debugBool = true; // print out to debug files and plotting files - use with low track #
		Tracker::instance()->setTrackNumber(tracksInput);
		Logger::Instance()->write(Logger::WARNING,  "******PLOTTING MODE*****");
		Tracker::instance()->setOffset1(offset1);
		Tracker::instance()->setOffset2(offset2);
	}
	else if (compareStr == "n" || compareStr == "a") {
		debugBool = false; // print out to debug files
		plotBool = false;  // print out to plotting files
		Tracker::instance()->setTrackNumber(tracksInput);
		Tracker::instance()->setOffset1(offset1);
		Tracker::instance()->setOffset2(offset2);
	}
	else {
		Logger::Instance()->write(Logger::ERROR, "Please specify verbosity flag. (e.g. debug [d], plot[p] or align/normal [a/n])");
	}

	Logger::Instance()->setUseColor(false); // will be re-enabled below [to use custom colour output to terminal]
	stringstream msg0, msg01, msg02, msg1;
	Logger::Instance()->write(Logger::NOTE, "");
	msg0 << Logger::blue() <<  "*********************************************************************" << Logger::def();
	Logger::Instance()->write(Logger::NOTE, msg0.str());
	msg01 << Logger::yellow() << "  g-2 Tracker Alignment (v0.5) - Gleb Lukicov (UCL) - November 2017         " << Logger::def();
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

	// See https://github.com/glukicov/alignTrack for instructions to generate random numbers
	try {
		Tracker::instance()->set_uniform_file("uniform_ran.txt");
		Tracker::instance()->set_gaussian_file("gaussian_ran.txt");
	}

	catch (ios_base::failure& e) {
		cerr << "Filestream exception caught: " << e.what() << endl;
		cerr << "Please ensure valid filenames are specified!" << endl;
		cerr << "Have random numbers been generated with randomIntGenerator.py?" << endl;
		return 1;
	}
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
	ofstream debug_mp2("Tracker_d_mille.txt");  debug_mp2 << fixed << setprecision(setPrecision);  //Inputs into binary files
	ofstream debug_calc("Tracker_d_calc.txt");  debug_calc << fixed << setprecision(setPrecision);     // intermediate calculation
	ofstream debug_mis("Tracker_d_mis.txt");    debug_mis << fixed << setprecision(setPrecision);    // Misalignment
	ofstream debug_geom("Tracker_d_geom.txt");  debug_geom << fixed << setprecision(setPrecision);    // Geometry
	ofstream debug_off("Tracker_d_off.txt");    debug_off << fixed << setprecision(setPrecision); // Offsets/Missed hits
	ofstream debug_mc("Tracker_d_MC.txt");      debug_mc << fixed << setprecision(setPrecision);  // Final results from MC
	ofstream debug_con("Tracker_d_con.txt");    debug_con << fixed << setprecision(setPrecision);   //Constraints
	ofstream plot_gen("Tracker_p_gen.txt");     plot_gen << fixed << setprecision(setPrecision);  //Truth Track points
	ofstream plot_fit("Tracker_p_fit.txt");     plot_fit << fixed << setprecision(setPrecision);   //Reconstructed Track points
	ofstream contsants_plot("Tracker_p_constants.txt");  contsants_plot << fixed << setprecision(setPrecision);   // passing constants (e.g. strawN to python script)
	ofstream plot_hits_gen("Tracker_p_hits_gen.txt");  plot_hits_gen << fixed << setprecision(setPrecision); // Truth Hits points
	ofstream plot_hits_fit("Tracker_p_hits_fit.txt");   plot_hits_fit << fixed << setprecision(setPrecision); // Recon Hits points
	ofstream pede_mis("Tracker_pede_mis.txt");  pede_mis << fixed << setprecision(setPrecision);  // Misalignments
	ofstream timeFile("Tracker_time.txt", ios_base::app);  timeFile << fixed << setprecision(setPrecision);  // Misalignments
	ofstream metric("Tracker_metric.txt"); stringstream metricStr; metricStr.str(""); //metric << fixed << setprecision(setPrecision);  // Misalignments

	//ofstream debug_append("Tracker_d_append.txt", std::ios_base::app);  debug_append << fixed << setprecision(setPrecision);  // Misalignments

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
	TH1F* h_sigma = new TH1F("h_sigma", "Resolution (#sigma)",  49,  Tracker::instance()->getResolution() - 0.001,
	                         Tracker::instance()->getResolution() + 0.001); // F=float bins, name, title, nBins, Min, Max
	TH1F* h_res_MP2 = new TH1F("h_res_MP2", "Residuals: Recon",  199, -0.06, 0.06);
	TH1F* h_dca = new TH1F("h_dca", "DCA",  149,  -0.1, Tracker::instance()->getStrawRadius() + 0.25);
	TH1F* h_dca_unsmeared = new TH1F("h_dca_unsmeared", "Unsmeared DCA",  98,  -0.1, Tracker::instance()->getStrawRadius() + 0.25);
	TH1I* h_id_dca = new TH1I("h_id_dca", "Straw IDs", Tracker::instance()->getStrawN(), 0, Tracker::instance()->getStrawN());
	// Track-generation-based
	TH1F* h_slope = new TH1F("h_slope", "Slope: Truth",  80,  -0.02, 0.02);
	TH1F* h_intercept = new TH1F("h_intercept", "Intercept: Truth ",  88,  -1.3, 1.3);
	TH1F* h_recon_slope = new TH1F("h_recon_slope", "Slope: Recon", 80,  -0.02, 0.02);
	TH1F* h_recon_intercept = new TH1F("h_recon_intercept", "Intercept: Recon",  99,  -1.3, 1.3);
	TH1F* h_x0 = new TH1F("h_x0", "Truth #x_{0}",  99,  -2, 2);
	TH1F* h_x1 = new TH1F("h_x1", "Truth #x_{1}",  99,  -3, 3);
	//Track/Hits-based
	TH1F* h_track_true = new TH1F("h_track_true", "All track points: Truth",  49,  -3, 3);
	TH1F* h_track_recon = new TH1F("h_track_recon", "All track points: Recon",  49,  -3, 3);
	TH1F* h_track_TR_diff = new TH1F("h_track_TR_diff", "#Delta (Recon-True) track points",  149, -0.03, 0.03);
	TH1I* h_labels = new TH1I("h_labels", "Labels in PEDE", 8 , 0, 8);
	TH1F* h_residual_true = new TH1F("h_residual_true", "Residuals: Truth", 500, -0.06, 0.06);
	TH1F* h_chi2_true = new TH1F("h_chi2_true", "#Chi^{2}: Truth", 40, -1, 50);
	TH1F* h_residual_recon = new TH1F("h_residual_recon", "Residuals: Recon", 199, -0.2, 0.2);
	TH1F* h_chi2_recon = new TH1F("h_chi2_recon", "#Chi^{2}: Recon", 250, 0, 120);
	TH1I* h_hitCount = new TH1I("h_hitCount", "Hit count", 32 , 0, 32);
	TH1F* h_reconMinusTrue_track_slope = new TH1F("h_reconMinusTrue_track_slope", "#Delta (Recon - True) Slope",  119,  -0.002, 0.002);
	TH1F* h_reconMinusTrue_track_intercept = new TH1F("h_reconMinusTrue_track_intercept", " #Delta (Recon - True) Intercept",  119,  -0.06, 0.06);
	TH1F* h_pval = new TH1F("p_value", "p-value", 48, -0.1, 1.1);
	TH1F* h_chi2_circle = new TH1F("h_chi2_circle", "#Chi^{2}: circle-fit", 89, -0.1, 90);
	TH1F* h_chi2_circle_ndf = new TH1F("h_chi2_circle_ndf", "#Chi^{2}/ndf: circle-fit", 89, -0.1, 4);
	TH1F* h_driftRad = new TH1F("h_driftRad", "Drift Rad: circle fit",  149,  -0.1, Tracker::instance()->getStrawRadius() + 0.25);
	TH1F* h_DLC1 = new TH1F("h_DLC1", "DLC1: All Modules",  149,  -1.1, 1.1); h_DLC1->SetDirectory(cd_PEDE);
	TH1F* h_DLC2 = new TH1F("h_DLC2", "DLC2: All Modules",  879,  -65.0, 65.0); h_DLC2->SetDirectory(cd_PEDE);
	TH1F* h_DGL1 = new TH1F("h_DGL1", "DGL1: All Modules",  149,  -1.1, 1.1); h_DGL1->SetDirectory(cd_PEDE);

	// "special" histos
	THStack* hs_hits_recon = new THStack("hs_hits_recon", "");
	TH2F* h_res_x_z = new TH2F("h_res_x_z", "Residuals vs z", 600, 0, 18 * Tracker::instance()->getModuleN(), 79, -0.1, 0.1);
	h_res_x_z->SetDirectory(cd_All_Hits); h_res_x_z->GetXaxis()->SetTitle("cm");  h_res_x_z->GetYaxis()->SetTitle("cm");
	TH2F* h_SD_z_res_Recon = new TH2F("h_SD_z_res_Recon", "Residuals SD: Recon", 600, 0, 18 * Tracker::instance()->getModuleN(), 59, 80, 400);
	h_SD_z_res_Recon->SetDirectory(cd_All_Hits); h_SD_z_res_Recon->GetXaxis()->SetTitle("Module/Layer separation [cm]");  h_SD_z_res_Recon->GetYaxis()->SetTitle("Residual SD [um]");
	TH2F* h_SD_z_res_Est = new TH2F("h_SD_z_res_Est", "Residuals SD: Expect", 600, 0, 18 * Tracker::instance()->getModuleN(), 59, 120, 150);
	h_SD_z_res_Est->SetDirectory(cd_All_Hits); h_SD_z_res_Est->GetXaxis()->SetTitle("Module/Layer separation [cm]");  h_SD_z_res_Est->GetYaxis()->SetTitle("Residual SD [um]");
	TH2F* h_Pulls_z = new TH2F("h_Pulls_z", "Pulls", 600, 0, 18 * Tracker::instance()->getModuleN(), 59, -1, 10);
	h_Pulls_z->SetDirectory(cd_All_Hits); h_Pulls_z->GetXaxis()->SetTitle("Module/Layer separation [cm]");  h_Pulls_z->GetYaxis()->SetTitle("Measurement Pulls [cm]");

	//Use array of pointer of type TH1x to set axis titles and directories
	TH1F* cmTitle[] = {h_reconMinusTrue_track_intercept, h_sigma, h_res_MP2, h_dca, h_track_true, h_track_recon,
	                   h_intercept, h_x0, h_x1, h_recon_intercept, h_residual_true, h_residual_recon,
	                   h_driftRad, h_track_TR_diff, h_dca_unsmeared
	                  };
	for (int i = 0; i < (int) sizeof( cmTitle ) / sizeof( cmTitle[0] ); i++) {
		TH1F* temp = cmTitle[i];
		cmTitle[i]->SetXTitle("[cm]");
	}
	TH1F* cdAllHits_F[] = {h_sigma, h_res_MP2, h_dca, h_track_true, h_track_recon, h_residual_true, h_chi2_true, h_residual_recon,
	                       h_chi2_recon, h_driftRad, h_track_TR_diff, h_dca_unsmeared
	                      };
	TH1F* cdTracks_F[] = {h_intercept, h_slope, h_x0, h_x1, h_reconMinusTrue_track_slope, h_reconMinusTrue_track_intercept,
	                      h_recon_slope, h_recon_intercept, h_pval, h_chi2_circle, h_chi2_circle_ndf
	                     };
	TH1I* cdAllHits_I[] = {h_labels, h_hitCount, h_id_dca};
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
	float z_jump_bin = 5.0; // the increment in z for consecutive layers
	for (int i_module = 0; i_module < Tracker::instance()->getModuleN(); i_module++) {
		for (int i_view = 0; i_view < Tracker::instance()->getViewN(); i_view++) {
			for (int i_layer = 0; i_layer < Tracker::instance()->getLayerN(); i_layer++) {

				string UV = Tracker::instance()->getUVmapping(i_view, i_layer); // converting view/layer ID into conventional labels

				h_name.str(""); h_name << "h_dca_M_" << i_module << "_" << UV;
				h_title.str(""); h_title << "DCA M" << i_module << " " << UV ;
				auto hl1 = new TH1F(h_name.str().c_str(), h_title.str().c_str(), 99,  -0.1, Tracker::instance()->getStrawRadius() + 0.25);
				hl1->GetXaxis()->SetTitle("[cm]"); hl1->SetDirectory(cd_UV);

				h_name.str(""); h_name << "h_strawID_M_" << i_module << "_" << UV;
				h_title.str(""); h_title << "strawID M" << i_module << " " << UV ;
				auto hl2 = new TH1I(h_name.str().c_str(), h_title.str().c_str(), Tracker::instance()->getStrawN(), 0, Tracker::instance()->getStrawN());
				hl2->GetXaxis()->SetTitle("Straw ID [0-31]"); hl2->SetDirectory(cd_UV);

				h_name.str(""); h_name << "h_LR_M_" << i_module << "_" << UV;
				h_title.str(""); h_title << "LR in M" << i_module << " " << UV ;
				auto hl3 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  4, -2, 2);
				hl3->GetXaxis()->SetTitle("L= - 1.0; R = +1.0 "); hl3->SetDirectory(cd_UV);

				h_name.str(""); h_name << "h_residual_recon_M_" << i_module << "_" << UV;
				h_title.str(""); h_title << "Residuals Recon M" << i_module << " " << UV ;
				auto hl4 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  149, -0.2, 0.2);
				hl4->GetXaxis()->SetTitle("[cm]"); hl4->SetDirectory(cd_UV);

				h_name.str(""); h_name << "h_pull_M_" << i_module << "_" << UV;
				h_title.str(""); h_title << "Pull M" << i_module << " " << UV ;
				auto hl6 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  149, -15.0, 15.0); hl6->SetDirectory(cd_UV);

				// Parasitic booking with specific binning
				for (int i_LC = 0; i_LC < 2; i_LC++) { //loop over derivatives
					for (int i_LR = 0; i_LR < 2; i_LR++) {
						for (int i_RS = 0; i_RS < 2; i_RS++) {
							h_name.str(""); h_name << "h_" << nameLC[i_LC] << "_M" << i_module << "_" << UV << "_S3_" << valueLR[i_LR] << "_" << nameResSign[i_RS];
							h_title.str(""); h_title << nameLC[i_LC] << "_M" << i_module << UV << "_S3_" << nameLR[i_LR] << "_" << nameResSign[i_RS];
							if (i_LC == 0) {
								if ( (nameResSign[i_RS] == 'P' && nameLR[i_LR] == 'L') || (nameResSign[i_RS] == 'N' && nameLR[i_LR] == 'L')  ) {
									auto hl7 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  249, -1.00001, -0.99995); hl7->SetDirectory(cd_PEDE);
								}
								if ( (nameResSign[i_RS] == 'N' && nameLR[i_LR] == 'R') ||  (nameResSign[i_RS] == 'P' && nameLR[i_LR] == 'R') ) {
									auto hl7 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  249, 0.99995, 1.00001); hl7->SetDirectory(cd_PEDE);
								}
							} // lc1
							if (i_LC == 1) {
								if (nameLR[i_LR] == 'L') {
									auto hl7 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  89, -z_jump_bin - 0.01, -z_jump_bin + 0.01); hl7->SetDirectory(cd_PEDE);
								}
								if (nameLR[i_LR] == 'R') {
									auto hl7 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  89, z_jump_bin - 0.01, z_jump_bin + 0.01); hl7->SetDirectory(cd_PEDE);
								}
							} // lc2
						} // Res sign P/N
					} // LR
				} // lc1-2
				if (i_layer == 0) { z_jump_bin += 0.515; }
			} // layers
			if (i_view == 0) { z_jump_bin += 2.020; }
		} // views
		z_jump_bin += 13.735;
	} // modules

	//Modules
	for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++) {

		h_name.str(""); h_name << "h_DCA_Module_" << i_module;
		h_title.str(""); h_title << "DCA M" << i_module;
		auto hm3 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  99, -0.1, 0.4);
		hm3->GetXaxis()->SetTitle("[cm]"); hm3->SetDirectory(cd_Modules);

		h_name.str(""); h_name << "h_Residuals_Module_" << i_module;
		h_title.str(""); h_title << "Residuals Recon M" << i_module;
		auto hm4 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  199, -0.2, 0.2);
		hm4->GetXaxis()->SetTitle("[cm]"); hm4->SetDirectory(cd_Modules);

		h_name.str(""); h_name << "h_pull_M_" << i_module;
		h_title.str(""); h_title << "Pull M" << i_module;
		auto hm5 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  149, -15.0, 15.0);
		hm5->SetDirectory(cd_Modules);

		h_name.str(""); h_name << "h_DLC2_M" << i_module;
		h_title.str(""); h_title << "DLC2_M" << i_module;
		auto hm7 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  631, -65, 65);
		hm7->SetDirectory(cd_PEDE);

		for (int i_LR = 0; i_LR < 2; i_LR++) {
			h_name.str(""); h_name << "h_Residuals_Module_" << i_module << "_" <<  valueLR[i_LR];
			h_title.str(""); h_title << "Residuals Recon M" << i_module  << " " << nameLR[i_LR];
			auto hm6 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  199, -0.2, 0.2);
			hm6->GetXaxis()->SetTitle("[cm]"); hm6->SetDirectory(cd_Modules);

			h_name.str(""); h_name << "h_DLC1_M" << i_module << "_" << valueLR[i_LR];
			h_title.str(""); h_title << "DLC1_M" << i_module  << " " << nameLR[i_LR];
			if (i_LR == 0) {auto hm8 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  247, -1.00001, -0.99995); hm8->SetDirectory(cd_PEDE);}
			if (i_LR == 1) {auto hm8 = new TH1F(h_name.str().c_str(), h_title.str().c_str(),  247, 0.99995, 1.00001); hm8->SetDirectory(cd_PEDE);}


		}
	}

	// Modules and Straws ["combing 4 layers into 1"]
	for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++) {
		for (int i_straw = 0 ; i_straw < Tracker::instance()->getStrawN(); i_straw++) {
			h_name.str(""); h_name << "h" << i_module << "_straw" << i_straw;
			h_title.str(""); h_title << "DCA S" << i_straw << " M" << i_module;
			auto hmS3 = new TH1F(h_name.str().c_str(), h_title.str().c_str(), 49, -0.1, 0.4);
			hmS3->GetXaxis()->SetTitle("[cm]"); hmS3->SetDirectory(cd_Straws);
		}
	}

//------------------------------------------Mille Routines---------------------------------------------------------//

	// Creating .bin, steering, and constrain files
	Mille M (outFileName, asBinary, writeZero);  // call to Mille.cc to create a .bin file
	helper << "Generating test data for g-2 Tracker Alignment in PEDE:" << endl;

	helper << fixed << setprecision(4);
	// SETTING GEOMETRY
	Tracker::instance()->setGeometry(debug_geom, debugBool);
	helper << "Geometry is set!" << endl << endl;

	// XXX: definition of broken lines here in the future

	// MISALIGNMENT
	Tracker::instance()->misalign(debug_mis, pede_mis, debugBool, metric);
	helper << "Misalignment is complete!" << endl;
	helper << fixed << setprecision(setPrecision);

	// Write a constraint file, for use with pede
	Tracker::instance()->write_constraint_file(constraint_file, debug_con, debugBool, metric);
	helper << "Constraints are written! [see Tracker_con.txt]" << endl;

	//Now writing the steering file
	Tracker::instance()->write_steering_file(steering_file, metric);
	helper << "Steering file was generated! [see Tracker_con.txt]" << endl;

	// Now set pre-sigma for know global parameters
	Tracker::instance()->write_presigma_file(presigma_file, metric);
	helper << "Presigma Parameter file was generated! [see Tracker_par.txt]" << endl;

	helper << "Calculating residuals..." << endl;

//------------------------------------------Main Mille Track Loop---------------------------------------------------------//
	bool StrongDebugBool = false;
	//Generating tracks
	for (int trackCount = 0; trackCount < Tracker::instance()->getTrackNumber(); trackCount++) {
		//float p=pow(10.0, 1+Tracker::instance()->generate_uniform());
		//scatterError=sqrt(Tracker::instance()->getWidth())*0.014/p;
		scatterError = 0; // set no scatterError for now
		char tmpNameResSign; //to assign sign label for the residual [+/-]

		if (debugBool && StrongDebugBool) { helper << "Track: " << trackCount << endl; }

		//Generating tracks
		MCData generated_MC = Tracker::instance()->MC_launch(scatterError, debug_calc, debug_off, debug_mc, plot_fit, plot_gen, plot_hits_gen,
		                      plot_hits_fit, debugBool);

		//First of all, check if track has not failed the cut
		if (generated_MC.cut == false) {

			for (int hitCount = 0; hitCount < generated_MC.hit_count; hitCount++) { //counting only hits going though detector

				//******************PEDE INPUTS******************************
				float rMeas_mp2 =  generated_MC.residuals[hitCount]; //Residual
				float sigma_mp2 = generated_MC.hit_sigmas[hitCount]; //Resolution
				//Number of local and global parameters
				const int nalc = 2;
				const int nagl = 1;
				//label to associate hits within different layers with a correct module
				// same as ModuleN+1 [already converted]
				int label_mp2 = generated_MC.i_hits[hitCount];

				//Variables for derivative calculations:
				float z = generated_MC.z_straw[hitCount];
				float x = generated_MC.x_straw[hitCount];
				float m = generated_MC.slope_recon;
				float c = generated_MC.intercept_recon;

				//Local derivatives
				float dlc1 = ( c + m * z - x ) / ( sqrt(m * m + 1) * abs(c + m * z - x) ) ; // "DCA magnitude"
				float dlc2 = ( (m * m + 1) * z * (c + m * z - x) - m * pow(abs(c + m * z - x), 2) ) / ( pow(m * m + 1, 1.5) * abs(c + m * z - x)  ) ;
				float derlc[nalc] = {dlc1, dlc2};
				//Global derivatives
				float dgl1 = ( c + m * z - x ) / ( sqrt(m * m + 1) * abs(c + m * z - x) ); // =DLC1
				float dergl[nagl] = {dgl1};
				//Labels
				int l1 = label_mp2;
				int label[nagl] = {l1};

				//TODO multiple scattering errors (no correlations) (for imodel == 1)
				//add break points multiple scattering later XXX (for imodel == 2)
				//! add 'broken lines' offsets for multiple scattering XXX (for imodel == 3)

				M.mille(nalc, derlc, nagl, dergl, label, rMeas_mp2, sigma_mp2);
				if (debugBool) {
					debug_mp2  << nalc << " " << derlc[0] << " " << derlc[1] << " " << nagl << " " << dergl[0] << " "  << label[0]  << " "
					           << rMeas_mp2 << "  " << sigma_mp2 << endl;
				}

				//***********************************Sanity Plots******************************************//
				//Getting the "detector coordinates"
				int moduleN = generated_MC.Module_i[hitCount];
				int viewN = generated_MC.View_i[hitCount];
				int layerN = generated_MC.Layer_i[hitCount];
				int strawID = generated_MC.Straw_i[hitCount];
				string UV = Tracker::instance()->getUVmapping(viewN, layerN); // converting view/layer ID into conventional labels

				//Fill for all hits
				h_res_MP2 -> Fill (rMeas_mp2); // residuals
				h_sigma -> Fill(sigma_mp2); // errors
				h_dca->Fill(generated_MC.dca[hitCount]); // DCA
				h_dca_unsmeared->Fill(generated_MC.dca_unsmeared[hitCount]);
				h_labels->Fill(l1);
				h_id_dca ->Fill(strawID);
				h_driftRad->Fill(generated_MC.driftRad[hitCount]);
				if (generated_MC.driftRad[hitCount] < 0) negDCA++;
				h_DLC1->Fill(dlc1); h_DLC2->Fill(dlc2); h_DGL1->Fill(dgl1);

				//Track-based hit parameters
				h_res_x_z->Fill(generated_MC.z_hits[hitCount], generated_MC.residuals[hitCount]);
				h_track_true->Fill(generated_MC.x_track_true[hitCount]);
				h_track_recon->Fill(generated_MC.x_track_recon[hitCount]);
				h_track_TR_diff->Fill(generated_MC.x_track_recon[hitCount] - generated_MC.x_track_true[hitCount]);

				//Calculating Chi2 stats:
				float residual_gen = generated_MC.residuals_gen[hitCount];
				h_residual_true->Fill(residual_gen);
				residuals_true_sum_2 += pow(residual_gen / sigma_mp2, 2);
				h_residual_recon->Fill(rMeas_mp2); //already used as input to mille
				residuals_recon_sum_2 += pow(rMeas_mp2 / sigma_mp2, 2);

				//Fill for hits in modules/layers/straws
				h_name.str(""); h_name << "UV/h_dca_M_" << moduleN << "_" << UV;
				TH1F* h1 = (TH1F*)file->Get( h_name.str().c_str() );
				h1->Fill(generated_MC.dca[hitCount]);

				h_name.str(""); h_name << "UV/h_strawID_M_" << moduleN << "_" << UV;
				TH1I* h2 = (TH1I*)file->Get( h_name.str().c_str() );
				h2 ->Fill(strawID);

				h_name.str(""); h_name << "UV/h_LR_M_" << moduleN << "_" << UV;
				TH1F* h3 = (TH1F*)file->Get( h_name.str().c_str() );
				h3 ->Fill(generated_MC.LR[hitCount]);

				h_name.str(""); h_name << "Straws/h" << moduleN << "_straw" << strawID;
				TH1F* h4 = (TH1F*)file->Get( h_name.str().c_str() );
				h4->Fill(generated_MC.dca[hitCount]);
				h4-> SetFillColor(colourVector[strawID]);

				h_name.str(""); h_name << "Modules/h_DCA_Module_" << moduleN;
				TH1F* h7 = (TH1F*)file->Get( h_name.str().c_str() );
				h7->Fill(generated_MC.dca[hitCount]);

				h_name.str(""); h_name << "Modules/h_Residuals_Module_" << moduleN;
				TH1F* h8 = (TH1F*)file->Get( h_name.str().c_str() );
				h8->Fill(rMeas_mp2);

				h_name.str(""); h_name << "UV/h_residual_recon_M_" << moduleN << "_" << UV;
				TH1F* h9 = (TH1F*)file->Get( h_name.str().c_str() );
				h9->Fill(rMeas_mp2);

				h_name.str(""); h_name << "UV/h_pull_M_" << moduleN << "_" << UV;
				TH1F* h11 = (TH1F*)file->Get( h_name.str().c_str() );
				//h11->Fill( rMeas_mp2 / sqrt( residual_gen - Tracker::instance()->get_sigma_recon_estimated(moduleN, viewN, layerN) ) );
				h11->Fill( rMeas_mp2 / Tracker::instance()->getResolution());

				h_name.str(""); h_name << "Modules/h_pull_M_" << moduleN;
				TH1F* h12 = (TH1F*)file->Get( h_name.str().c_str() );
				h12->Fill( rMeas_mp2 / Tracker::instance()->getResolution());

				h_name.str(""); h_name << "Modules/h_Residuals_Module_" << moduleN << "_" << generated_MC.LR[hitCount];
				TH1F* h13 = (TH1F*)file->Get( h_name.str().c_str() );
				h13->Fill(rMeas_mp2);

				h_name.str(""); h_name << "PEDE/h_DLC2_M" << moduleN;
				TH1F* h14 = (TH1F*)file->Get( h_name.str().c_str() );
				h14->Fill(dlc2);

				h_name.str(""); h_name << "PEDE/h_DLC1_M" << moduleN << "_" << generated_MC.LR[hitCount];
				TH1F* h15 = (TH1F*)file->Get( h_name.str().c_str() );
				h15->Fill(dlc1);

				// Now fill DLC1, DLC2 plots for U0-V1 for each module but only for straw 3 (4th straw from the "top")
				if (strawID == 3) {
					if (rMeas_mp2 > 0) {tmpNameResSign = 'P';}// DCA > driftRad
					else {tmpNameResSign = 'N';} // DCA < driftRad

					//loop over booked histos and fill
					for (int i_LC = 0; i_LC < 2; i_LC++) {
						h_name.str(""); h_name << "PEDE/h_" << nameLC[i_LC] << "_M" << moduleN << "_" << UV << "_S3_" << generated_MC.LR[hitCount] << "_" << tmpNameResSign;
						TH1F* h16 = (TH1F*)file->Get( h_name.str().c_str() );
						if (i_LC == 0) {h16->Fill(dlc1);}
						if (i_LC == 1) {h16->Fill(dlc2);}
					}
				} // strawID=3

				hitsN++; //count hits
			} // end of hits loop

			//***********************************Sanity Plots: Once per Track******************************************/
			//For true tracks
			float chi2_true = residuals_true_sum_2;
			h_chi2_true->Fill(chi2_true);
			//For recon tracks
			float chi2_recon = residuals_recon_sum_2;
			h_chi2_recon->Fill(chi2_recon);
			//Resetting counters for next track
			residuals_true_sum_2 = 0;
			residuals_recon_sum_2 = 0;
			h_hitCount->Fill(generated_MC.hit_count);

			//Filling Track-based plots
			h_slope->Fill(generated_MC.slope_truth);
			h_intercept->Fill(generated_MC.intercept_truth);
			h_x0->Fill(generated_MC.x0);
			h_x1->Fill(generated_MC.x1);
			h_reconMinusTrue_track_intercept->Fill(generated_MC.intercept_truth - generated_MC.intercept_recon);
			h_reconMinusTrue_track_slope->Fill(generated_MC.slope_truth - generated_MC.slope_recon);
			h_recon_slope->Fill(generated_MC.slope_recon);
			h_recon_intercept->Fill(generated_MC.intercept_recon);
			h_pval->Fill(generated_MC.p_value);
			h_chi2_circle->Fill(generated_MC.chi2_circle);
			h_chi2_circle_ndf->Fill(generated_MC.chi2_circle / (generated_MC.hit_count - 2));

			// XXX additional measurements from MS IF (imodel == 2) THEN
			//IF (imodel >= 3) THEN

			M.end(); // Write buffer (set of derivatives with same local parameters) to file.
			recordN++; // count records (i.e. written tracks);

		} // cut on DCA check

	} // end of track count // End of Mille // End of collecting residual records
	helper << "Mille residual-accumulation routine completed! [see Tracker_data.bin]" << endl;

	//Passing constants to plotting script

	contsants_plot << Tracker::instance()->getModuleN() << " " << Tracker::instance()->getViewN() << " "
	               << Tracker::instance()->getLayerN() << " " << Tracker::instance()->getStrawN() << " " << recordN << " "
	               << Tracker::instance()->getBeamOffset()   << " " << Tracker::instance()->getBeamStart() << " " <<  Tracker::instance()->getBeamPositionLength()
	               << "  " << Tracker::instance()->getBeamStop() <<  endl;



	//------------------------------------------ROOT: Fitting Functions---------------------------------------------------------//

	helper << endl;
	helper << "-------------------------------------------------------------------------" << endl;
	helper << "ROOT fitting parameters and output:" << endl;

	helper << fixed << setprecision(2);
	float LRSkewness = 0.0;
	// PEDE Plots calculations for U0 only
	for (int i_module = 0; i_module < Tracker::instance()->getModuleN(); i_module++) {
		for (int i_view = 0; i_view < 1; i_view++) {
			for (int i_layer = 0; i_layer < 1; i_layer++) {
				string UV = Tracker::instance()->getUVmapping(i_view, i_layer);
				for (int i_LR = 0; i_LR < 2; i_LR++) {
					for (int i_RS = 0; i_RS < 2; i_RS++) {
						h_name.str(""); h_name << "PEDE/h_DLC2_M" << i_module << "_" << UV << "_S3_" << valueLR[i_LR] << "_" << nameResSign[i_RS];
						TH1F* hp1 = (TH1F*)file->Get( h_name.str().c_str() );
						LRSkewness += hp1->GetSkewness();
					}
					//helper << "M" << i_module << UV << "_S3_" << nameLR[i_LR] << " :: MeanSkewness= " << LRSkewness / 2.0 << endl;
					LRSkewness = 0.0;
				}
			}
		}
	}
	helper << fixed << setprecision(setPrecision);

	// Store alignment parameters from measurements
	vector<float> sigma_recon_actual;
	vector<float> sigmaError_recon_actual;
	vector<float> pull_actual;
	vector<float> pull_actual_SD;
	vector<float> Res_mean;
	vector<float> Res_mean_SD;

	//Filling TH2 for residual SD and pulls
	int z_counter = 0;
	for (int i_module = 0; i_module < Tracker::instance()->getModuleN(); i_module++) {
		for (int i_view = 0; i_view < Tracker::instance()->getViewN(); i_view++) {
			for (int i_layer = 0; i_layer < Tracker::instance()->getLayerN(); i_layer++) {
				string UV = Tracker::instance()->getUVmapping(i_view, i_layer);
				// Residuals
				h_name.str(""); h_name << "UV/h_residual_recon_M_" << i_module << "_" << UV;
				TH1F* hRes_actual = (TH1F*)file->Get( h_name.str().c_str() );
				sigma_recon_actual.push_back(hRes_actual->GetStdDev() * 10000); // um ->cm
				sigmaError_recon_actual.push_back(hRes_actual->GetStdDevError() * 10000); // um ->cm
				h_SD_z_res_Recon->Fill(Tracker::instance()->getZDistance(z_counter), hRes_actual->GetStdDev() * 10000); // um ->cm
				h_SD_z_res_Recon->SetBinError(Tracker::instance()->getZDistance(z_counter), hRes_actual->GetStdDev() * 10000, hRes_actual->GetStdDevError() * 10000); // um ->cm
				h_SD_z_res_Recon->SetMarkerStyle(33); h_SD_z_res_Recon->SetMarkerColor(kRed);
				Res_mean.push_back(hRes_actual->GetMean());
				Res_mean_SD.push_back(hRes_actual->GetStdDev());
				h_SD_z_res_Est->Fill(Tracker::instance()->getZDistance(z_counter) , Tracker::instance()->get_sigma_recon_estimated(i_module, i_view, i_layer) * 10000); // um ->cm
				h_SD_z_res_Est->SetMarkerStyle(33);
				h_SD_z_res_Est->SetMarkerColor(kBlue);
				// Pulls
				h_name.str(""); h_name << "UV/h_pull_M_" << i_module << "_" << UV;
				TH1F* h_pull = (TH1F*)file->Get( h_name.str().c_str() );
				pull_actual.push_back(h_pull->GetMean());
				pull_actual_SD.push_back(h_pull->GetStdDev());
				h_Pulls_z->Fill(Tracker::instance()->getZDistance(z_counter), h_pull->GetStdDev());
				h_Pulls_z->SetBinError(Tracker::instance()->getZDistance(z_counter), h_pull->GetStdDev(), h_pull->GetStdDevError());
				h_Pulls_z->SetMarkerStyle(33); h_Pulls_z->SetMarkerColor(kRed);
				z_counter++;
			}// layer
		} // view
	} // modules

	//Dealing with Chi2 and residual [true] fits
	//Use Pearson chi-square method, using expected errors instead of the observed one given by TH1::GetBinError (default case).
	//The expected error is instead estimated from the the square-root of the bin function value.
	TF1* chi2pdf = new TF1("chi2pdf", "[2]*ROOT::Math::chisquared_pdf(x,[0],[1])", 0, 40);
	chi2pdf->SetParameters(Tracker::instance()->get_Chi2_recon_estimated(), 0., h_chi2_true->Integral("WIDTH"));
	h_chi2_true->Fit("chi2pdf", "Q");
	h_residual_true->Fit("gaus", "Q");
	Chi2_recon_actual = h_chi2_recon->GetMean();

	//Residuals SD per layer
	vector<float> zDistance = Tracker::instance()->getZDistanceVector();
	vector<float> sigma_recon_estimated = Tracker::instance()->get_sigma_recon_estimatedVector();
	TCanvas *cRes = new TCanvas("cRes", "cRes", 700, 500);
	const Int_t n = Tracker::instance()->getLayerTotalN();
	Float_t* z_distance  = &zDistance[0];
	Float_t* Res_Recon_SD  = &sigma_recon_actual[0];
	Float_t* Res_Recon_SD_error = &sigmaError_recon_actual[0];
	auto gr = new TGraphErrors(n, z_distance, Res_Recon_SD, 0, Res_Recon_SD_error);
	gr->SetTitle("Residuals SD");
	gr->SetMarkerColor(kWhite);
	gr->SetLineColor(kRed);
	gr->SetMarkerStyle(1);
	gr->Draw("A*");
	Float_t* Res_Est_SD  = &sigma_recon_estimated[0];
	for (int i = 0; i < n; i++) {
		TMarker *m1 = new TMarker(z_distance[i], Res_Est_SD[i] * 10000, 20);
		m1->SetMarkerColor(kBlue); m1->Draw();
	}
	gr->GetXaxis()->SetTitle("Module/Layer separation [cm]");
	gr->GetYaxis()->SetTitle("Residual SD [um]");
	auto axis = gr->GetXaxis();
	axis->SetLimits(0., 18 * Tracker::instance()->getModuleN());              // along X
	gr->GetHistogram()->SetMaximum(148.);   // along
	gr->GetHistogram()->SetMinimum(130.);  //   Y
	//cRes->Print("Print/cRes.png");
	cRes->Write();

	//Pulls per layer
	TCanvas *cPull = new TCanvas("cPull", "cPull", 200, 10, 700, 500);
	Float_t* Pull_Recon  = &pull_actual[0];
	Float_t* Pull_Recon_SD = &pull_actual_SD[0];
	auto gr2 = new TGraphErrors(n, z_distance, Pull_Recon, 0, Pull_Recon_SD);
	gr2->SetTitle("Pulls [Error = SD]");
	gr2->SetMarkerColor(kWhite);
	gr2->SetLineColor(kRed);
	gr2->SetMarkerStyle(1);
	gr2->Draw("A*");
	gr2->GetXaxis()->SetTitle("Module/Layer separation [cm]");
	gr2->GetYaxis()->SetTitle("Pulls per layer [Error = SD]");
	auto axis2 = gr2->GetXaxis();
	axis2->SetLimits(0., 18 * Tracker::instance()->getModuleN());              // along X
	gr2->GetHistogram()->SetMaximum(5.0);   // along
	gr2->GetHistogram()->SetMinimum(-3.0);  //   Y
	//cPull->Print("Print/cPull.png");
	cPull->Write();

	//Residual means per layer: reveal the shear affect of misalignment
	TCanvas *cResMean = new TCanvas("cResMean", "cResMean", 200, 10, 700, 500);
	Float_t* Res_meanR  = &Res_mean[0];
	Float_t* Res_mean_SDR = &Res_mean_SD[0];
	auto gr3 = new TGraphErrors(n, z_distance, Res_meanR, 0, Res_mean_SDR);
	gr3->SetTitle("Residual means [Error = SD]");
	gr3->SetMarkerColor(kWhite);
	gr3->SetLineColor(kRed);
	gr3->SetMarkerStyle(1);
	gr3->Draw("A*");
	gr3->GetXaxis()->SetTitle("Module/Layer separation [cm]");
	gr3->GetYaxis()->SetTitle("Residual means per layer [Error = SD] [cm]");
	for (int i = 0; i < n; i++) {
		TMarker *m2 = new TMarker(z_distance[i], Tracker::instance()->get_shearMis(i), 20);
		m2->SetMarkerColor(kBlue);
		m2->Draw();
	}
	auto axis3 = gr3->GetXaxis();
	axis3->SetLimits(0., 18 * Tracker::instance()->getModuleN());              // along X
	gr3->GetHistogram()->SetMaximum(0.2);   // along
	gr3->GetHistogram()->SetMinimum(-0.2);  //   Y
	//cResMean->Print("Print/cResMean.png");
	cResMean->Write();

	// if (recordN >= 100) { // below 100 tracks these plots are useless
	// 	TCanvas *cChi2 = new TCanvas("cChi2", "cChi2", 700, 700);
	// 	gStyle->SetOptStat("ourRmMe");
	// 	gStyle->SetOptFit(1111);
	// 	chi2pdf->SetParameters(Tracker::instance()->get_Chi2_recon_estimated(), 0., h_chi2_recon->Integral("WIDTH"));
	// 	h_chi2_recon->SetBinErrorOption(TH1::kPoisson); // errors from Poisson interval at 68.3% (1 sigma)
	// 	h_chi2_recon->Fit("chi2pdf", "Q");
	// 	cChi2->Clear(); // Fit does not draw into correct pad
	// 	auto rp1 = new TRatioPlot(h_chi2_recon, "errasym");
	// 	rp1->SetGraphDrawOpt("P");
	// 	rp1->SetSeparationMargin(0.0);
	// 	cChi2->SetTicks(0, 1);
	// 	rp1->Draw("noconfint");
	// 	cChi2->Update();
	// 	rp1->GetLowerRefYaxis()->SetTitle("Frac. Error");
	// 	//cChi2->Print("FoM_Chi2_recon.C");
	// 	//cChi2->Print("Print/FoM_Chi2_recon.png");
	// 	cChi2->Write();
	// 	//TODO fix malloc problem closing the canvas in ROOT Browser
	// 	delete chi2pdf;
	// }

	// Debug-style plots:
	if (1 == -1) {
		//Residuals per module
		TCanvas *cResAllM = new TCanvas("cResAllM", "cResAllM", 700, 900);
		TText T; T.SetTextFont(42); T.SetTextAlign(21);
		for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++) {
			h_name.str(""); h_name << "Modules/h_Residuals_Module_" << i_module;
			TH1F* hd1 = (TH1F*)file->Get( h_name.str().c_str() );
			hd1->SetFillColor(colourVector[i_module]);
			hd1->Draw("same");
			gStyle->SetOptStat("");
			hd1->SetTitle("");
		}
		T.DrawTextNDC(.5, .95, "Residuals");
		cResAllM->Write();


		//Stacked DCA per straw in each module
		vector<THStack*> stackModule;
		for (int i_module = 0; i_module < Tracker::instance()->getModuleN(); i_module++) {
			h_name.str(""); h_name << "DCA_Module_" << i_module;
			THStack* stack = new THStack(h_name.str().c_str(), " ");
			stackModule.push_back(stack);
		}
		//THStack* stackModule[] = {hs_DCA_Module_0, hs_DCA_Module_1, hs_DCA_Module_2, hs_DCA_Module_3, hs_DCA_Module_4, hs_DCA_Module_5, hs_DCA_Module_6, hs_DCA_Module_7};
		TCanvas *cDCA = new TCanvas("cDCA", "cDCA", 700, 900);
		T.SetTextFont(42); T.SetTextAlign(21);
		cDCA->Divide((Tracker::instance()->getModuleN() + 1) / 2, (Tracker::instance()->getModuleN() + 1) / 2);
		for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++) {
			for (int i_straw = 0 ; i_straw < Tracker::instance()->getStrawN(); i_straw++) {
				h_name.str(""); h_name << "Straws/h" << i_module << "_straw" << i_straw;
				TH1F* hs3 = (TH1F*)file->Get( h_name.str().c_str() );
				stackModule[i_module]->Add(hs3);
			}
		}
		for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++) {
			h_title.str(""); h_title << "Module " << i_module << " DCA per straw";
			cDCA->cd(i_module + 1); stackModule[i_module]->Draw(); T.DrawTextNDC(.5, .95, h_title.str().c_str());
		}
		cDCA->Write();

		//Stacked reconstructed hits
		TCanvas *cStack = new TCanvas("cStack", "cStack", 700, 900);
		T.SetTextFont(42); T.SetTextAlign(21);
		for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++) {
			h_name.str(""); h_name << "Modules/hs_hits_recon_Module" << i_module;
			TH1F* hs4 = (TH1F*)file->Get( h_name.str().c_str() );
			hs_hits_recon->Add(hs4);
		}
		cStack->Divide(1, 1);
		cStack->cd(1);  hs_hits_recon->Draw(); T.DrawTextNDC(.5, .95, "Recon Hits per Module"); hs_hits_recon->GetXaxis()->SetTitle("[cm]");
		cStack->Write();
	} // strong plotting

	helper << "-------------------------------------------------------------------------" << endl;
	helper << " " << endl;
	helper << Tracker::instance()->getTrackNumber() << " tracks requested; " << recordN << " generated with " << hitsN << " hits." << endl;
	float rejectedTracks = (Tracker::instance()->getTrackNumber() - recordN) / float(Tracker::instance()->getTrackNumber()) * 100.0;
	helper << Tracker::instance()->getTrackNumber() - recordN << " records rejected (" << rejectedTracks << " %)." << endl;
	helper << "Number of DCA (==drift radii) smeared below 0 is " << negDCA << endl;
	stringstream out1, out2, out3, out4, out5;
	out1 << "Expected Mean Chi2 (for a general straight line fit) " << Tracker::instance()->get_Chi2_recon_estimated();
	out2 << "Measured Mean Chi2 (circle fit) " << Chi2_recon_actual;
	Logger::Instance()->write(Logger::WARNING, out1.str());
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
	helper << "Normal random numbers were used " << RandomBuffer::instance()->getNormTotal() << " times" << endl;
	helper << "Gaussian random numbers were used " << RandomBuffer::instance()->getGausTotal() << " times" << endl;
	if (debugBool) {
		Logger::Instance()->write(Logger::WARNING, "Text debug files were produced: ls Tracker_d_*.txt");
	}
	if (plotBool) {
		Logger::Instance()->write(Logger::WARNING, "./PlotMP2.py to draw geometry with hits");
	}

	// Close text files
	constraint_file.close();
	steering_file.close();
	debug_mp2.close();
	debug_calc.close();
	debug_mis.close();
	debug_geom.close();
	debug_off.close();
	debug_mc.close();
	debug_con.close();

	file->Write();
	file->Close();
	delete file;

	metric << "| R: " << Tracker::instance()->getResolution() * 1e4 << " um "
	       << "| DCA Cut of " << Tracker::instance()->getTrackCut() * 1e4 << " um : " << boolYN[Tracker::instance()->getTrackCutBool()]
	       << "| Hit rej.: " << boolYN[Tracker::instance()->getHitCutStatus()]
	       << "| Truth LR : " << boolYN[Tracker::instance()->getLRStatus()]
	       << "| p-value cut (<): " << Tracker::instance()->getPValCut() ;

	helper << endl;
	helper << "Programme log written to: Tracker_log.txt" << endl;
	helper << fixed << setprecision(4);
	t_cpu = clock() - t_cpu;
	auto t_end = chrono::high_resolution_clock::now();
	helper << "Programme execution took " <<  t_cpu << " CPU clicks (" << ((float)t_cpu) / CLOCKS_PER_SEC << " s)." << " Wall clock time passed: "
	       << chrono::duration<double>(t_end - t_start).count() << " s." << endl;
	timeFile << chrono::duration<double>(t_end - t_start).count() << endl;
	time_t now = time(0);
	char* dt = ctime(&now);
	//helper << "Peak RAM use: " << Tracker::instance()->getPeakRSS( ) / 1e9 << " GB" << endl;
	helper << "The C++ compiler used: " << true_cxx << " " << true_cxx_ver
	       << " Job finished on: " << dt << endl;
	return 0;
} //end of main