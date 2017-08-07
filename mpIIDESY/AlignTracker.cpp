/*
*   Gleb Lukicov (g.lukicov@ucl.ac.uk) @ Fermilab
*   Created: 17 April 2017  
*   Modified: 4 August 2017 
----------------------------------------------------------------
This programme uses MC methods to produce a .bin data file for the 
PEDE routine, to align the tracking detector for the g-2 
experiment.
Methods and functions are contained in AlignTracker_methods.cpp (Tracker class)

============================= Test Model v0.3 ======================================= 
*Simple 2D case (1D alignment along x).
* No B-field, straight tracks, no MS, 100% efficiency. 
*(!) NATURAL UNITS (same as PEDE): cm, rad, GeV [natural units will be omitted] 

(!) Module, Straw, View, Layer labelling starts from 0 [+1 is put by hand as MP2 doesn't accept 0 as a label]

Terminology: detector=module (==the alignment object)

Gaussian (gaus) Smearing: Mean=0, RMS=1  (on hits, resolution, etc.)
Beam originates (in z direction) from z=0, straight tracks, no scattering, only resolution (gaus) smearing,
    tracks are lines between z=beamStart and z=beamStop, with x1 and x2 = Constant * rand[0,1]
Resolution (i.e. tracker systematic) is implemented. 
Misalignment is implemented "manually" for each module in +/- x-direction [all layers in module are dependent (i.e. move together)]

Hit rejection:  -not implemented 

Constrains TODO 
Steering options TODO 
                       			--Tracker Geometry--:
Spacing (in x) between straws in layer is 0.606 [strawSpacing]
Spacing between layers in same view is 0.515 [layerSpacing]
Spacing between U and V views (last and first layers in views) is 2.020 [viewSpaing]
Spacing between modules (last and first layers in modules) 13.735 [moduleSpacing]
Layers in a view are (e.g. U0 vs. U1) have an extra relative x displacement of 0.303  [layerDisplacement]

For more information see: 
http://gm2-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=4375&filename=Orientation.pdf&version=1 

=======================================================================================

#TODO Model type is implemented via imodel = 0, 1 .....  

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
  std::ostringstream ss;
  ss << a << '.' << b << '.' << c;
  return ss.str();
}

//***************MAIN****************//
int main(int argc, char* argv[]){
    clock_t t_cpu; // CPU ticks for programme execution
    t_cpu = clock(); 
    auto t_start = std::chrono::high_resolution_clock::now(); // Wall clock ticks

    //Determining compiler used:
    std::string true_cxx = 
    #ifdef __clang__
        "clang++";
    #else
        "g++";
     #endif

    std::string true_cxx_ver =
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
    bool debugBool = false; // './AlignTracker n' - for normal, of ./AlignTracker d' - for verbose debug output
    bool plotBool = false; // './AlignTracker p' - for plotting with PlotGen.py 
    bool strongPlotting = false;  // For canvas to be drawn and saved as .png   // XXX HACK (set by hand here)
    //Set up counters for hits and records (tracks)
    int hitsN = 0; // actually recorded (i.e. non-rejected hits)
    int recordN=0; //records = tracks
    float scatterError; // multiple scattering error [calculated here and passed back to the Tracker class]
    float residuals_track_sum_2=0.0;
    float residuals_fit_sum_2=0.0;
    float sigma_fit_calc;  // estimated value for RMS on the residuals for the fit [need global scope for printout]
    const Color_t colourVector[]={kMagenta, kOrange, kBlue, kGreen, kYellow, kRed, kGray, kBlack};
       
    //Tell the logger to only show message at INFO level or above
    Logger::Instance()->setLogLevel(Logger::NOTE); 
    //Tell the logger to throw exceptions when ERROR messages are received
    Logger::Instance()->enableCriticalErrorThrow();

    //TApplication theApp("App",&argc, argv); //for Canvas display 

    // Check if correct number of arguments specified, exiting if not
    if (argc > 3) { Logger::Instance()->write(Logger::ERROR, "Too many arguments -  please specify verbosity flag. (e.g. debug [d], plot[p] or align/normal [n])");} 
        else if (argc < 3) {Logger::Instance()->write(Logger::ERROR, "Too few arguments - please specify verbosity flag. (e.g. e.g. debug [d], plot[p] or align/normal [n])");} 
        else { // Set filenames to read random numbers from, using arguments. Catch exception if these files do not exist.
            try {
                compareStr = argv[1];
                tracksInput = stoi(argv[2]);
            } 
        catch (ios_base::failure& e) {
            Logger::Instance()->write(Logger::ERROR, "Exception caught: " + string(e.what()) + "\nPlease ensure valid verbosity level specified!");
        }
    } // end of 2nd else [correct # arguments]

    //this is also passed to Tracker functions, with debug file names
    if (compareStr=="d"){
    debugBool = true; // print out to debug files [and verbose cout output]
    Tracker::instance()->setTrackNumber(tracksInput);
    Logger::Instance()->write(Logger::WARNING,  "DEBUG MODE"); }
    else if (compareStr=="p"){
    plotBool = true; 
    debugBool = true; // print out to debug files [and verbose cout output]
    Tracker::instance()->setTrackNumber(tracksInput); 
    Logger::Instance()->write(Logger::WARNING,  "PLOTTING MODE"); //setting small statistics for visual plotting of tracks
    } 
    else if (compareStr=="n" || compareStr=="a"){
    debugBool = false; // print out to debug files
    plotBool = false;  
    Tracker::instance()->setTrackNumber(tracksInput); 
    }
    else{
        Logger::Instance()->write(Logger::ERROR, "Please specify verbosity flag. (e.g. debug [d], plot[p] or align/normal [n])");
    }
    
    Logger::Instance()->setUseColor(false); // will be re-enabled below [to use custom colour output to terminal]
    std::stringstream msg0, msg01, msg02, msg1;
    Logger::Instance()->write(Logger::NOTE, "");
    msg0 << Logger::blue() <<  "*****************************************************************" << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg0.str());
    msg01 << Logger::yellow() << "  g-2 Tracker Alignment (v0.3) - Gleb Lukicov (UCL) - August 2017            " << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg01.str());
    msg1 << Logger::blue() <<  "*****************************************************************" << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg1.str());
    Logger::Instance()->setUseColor(true); // back to default colours 

    cout << "Simple Alignment Model with " << Tracker::instance()->getModuleN() << " tracker modules, having " << Tracker::instance()->getStrawN() 
    	<< " straws per layer." << endl;
    cout << "["<< Tracker::instance()->getLayerN() << " layers per module; " << Tracker::instance()->getViewN() << " views per module]." << endl;
    cout << "Total of " << Tracker::instance()->getLayerTotalN() << " measurement layers." << endl; 
    cout << "No B-field, Straight Tracks (general lines), 100% efficiency." << endl;
    cout << "No Hit rejection:" << endl;
    // DCA > StrawSpacing [" << Tracker::instance()->getStrawSpacing() << " cm]." << endl;
    cout << "Parallel Tracks: single hit per layer allowed [shortest DCA is chosen as the hit]." << endl;
    cout << "Resolution is " << Tracker::instance()->getResolution() << " cm  [hit smearing]." << endl;  


    // See https://github.com/glukicov/alignTrack for instructions to generate random numbers
    try {
        Tracker::instance()->set_uniform_file("uniform_ran.txt"); 
        Tracker::instance()->set_gaussian_file("gaussian_ran.txt");
        } 
        
    catch (ios_base::failure& e) {
        cerr << "Filestream exception caught: " << e.what() << endl;
        cerr << "Please ensure valid filenames are specified!" << endl; 
        return 1;
        } 
 //------------------------------------------Setting Output Files---------------------------------------------------------//     
 
    //arguments for Mille constructor:
    const char* outFileName = "Tracker_data.bin"; 
    bool asBinary = true; // set false for debugging
    bool writeZero = false; // to write 0 LC/DLC labels - not accepted by pede 
    
    //Constraints and Steering files are the inputs to Pede (together with binary file). 
    string conFileName = "Tracker_con.txt";
    string strFileName = "Tracker_str.txt";
    ofstream constraint_file(conFileName);
    ofstream steering_file(strFileName);
    cout << fixed << setprecision(setPrecision); // set precision of standard screen output

    //Debug files [only filled with "d" option]: file streams for debug files 
    ofstream debug_mp2("Tracker_d_mille.txt");               //Inputs into binary files
    ofstream debug_calc("Tracker_d_calc.txt");    // intermediate calculation 
    ofstream debug_mis("Tracker_d_mis.txt");    // Misalignment 
    ofstream debug_geom("Tracker_d_geom.txt");    // Geometry 
    ofstream debug_off("Tracker_d_off.txt");    // Offsets/Missed hits
    ofstream debug_mc("Tracker_d_MC.txt");     // Final results from MC
    ofstream debug_con("Tracker_d_con.txt");    //Constraints
    ofstream plot_gen("Tracker_p_gen.txt");    //Truth Track points
    ofstream plot_fit("Tracker_p_fit.txt");    //Reconstructed Track points
    ofstream contsants_plot("Tracker_p_constants.txt");    // passing constants (e.g. strawN to python script)
    ofstream plot_hits_gen("Tracker_p_hits_gen.txt");    // Truth Hits points
    ofstream plot_hits_fit("Tracker_p_hits_fit.txt");    // Truth Hits points
    ofstream pede_mis("Tracker_pede_mis.txt");    // Misalignments 
    ofstream bias; bias.open("bias.txt", ios_base::app);
    
    // Setting fixed precision for floating point values
    debug_mp2 << fixed << setprecision(setPrecision); 
    debug_calc << fixed << setprecision(setPrecision); 
    debug_mis << fixed << setprecision(setPrecision); 
    debug_geom << fixed << setprecision(setPrecision); 
    debug_off << fixed << setprecision(setPrecision); 
    debug_mc << fixed << setprecision(setPrecision); 
    debug_con << fixed << setprecision(setPrecision);
    plot_gen << fixed << setprecision(setPrecision);
    plot_fit << fixed << setprecision(setPrecision);
    contsants_plot << fixed << setprecision(setPrecision);
    pede_mis << fixed << setprecision(setPrecision);
  
   //------------------------------------------ROOT: Booking etc.---------------------------------------------------------//   

    //output ROOT file
    TFile* file = new TFile("Tracker.root", "recreate");  // recreate = overwrite if already exists
    //create a subdirectories 
    TDirectory* cd_All_Hits = file->mkdir("All_Hits");
    TDirectory* cd_UV = file->mkdir("UV");
    TDirectory* cd_Modules = file->mkdir("Modules");
    TDirectory* cd_Tracks = file->mkdir("Tracks");
    TDirectory* cd_Straws = file->mkdir("Straws");
    // Book histograms [once only]
    // Key quantities 
    TH1F* h_sigma = new TH1F("h_sigma", "MP2 Input: Detector Resolution (sigma) [cm]",  49,  Tracker::instance()->getResolution()-0.001, 
    	Tracker::instance()->getResolution()+0.001); // F=float bins, name, title, nBins, Min, Max
    TH1F* h_hits_MP2 = new TH1F("h_hits_MP2", "MP2 Input: Residuals from fitted line to ideal geometry [cm]",  99, -0.1, 0.1);
    TH1F* h_dca = new TH1F("h_dca", "DCA (to misaligned detector from generated track)",  149,  -0.05, Tracker::instance()->getStrawRadius()+0.25);
    TH1I* h_id_dca = new TH1I("h_id_dca", "ID for hit straws", Tracker::instance()->getStrawN(), 0, Tracker::instance()->getStrawN());
    // Track-generation-based
    TH1F* h_slope = new TH1F("h_slope", "Slope ",  80,  -0.05, 0.05);
    TH1F* h_intercept = new TH1F("h_intercept", "Intercept ",  99,  -3, 3);
    TH1F* h_x0 = new TH1F("h_x0", "Generated x0 of the track",  99,  -3, 3);
    TH1F* h_x1 = new TH1F("h_x1", "Generated x1 of the track",  99,  -3, 3);
    //Track/Hits-based
    TH1F* h_track_true = new TH1F("h_track_true", "True track position (the x of the generated track in-line with a layer)",  49,  -3, 3);
    TH1F* h_hits_true = new TH1F("h_hits_true", "True hit position (the x of the generated and smeared hit)",  149,  -3, 3);
    TH1F* h_hits_recon = new TH1F("h_hits_recon", "Reconstructed x position of the hits in ideal detector",  149,  -3, 3);
    TH1F* h_track_recon = new TH1F("h_track_recon", "Reconstructed x of the fitted track (to ideal geometry)",  149,  -(Tracker::instance()->getBeamOffset()+3), 
    	Tracker::instance()->getBeamPositionLength()+1);
    TH1I* h_labels = new TH1I("h_labels", "Labels in PEDE", 8 , 0, 8);
    TH1F* h_residual_track = new TH1F("h_residual_track", "Residuals for generated tracks", 500, -0.4, 0.4);
    TH1F* h_chi2_track = new TH1F("h_chi2_track", "Chi2 for generated tracks", 40, -1, 100);
    TH1F* h_residual_fit = new TH1F("h_residual_fit", "Residuals for fitted tracks", 500, -0.2, 0.2);
    TH1F* h_chi2_fit = new TH1F("h_chi2_fit", "Chi2 for fitted tracks", 59, -1, 250);
    TH1I* h_hitCount = new TH1I("h_hitCount", "Total Hit count per track", 32 , 0, 32);
    TH1F* h_reconMinusTrue_track = new TH1F("h_reconMinusTrue_line", "Reconstructed - True X position of the lines",  149,  -0.1, 0.1);
    TH1F* h_reconMinusTrue_hits = new TH1F("h_reconMinusTrue_hits", "Reconstructed - True X position of the hits",  169,  -0.1, 0.2);
    TH1F* h_reconMinusTrue_track_slope = new TH1F("h_reconMinusTrue_track_slope", "Reconstructed - True Track Slope",  99,  -0.002, 0.002);
    TH1F* h_reconMinusTrue_track_intercept = new TH1F("h_reconMinusTrue_track_intercept", "Reconstructed - True Track Intercept",  39,  -0.06, 0.07);
    TH1F* h_frac_Dslope = new TH1F("h_frac_Dslope", "(True-Recon)/True Track slope",  199,  -1.1, 1.1);
    TH1F* h_frac_Dintercept = new TH1F("h_frac_Dintercept", "(True-Recon)/True Track intercept",  199,  -1.1, 1.1);
    TH1F* h_meanXRecon = new TH1F("h_meanXRecon", "Mean X of recon track", 39, -2.2, 2.2);

    // "special" histos
    TH2F* h_res_x_z = new TH2F("h_res_x_z", "Residuals vs z", 600, 0, 60, 49, -0.08, 0.08);
    h_res_x_z->SetDirectory(cd_All_Hits); h_res_x_z->GetXaxis()->SetTitle("cm");  h_res_x_z->GetYaxis()->SetTitle("cm");
    THStack* hs_hits_recon = new THStack("hs_hits_recon", "");
    
    std::stringstream h_name;
    std::stringstream h_title;
    
    // Modules, views, layers
	for (int i_module=0; i_module< Tracker::instance()->getModuleN(); i_module++){
	    for (int i_view=0; i_view< Tracker::instance()->getViewN(); i_view++){
		    for (int i_layer=0; i_layer< Tracker::instance()->getLayerN(); i_layer++){

		    	string UV = Tracker::instance()->getUVmapping(i_view, i_layer); // converting view/layer ID into conventional labels

		    	h_name.str(""); h_name << "h_dca_M_" << i_module << "_" << UV;
		        h_title.str(""); h_title << "DCA in Module " << i_module << " " << UV ;
		        auto hl1 = new TH1F(h_name.str().c_str(),h_title.str().c_str(), 100,  -0.05, Tracker::instance()->getStrawRadius()+0.25);
		        hl1->GetXaxis()->SetTitle("[cm]"); hl1->SetDirectory(cd_UV);
		        
		        h_name.str(""); h_name << "h_strawID_M_" << i_module << "_" << UV;
		        h_title.str(""); h_title << "strawID in Module " << i_module << " " << UV ;
		        auto hl2 = new TH1I(h_name.str().c_str(),h_title.str().c_str(), Tracker::instance()->getStrawN(), 0, Tracker::instance()->getStrawN());
		        hl2->GetXaxis()->SetTitle("Straw ID [0-31]"); hl2->SetDirectory(cd_UV);
		        
		        h_name.str(""); h_name << "h_LR_M_" << i_module << "_" << UV;
		        h_title.str(""); h_title << "Left-Right hit distribution in Module " << i_module << " " << UV ;
		        auto hl3 = new TH1F(h_name.str().c_str(),h_title.str().c_str(),  4, -2, 2);
		        hl3->GetXaxis()->SetTitle("L= - 1.0; R = +1.0 "); hl3->SetDirectory(cd_UV);

		    	h_name.str(""); h_name << "h_residual_fit_M_" << i_module << "_" <<UV;
		        h_title.str(""); h_title << "Residuals to recon line for Module " << i_module << " " << UV ;
		        auto hl4 = new TH1F(h_name.str().c_str(),h_title.str().c_str(),  149, -0.08, 0.08);
		        hl4->GetXaxis()->SetTitle("[cm]"); hl4->SetDirectory(cd_UV);

		        h_name.str(""); h_name << "h_line_jitter_M_" << i_module << "_" << UV;
		        h_title.str(""); h_title << "Line Jitter for Module " << i_module << " " << UV ;
		        auto hl5 = new TH1F(h_name.str().c_str(),h_title.str().c_str(),  149, -0.03, 0.03);
		        hl5->GetXaxis()->SetTitle("[cm]"); hl5->SetDirectory(cd_UV);

		    } // layers
		} // views
	} // modules

       
  	//Modules
    for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++) {
        h_name.str(""); h_name << "hs_hits_recon_Module" << i_module;
        h_title.str(""); h_title << " ";
        auto hm1 = new TH1F(h_name.str().c_str(),h_title.str().c_str(), 49, -2.5, 2.5);
        hm1->GetXaxis()->SetTitle("[cm]"); hm1->SetDirectory(cd_Modules);

        h_name.str(""); h_name << "h_reconMinusTrue_line_Module_" << i_module;
    	h_title.str(""); h_title << "Recon vs True Track in Module " << i_module;
    	auto hm2 = new TH1F(h_name.str().c_str(),h_title.str().c_str(),  49, -0.01, 0.01);
    	hm2->GetXaxis()->SetTitle("[cm]"); hm2->SetDirectory(cd_Modules);

    	h_name.str(""); h_name << "h_DCA_Module_" << i_module;
    	h_title.str(""); h_title << "DCA in Module " << i_module;
    	auto hm3 = new TH1F(h_name.str().c_str(),h_title.str().c_str(),  49, -0.1, 0.4);
    	hm3->GetXaxis()->SetTitle("[cm]"); hm3->SetDirectory(cd_Modules);

    	h_name.str(""); h_name << "h_Residuals_Module_" << i_module;
    	h_title.str(""); h_title << "Residuals in Module " << i_module;
    	auto hm4 = new TH1F(h_name.str().c_str(),h_title.str().c_str(),  199, -0.2, 0.2);
    	hm4->GetXaxis()->SetTitle("[cm]"); hm4->SetDirectory(cd_Modules);
    }
       
    // Modules and Straws ["combing 4 layers into 1"]
    for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++) {
        for (int i_straw = 0 ; i_straw < Tracker::instance()->getStrawN(); i_straw++) {
            h_name.str(""); h_name << "h" << i_module << "_straw" << i_straw;
            h_title.str(""); h_title << "DCA Straw " << i_straw << " M" << i_module;
            auto hms1 = new TH1F(h_name.str().c_str(),h_title.str().c_str(), 49, -0.1, 0.4);
            hms1->GetXaxis()->SetTitle("[cm]"); hms1->SetDirectory(cd_Straws);
        }
    }
    
    //Use array of pointer of type TH1x to set axis titles and directories 
    TH1F* cmTitle[] = {h_reconMinusTrue_track_intercept, h_sigma, h_hits_MP2, h_dca, h_track_true, h_track_recon,
    	h_intercept, h_x0, h_x1, h_residual_track, h_hits_true, h_hits_recon, h_residual_fit, h_reconMinusTrue_hits, h_reconMinusTrue_track, h_frac_Dintercept, h_meanXRecon};
    for (int i=0; i<(int) sizeof( cmTitle ) / sizeof( cmTitle[0] ); i++){
        TH1F* temp = cmTitle[i];
        cmTitle[i]->SetXTitle("[cm]");
    }
    TH1F* cdAllHits_F[] = {h_sigma, h_hits_MP2, h_dca, h_track_true, h_track_recon, h_hits_true, h_hits_recon, h_residual_track, h_chi2_track, h_residual_fit, 
    	h_chi2_fit, h_reconMinusTrue_hits, h_reconMinusTrue_track}; 
    TH1F* cdTracks_F[] = {h_frac_Dintercept, h_frac_Dslope, h_intercept, h_slope, h_x0, h_x1, h_reconMinusTrue_track_slope, h_reconMinusTrue_track_intercept, h_meanXRecon}; 
    TH1I* cdAllHits_I[] = {h_labels, h_hitCount, h_id_dca};
    for (int i=0; i<(int) sizeof( cdAllHits_F ) / sizeof( cdAllHits_F[0] ); i++){
        cdAllHits_F[i]->SetDirectory(cd_All_Hits);
    }
    for (int i=0; i<(int) sizeof( cdTracks_F ) / sizeof( cdTracks_F[0] ); i++){
        cdTracks_F[i]->SetDirectory(cd_Tracks);
    }
    for (int i=0; i<(int) sizeof( cdAllHits_I ) / sizeof( cdAllHits_I[0] ); i++){
        cdAllHits_I[i]->SetDirectory(cd_All_Hits);
    }

//------------------------------------------Mille Routines---------------------------------------------------------// 

    // Creating .bin, steering, and constrain files
    Mille m (outFileName, asBinary, writeZero);  // call to Mille.cc to create a .bin file
    cout << "Generating test data for g-2 Tracker Alignment in PEDE:" << endl;
       
    //Passing constants to plotting script
    if (plotBool || debugBool){
        contsants_plot << Tracker::instance()->getModuleN() << " " << Tracker::instance()->getViewN() << " " 
                      << Tracker::instance()->getLayerN() << " " << Tracker::instance()->getStrawN() << " " << Tracker::instance()->getTrackNumber() << " " 
                      << Tracker::instance()->getBeamOffset()   << " " << Tracker::instance()->getBeamStart() << " " <<  Tracker::instance()->getBeamPositionLength()  
                      << "  " << Tracker::instance()->getBeamStop() <<  endl;
    }

    // SETTING GEOMETRY
    Tracker::instance()->setGeometry(debug_geom, debugBool);
    cout<< "Geometry is set!" << endl << endl;
    
    // XXX: definition of broken lines here in the future
   
    // MISALIGNMENT
    Tracker::instance()->misalign(debug_mis, debugBool); 
    cout<< "Misalignment is complete!" << endl << endl;
    float sigma_calc = sqrt(pow(Tracker::instance()->getResolution(),2)-pow(Tracker::instance()->getResolution(),2)/float(Tracker::instance()->getLayerTotalN()));
    float Chi2_calc=float(Tracker::instance()->getLayerTotalN());
    float layers_per_modules = float(Tracker::instance()->getViewN())*Tracker::instance()->getLayerN();
    cout << "Manual Misalignment: " << endl;
    for (int i_module=0; i_module<Tracker::instance()->getModuleN(); i_module++){
        cout << showpos << "Module " << i_module <<" :: Characteristic:  " << Tracker::instance()->getSdevX(i_module) << " cm. ";
        if (plotBool || debugBool){pede_mis << Tracker::instance()->getSdevX(i_module) << " "; }
        float indivMis = Tracker::instance()->getSdevX(i_module)-Tracker::instance()->getOverallMis();
        cout << showpos << "Relative: " << indivMis << " cm." << endl;
        //cout << " Chi2_calc" << Chi2_calc << endl;
        Chi2_calc += layers_per_modules * (pow(indivMis,2) - 2*sigma_calc*indivMis)/(pow(sigma_calc,2));
    }
    
    cout << noshowpos; 
    cout << "The overall misalignment was " << Tracker::instance()->getOverallMis() << " cm" <<  endl;
    cout << "The estimated mean residual fit Chi2 = " << Chi2_calc <<  endl;
    cout << "The estimated residual resolution is = " << sigma_calc << " cm" <<  endl;
    cout << "The estimated track 'jitter' is = " << Tracker::instance()->getResolution()/sqrt(Tracker::instance()->getLayerTotalN())  << " cm" <<  endl;
   
    // Write a constraint file, for use with pede
    Tracker::instance()->write_constraint_file(constraint_file, debug_con, debugBool);
    cout<< "Constraints are written! [see Tracker_con.txt]" << endl;

    //Now writing the steering file
    if(steering_file.is_open()){
        steering_file <<  "* g-2 Tracker Alignment: PEDE Steering File" << endl
        << " "  << endl
        << "Tracker_con.txt   ! constraints text file " << endl
        << "Cfiles ! following bin files are Cfiles" << endl 
        << "Tracker_data.bin   ! binary data file" << endl
        << "method inversion 1 0.01" << endl
        << "printrecord  -1 -1      ! debug printout for bad data records" << endl
        << "printrecord 1 -1 ! produces mpdebug.txt"<< endl     //  
        << " "  << endl
        << "end ! optional for end-of-data"<< endl;
    } // end of str file 
    cout<< "Steering file was generated! [see Tracker_con.txt]" << endl;

    cout<< "Calculating residuals..." << endl;
    
 //------------------------------------------Main Mille Track Loop---------------------------------------------------------//    
    bool StrongDebugBool = false;
    //Generating tracks 
    for (int trackCount=0; trackCount<Tracker::instance()->getTrackNumber(); trackCount++){  
        //float p=pow(10.0, 1+Tracker::instance()->generate_uniform());
        //scatterError=sqrt(Tracker::instance()->getWidth())*0.014/p;
        scatterError = 0; // set no scatterError for now 
        
        if (debugBool && StrongDebugBool){
            cout << "Track: " << trackCount << endl;
        }

        //Generating tracks 
        MCData generated_MC = Tracker::instance()->MC_launch(scatterError, debug_calc, debug_off, debug_mc, plot_fit, plot_gen, plot_hits_gen, 
        	plot_hits_fit, debugBool);
            
        for (int hitCount=0; hitCount<generated_MC.hit_count; hitCount++){  //counting only hits going though detector
            //calculating the layer and pixel from the hit number 
            
            //Number of local and global parameters 
            const int nalc = 2;
            const int nagl = 1;  

            int label_mp2 = generated_MC.i_hits[hitCount]; //label to associate hits within different layers with a correct module

            //Local derivatives
            float dlc1=Tracker::instance()->getProjectionX(label_mp2);
            float dlc2=generated_MC.z_hits[hitCount]*Tracker::instance()->getProjectionX(label_mp2);
            float derlc[nalc] = {dlc1, dlc2};
            //Global derivatives
            float dgl1 = Tracker::instance()->getProjectionX(label_mp2);
            float dergl[nagl] = {dgl1};  
            //Labels
            int l1 = label_mp2; 
            int label[nagl] = {l1}; 
            
             
            //TODO multiple scattering errors (no correlations) (for imodel == 1)
            //add break points multiple scattering later XXX (for imodel == 2)
            //! add 'broken lines' offsets for multiple scattering XXX (for imodel == 3)
           
            float rMeas_mp2 =  generated_MC.x_residuals[hitCount]; 
            float sigma_mp2 = generated_MC.hit_sigmas[hitCount]; 

            m.mille(nalc, derlc, nagl, dergl, label, rMeas_mp2, sigma_mp2);
            
            //***********************************Sanity Plots******************************************//
            //Fill for all hits
            h_hits_MP2 -> Fill (rMeas_mp2); // residuals 
            h_sigma -> Fill(sigma_mp2); // errors 
            h_dca->Fill(generated_MC.x_mis_dca[hitCount]); // DCA
            h_hits_true->Fill(generated_MC.x_hit_true[hitCount]); // True (smeared) hit position
            h_hits_recon->Fill(generated_MC.x_hit_recon[hitCount]); // Reconstructed hit position
            h_labels->Fill(l1);
            h_reconMinusTrue_hits->Fill(generated_MC.x_hit_true[hitCount]-generated_MC.x_hit_recon[hitCount]);
 			h_id_dca ->Fill(generated_MC.strawID[hitCount]);
            
            //Track-based hit parameters
            h_track_true->Fill(generated_MC.x_track_true[hitCount]); // True (generated) track position
            h_track_recon->Fill(generated_MC.x_track_recon[hitCount]); // Reconstructed (fitted) track position
            h_reconMinusTrue_track->Fill(generated_MC.x_track_true[hitCount]-generated_MC.x_track_recon[hitCount]);
                                              
            //Calculating Chi2 stats:
            float residual_gen = generated_MC.residuals_gen[hitCount]; 
            h_residual_track->Fill(residual_gen);
            residuals_track_sum_2+=pow(residual_gen/sigma_mp2,2);
            h_residual_fit->Fill(rMeas_mp2); //already used as input to mille
            sigma_fit_calc = sqrt(pow(sigma_mp2,2)-(pow(sigma_mp2,2)/pow(generated_MC.hit_count,1)));
            sigma_fit_calc = 0.0140;        // TODO get it from calculation 
            residuals_fit_sum_2+=pow(rMeas_mp2/sigma_fit_calc,2);
            
            //Fill for hits in modules/layers/straws
            string UV = Tracker::instance()->getUVmapping(generated_MC.View_i[hitCount], generated_MC.Layer_i[hitCount]); // converting view/layer ID into conventional labels
            
            h_name.str(""); h_name << "UV/h_dca_M_" << generated_MC.Module_i[hitCount] << "_" << UV;
            TH1F* h1 = (TH1F*)file->Get( h_name.str().c_str() );
            h1->Fill(generated_MC.x_mis_dca[hitCount]);
            
            h_name.str(""); h_name << "UV/h_strawID_M_" << generated_MC.Module_i[hitCount] << "_" << UV;
            TH1I* h2 = (TH1I*)file->Get( h_name.str().c_str() );
            h2 ->Fill(generated_MC.strawID[hitCount]);
            
            h_name.str(""); h_name << "UV/h_LR_M_" << generated_MC.Module_i[hitCount] << "_" << UV;
            TH1F* h3 = (TH1F*)file->Get( h_name.str().c_str() );
            h3 ->Fill(generated_MC.LR[hitCount]);

            h_name.str(""); h_name << "Straws/h" << generated_MC.Module_i[hitCount] << "_straw" << generated_MC.Straw_i[hitCount];
            TH1F* h4 = (TH1F*)file->Get( h_name.str().c_str() );
            h4->Fill(generated_MC.x_mis_dca[hitCount]);
            h4-> SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);

            h_name.str(""); h_name << "Modules/hs_hits_recon_Module" << generated_MC.Module_i[hitCount];
            TH1F* h5 = (TH1F*)file->Get( h_name.str().c_str() );
            h5->Fill(generated_MC.x_hit_recon[hitCount]);
            h5 -> SetFillColor(colourVector[generated_MC.Module_i[hitCount]]);

            h_name.str(""); h_name << "Modules/h_reconMinusTrue_line_Module_" << generated_MC.Module_i[hitCount];
            TH1F* h6 = (TH1F*)file->Get( h_name.str().c_str() );
            h6->Fill(generated_MC.x_track_true[hitCount]-generated_MC.x_track_recon[hitCount]);

	        h_name.str(""); h_name << "Modules/h_DCA_Module_" << generated_MC.Module_i[hitCount];
	        TH1F* h7 = (TH1F*)file->Get( h_name.str().c_str() );
	        h7->Fill(generated_MC.x_mis_dca[hitCount]);

	        h_name.str(""); h_name << "Modules/h_Residuals_Module_" << generated_MC.Module_i[hitCount];
	        TH1F* h8 = (TH1F*)file->Get( h_name.str().c_str() );
	        h8->Fill(rMeas_mp2);

	        h_name.str(""); h_name << "UV/h_residual_fit_M_" << generated_MC.Module_i[hitCount] << "_" << UV;
	        TH1F* h9 = (TH1F*)file->Get( h_name.str().c_str() );
	        h9->Fill(rMeas_mp2);

	        h_name.str(""); h_name << "UV/h_line_jitter_M_" << generated_MC.Module_i[hitCount] << "_" << UV;
	        TH1F* h10 = (TH1F*)file->Get( h_name.str().c_str() );
	        h10->Fill(generated_MC.x_track_true[hitCount]-generated_MC.x_track_recon[hitCount]);
 
        	if (debugBool){ debug_mp2  << nalc << " " << derlc[0] << " " << derlc[1] << " " << nagl << " " << dergl[0] << " "  << label[0]  << " " 
        	<< rMeas_mp2 << "  " << sigma_mp2 << endl;}
        	hitsN++; //count hits
        } // end of hits loop
         
         //***********************************Sanity Plots: Once per Track******************************************//
        
        
        bool trackCkeck = true;
        float reconHitCheck = 0.0;  

        for (int i=0; i<generated_MC.hit_count; i++){
            reconHitCheck = generated_MC.x_hit_recon[i];
            if (abs(reconHitCheck)>=1.0){
                trackCkeck = false;
            }
        }

        if (trackCkeck){
             for (int i=0; i<generated_MC.hit_count; i++){
                h_res_x_z->Fill(generated_MC.z_hits[i], generated_MC.x_residuals[i]);
            }
        }


        //For generated tracks
        float chi2_track=residuals_track_sum_2;
        h_chi2_track->Fill(chi2_track);
        //h_chi2_ndf_track->Fill(chi2_track/generated_MC.hit_count);  //ndf=#points - constraint [1 = linear] hit count goes from 0 to N points -1 
        //For fitted tracks
        float chi2_fit=residuals_fit_sum_2;
        h_chi2_fit->Fill(chi2_fit);
        //h_chi2_ndf_fit->Fill(chi2_fit/generated_MC.hit_count);  //ndf=#points - constraint [1 = linear] 
        //Resetting counters for next track
        residuals_track_sum_2=0;
        residuals_fit_sum_2=0;
        h_hitCount->Fill(generated_MC.hit_count);
        h_meanXRecon->Fill(generated_MC.meanXReconTrack);


        //Filling Track-based plots 
        h_slope->Fill(generated_MC.slope_truth);
        h_intercept->Fill(generated_MC.intercept_truth);
        h_x0->Fill(generated_MC.x0);
        h_x1->Fill(generated_MC.x1);
        h_reconMinusTrue_track_intercept->Fill(generated_MC.intercept_truth-generated_MC.intercept_recon);
        h_reconMinusTrue_track_slope->Fill(generated_MC.slope_truth-generated_MC.slope_recon);
        float frac_c = (generated_MC.intercept_truth-generated_MC.intercept_recon)/generated_MC.intercept_truth;
        h_frac_Dintercept->Fill(frac_c);
        float frac_m = (generated_MC.slope_truth-generated_MC.slope_recon)/generated_MC.slope_truth;
        h_frac_Dslope->Fill(frac_m);
        // XXX additional measurements from MS IF (imodel == 2) THEN
        //IF (imodel >= 3) THEN

        m.end(); // Write buffer (set of derivatives with same local parameters) to file.
        recordN++; // count records (i.e. tracks);
       
    } // end of track count
    
    cout << " " << endl;
    cout << "****** Calculated Sigma_fit= " << sigma_fit_calc << " *****" << endl; 
    cout << Tracker::instance()->getTrackNumber() << " tracks generated with " << hitsN << " hits." << endl;
    cout << recordN << " records written." << endl;
    float rejectsFrac=Tracker::instance()->getRejectedHitsDCA();
    rejectsFrac = rejectsFrac/(Tracker::instance()->getLayerTotalN()*recordN);
    cout << fixed << setprecision(1);
    cout << "Hits that missed a straw (DCA rejection for all layers): " << Tracker::instance()->getRejectedHitsDCA() << ". [" << rejectsFrac*100 << "%]" << endl;
    cout << "Multiple hits (for all layers): " << Tracker::instance()->getMultipleHitsLayer() << "." << endl;
    cout << "Ambiguity Hits: " << Tracker::instance()->getAmbiguityHit() << "." << endl;
    cout << " " << endl;
    cout << "Ready for PEDE algorithm: ./pede Tracker_str.txt" << endl; 
    cout << "Sanity Plots: root Tracker.root" << endl;
    Logger::Instance()->setUseColor(false); // will be re-enabled below
    // Millepede courtesy of John 
    std::stringstream msg2, msg3, msg4;
    msg2 << Logger::green() << "    _____________________________  \\  /" << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg2.str());
    msg3 << Logger::yellow() << "   {_|_|_|_|_|_|_|_|_|_|_|_|_|_|_( ͡° ͜ʖ ͡°) " << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg3.str());
    msg4 << Logger::red() << "    /\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/" << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg4.str());
    Logger::Instance()->setUseColor(true); // back to default colours 
    if (debugBool) {
    cout << "Normal random numbers were used " << RandomBuffer::instance()->getNormTotal() << " times" <<endl;
    cout << "Gaussian random numbers were used " << RandomBuffer::instance()->getGausTotal() << " times" <<endl;
    Logger::Instance()->write(Logger::WARNING, "Text debug files were produced: ls Tracker_d_*.txt");
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
    
    //------------------------------------------ROOT: Fitting Functions---------------------------------------------------------// 

    cout << endl;
    cout << "-------------------------------------------------------------------------"<< endl; 
    cout << "ROOT fitting parameters and output:" << endl; 
    
    //h_reconMinusTrue_track->Fit("gaus", "Q");
    TF1* chi2pdf = new TF1("chi2pdf","[2]*ROOT::Math::chisquared_pdf(x,[0],[1])",0,40);
    chi2pdf->SetParameters(15, 0., h_chi2_track->Integral("WIDTH")); 
    h_chi2_track->Fit("chi2pdf", "Q"); //Use Pearson chi-square method, using expected errors instead of the observed one given by TH1::GetBinError (default case). 
    //The expected error is instead estimated from the the square-root of the bin function value.
    //h_chi2_fit->Fit("chi2pdf");
	TF1* gausFit = new TF1("gausFit","[2]*ROOT::Math::gaussian_pdf(x,[0],[1])", -0.06, 0.06);
    gausFit->SetParameters(0.01405, 0.0, h_residual_track->Integral("WIDTH"));     
    h_residual_track->Fit("gaus", "Q");
    //h_residual_fit->Fit("gausFit"); 

    float biasMean = h_reconMinusTrue_track->GetMean();
    float biasMeanError = h_reconMinusTrue_track->GetMeanError();

    bias << biasMean << " " << biasMeanError << endl;

if (strongPlotting){
   
    //Residuals per module 
    TCanvas *csg = new TCanvas("csg","csg",700,900);
    TText T; T.SetTextFont(42); T.SetTextAlign(21);
    for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++){
    	h_name.str(""); h_name << "Modules/h_Residuals_Module_" << i_module;
	    TH1F* hd1 = (TH1F*)file->Get( h_name.str().c_str() );
	    hd1->SetFillColor(colourVector[i_module]);
	    hd1->Draw("same");
	    gStyle->SetOptStat("");
	    hd1->SetTitle("");
    }
    T.DrawTextNDC(.5,.95,"Residuals in all Modules");
    csg->Print("residuals_func.png");
     
	//Stacked DCA per straw in each module
	THStack* hs_DCA_Module_0 = new THStack("hs_DCA_Module_0", "");
    THStack* hs_DCA_Module_1 = new THStack("hs_DCA_Module_1", "");
    THStack* hs_DCA_Module_2 = new THStack("hs_DCA_Module_2", "");
    THStack* hs_DCA_Module_3 = new THStack("hs_DCA_Module_3", "");
    THStack* hs_DCA_Module_4 = new THStack("hs_DCA_Module_0", "");
    THStack* hs_DCA_Module_5 = new THStack("hs_DCA_Module_1", "");
    THStack* hs_DCA_Module_6 = new THStack("hs_DCA_Module_2", "");
    THStack* hs_DCA_Module_7 = new THStack("hs_DCA_Module_3", "");
    THStack* stackModule[] = {hs_DCA_Module_0, hs_DCA_Module_1, hs_DCA_Module_2, hs_DCA_Module_3, hs_DCA_Module_4, hs_DCA_Module_5, hs_DCA_Module_6, hs_DCA_Module_7};
    TCanvas *cs = new TCanvas("cs","cs",700,900);
    T.SetTextFont(42); T.SetTextAlign(21);
    cs->Divide((Tracker::instance()->getModuleN()+1)/2,(Tracker::instance()->getModuleN()+1)/2);
    for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++) {
    	for (int i_straw = 0 ; i_straw < Tracker::instance()->getStrawN(); i_straw++) {
            h_name.str(""); h_name << "Straws/h" << i_module << "_straw" << i_straw;
	        TH1F* hs3 = (TH1F*)file->Get( h_name.str().c_str() );
    		stackModule[i_module]->Add(hs3);
        }
    }    
    for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++) {
    	h_title.str(""); h_title << "Module " << i_module << " DCA per straw";
    	cs->cd(i_module+1); stackModule[i_module]->Draw(); T.DrawTextNDC(.5,.95,h_title.str().c_str());
    }
    cs->Print("stack_dca.png");

    //Stacked reconstructed hits
    TCanvas *csr = new TCanvas("csr","csr",700,900);
    T.SetTextFont(42); T.SetTextAlign(21);
    for (int i_module = 0 ; i_module < Tracker::instance()->getModuleN(); i_module++) {
            h_name.str(""); h_name << "Modules/hs_hits_recon_Module" << i_module;
            TH1F* hs4 = (TH1F*)file->Get( h_name.str().c_str() );
            hs_hits_recon->Add(hs4);
    }    
    csr->Divide(1,1);
    csr->cd(1);  hs_hits_recon->Draw(); T.DrawTextNDC(.5,.95,"Recon Hits per Module"); hs_hits_recon->GetXaxis()->SetTitle("[cm]");
    csr->Print("stack_recon_hits.png");

} // strong plotting 

    file->Write();
    file->Close(); //good habit!
    cout << "-------------------------------------------------------------------------"<< endl; 
    cout << endl;
    
    cout << fixed << setprecision(4);
    t_cpu = clock() - t_cpu;
    auto t_end = std::chrono::high_resolution_clock::now();
    cout << "Programme execution took " <<  t_cpu << " CPU clicks (" << ((float)t_cpu)/CLOCKS_PER_SEC << " s)." << " Wall clock time passed: " 
    	<< std::chrono::duration<double>(t_end-t_start).count() << " s." << endl;
    time_t now = time(0);
    char* dt = ctime(&now);
    cout << "Peak RAM use: " << Tracker::instance()->getPeakRSS( )/1e9 << " GB. The C++ compiler used: " << true_cxx << " " << true_cxx_ver
    	<<" Job finished on: " << dt << endl;
    //if(plotBool){theApp.Run();} //ctrl+c to exit //for Canvas in ROOT
    return 0; 
} //end of main
