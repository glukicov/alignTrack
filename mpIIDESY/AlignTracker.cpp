/*
*   Gleb Lukicov (g.lukicov@ucl.ac.uk) @ Fermilab
*   Created: 17 April 2017  
*   Modified: 3 May 2017 
----------------------------------------------------------------
This programme uses MC methods to produce a .bin data file for the 
PEDE routine, to align the tracking detector for the g-2 
experiment.
Methods and functions are contained in AlignTracker_methods.cpp (Tracker class)

============================= Test Model v0.2 =======================================
*Simple 2D case (1D alignment along x).
* No B-field, straight tracks, no MS, 100% efficiency. 
*(!) NATURAL UNITS (same as PEDE): cm, rad, GeV [natural units will be omitted] 

(!) Module, Straw, View, Layer labelling starts from 0 [+1 is put by hand as MP2 doesn't accept 0 as a label]

Terminology: detector=module (==the alignment object)

Gaussian (gaus) Smearing: Mean=0, RMS=1  (on hits, resolution, etc.)
Beam originates (in z direction) from z=0, straight tracks, no scattering, only resolution (gaus) smearing,
    tracks are lines between z=beamStart and z=beamStop, with x1 and x2 = Constant * rand[0,1]
Resolution (i.e. tracker systematic) is implemented. 
Misalignment manually by 0.15 for each module in +/- x-direction [all layers in module are dependent (i.e. move together)]

Hit rejection:  1) outside of layer
                2) more than 0.25 cm away from straw centre (i.e straw radius) 

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

    //////----Variable Initialisation-------///////////
    int imodel = 0;  //Model type (see above) TODO implement this as an argument to main [for broken lines, MF, etc.]
    int setPrecision = 6; // precision (# decimal points) of printout for debug text files and cout
    string compareStr; //for debug vs. normal output as specified by the user
    int tracksInput; // number of tracks to generate as specified by the user 
    bool debugBool = false; // './AlignTracker n' - for normal, of ./AlignTracker d' - for verbose debug output
    bool plotBool = false; // './AlignTracker p' - for plotting with PlotGen.py 
    //Set up counters for hits and records (tracks)
    int hitsN = 0; // actually recorded (i.e. non-rejected hits)
    int recordN=0; //records = tracks
    float scatterError; // multiple scattering error [calculated here and passed back to the Tracker class]
    float residuals_track_sum_2=0.0;
    float residuals_fit_sum_2=0.0;
    const Color_t colourVector[]={kWhite, kBlack, kBlue, kGreen, kYellow, kRed, kGray, kMagenta};
       
    //Tell the logger to only show message at INFO level or above
    Logger::Instance()->setLogLevel(Logger::NOTE); 
    //Tell the logger to throw exceptions when ERROR messages are received
    Logger::Instance()->enableCriticalErrorThrow();

    //TApplication theApp("App",&argc, argv); //for Canvas in ROOT

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
        Logger::Instance()->write(Logger::ERROR, "Please specify verbosity flag. (e.g. e.g. debug [d], plot[p] or align/normal [n])");
    }
    
    Logger::Instance()->setUseColor(false); // will be re-enabled below [to use custom colour output to terminal]
    std::stringstream msg0, msg01, msg02, msg1;
    Logger::Instance()->write(Logger::NOTE, "");
    msg0 << Logger::blue() <<  "*****************************************************************" << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg0.str());
    msg01 << Logger::yellow() << "  g-2 Tracker Alignment (v0.2) - Gleb Lukicov (UCL) - June 2017            " << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg01.str());
    msg1 << Logger::blue() <<  "*****************************************************************" << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg1.str());
    Logger::Instance()->setUseColor(true); // back to default colours 

    cout << "Simple Alignment Model with " << Tracker::instance()->getModuleN() << " tracker modules, having " << Tracker::instance()->getStrawN() << " straws per layer." << endl;
    cout << "["<< Tracker::instance()->getLayerN() << " layers per module; " << Tracker::instance()->getViewN() << " views per module]." << endl;
    cout << "Total of " << Tracker::instance()->getLayerTotalN() << " measurement layers." << endl; 
    cout << "No B-field, Straight (parallel) Tracks, 100% efficiency." << endl;
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

    //Debug files [only filled with "d" option]
    string mp2_debugFileName = "Tracker_d_mille.txt"; // Inputs into binary files
    string cacl_debugFileName = "Tracker_d_calc.txt"; // intermediate calculation 
    string mis_debugFileName = "Tracker_d_mis.txt";  // Misalignment 
    string geom_debugFileName = "Tracker_d_geom.txt"; // Geometry 
    string off_debugFileName = "Tracker_d_off.txt";  // Offsets/Missed hits
    string MC_debugFileName = "Tracker_d_MC.txt";  // Final results from MC
    string con_debugFileName = "Tracker_d_con.txt"; //Constraints
    string gen_plotFileName = "Tracker_p_gen.txt"; // Truth Track points
    string fit_plotFileName = "Tracker_p_fit.txt"; // Reconstructed Track points
    string contsants_plotFileName = "Tracker_p_constants.txt"; // passing constants (e.g. strawN to python script)
    string hits_gen_plotFileName = "Tracker_p_hits_gen.txt"; // Truth Hits points

    // file streams for debug files 
    ofstream debug_mp2(mp2_debugFileName);
    ofstream debug_calc(cacl_debugFileName);
    ofstream debug_mis(mis_debugFileName);
    ofstream debug_geom(geom_debugFileName);
    ofstream debug_off(off_debugFileName);
    ofstream debug_mc(MC_debugFileName); 
    ofstream debug_con(con_debugFileName);
    ofstream plot_gen(gen_plotFileName);
    ofstream plot_fit(fit_plotFileName);
    ofstream contsants_plot(contsants_plotFileName);
    ofstream plot_hits_gen(hits_gen_plotFileName);
    
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
    
    
    //output ROOT file
    TFile* file = new TFile("Tracker.root", "recreate");  // recreate = overwrite if already exists
    //create a subdirectories 
    TDirectory* cd_All_Hits = file->mkdir("All_Hits");
    TDirectory* cd_Layers = file->mkdir("Layers");
    TDirectory* cd_Modules = file->mkdir("Modules");
    TDirectory* cd_Views = file->mkdir("Views");
    TDirectory* cd_Tracks = file->mkdir("Tracks");
    TDirectory* cd_Debug = file->mkdir("Debug");
    TDirectory* cd_Straws = file->mkdir("Straws");
     // Book histograms [once only]
    TH1F* h_sigma = new TH1F("h_sigma", "MP2 Input: Detector Resolution (sigma) [cm]",  500,  Tracker::instance()->getResolution()-0.001, Tracker::instance()->getResolution()+0.001); // F=float bins, name, title, nBins, Min, Max
    TH1F* h_hits_MP2 = new TH1F("h_hits_MP2", "MP2 Input: Residuals from fitted line to ideal geometry [cm]",  1000, -0.2, 0.2);
    TH1F* h_det = new TH1F("h_det", "DCA (to misaligned detector from generated track)",  500,  -0.05, Tracker::instance()->getStrawRadius()+0.25);
    TH1F* h_true = new TH1F("h_true", "True track position (the x of the generated track in-line with a layer)",  300,  -1.5, 1.5);
    TH1F* h_hits_true = new TH1F("h_hits_true", "True hit position (the x of the generated and smeared hit)",  300,  -1.5, 1.5);
    TH1F* h_recon = new TH1F("h_recon", "Reconstructed X position of the hits in ideal detector",  500,  -1.5, 1.5);
    TH1F* h_fit = new TH1F("h_fit", "Reconstructed x of the fitted line (to ideal geometry)",  500,  -(Tracker::instance()->getBeamOffset()+3), Tracker::instance()->getBeamPositionLength()+1);
    TH1I* h_labels = new TH1I("h_labels", "Labels in PEDE", 8 , 0, 8);
    TH1F* h_resiudal_track = new TH1F("h_resiudal_track", "Residuals for generated tracks", 500, -0.4, 0.4);
    TH1F* h_chi2_track = new TH1F("h_chi2_track", "Chi2 for generated tracks", 40, -1, 50);
    //TH1F* h_chi2_ndf_track = new TH1F("h_chi2_ndf_track", "Chi2/ndf for generated tracks", 60, -1, 5);
    TH1F* h_resiudal_fit = new TH1F("h_resiudal_fit", "Residuals for fitted tracks", 500, -0.6, 0.6);
    TH1F* h_chi2_fit = new TH1F("h_chi2_fit", "Chi2 for fitted tracks", 1500, -1, 1000);
    //TH1F* h_chi2_ndf_fit = new TH1F("h_chi2_ndf_fit", "Chi2/ndf for fitted tracks", 200, -1, 30);
    TH1I* h_hitCount = new TH1I("h_hitCount", "Total Hit count per track", 32 , 0, 32);
    TH1F* h_reconMinusTrue_track = new TH1F("h_reconMinusTrue_line", "Reconstructed - True X position of the lines",  500,  -0.1, 0.1);
    TH1F* h_reconMinusTrue_hits = new TH1F("h_reconMinusTrue_hits", "Reconstructed - True X position of the hits",  499,  -0.1, 0.1);
    
    TH1I* h_id_dca = new TH1I("h_id_dca", "ID for hit straws", Tracker::instance()->getStrawN(), 0, Tracker::instance()->getStrawN());
    TH1F* h_rightTail = new TH1F("h_rightTail", "Reconstructed X hits > 1",  50, 0.9, 1.1);
    TH1F* h_leftTail = new TH1F("h_leftTail", "Reconstructed X hits < -1",  50,  -0.9,-1.1);
    TH1F* h_rightTail_true = new TH1F("h_rightTail_true", "True X hits > 1",  50, 0.9, 1.1);
    TH1F* h_leftTail_true = new TH1F("h_leftTail_true", "True X hits < -1",  50,  -0.9,-1.1);
    
    //Use array of pointer of type TH1x to set axis titles and directories 
    TH1F* cmTitle[] = {h_sigma, h_hits_MP2, h_det, h_true, h_recon, h_fit, h_resiudal_track, h_hits_true, h_resiudal_fit, h_reconMinusTrue_hits, h_reconMinusTrue_track};
    for (int i=0; i<(int) sizeof( cmTitle ) / sizeof( cmTitle[0] ); i++){
    	TH1F* temp = cmTitle[i];
    	cmTitle[i]->SetXTitle("[cm]");
    }
    TH1F* cdAllHits_F[] = {h_sigma, h_hits_MP2, h_det, h_true, h_recon, h_fit, h_hits_true, h_resiudal_track, h_chi2_track, h_resiudal_fit, h_chi2_fit, h_reconMinusTrue_hits, h_reconMinusTrue_track}; 
    TH1I* cdAllHits_I[] = {h_labels, h_hitCount, h_id_dca};
    TH1F* cdDebug[] = {h_rightTail, h_leftTail, h_rightTail_true, h_leftTail_true};
    for (int i=0; i<(int) sizeof( cdAllHits_F ) / sizeof( cdAllHits_F[0] ); i++){
    	cdAllHits_F[i]->SetDirectory(cd_All_Hits);
    }
    for (int i=0; i<(int) sizeof( cdAllHits_I ) / sizeof( cdAllHits_I[0] ); i++){
    	cdAllHits_I[i]->SetDirectory(cd_All_Hits);
    }
    for (int i=0; i<(int) sizeof( cdDebug ) / sizeof( cdDebug[0] ); i++){
    	cdDebug[i]->SetDirectory(cd_Debug);
    }

    std::stringstream h_name;
    std::stringstream h_title;
    //Booking histograms for TOTAL # layers
    for (int i = 0 ; i < Tracker::instance()->getLayerTotalN(); i++) {
        h_name.str(""); h_name << "h_det_layer_" << i;
        h_title.str(""); h_title << "DCA in layer " << i;
        auto h1 = new TH1F(h_name.str().c_str(),h_title.str().c_str(), 100,  -0.05, Tracker::instance()->getStrawRadius()+0.25);
        h1->GetXaxis()->SetTitle("[cm]");
        
        
        h_name.str(""); h_name << "h_strawID_layer_" << i;
        h_title.str(""); h_title << "strawID in layer " << i;
        auto h2 = new TH1I(h_name.str().c_str(),h_title.str().c_str(), Tracker::instance()->getStrawN(), 0, Tracker::instance()->getStrawN());
        h2->GetXaxis()->SetTitle("Straw ID [0-31]");
        

        h_name.str(""); h_name << "h_LR_layer_" << i;
        h_title.str(""); h_title << "Left-Right hit distribution in layer" << i;
        auto h3 = new TH1F(h_name.str().c_str(),h_title.str().c_str(),  4, -2, 2);
        h3->GetXaxis()->SetTitle("L= - 1.0; R = +1.0 ");

        h1->SetDirectory(cd_Layers);
        h2->SetDirectory(cd_Layers);
        h3->SetDirectory(cd_Layers);
    }

   
    
    TH1F* hres_0 = new TH1F("h_Residuals_module_0", "", 200, -0.2, 0.2);
    TH1F* hres_1 = new TH1F("h_Residuals_module_1","Residuals in Module 1", 200, -0.2, 0.2);
    TH1F* hres_2 = new TH1F("h_Residuals_module_2","Residuals in Module 2", 200, -0.2, 0.2);
    TH1F* hres_3 = new TH1F("h_Residuals_module_3","Residuals in Module 3", 200, -0.2, 0.2);
    TH1F* hres_4 = new TH1F("h_Residuals_module_4","Residuals in Module 4", 200, -0.2, 0.2);
    TH1F* hres_5 = new TH1F("h_Residuals_module_5","Residuals in Module 5", 200, -0.2, 0.2);
    
    //Booking stack plots for canvas 
    THStack* hs_DCA_Module_0 = new THStack("hs_DCA_Module_0", "");
    THStack* hs_DCA_Module_1 = new THStack("hs_DCA_Module_1", "");
    THStack* hs_DCA_Module_2 = new THStack("hs_DCA_Module_2", "");
    THStack* hs_DCA_Module_3 = new THStack("hs_DCA_Module_3", "");

    TH1F* h0_straw1 = new TH1F("h0_straw1", "DCA Straw 1-M0", 99, -0.1, 0.4);
    TH1F* h0_straw2 = new TH1F("h0_straw2", "DCA Straw 2-M0", 99, -0.1, 0.4);
    TH1F* h0_straw3 = new TH1F("h0_straw3", "DCA Straw 3-M0", 99, -0.1, 0.4);
    TH1F* h0_straw4 = new TH1F("h0_straw4", "DCA Straw 4-M0", 99, -0.1, 0.4);
    TH1F* h0_straw5 = new TH1F("h0_straw5", "DCA Straw 5-M0", 99, -0.1, 0.4);

    TH1F* h1_straw1 = new TH1F("h2_straw1", "DCA Straw 1-M1", 99, -0.1, 0.4);
    TH1F* h1_straw2 = new TH1F("h2_straw2", "DCA Straw 2-M1", 99, -0.1, 0.4);
    TH1F* h1_straw3 = new TH1F("h2_straw3", "DCA Straw 3-M1", 99, -0.1, 0.4);
    TH1F* h1_straw4 = new TH1F("h2_straw4", "DCA Straw 4-M1", 99, -0.1, 0.4);
    TH1F* h1_straw5 = new TH1F("h2_straw5", "DCA Straw 5-M1", 99, -0.1, 0.4);

    TH1F* h2_straw1 = new TH1F("h3_straw1", "DCA Straw 1-M2", 99, -0.1, 0.4);
    TH1F* h2_straw2 = new TH1F("h3_straw2", "DCA Straw 2-M2", 99, -0.1, 0.4);
    TH1F* h2_straw3 = new TH1F("h3_straw3", "DCA Straw 3-M2", 99, -0.1, 0.4);
    TH1F* h2_straw4 = new TH1F("h3_straw4", "DCA Straw 4-M2", 99, -0.1, 0.4);
    TH1F* h2_straw5 = new TH1F("h1_straw5", "DCA Straw 5-M2", 99, -0.1, 0.4);

    TH1F* h3_straw1 = new TH1F("h4_straw1", "DCA Straw 1-M3", 99, -0.1, 0.4);
    TH1F* h3_straw2 = new TH1F("h4_straw2", "DCA Straw 2-M3", 99, -0.1, 0.4);
    TH1F* h3_straw3 = new TH1F("h4_straw3", "DCA Straw 3-M3", 99, -0.1, 0.4);
    TH1F* h3_straw4 = new TH1F("h4_straw4", "DCA Straw 4-M3", 99, -0.1, 0.4);
    TH1F* h3_straw5 = new TH1F("h4_straw5", "DCA Straw 5-M3", 99, -0.1, 0.4);


    TH1F* h_DCA_Module_0 = new TH1F ("h_DCA_Module_0", "DCA in Module 0", 99, -0.1, 0.4);
    TH1F* h_DCA_Module_1 = new TH1F ("h_DCA_Module_1", "DCA in Module 1", 99, -0.1, 0.4);
    TH1F* h_DCA_Module_2 = new TH1F ("h_DCA_Module_2", "DCA in Module 2", 99, -0.1, 0.4);
    TH1F* h_DCA_Module_3 = new TH1F ("h_DCA_Module_3", "DCA in Module 3", 99, -0.1, 0.4);
    h_DCA_Module_0 -> SetDirectory(cd_Modules);
    h_DCA_Module_1 -> SetDirectory(cd_Modules);
    h_DCA_Module_2 -> SetDirectory(cd_Modules);
    h_DCA_Module_3 -> SetDirectory(cd_Modules);

    for (int i_module=0; i_module< Tracker::instance()->getModuleN(); i_module++){
         for (int i_view=0; i_view< Tracker::instance()->getViewN(); i_view++){
            for (int i_layer=0; i_layer< Tracker::instance()->getLayerN(); i_layer++){ 
                for (int i_straw=0; i_straw< Tracker::instance()->getStrawN(); i_straw++){
                        h_name.str(""); h_name << "h_DCA_Module_" << i_module <<"_View_" << i_view << "_Layer_"<<i_layer<<"_Straw_"<<i_straw;
                        h_title.str(""); h_title<< "DCA in Module_" << i_module <<"_View_" << i_view << "_Layer_"<<i_layer<<"_Straw_"<<i_straw;
                        auto hs = new TH1F(h_name.str().c_str(),h_title.str().c_str(), 50, -0.1, 0.4);
                        hs->GetXaxis()->SetTitle("[cm]");
                        hs->SetDirectory(cd_Straws);
                }
            }
        }
    }


    // Creating .bin, steering, and constrain files
    Mille m (outFileName, asBinary, writeZero);  // call to Mille.cc to create a .bin file
    cout << "Generating test data for g-2 Tracker Alignment in PEDE:" << endl;
       
    //Passing constants to plotting script
    if (plotBool){
        contsants_plot << Tracker::instance()->getModuleN() << " " << Tracker::instance()->getViewN() << " " 
                      << Tracker::instance()->getLayerN() << " " << Tracker::instance()->getStrawN() << " " << Tracker::instance()->getTrackNumber() << " " 
                      << Tracker::instance()->getBeamOffset()   << " " << Tracker::instance()->getBeamStart() << " " <<  Tracker::instance()->getBeamPositionLength()  << "  " << Tracker::instance()->getBeamStop() <<  endl;
    }

    // SETTING GEOMETRY
    Tracker::instance()->setGeometry(debug_geom, debugBool);
    cout<< "Geometry is set!" << endl;
    
    // XXX: definition of broken lines here in the future
   
    // MISALIGNMENT
    
    Tracker::instance()->misalign(debug_mis, debugBool); 
    cout<< "Misalignment is complete!" << endl;
    float sigma_calc = sqrt(pow(Tracker::instance()->getResolution(),2)-pow(Tracker::instance()->getResolution(),2)/float(Tracker::instance()->getLayerTotalN()));
    float Chi2_calc=float(Tracker::instance()->getLayerTotalN());
    float layers_per_modules = float(Tracker::instance()->getViewN())*Tracker::instance()->getLayerN();
    cout << "Manual Misalignment: " << endl;
    for (int i_module=0; i_module<Tracker::instance()->getModuleN(); i_module++){
        cout << "Module " << i_module <<" :: Relative:  " << Tracker::instance()->getSdevX(i_module) << " cm. ";
        float indivMis = Tracker::instance()->getSdevX(i_module)-Tracker::instance()->getOverallMis();
        cout << "Individual: " << indivMis << " cm." << endl;
        //cout << " Chi2_calc" << Chi2_calc << endl;
        Chi2_calc += layers_per_modules * (pow(indivMis,2) - 2*sigma_calc*indivMis)/(pow(sigma_calc,2));
    }
    
    cout << "The overall misalignment was " << Tracker::instance()->getOverallMis() << " cm" <<  endl;
    cout << "The catapulted Chi2 = " << Chi2_calc <<  endl;
   
    // Write a constraint file, for use with pede
    Tracker::instance()->write_constraint_file(constraint_file, debug_con, debugBool);
    cout<< "Constraints are written! [see Tracker_con.txt]" << endl;

    //Now writing the steering file
    //TODO what steering parameters does the tracker require? 
    if(steering_file.is_open()){
        
        steering_file <<  "*            Default test steering file" << endl
        << "Cfiles ! following bin files are Cfiles" << endl 
        << "Tracker_con.txt   ! constraints text file " << endl
        << "Tracker_data.bin   ! binary data file" << endl
         << "fortranfiles ! following bin files are fortran" << endl
        //<< "*outlierrejection 100.0 ! reject if Chi^2/Ndf >" << endl
        //<< "*outliersuppression 3   ! 3 local_fit iterations" << endl
        << "*hugecut 50.0     !cut factor in iteration 0" << endl
        << "*chisqcut 1.0 1.0 ! cut factor in iterations 1 and 2" << endl
        << "*entries  10 ! lower limit on number of entries/parameter" << endl
        //<< "" << endl 
        << "*pairentries 10 ! lower limit on number of parameter pairs"  << endl
        << "                ! (not yet!)"            << endl
        << "*printrecord   1  2      ! debug printout for records" << endl
       // << "" << endl
        << "printrecord  -1 -1      ! debug printout for bad data records" << endl
        //<< "" << endl
        << "*outlierdownweighting  2 ! number of internal iterations (> 1)"<< endl
        << "*dwfractioncut      0.2  ! 0 < value < 0.5"<< endl
        << "*presigma           0.01 ! default value for presigma"<< endl
        << "*regularisation 1.0      ! regularisation factor"<< endl
        << "*regularisation 1.0 0.01 ! regularisation factor, pre-sigma"<< endl
        << " " << endl
        << "*bandwidth 0         ! width of precond. band matrix"<< endl
        << "*method HIP 3 0.001 ! diagonalization      "<< endl   // XXX this method ignores constraints
        << "*method bandcholesky 3 0.001 ! diagonalization      "<< endl 
        << "*method cholesky 3 0.001 ! diagonalization      "<< endl 
        << "*method diagonalization 3 0.001 ! diagonalization      "<< endl
        << "*method fullMINRES       3 0.01 ! minimal residual     "<< endl
        << "*method sparseMINRES     3 0.01 ! minimal residual     "<< endl
        << "*mrestol      1.0D-8          ! epsilon for MINRES"<< endl
        << "*method inversion       3 0.001 ! Gauss matrix inversion"<< endl
        << "* last method is applied"<< endl
        << "*matiter      3  ! recalculate matrix in iterations" << endl
        << "printrecord 1 -1" << endl     // XXX produces mpdebug.txt
        << " "  << endl
        << "end ! optional for end-of-data"<< endl;
    } // end of str file 
    cout<< "Steering file was generated! [see Tracker_con.txt]" << endl;

    cout<< "Calculating residuals..." << endl;
    

    bool StrongDebugBool = false;
    //Generating tracks 
    for (int trackCount=0; trackCount<Tracker::instance()->getTrackNumber(); trackCount++){  //Track number is set in methods XXX is it the best place for it? 
        //float p=pow(10.0, 1+Tracker::instance()->generate_uniform());
        //scatterError=sqrt(Tracker::instance()->getWidth())*0.014/p;
        scatterError = 0; // set no scatterError for now 
        
        if (debugBool && StrongDebugBool){
            cout << "Track: " << trackCount << endl;
        }

        //Generating tracks 
        MCData generated_MC = Tracker::instance()->MC_launch(scatterError, debug_calc, debug_off, debug_mc, plot_fit, plot_gen, plot_hits_gen, debugBool);

            
        for (int hitCount=0; hitCount<generated_MC.hit_count; hitCount++){  //counting only hits going though detector
            //calculating the layer and pixel from the hit number 
            
            //Number of local and global parameters 
            //const int nalc = 2;
            const int nalc = 1; 
            const int nagl = 1;  

            int label_mp2 = generated_MC.i_hits[hitCount]; //label to associate hits within different layers with a correct module

            //Local derivatives
            float dlc1=Tracker::instance()->getProjectionX(label_mp2);
            //float dlc2=generated_MC.z_hits[hitCount]*Tracker::instance()->getProjectionX(label_mp2);
            //float derlc[nalc] = {dlc1, dlc2};
            float derlc[nalc] = {dlc1};
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
            
            //Sanity Plots
            //Fill for all hits
            h_hits_MP2 -> Fill (rMeas_mp2); // residuals 
            h_sigma -> Fill(sigma_mp2); // errors 
            h_det->Fill(generated_MC.x_mis_dca[hitCount]); // DCA
            h_hits_true->Fill(generated_MC.x_hit_true[hitCount]); // True (smeared) hit position
            h_true->Fill(generated_MC.x_track_true[hitCount]); // True (generated) track position
            h_recon->Fill(generated_MC.x_hit_recon[hitCount]); // Reconstructed hit position
            h_fit->Fill(generated_MC.x_track_recon[hitCount]); // Reconstructed (fitted) track position
            h_labels->Fill(l1);
            h_reconMinusTrue_track->Fill(generated_MC.x_track_true[hitCount]-generated_MC.x_track_recon[hitCount]);
            h_reconMinusTrue_hits->Fill(generated_MC.x_hit_true[hitCount]-generated_MC.x_hit_recon[hitCount]);
 			h_id_dca ->Fill(generated_MC.strawID[hitCount]);

            //Debug Plots
            if (generated_MC.x_hit_recon[hitCount] > 1.0){
                h_rightTail -> Fill(generated_MC.x_hit_recon[hitCount]);
                h_rightTail -> Fill(2-generated_MC.x_hit_recon[hitCount]);
            }

            if ( generated_MC.x_hit_recon[hitCount] < -1.0){
                h_leftTail -> Fill(generated_MC.x_hit_recon[hitCount]);
                h_leftTail -> Fill(-2+abs(generated_MC.x_hit_recon[hitCount]));
            }

            if (generated_MC.x_hit_true[hitCount] > 1.0){
                h_rightTail_true -> Fill(generated_MC.x_hit_true[hitCount]);
                h_rightTail_true -> Fill(2-generated_MC.x_hit_true[hitCount]);
            }

            if ( generated_MC.x_hit_true[hitCount] < -1.0){
                h_leftTail_true -> Fill(generated_MC.x_hit_true[hitCount]);
                h_leftTail_true -> Fill(-2+abs(generated_MC.x_hit_true[hitCount]));
            }

                        
            //Calculating Chi2 stats:
            float residual_gen = generated_MC.residuals_gen[hitCount]; 
            h_resiudal_track->Fill(residual_gen);
            residuals_track_sum_2+=pow(residual_gen/sigma_mp2,2);
            h_resiudal_fit->Fill(rMeas_mp2); //already used as input to mille
            float sigma_fit = sqrt(pow(sigma_mp2,2)-(pow(sigma_mp2,2)/pow(generated_MC.hit_count,1)));
            //cout << "sigma_fit= " << sigma_fit ; 
            residuals_fit_sum_2+=pow(rMeas_mp2/sigma_fit,2);
            
            //Fill for hits in layers
            h_name.str(""); h_name << "Layers/h_det_layer_" << hitCount;
            TH1F* h1 = (TH1F*)file->Get( h_name.str().c_str() );
            h1->Fill(generated_MC.x_mis_dca[hitCount]);
            h_name.str(""); h_name << "Layers/h_strawID_layer_" << hitCount;
            TH1I* h2 = (TH1I*)file->Get( h_name.str().c_str() );
            h2 ->Fill(generated_MC.strawID[hitCount]);
            h_name.str(""); h_name << "Layers/h_LR_layer_" << hitCount;
            TH1F* h3 = (TH1F*)file->Get( h_name.str().c_str() );
            h3 ->Fill(generated_MC.LR[hitCount]);


            // h_name.str(""); h_name << "Modules/h_Residuals_module_" << generated_MC.Module_i[hitCount];
            // TH1F* h4 = (TH1F*)file->Get(h_name.str().c_str() );
            // h4->Fill(rMeas_mp2);
           
            h_name.str(""); h_name << "Straws/h_DCA_Module_" << generated_MC.Module_i[hitCount] <<"_View_" <<  generated_MC.View_i[hitCount] << "_Layer_"<< generated_MC.Layer_i[hitCount]<<"_Straw_"<< generated_MC.Straw_i[hitCount];
            TH1F* h5 = (TH1F*)file->Get(h_name.str().c_str() );       
            h5->Fill(generated_MC.x_mis_dca[hitCount]);
            h5->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);
            
       
            if (generated_MC.Module_i[hitCount] == 0){ hres_0 -> Fill(rMeas_mp2);}
            if (generated_MC.Module_i[hitCount] == 1){ hres_1 -> Fill(rMeas_mp2);}
            if (generated_MC.Module_i[hitCount] == 2){ hres_2 -> Fill(rMeas_mp2);}
            if (generated_MC.Module_i[hitCount] == 3){ hres_3 -> Fill(rMeas_mp2);}
            if (generated_MC.Module_i[hitCount] == 4){ hres_4 -> Fill(rMeas_mp2);}
            if (generated_MC.Module_i[hitCount] == 5){ hres_5 -> Fill(rMeas_mp2);}

        #if 0
            if (generated_MC.Layer_i[hitCount] ==1 && generated_MC.View_i[hitCount] ==0 ){
            //if (1==1){
            if (generated_MC.Module_i[hitCount] == 0){
                //hs_DCA_Module_0->Add(h5);
                h_DCA_Module_0->Fill(generated_MC.x_mis_dca[hitCount]);

                if (generated_MC.Straw_i[hitCount]== 1){ h0_straw1 -> Fill(generated_MC.x_mis_dca[hitCount]); h0_straw1->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
                if (generated_MC.Straw_i[hitCount]== 2){ h0_straw2 -> Fill(generated_MC.x_mis_dca[hitCount]); h0_straw2->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
                if (generated_MC.Straw_i[hitCount]== 3){ h0_straw3 -> Fill(generated_MC.x_mis_dca[hitCount]); h0_straw3->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
                if (generated_MC.Straw_i[hitCount]== 4){ h0_straw4 -> Fill(generated_MC.x_mis_dca[hitCount]); h0_straw4->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
                if (generated_MC.Straw_i[hitCount]== 5){ h0_straw5 -> Fill(generated_MC.x_mis_dca[hitCount]); h0_straw5->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            
            }
            if (generated_MC.Module_i[hitCount] == 1){
            //hs_DCA_Module_1->Add(h5);
            h_DCA_Module_1->Fill(generated_MC.x_mis_dca[hitCount]);
            
            if (generated_MC.Straw_i[hitCount]== 1){ h1_straw1 -> Fill(generated_MC.x_mis_dca[hitCount]); h1_straw1->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            if (generated_MC.Straw_i[hitCount]== 2){ h1_straw2 -> Fill(generated_MC.x_mis_dca[hitCount]); h1_straw2->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            if (generated_MC.Straw_i[hitCount]== 3){ h1_straw3 -> Fill(generated_MC.x_mis_dca[hitCount]); h1_straw3->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            if (generated_MC.Straw_i[hitCount]== 4){ h1_straw4 -> Fill(generated_MC.x_mis_dca[hitCount]); h1_straw4->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            if (generated_MC.Straw_i[hitCount]== 5){ h1_straw5 -> Fill(generated_MC.x_mis_dca[hitCount]); h1_straw5->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            

            }
            
            if (generated_MC.Module_i[hitCount] == 2){
            //hs_DCA_Module_2->Add(h5);
            h_DCA_Module_2->Fill(generated_MC.x_mis_dca[hitCount]);
            if (generated_MC.Straw_i[hitCount]== 1){ h2_straw1 -> Fill(generated_MC.x_mis_dca[hitCount]); h2_straw1->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            if (generated_MC.Straw_i[hitCount]== 2){ h2_straw2 -> Fill(generated_MC.x_mis_dca[hitCount]); h2_straw2->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            if (generated_MC.Straw_i[hitCount]== 3){ h2_straw3 -> Fill(generated_MC.x_mis_dca[hitCount]); h2_straw3->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            if (generated_MC.Straw_i[hitCount]== 4){ h2_straw4 -> Fill(generated_MC.x_mis_dca[hitCount]); h2_straw4->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            if (generated_MC.Straw_i[hitCount]== 5){ h2_straw5 -> Fill(generated_MC.x_mis_dca[hitCount]); h2_straw5->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            
            }
            
            if (generated_MC.Module_i[hitCount] == 3){
            //hs_DCA_Module_3->Add(h5);
            h_DCA_Module_3->Fill(generated_MC.x_mis_dca[hitCount]);
            if (generated_MC.Straw_i[hitCount]== 1){ h3_straw1 -> Fill(generated_MC.x_mis_dca[hitCount]); h3_straw1->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            if (generated_MC.Straw_i[hitCount]== 2){ h3_straw2 -> Fill(generated_MC.x_mis_dca[hitCount]); h3_straw2->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            if (generated_MC.Straw_i[hitCount]== 3){ h3_straw3 -> Fill(generated_MC.x_mis_dca[hitCount]); h3_straw3->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            if (generated_MC.Straw_i[hitCount]== 4){ h3_straw4 -> Fill(generated_MC.x_mis_dca[hitCount]); h3_straw4->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
            if (generated_MC.Straw_i[hitCount]== 5){ h3_straw5 -> Fill(generated_MC.x_mis_dca[hitCount]); h3_straw5->SetFillColor(colourVector[generated_MC.Straw_i[hitCount]]);}
        
            }
        }
        #endif

            if (debugBool){ debug_mp2  << nalc << " " << derlc[0] <<  " " << nagl << " " << dergl[0] << " "  << label[0]  << " " << rMeas_mp2 << "  " << sigma_mp2 << endl;}
            hitsN++; //count hits
        } // end of hits loop
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
        
        // XXX additional measurements from MS IF (imodel == 2) THEN
        //IF (imodel >= 3) THEN

        m.end(); // Write buffer (set of derivatives with same local parameters) to file.
        recordN++; // count records (i.e. tracks);
    
       
    } // end of track count

    
    cout << " " << endl;
    cout << Tracker::instance()->getTrackNumber() << " tracks generated with " << hitsN << " hits." << endl;
    cout << recordN << " records written." << endl;
    float rejectsFrac=Tracker::instance()->getRejectedHitsDCA();
    rejectsFrac = rejectsFrac/(Tracker::instance()->getLayerTotalN()*recordN);
    cout << fixed << setprecision(1);
    cout << "Hits that missed a straw (DCA rejection for all layers): " << Tracker::instance()->getRejectedHitsDCA() << ". [" << rejectsFrac*100 << "%]" << endl;
    cout << "Multiple hits (for all layers): " << Tracker::instance()->getMultipleHitsLayer() << "." << endl;
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
    
    //ROOT stuff
    cout << endl;
    cout << "-------------------------------------------------------------------------"<< endl; 
    cout << "ROOT fitting parameters:" << endl; 
    
    h_reconMinusTrue_track->Fit("gaus");
    TF1* chi2pdf = new TF1("chi2pdf","[2]*ROOT::Math::chisquared_pdf(x,[0],[1])",0,40);
    chi2pdf->SetParameters(15, 0., h_chi2_track->Integral("WIDTH")); 
    h_chi2_track->Fit("chi2pdf"); //Use Pearson chi-square method, using expected errors instead of the observed one given by TH1::GetBinError (default case). The expected error is instead estimated from the the square-root of the bin function value.
    //h_chi2_fit->Fit("chi2pdf");
    h_resiudal_track->Fit("gaus");
    //h_resiudal_fit->Fit("gaus"); 


    //hres_0->Draw(); 
    TCanvas *csg = new TCanvas("cs","cs",700,900);
    TText Tg; Tg.SetTextFont(42); Tg.SetTextAlign(21);
    hres_0->Add(hres_5);
    hres_0->Draw();
    
    Tg.DrawTextNDC(.5,.95,"Residuals in all Modules");
    hres_1->Draw("same");
    hres_2->Draw("same");
    hres_3->Draw("same");
    hres_4->Draw("same");
    TF1* gaussian0 = new TF1("gaussian","[0]*TMath::Gaus(x,[1],[2])", -0.08-Tracker::instance()->getOverallMis() , 0.08-Tracker::instance()->getOverallMis());
    gaussian0->SetParameters(4480, 0-Tracker::instance()->getOverallMis(), 0.01452);
    TF1* gaussian1 = new TF1("gaussian1","[0]*TMath::Gaus(x,[1],[2])", -0.16 -Tracker::instance()->getOverallMis() , 0.08 -Tracker::instance()->getOverallMis());
    gaussian1->SetParameters(2230, -0.03-Tracker::instance()->getOverallMis(), 0.01452);
    TF1* gaussian2 = new TF1("gaussian3","[0]*TMath::Gaus(x,[1],[2])", -0.08-Tracker::instance()->getOverallMis(), 0.16-Tracker::instance()->getOverallMis());
    gaussian2->SetParameters(2230, 0.1-Tracker::instance()->getOverallMis(), 0.01452);
    TF1* gaussian3 = new TF1("gaussian3","[0]*TMath::Gaus(x,[1],[2])", -0.12 -Tracker::instance()->getOverallMis() , 0.04 -Tracker::instance()->getOverallMis());
     gaussian3->SetParameters(2230, -0.045-Tracker::instance()->getOverallMis(), 0.01452);
    TF1* gaussian4 = new TF1("gaussian4","[0]*TMath::Gaus(x,[1],[2])", -0.04-Tracker::instance()->getOverallMis(), 0.12-Tracker::instance()->getOverallMis());
    gaussian4->SetParameters(2230, 0.05-Tracker::instance()->getOverallMis(), 0.01452);


    gaussian0->SetLineColor(kRed);
    gaussian1->SetLineColor(kBlack);
    gaussian2->SetLineColor(kGreen);
    gaussian3->SetLineColor(kOrange);
    gaussian4->SetLineColor(kBlue);
    gaussian0->Draw("same");
    gaussian1->Draw("same");
    gaussian2->Draw("same");
    gaussian3->Draw("same");
    gaussian4->Draw("same");
    gStyle->SetOptStat("oue");
    csg->Print("residuals_func.png");
    //gROOT->ForceStyle();

   #if 0
    TCanvas *cs = new TCanvas("cs","cs",700,900);
    TText T; T.SetTextFont(42); T.SetTextAlign(21);
    hs_DCA_Module_0->Add(h0_straw1);
    hs_DCA_Module_0->Add(h0_straw2);
    hs_DCA_Module_0->Add(h0_straw3);
    hs_DCA_Module_0->Add(h0_straw4);
    hs_DCA_Module_0->Add(h0_straw5);

    hs_DCA_Module_1->Add(h1_straw1);
    hs_DCA_Module_1->Add(h1_straw2);
    hs_DCA_Module_1->Add(h1_straw3);
    hs_DCA_Module_1->Add(h1_straw4);
    hs_DCA_Module_1->Add(h1_straw5);

    hs_DCA_Module_2->Add(h2_straw1);
    hs_DCA_Module_2->Add(h2_straw2);
    hs_DCA_Module_2->Add(h2_straw3);
    hs_DCA_Module_2->Add(h2_straw4);
    hs_DCA_Module_2->Add(h2_straw5);

    hs_DCA_Module_3->Add(h3_straw1);
    hs_DCA_Module_3->Add(h3_straw2);
    hs_DCA_Module_3->Add(h3_straw3);
    hs_DCA_Module_3->Add(h3_straw4);
    hs_DCA_Module_3->Add(h3_straw5);
    cs->Divide(2,2);
    cs->cd(1); hs_DCA_Module_0->Draw(); T.DrawTextNDC(.5,.95,"Module 0 DCA per straw");
    cs->cd(2); hs_DCA_Module_1->Draw(); T.DrawTextNDC(.5,.95,"Module 1 DCA per straw");
    cs->cd(3); hs_DCA_Module_2->Draw(); T.DrawTextNDC(.5,.95,"Module 2 DCA per straw");
    cs->cd(4); hs_DCA_Module_3->Draw(); T.DrawTextNDC(.5,.95,"Module 3 DCA per straw");
    
    cs->Print("stack_4.png");
    #endif

    // TCanvas *csh = new TCanvas("cs","cs",700,900);
    // TText Th; Th.SetTextFont(42); Th.SetTextAlign(21);
    // cs->Divide(2,2);
    // cs->cd(1); hs_DCA_Module_0->Draw("nostack"); Th.DrawTextNDC(.5,.95,"Module 0 DCA per straw");
    // cs->cd(2); hs_DCA_Module_1->Draw("nostack"); Th.DrawTextNDC(.5,.95,"Module 1 DCA per straw");
    // cs->cd(3); hs_DCA_Module_2->Draw("nostack"); Th.DrawTextNDC(.5,.95,"Module 2 DCA per straw");
    // cs->cd(4); hs_DCA_Module_3->Draw("nostack"); Th.DrawTextNDC(.5,.95,"Module 3 DCA per straw");
    // csh->Print("stack.png");
    // TCanvas *css = new TCanvas("css","css",700,900);
    // TText Ts; Ts.SetTextFont(42); Ts.SetTextAlign(21);
    // css->Divide(2,2);
    // css->cd(1); hs_DCA_Module_0->Draw(""); Ts.DrawTextNDC(.5,.95,"Module 0 DCA per straw: stacked");
    // css->cd(2); hs_DCA_Module_1->Draw(""); Ts.DrawTextNDC(.5,.95,"Module 1 DCA per straw: stacked");
    // css->cd(3); hs_DCA_Module_2->Draw(""); Ts.DrawTextNDC(.5,.95,"Module 2 DCA per straw: stacked");
    // css->cd(4); hs_DCA_Module_3->Draw(""); Ts.DrawTextNDC(.5,.95,"Module 3 DCA per straw: stacked");
    // css->Print("stack_s.png");


    file->Write();
    file->Close(); //good habit!
    cout << "-------------------------------------------------------------------------"<< endl; 
    cout << endl;
    
    cout << fixed << setprecision(4);
    t_cpu = clock() - t_cpu;
    auto t_end = std::chrono::high_resolution_clock::now();
    cout << "Programme execution took " <<  t_cpu << " CPU clicks (" << ((float)t_cpu)/CLOCKS_PER_SEC << " s)." << " Wall clock time passed: " << std::chrono::duration<double>(t_end-t_start).count() << " s." << endl;
    time_t now = time(0);
    char* dt = ctime(&now);
    cout << "Peak RAM use: " << Tracker::instance()->getPeakRSS( )/1e9 << " GB. The C++ compiler used: " << true_cxx << " " << true_cxx_ver<<" Job finished on: " << dt << endl;
    //if(plotBool){theApp.Run();} //ctrl+c to exit //for Canvas in ROOT
    return 0; 
} //end of main
