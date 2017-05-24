    /*
*  TODO: 
    -TTree,Ntuples, etc for ROOT plotting - separate macro: use more advances ROOT methods properly 
    -Graphically plot what's going on pyROOT/Python
    -Use git tag for new version: small increments 0.x, large changes x.0, additional features
        (e.g MF) use a, b, c i.e. v1.3a. Keep version history wiki on git (read-up on that...)
    -End of programme print out more used const and non-const par. 
*
*   Gleb Lukicov (g.lukicov@ucl.ac.uk) @ Fermilab
*   Created: 17 April 2017  
*   Modified: 3 May 2017 
----------------------------------------------------------------
This programme uses MC methods to produce a .bin data file for the 
PEDE routine, to align the tracking detector for the g-2 
experiment.
Methods and functions are contained in AlignTracker_methods.cpp (Tracker class)

============================= Test Model v0.1 =======================================
*Simple 2D case (1D alignment along x) of simplified tracker station of 2 small trackers.
* No B-field, straight tracks, no MS, 100% efficiency. 
*(!) NATURAL UNITS (same as PEDE): cm, rad, GeV [natural units will be omitted] 

(!) Module, Straw, View, etc. labelling starts from 0 [+1 is put by hand as MP2 doesn't seem to accept 0 as a label]

Terminology: for now straw=central wire, detector=module (==the alignment object)

2 trackers with 16 straws in length (instead of 32), 4 layers per tracker: 64 straws per tracker.
    total of 128 straws in 8 layers. 
Gaussian (gaus) Smearing: Mean=0, RMS=1  (on hits, resolution, etc.)
Beam originates (in z direction) from z=0, straight tracks, no scattering, only resolution (gaus) smearing,
    tracks are lines between z=0 and z=25, with x1 and x2 = 10 * rand[0,1] (=0-10)
Resolution (i.e. tracker systematic) is   100um = 0.01 cm 
Misalignment manually by 0.15 for each module in +/- x-direction [all layers in module are dependent (i.e. move together)]

Hit rejection:  1) outside of layer
                2) more than 0.25 cm away from straw centre (i.e straw radius) 

Labels Module 0, Module 1. 
Constrains TODO ??? 
Steering options TODO  ??? do some reading
                                    Tracker Geometry:
First layer fist straw (A1) is at z=50, x=0; spacing (in x) between straws in layer is 0.606 [strawSpacing]
Spacing between layers in same view is 0.515 [layerSpacing]
Spacing between U and V views (last and first layers in views) is 2.020 [viewSpaing]
Spacing between modules (last and first layers in modules) 13.735 [moduleSpacing]
Layers in a view are (e.g. U0 vs. U1) have an extra relative x displacement of 0.303  [layerDisplacement]

For more information see: 
http://gm2-docdb.fnal.gov:8080/cgi-bin/RetrieveFile?docid=4375&filename=Orientation.pdf&version=1 

=======================================================================================

#TODO Model type is implemented via imodel = 0, 1 .....  

#Verbosity level implemented via command line command: n = normal, d = debug [e.g. './AlignTracker d']

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
    string compareStr; //for debug vs. normal output as specified by the user
    int tracksInput; // number of tracks to generate as specified by the user 
    bool debugBool = false; // './AlignTracker n' - for normal, of ./AlignTracker d' - for verbose debug output
    bool plotBool = false; // './AlignTracker p' - for plotting small stats
    bool debugBoolStrong = true; //XXX hack
    //float rand_num; //to store random (int number from buffer + max )/(2.0 * max) [0, 1]   
    //float twoR = 2.0;  //For normalisation of uniform random numbers [0,1]
    //Set up counters for hits and records (tracks)
    int hitsN = 0; // actually recorded (i.e. non-rejected hits)
    int recordN=0; //records = tracks
    float scatterError; // multiple scattering error [calculated here and passed back to the Tracker class]
    
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
    msg01 << Logger::yellow() << "  g-2 Tracker Alignment (v0.1) - Gleb Lukicov (UCL) - May 2017            " << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg01.str());
    msg1 << Logger::blue() <<  "*****************************************************************" << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg1.str());
    Logger::Instance()->setUseColor(true); // back to default colours 

    cout << "Simple Alignment Model with " << Tracker::instance()->getModuleN() << " tracker modules, having " << Tracker::instance()->getStrawN() << " straws per layer " << endl;
    cout << "["<< Tracker::instance()->getLayerN() << " layers per module; " << Tracker::instance()->getViewN() << " views per module]." << endl;
    cout << "No B-field, Straight (parallel) Tracks, 100\% efficiency" << endl;
    cout << "Hit rejection: DCA > StrawRadius [" << Tracker::instance()->getStrawRadius() << " cm]" << endl;
    cout << "Parallel Tracks: single hit per layer allowed [shortest DCA is chosen as the hit]" << endl;


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
    //Debug files [only filled with "d" option]
    string mp2_debugFileName = "Tracker_d_mille.txt"; // Inputs into binary files
    string cacl_debugFileName = "Tracker_d_calc.txt"; // intermediate calculation 
    string mis_debugFileName = "Tracker_d_mis.txt";  // Misalignment 
    string geom_debugFileName = "Tracker_d_geom.txt"; // Geometry 
    string off_debugFileName = "Tracker_d_off.txt";  // Offsets/Missed hits
    string MC_debugFileName = "Tracker_d_MC.txt";  // Final results from MC
    string con_debugFileName = "Tracker_d_con.txt"; //Constraints
    string gen_plotFileName = "Tracker_p_gen.txt"; //
    string fit_plotFileName = "Tracker_p_fit.txt"; //
  
        
    // TODO TTree -> separate Macro for plotting [see Mark's suggested code: check correct implementation for future] 
    //output ROOT file
    TFile* file = new TFile("Tracker.root", "recreate");  // recreate = overwrite if already exists
     // Book histograms
    TH1F* h_sigma = new TH1F("h_sigma", "Sigma [cm]",  100,  0.100, 0.015); // F=float bins, name, title, nBins, Min, Max
    TH1F* h_hits_MP2 = new TH1F("h_hits_MP2", "MP2 Hits: Residuals from fitted line to ideal geometry (with DCAs from misaligned geom.) [cm]",  400, 0.0, 0.2);
    //TH1F* h_slope = new TH1F("h_slope", "Slope ",  500,  -300, 300);
    //TH1F* h_c = new TH1F("h_c", "Intercept ",  500,  -300, 300);
    TH1F* h_det = new TH1F("h_det", "DCA (to misaligned detector from generated track)",  500,  -0.05, 0.4);
    TH1F* h_true = new TH1F("h_true", "True hits (the x of the track in-line with a layer)",  500,  -1, 3);
    TH1F* h_ideal = new TH1F("h_x_ideal", "Ideal X position of the hits: dca (from mis.) + ideal position of a straw",  500,  -1, 3);
    TH1F* h_fit = new TH1F("h_fit", "Reconstructed x of the fitted line (to ideal geometry)",  500,  -1, 3);

    //std::unique_ptr<RootManager> rm_;  // gm2util/common/RootManager.hh" TODO  
    // std::string name;
    // std::stringstream sh; std::stringstream st;
    // for (int i_layer=0; i_layer < Tracker::instance()->getLayerN(); i_layer++){
    //     sh.str(""); sh << "h_det_layer_" << i_layer;
    //     auto hx = new TH1F(sh.str().c_str(),st.str().c_str(), 500,  -0.05, 0.4);
    //     hx->GetXaxis()->SetTitle("DCA [cm]");
    //     //rm_->Add(currentDirName.str(), hx);
    // }

    TH1F* h_det_layer_0 = new TH1F("h_det_layer_0", "DCA_layer_0",  500,  -0.05, 0.4);
    TH1F* h_det_layer_1 = new TH1F("h_det_layer_1", "DCA_layer_1",  500,  -0.05, 0.4);
    TH1F* h_det_layer_2 = new TH1F("h_det_layer_2", "DCA_layer_2",  500,  -0.05, 0.4);
    TH1F* h_det_layer_3 = new TH1F("h_det_layer_3", "DCA_layer_3",  500,  -0.05, 0.4);

    TH1I* h_strawID_layer_0 = new TH1I("h_strawID_layer_0", "h_strawID_layer_0", 5, -1, 4);
    TH1I* h_strawID_layer_1 = new TH1I("h_strawID_layer_1", "h_strawID_layer_1", 5, -1, 4);
    TH1I* h_strawID_layer_2 = new TH1I("h_strawID_layer_2", "h_strawID_layer_2", 5, -1, 4);
    TH1I* h_strawID_layer_3 = new TH1I("h_strawID_layer_3", "h_strawID_layer_3", 5, -1, 4);

    TH1F* h_LR_layer_0 = new TH1F("h_LR_layer_0", "h_LR_layer_0", 4, -2, 2);
    TH1F* h_LR_layer_1 = new TH1F("h_LR_layer_1", "h_LR_layer_1", 4, -2, 2);
    TH1F* h_LR_layer_2 = new TH1F("h_LR_layer_2", "h_LR_layer_2", 4, -2, 2);
    TH1F* h_LR_layer_3 = new TH1F("h_LR_layer_3", "h_LR_layer_3", 4, -2, 2);
    
    // Creating .bin, steering, and constrain files
    Mille m (outFileName, asBinary, writeZero);  // call to Mille.cc to create a .bin file
    cout << "Generating test data for g-2 Tracker Alignment in PEDE:" << endl;
   
    // file streams for constrain, steering, and debug files 
    ofstream constraint_file(conFileName);
    ofstream steering_file(strFileName);
    ofstream debug_mp2(mp2_debugFileName);
    ofstream debug_calc(cacl_debugFileName);
    ofstream debug_mis(mis_debugFileName);
    ofstream debug_geom(geom_debugFileName);
    ofstream debug_off(off_debugFileName);
    ofstream debug_mc(MC_debugFileName); 
    ofstream debug_con(con_debugFileName);
    ofstream plot_gen(gen_plotFileName);
    ofstream plot_fit(fit_plotFileName);
    
    // Setting fixed precision for floating point values
    cout << fixed << setprecision(6); 
    debug_mp2 << fixed << setprecision(6); 
    debug_calc << fixed << setprecision(6); 
    debug_mis << fixed << setprecision(6); 
    debug_geom << fixed << setprecision(6); 
    debug_off << fixed << setprecision(6); 
    debug_mc << fixed << setprecision(6); 
    debug_con << fixed << setprecision(6);
    plot_gen << fixed << setprecision(9);
    plot_fit << fixed << setprecision(9);
  

    // SETTING GEOMETRY
    Tracker::instance()->setGeometry(debug_geom, debugBool);
    cout<< "Geometry is set!" << endl;
    
    // XXX: definition of broken lines here in the future
   
    // MISALIGNMENT
    Tracker::instance()->misalign(debug_mis, debugBool); 
    cout<< "Misalignment is complete!" << endl;
    
    // Write a constraint file, for use with pede
    Tracker::instance()->write_constraint_file(constraint_file, debug_con, debugBool);
    cout<< "Constraints are written!" << endl;

    //Now writing the steering file
    //TODO what steering parameters does tracker require? 
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
        << "*printrecord  -1 -1      ! debug printout for bad data records" << endl
        //<< "" << endl
        << "*outlierdownweighting  2 ! number of internal iterations (> 1)"<< endl
        << "*dwfractioncut      0.2  ! 0 < value < 0.5"<< endl
        << "*presigma           0.01 ! default value for presigma"<< endl
        << "*regularisation 1.0      ! regularisation factor"<< endl
        << "*regularisation 1.0 0.01 ! regularisation factor, pre-sigma"<< endl
        << " " << endl
        << "*bandwidth 0         ! width of precond. band matrix"<< endl
        << "method diagonalization 3 0.001 ! diagonalization      "<< endl
        << "method fullMINRES       3 0.01 ! minimal residual     "<< endl
        << "method sparseMINRES     3 0.01 ! minimal residual     "<< endl
        << "*mrestol      1.0D-8          ! epsilon for MINRES"<< endl
        << "method inversion       3 0.001 ! Gauss matrix inversion"<< endl
        << "* last method is applied"<< endl
        << "*matiter      3  ! recalculate matrix in iterations" << endl
        << " "  << endl
        << "end ! optional for end-of-data"<< endl;
    } // end of str file 
    cout<< "Steering file was generated!" << endl;

    cout<< "Generating data..." << endl;
    //Generating particles with energies: 10-100 [GeV]
    for (int trackCount=0; trackCount<Tracker::instance()->getTrackNumber(); trackCount++){  //Track number is set in methods XXX is it the best place for it? 
        //rand_num = (RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max());
        //float p=pow(10.0, 1+rand_num);
        //scatterError=sqrt(Tracker::instance()->getWidth())*0.014/p;
        scatterError = 0; // set no scatterError for now 
        
        if (debugBool){
            cout << "Track: " << trackCount << endl;
        }

        //Generating tracks 
        MCData generated_MC = Tracker::instance()->MC(scatterError, debug_calc, debug_off, debug_mc, plot_fit, debugBool);

            
        for (int hitCount=0; hitCount<generated_MC.hit_count; hitCount++){  //counting only hits going though detector
            //calculating the layer and pixel from the hit number 
            
            //Number of local and global parameters 
            const int nalc = 2;  
            const int nagl = 1;  

            int label_mp2 = generated_MC.i_hits[hitCount]; //label to associate hits within different layers with a correct module

            //Local derivatives
            // TODO check the maths...
            float dlc1=Tracker::instance()->getProjectionX(label_mp2);
            float dlc2=generated_MC.z_hits[hitCount]*Tracker::instance()->getProjectionX(label_mp2);
            float derlc[nalc] = {dlc1, dlc2};
            //Global derivatives
            float dgl1 = Tracker::instance()->getProjectionX(label_mp2);
            float dergl[nagl] = {dgl1};  
            //Labels 
            /// TODO check that properly
            int l1 = label_mp2+1; //Millepede doesn't like 0 as a label XXX 
            int label[nagl] = {l1}; 
             
            //TODO multiple scattering errors (no correlations) (for imodel == 1)
            //add break points multiple scattering later XXX (for imodel == 2)
            //! add 'broken lines' offsets for multiple scattering XXX (for imodel == 3)
           
            float rMeas_mp2 =  generated_MC.x_hits[hitCount]; 
            float sigma_mp2 = generated_MC.hit_sigmas[hitCount]; 

            //XXX passing by reference/pointer vs as variable (?) - makes any difference
            m.mille(nalc, derlc, nagl, dergl, label, rMeas_mp2, sigma_mp2);
            
            //Sanity Plots

            //Fill for hits
            h_hits_MP2 -> Fill (rMeas_mp2); 
            h_sigma -> Fill(sigma_mp2);
            h_det->Fill(generated_MC.x_mis_dca[hitCount]);
            h_true->Fill(generated_MC.x_track[hitCount]);
            h_ideal->Fill(generated_MC.x_ideal[hitCount]);
            h_fit->Fill(generated_MC.x_fitted[hitCount]);

            if (hitCount ==0){h_det_layer_0->Fill(generated_MC.x_mis_dca[hitCount]); h_strawID_layer_0->Fill(generated_MC.strawID[hitCount]); h_LR_layer_0->Fill(generated_MC.LR[hitCount]);} 
            if (hitCount ==1){h_det_layer_1->Fill(generated_MC.x_mis_dca[hitCount]); h_strawID_layer_1->Fill(generated_MC.strawID[hitCount]); h_LR_layer_1->Fill(generated_MC.LR[hitCount]);} 
            if (hitCount ==2){h_det_layer_2->Fill(generated_MC.x_mis_dca[hitCount]); h_strawID_layer_2->Fill(generated_MC.strawID[hitCount]); h_LR_layer_2->Fill(generated_MC.LR[hitCount]);} 
            if (hitCount ==3){h_det_layer_3->Fill(generated_MC.x_mis_dca[hitCount]); h_strawID_layer_3->Fill(generated_MC.strawID[hitCount]); h_LR_layer_3->Fill(generated_MC.LR[hitCount]);} 

            //Fill for tracks
            if (hitCount ==0){
                
                // h_slope->Fill(generated_MC.x_m[hitCount]);
                // h_c->Fill(generated_MC.x_c[hitCount]);
                plot_gen << generated_MC.x0_gen[hitCount] << " " << generated_MC.z0_gen[hitCount] << " " << generated_MC.x1_gen[hitCount] << " " << generated_MC.z1_gen[hitCount] << endl;
            }
            debug_mp2  << nalc << " " << derlc[0] << " " << derlc[1] << " " << nagl << " " << dergl[0] << " "  << label[0]  << " " << rMeas_mp2 << "  " << sigma_mp2 << endl;
            hitsN++; //count hits
        } // end of hits loop
        
     
        //TODO impliment via ROOT Manager!  
        // for (int i_layer=0; i_layer < Tracker::instance()->getLayerN(); i_layer++){
        //     sh.str(""); sh << "h_det_layer_" << i_layer;
        //     auto h = rm_->Get<TH1F*>(currentDirName.str(),sh.str());
        //     h->Fill(generated_MC.x_mis_dca[i_layer]);
        // }

        // XXX additional measurements from MS IF (imodel == 2) THEN
        //IF (imodel >= 3) THEN

        m.end(); // Write buffer (set of derivatives with same local parameters) to file.
        recordN++; // count records (i.e. tracks);
    
        
    } // end of track count
    
    cout << " " << endl;
    cout << Tracker::instance()->getTrackNumber() << " tracks generated with " << hitsN << " hits." << endl;
    cout << recordN << " records written." << endl;
    cout << "There was " << Tracker::instance()->getRejectedHitsDCA() << " hits that missed a straw in TOTAL." << endl;
    cout << "Additionally, " << Tracker::instance()->getMultipleHitsLayer() << " multiple hits in layers in TOTAL." << endl;
    cout << " " << endl;
    cout << "Ready for PEDE algorithm: ./pede Tracker_str.txt" << endl; 
    cout << "Sanity Plots: root Tracker.root" << endl;
    cout << "Manual Misalignment was " << Tracker::instance()->getSdevX(0) << " cm for Module 0, and " << Tracker::instance()->getSdevX(1) << " cm for Module 1." << endl; 
    cout << "Resolution was " << Tracker::instance()->getResolution() << " cm" << endl;  
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
    h_sigma->SetXTitle( "[cm]");
    h_sigma->SetXTitle( "[cm]");
    h_det->SetXTitle( "DCA [cm]");
    h_true->SetXTitle( "[cm]");
    h_ideal->SetXTitle( "[cm]");
    h_fit->SetXTitle( "[cm]");
    file->Write();
    file->Close(); //good habit!
    
    t_cpu = clock() - t_cpu;
    auto t_end = std::chrono::high_resolution_clock::now();
    cout << "Programme execution took " <<  t_cpu << " CPU clicks (" << ((float)t_cpu)/CLOCKS_PER_SEC << " s)." << " Wall clock time passed: " << std::chrono::duration<double>(t_end-t_start).count() << " s." << endl;
    time_t now = time(0);
    char* dt = ctime(&now);
    cout << "Peak RAM use: " << Tracker::instance()->getPeakRSS( )/1e9 << " GB. The C++ compiler used: " << true_cxx << " " << true_cxx_ver<<" Job finished on: " << dt << endl;
    //if(plotBool){theApp.Run();} //ctrl+c to exit //for Canvas in ROOT
    return 0; 
} //end of main 
