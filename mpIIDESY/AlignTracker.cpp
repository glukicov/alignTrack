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
First layer fist straw (A1) is at z=50, x=0; spacing (in x) between straws in layer is ... XXX [strawSpacing]
Spacing between layers in same view is ... XXX [layerSpacing]
Spacing between U and V views (last and first layers in views) is ... XXX [viewSpaing]
Spacing between modules (last and first layers in modules) ... XXX [moduleSpacing]
Layers in a view are (e.g. U0 vs. U1) have an extra relative x displacement of ... XXX  [layerDisplacement]

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

//***************MAIN****************//
int main(int argc, char* argv[]){

    //////----Variable Initialisation-------///////////
    int imodel = 0;  //Model type (see above) TODO implement this as an argument to main [for broken lines, MF, etc.] 
    string compareStr; //for debug vs. normal output as specified by the user
    bool debugBool = false; // './AlignTracker n' - for normal, of ./AlignTracker d' - for verbose debug output
    //float rand_num; //to store random (int number from buffer + max )/(2.0 * max) [0, 1]   
    //float twoR = 2.0;  //For normalisation of uniform random numbers [0,1]
    //Set up counters for hits and records (tracks)
    int hitsN = 0; // actually recorded (i.e. non-rejected hits)
    int recordN=0; //records = tracks 
    float scatterError; // multiple scattering error [calculated here and passed back to Tracker class]
    
    //Tell the logger to only show message at INFO level or above
    Logger::Instance()->setLogLevel(Logger::NOTE); 
    //Tell the logger to throw exceptions when ERROR messages are received
    Logger::Instance()->enableCriticalErrorThrow();

    // Check if correct number of arguments specified, exiting if not
    if (argc > 2) { Logger::Instance()->write(Logger::ERROR, "Too many arguments -  please specify verbosity flag. (e.g. debug or align)");} 
        else if (argc < 2) {Logger::Instance()->write(Logger::ERROR, "Too few arguments - please specify verbosity flag. (e.g. debug or align)");} 
        else { // Set filenames to read random numbers from, using arguments. Catch exception if these files do not exist.
            try {compareStr = argv[1];} 
        catch (ios_base::failure& e) {
            Logger::Instance()->write(Logger::ERROR, "Filestream exception caught: " + string(e.what()) + "\nPlease ensure valid filenames are specified!");
        }
    } // end of 2nd else [correct # arguments]

    //this is also passed to Tracker functions, with debug file names
    if (compareStr=="d"){debugBool = true; // print out to debug files
    Logger::Instance()->write(Logger::WARNING,  "DEBUG MODE"); }
    else{ debugBool = false; // print out to debug files
    }
    
    Logger::Instance()->setUseColor(false); // will be re-enabled below [to use custom colour output to terminal]
    std::stringstream msg0, msg01, msg02, msg1;
    Logger::Instance()->write(Logger::NOTE, "");
    msg0 << Logger::blue() <<  "*************************************************************" << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg0.str());
    msg01 << Logger::yellow() << "   g-2  Tracker Alignment (v0.1) - Gleb Lukicov (UCL) - May 2017            " << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg01.str());
    msg1 << Logger::blue() <<  "*************************************************************" << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg1.str());
    Logger::Instance()->setUseColor(true); // back to default colours 

    cout << "Simple Alignment Model with " << Tracker::instance()->getModuleN() << " traker modules, having " << Tracker::instance()->getStrawN() << " straws per layer (out of 4) layers." << endl;
    cout << "No B-field, Straight Tracks, 100\% efficiency" << endl; 


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
   
        
    // TODO TTree -> separate Macro for plotting [see Mark's suggested code: check correct implementation for future] 
    //output ROOT file
    TFile* file = new TFile("Tracker.root", "recreate");  // recreate = overwrite if already exists
     // Book histograms
    TH1F* h_sigma = new TH1F("h_sigma", "Sigma [cm]",  100,  0, 0.02); // F=float bins, name, title, nBins, Min, Max
    TH1F* h_hits_MP2 = new TH1F("h_hits_MP2", "MP2 Hits [cm]",  400, -0.1, 0.1);
    TH2F* h_gen = new TH2F("h_gen", "Generated track points", 15, -3, 12, 30, -3, 27);
    TH1F* h_slope = new TH1F("h_slope", "Slope ",  500,  -300, 300);
    TH1F* h_c = new TH1F("h_c", "Intercept ",  500,  -300, 300);
    TH1F* h_det = new TH1F("h_det", "DCA (det hits)",  500,  -0.1, 0.4);
    TH1F* h_true = new TH1F("h_true", "True hits",  500,  -0.2, 0.8);
    TH1F* h_mis = new TH1F("h_mis", "Misal. hits",  500,  -0.2, 0.8);
    
    // Creating .bin, steering, and constrain files
    Mille m (outFileName, asBinary, writeZero);  // call to Mille.cc to create a .bin file
    cout << "Generating test data for g-2 Tracker Alignment in PEDE." << endl;
   
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

    // Setting fixed precision for floating point values
    cout << fixed << setprecision(6); 
    debug_mp2 << fixed << setprecision(6); 
    debug_calc << fixed << setprecision(6); 
    debug_mis << fixed << setprecision(6); 
    debug_geom << fixed << setprecision(6); 
    debug_off << fixed << setprecision(6); 
    debug_mc << fixed << setprecision(6); 
    debug_con << fixed << setprecision(6);

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


    //Generating particles with energies: 10-100 [GeV]
    for (int trackCount=0; trackCount<Tracker::instance()->getTrackCount(); trackCount++){  //Track count is set in methods XXX is it the best place for it? 
        //rand_num = (RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max());
        //float p=pow(10.0, 1+rand_num);
        //scatterError=sqrt(Tracker::instance()->getWidth())*0.014/p;
        scatterError = 0; // set no scatterError for now 
        
        //Generating tracks 
        LineData generated_line = Tracker::instance()->MC(scatterError, debug_calc, debug_off, debug_mc, debugBool);

        //Fill for tracks
        h_gen->Fill(generated_line.x0_gen[trackCount],generated_line.z0_gen[trackCount]);
        h_gen->Fill(generated_line.x1_gen[trackCount],generated_line.z1_gen[trackCount]);
        h_slope->Fill(generated_line.x_m[trackCount]);
        h_c->Fill(generated_line.x_c[trackCount]);
            
        for (int hitCount=0; hitCount<generated_line.hit_count; hitCount++){  //counting only hits going though detector
            //calculating the layer and pixel from the hit number 
            
            //Number of local and global parameters 
            const int nalc = 2;  
            const int nagl = 1;  

            int label_mp2 = generated_line.i_hits[hitCount]; //label to associate hits within different layers with a correct module

            //Local derivatives
            // TODO do the maths...
            float dlc1=Tracker::instance()->getProjectionX(label_mp2);
            float dlc2=generated_line.z_hits[hitCount]*Tracker::instance()->getProjectionX(label_mp2);
            float derlc[nalc] = {dlc1, dlc2};
            //Global derivatives
            float dgl1 = Tracker::instance()->getProjectionX(label_mp2);
            float dergl[nagl] = {dgl1};  
            //Labels 
            /// TODO check that properly
            int l1 = label_mp2+1; //Millepede doesn't like 0 as a label apparently ... XXX 
            int label[nagl] = {l1}; 
             
            //TODO multiple scattering errors (no correlations) (for imodel == 1)
            //add break points multiple scattering later XXX (for imodel == 2)
            //! add 'broken lines' offsets for multiple scattering XXX (for imodel == 3)
           
            float rMeas_mp2 =  generated_line.x_hits[hitCount]; 
            float sigma_mp2 = generated_line.hit_sigmas[hitCount]; 

            //XXX passing by reference/pointer vs as variable (?) - makes any difference
            m.mille(nalc, derlc, nagl, dergl, label, rMeas_mp2, sigma_mp2);
            
            //Sanity Plots

            //Fill for hits
             h_hits_MP2 -> Fill (rMeas_mp2); 
             h_sigma -> Fill(sigma_mp2);
             h_det->Fill(generated_line.x_det[hitCount]);
             h_mis->Fill(generated_line.x_mis[hitCount]);
             h_true->Fill(generated_line.x_true[hitCount]);

                       
            debug_mp2  << nalc << " " << derlc[0] << " " << derlc[1] << " " << nagl << " " << dergl[0] << " "  << label[0]  << " " << rMeas_mp2 << "  " << sigma_mp2 << endl;


            hitsN++; //count hits
        } // end of hits loop
         
        // XXX additional measurements from MS IF (imodel == 2) THEN
        //IF (imodel >= 3) THEN

        m.end(); // Write buffer (set of derivatives with same local parameters) to file.
        recordN++; // count records (i.e. tracks);
       
    } // end of track count
    
    
    cout << " " << endl;
    cout << Tracker::instance()->getTrackCount() << " tracks generated with " << hitsN << " hits." << endl;
    cout << recordN << " records written." << endl;
    cout << " " << endl;
    cout << "Ready for PEDE algorithm: ./pede Tracker_str.txt" << endl; 
    cout << "Sanity Plots: root Tracker.root" << endl;
    cout << "Manual Misalignment was " << Tracker::instance()->getDispX() << " cm: " << Tracker::instance()->getSdevX(0) << " for Module 0 and " << Tracker::instance()->getSdevX(1) << " for Module 1." << endl; 
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
    gStyle->SetOptStat(1111111); // TODO make this work with TFile? or just separate ROOT macro?
    h_gen->SetMarkerStyle(30);
    h_gen->SetMarkerColor(kBlue);
    h_gen->Draw();
    h_sigma ->Draw();
    h_hits_MP2 ->Draw();
    h_slope ->Draw(); 
    h_c ->Draw();
    h_det ->Draw();
    h_true ->Draw();
    h_mis ->Draw();
    
    file->Write();
    file->Close(); //good habit!
    
   
    return 0; 
} //end of main 
