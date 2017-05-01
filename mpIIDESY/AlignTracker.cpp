/*
*  TODOs: 
* TTree,Ntuples, etc for ROOT plotting - separate macro? 
*
*   Gleb Lukicov (g.lukicov@ucl.ac.uk) @ Fermilab
*   Created: 17 April 2017
*   Modified: 1 March 2017 
----------------------------------------------------------------
This programme uses MC methods to produce a .bin data file for the 
PEDE routine, to align the tracking detector for the g-2 
experiment.
Methods and functions are contained in AlignTracker_methods.cpp (Tracker class)

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
    float rand_num; //to store random (int number from buffer + max )/(2.0 * max) [0, 1]   
    float twoR = 2.0;  //For normalisation of uniform random numbers [0,1]
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
    msg01 << Logger::yellow() << "   g-2  Tracker Alignment - Gleb Lukicov (UCL) - April 2017            " << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg01.str());
    msg1 << Logger::blue() <<  "*************************************************************" << Logger::def();
    Logger::Instance()->write(Logger::NOTE,msg1.str());
    Logger::Instance()->setUseColor(true); // back to default colours 

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
    TH1F* h_sigma = new TH1F("h_sigma", "Sigma [cm]",  100,  0, 0.004); // F=float bins, name, title, nBins, Min, Max
    TH1F* h_hits_MP2 = new TH1F("h_hits_MP2", "MP2 Hits [cm]",  400, -0.05, 0.06); 
    TH2F* h_true_hits = new TH2F("h_true_hits", "True X vs Y hits [cm]", 100,  -10, 10, 100, -10, 10);
    TH2F* h_det_hits = new TH2F("h_det_hits", "Tracker X vs Y hits [cm]", 100,  -10, 10, 100, -10, 10);
    TH2F* h_mis_hits = new TH2F("h_mis_hits", "Misaligned X vs Y hits [cm]", 100,  -10, 10, 100, -10, 10);

 
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
    cout << fixed << setprecision(9); 
    debug_mp2 << fixed << setprecision(9); 
    debug_calc << fixed << setprecision(9); 
    debug_mis << fixed << setprecision(9); 
    debug_geom << fixed << setprecision(9); 
    debug_off << fixed << setprecision(9); 
    debug_mc << fixed << setprecision(9); 
    debug_con << fixed << setprecision(9);

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
    for (int icount=0; icount<Tracker::instance()->getTrackCount(); icount++){  //Track count is set in methods XXX is it the best place for it? 
        rand_num = (RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max());
        float p=pow(10.0, 1+rand_num);
        scatterError=sqrt(Tracker::instance()->getWidth())*0.014/p;  
        
        //Generating tracks 
        LineData generated_line = Tracker::instance()->MC(scatterError, debug_calc, debug_off, debug_mc, debugBool);
       
        for (int i=0; i<generated_line.hit_count; i++){  //counting only hits going though detector
            //calculating the layer and pixel from the hit number 
            int lyr = (generated_line.i_hits[i]/Tracker::instance()->getPixelXYN());  // [1-14] //This the layer id
            int im = generated_line.i_hits[i]%Tracker::instance()->getPixelXYN();  // [0-49] //This the pixel id
            //  computes the remainder of the division of plane (e.g. MOD(693, 50) = 693 - 50*13 = 43) 

            //Number of local and global parameters 
            const int nalc = 4; 
            const int nagl = 2;  
            
            //Local derivatives
            float dlc1=Tracker::instance()->getProjectionX()[lyr];
            float dlc2=Tracker::instance()->getProjectionY()[lyr];
            float dlc3=generated_line.x_hits[i]*Tracker::instance()->getProjectionX()[lyr];
            float dlc4=generated_line.x_hits[i]*Tracker::instance()->getProjectionY()[lyr];  
            float derlc[nalc] = {dlc1, dlc2, dlc3, dlc4};
            //Global derivatives
            float dgl1 = Tracker::instance()->getProjectionX()[lyr];
            float dgl2 = Tracker::instance()->getProjectionY()[lyr];
            float dergl[nagl] = {dgl1, dgl2};  
            //Labels 
            int l1 = im+Tracker::instance()->getPixelXYN()*Tracker::instance()->getLayer()[lyr];  
            int l2 = im+Tracker::instance()->getPixelXYN()*Tracker::instance()->getLayer()[lyr]+1000; 
            int label[nalc] = {l1, l2}; 
            
            //TODO multiple scattering errors (no correlations) (for imodel == 1)
            //add break points multiple scattering later XXX (for imodel == 2)
            //! add 'broken lines' offsets for multiple scattering XXX (for imodel == 3)
           
            float rMeas_mp2 =  generated_line.y_hits[i]; 
            float sigma_mp2 = generated_line.hit_sigmas[i]; 

            //Sanity Plots 
            h_sigma -> Fill(sigma_mp2);
            h_hits_MP2 -> Fill (rMeas_mp2); 
            h_true_hits ->  Fill(generated_line.x_true[i], generated_line.y_true[i]);
            //h_true_hits->SetContour(ncol);
            //gStyle->SetPalette(100,colors);
            h_true_hits->Draw("colz");
            h_det_hits ->  Fill(generated_line.x_det[i], generated_line.y_det[i]);
            h_det_hits->Draw("colz");
            h_mis_hits ->  Fill(generated_line.x_mis[i], generated_line.y_mis[i]);
            h_mis_hits->Draw("colz");
                      
            //XXX passing by reference/pointer vs as variable (?) - makes any difference
            m.mille(nalc, derlc, nagl, dergl, label, rMeas_mp2, sigma_mp2);

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
    file->Write();
    file->Close(); //good habit!
   
    return 0; 
} //end of main 
