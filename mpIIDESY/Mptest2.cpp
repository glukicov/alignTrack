// TODOs: 


// B) See example code for improvements: encasulation, C++ design, classess [with private and public vars/methods], 
// Maps, Histogramm class, [singelton] Parametrs::Instance();   
// 1) add logger from gm2trackerdaq  
// 2) Plot hits in ROOT - sanity plots  [check for loops for <= vs <].
// 3) add [iModel=0], generate .bin file from Fortran for simplest case and compare.
// 4) extend to broken-lines, scattering etc.

/*
* 
* Translation into C++11 of mptest2.f90 (original description below) 
* pass data as number of local and global parameters, and their derivatives, and residual and their sigma (smearing/accuracy)
* Mille will then pack this in .bin to be used by PEDE routine.
*
* Gleb Lukicov 11 Jan 2017
-----------------------------------------------------
! Code converted using TO_F90 by Alan Miller
! Date: 2012-03-16  Time: 11:08:55

!> \file
!! MC for simple 10 layer silicon strip tracker.
!!
!! \author Claus Kleinwort, DESY, 2009

                              

!! No B-field, straight tracks. Selected with command line option '-t=track-model'
!! The \a track-models differ in the implementation of multiple scattering (errors):
!! - \c SL0: Ignore multiple scattering. Fit 4 track parameters.
!! - \c SLE: Ignore correlations due to multiple scattering, use only diagonal of
!!   m.s. covariance matrix. Fit 4 track parameters.
!! - \c BP: Intoduce explicit scattering angles at each scatterer.
!!   Fit 4+2*(\ref mptest2::layerN "layerN"-2) parameters.
!!   Matrix of corresponding linear equation system is full and solution
!!   is obtained by inversion (time ~ parameters^3).
!! - \c BRLF: Use (fine) broken lines (see \ref ref_sec). Multiple scattering kinks
!!   are described by triplets of offsets at scatterers as track parameters.
!!   Fit 4+2*(\ref mptest2::layerN "layerN"-2) parameters. Matrix of corresponding
!!   linear equation system has band structure and solution
!!   is obtained by root-free Cholesky decomposition (time ~ parameters).
!! - \c BRLC: Use (coarse) broken lines. Similar to \c BRLF, but with stereo layers
!!   combined into single layer/scatterer. Fit 4+2*(\ref mptest2::detectorN "detectorN"-2) parameters.
!!
!! MC for simple silicon strip tracker:
!! - 10 silicon detector layers
!! - 50 pixels per layer (1*2cm)
!! - 10 cm spacing, no B-field
!! - layers 1,4,7,10 have additional +/-5deg stereo modules   
!! - intrinsic resolution 20mu, 2% X0 per strip module
!! - uniform track offsets/slopes
!! - momentum: log10(p) 10..100 GeV uniform
!!
!! Global parameters:
!! - Position offsets (2D) in measurement plane per module (alignment).

// Generate test files.
//
// Create text and binary files.
//
//      unit  8: textfile Mp2str.txt   = steering file
//      unit  9: textfile Mp2con.txt   = constraint file
//      unit 51: binary file Mp2test.bin, written using CALL MILLE(.)
//      existing file are removed
//
// \param [in] imodel track model
//
//           0: 'straight line', ignoring multiple scattering 
//           1: 'straight line', using diagonal of m.s. error matrix
//           2: 'break points'
//           3: 'broken lines', fine
//           4: 'broken lines', coarse (stereo layers combined)
//

!!
**/

#include "Mptest2.h"
#include "Logger.hh"  

using namespace std; 


//// -- Initialising logger staff -- /// 
// TODO add logger methods here
//const unsigned int logLevel = 4; // DEBUG



/////************MAIN***************/////////////
//TODO add arguments option 
//int main(int argc, int* argv[]){
int main(){

    //////////////////////////////LOGGER EXPERIMENTING /////////////
    try {

        //Tell the logger to only show message at INFO level or above
        Logger::Instance()->setLogLevel(Logger::INFO); 
    
        //Tell the logger to throw exceptions when ERROR messages are received
        Logger::Instance()->enableCriticalErrorThrow();

        //Send an INFO message (messages are just strings)
        Logger::Instance()->write(Logger::INFO,"Hello from Logger");

        //Send an INFO message with a number in it (make the string using a stringstream)
        std::stringstream msg;
        msg << "5.0 + 5.0 = " << (5.0+5.0);
        Logger::Instance()->write(Logger::INFO,msg.str());

        //Another way to send an INFO message with a number, using std::to_string to turn a number into a string
        long double a = 1.0;
        long double b = 6.0;
        Logger::Instance()->write(Logger::INFO,"1.0 + 6.0 = " + std::to_string(a+b));

        //Send an ERROR message, this will terminate the program
        Logger::Instance()->write(Logger::ERROR,"Something terrible happened");

    }

    //Catch Logger exceptions
    catch (CriticalError& e) {
        std::cerr << "A critical error occurred, exiting" << std::endl;
        //return -1; //Exit program wth an error code
    }



    //////////////////////////////////////////////////////////////////




    // Millepede courtesy of John 
    cout << endl;
    cout << "********************************************" << endl;
    cout << "*                 MPTEST 2                 *" << endl;
    cout << "********************************************" << endl;
    cout << endl;
    cout << "    _____________________________  \\  /" << endl;
    cout << "   {_|_|_|_|_|_|_|_|_|_|_|_|_|_|_( ͡° ͜ʖ ͡°) " << endl;
    cout << "    /\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/" << endl;
    cout << endl;

    try {
        Detector::instance()->set_uniform_file("uniform_ran.txt");
   
        Detector::instance()->set_gaussian_file("gaussian_ran.txt");
        
        } 
        
    catch (ios_base::failure& e) {
        cerr << "Filestream exception caught: " << e.what() << endl;
        cerr << "Please ensure valid filenames are specified!" << endl;
        return 1;
        } 


    //////----Variable Intialisation-------///////////
    //TODO rewrite as arguments to main
    int imodel = 0;  //XXX Model type (see above)
    int ip = 0;  // verbosity level of genlin2 [0= none, 1=verbose output] XXX

    //arguments for Mille constructor:
    const char* outFileName = "Mptest2.bin";
    bool asBinary = true; 
    bool writeZero = false;

    bool debugBool = true; // print out to debug files 

    string conFileName = "Mp2con.txt";
    string strFileName = "Mp2str.txt";
    string debugFileName = "Mp2debug.txt"; 
    string mp2_debugFileName = "Mp2debug_mp2.txt";  // for looking at parameters going to CALL MILLE
    //output ROOT file
    TFile* file = new TFile("Mptest2.root", "recreate");  // recreate = owerwrite if already exisists
     // Book histograms
    TH1F* h_sigma = new TH1F("h_sigma", "Sigma",  100,  0, 0.004); // D=double bins, name, title, nBins, Min, Max
    TH1F* h_Xhits = new TH1F("h_Xhits", "X Hits",  100,  -20, 20); // D=double bins, name, title, nBins, Min, Max
    TH1F* h_Yhits = new TH1F("h_Yhits", "Y Hits",  100,  -20, 20); // D=double bins, name, title, nBins, Min, Max
    TH1F* h_hits_MP2 = new TH1F("h_hits_MP2", "MP2 Hits",  100, -0.05, 0.05); // D=double bins, name, title, nBins, Min, Max
            
    
/*
    // Check if correct number of arguments specified, exiting if not
    if (argc > 3) {
        cout << "Too many arguments - please specify model and verbosity flag. (e.g. 0 0 by default)" << endl << endl;
        return 1;
    } else if (argc < 3) {
        cout << "Too few arguments - please specify model and verbosity flag. (e.g. 0 0 by default) " << endl << endl;
        return 1;
    } else {

        // Set filenames to read random numbers from, using arguments. Catch exception if these files do not exist.
        try {
            imodel = argv[1];
            ip = argv[2];
        } catch  {
            cerr << "Arguments to Mptest2.cpp failed" << e.what() << endl;
            return 1;
        }
    
    }
*/
    
   // Creating .bin, steering, and constrain files
    Mille m (outFileName, asBinary, writeZero);  // call to Mille.cc to create a .bin file
    cout << "Generating test data for Mp II.." << endl;
   
    // fstreams for str and cons files 
    ofstream constraint_file(conFileName);
    ofstream steering_file(strFileName);
    ofstream debug(debugFileName);
    ofstream debug_mp2(mp2_debugFileName);
   
    // GEOMETRY
    Detector::instance()->setGeometry();
    
    // XXX: definition of broken lines here in the future
   
    // MISALIGNMENT
    Detector::instance()->misalign(); 

    // Write constraint file, for use with pede TODO fix this 
    Detector::instance()->write_constraint_file(constraint_file);

    //Now writing steering and constraint files
    if(steering_file.is_open()){

        cout<< "Writing the steering file" << endl;
        
        steering_file <<  "*            Default test steering file" << endl
        << "Cfiles ! following bin files are Cfiles" << endl  
        << "Mp2con.txt   ! constraints text file " << endl
        << "Mptest2.bin   ! binary data file" << endl
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

    //Set up counters for hits and records (tracks)
    int hitsN = 0;
    int recordN=0;
    
    
    float scatterError = Detector::instance()->getScatterError(); // multiple scattering error

    //Generating particles with energies: 10..100 Gev
       for (int icount=0; icount<Detector::instance()->getTrackCount(); icount++){
        float p=pow(10.0, 1+RandomBuffer::instance()->get_uniform_number());
        scatterError=sqrt(Detector::instance()->getWidth())*0.014/p;

        //Generating tracks 
        LineData generated_line = Detector::instance()->genlin2();
        if (debugBool){
            debug << endl; 
            debug_mp2 << endl; 
            debug << "–––––––––––––––––––––––––––––––––––––––––––––––" <<  endl;
            debug_mp2 << "–––––––––––––––––––––––––––––––––––––––––––––––" <<  endl;
            debug << "Track # (C)        " << icount+1 << endl;
            debug_mp2 << "Track # (C)        " << icount+1 << endl;
        }
        for (int i=0; i<generated_line.hit_count; i++){
            //calculating the layer and pixel from the hit number - TODO make this more readable by adding extra variables/containers 
            int lyr = generated_line.i_hits[i]/Detector::instance()->getPixelXYN();  // [1-14]
            int im = generated_line.i_hits[i]%Detector::instance()->getPixelXYN();  // [0-49]
            
            //cout << "lyr= " << lyr << " im= " << im << endl;

            const int nalc = 4; 
            const int nagl = 2;  
            
            //Local derivatives
            float dlc1=Detector::instance()->getProjectionX()[lyr];
            float dlc2=Detector::instance()->getProjectionY()[lyr];
            float dlc3=generated_line.x_hits[i]*Detector::instance()->getProjectionX()[lyr];
            float dlc4=generated_line.x_hits[i]*Detector::instance()->getProjectionY()[lyr];  //XXX xhits again? 
            float derlc[nalc] = {dlc1, dlc2, dlc3, dlc4};
            //Global derivatives
            float dgl1 = Detector::instance()->getProjectionX()[lyr];
            float dgl2 = Detector::instance()->getProjectionY()[lyr];
            float dergl[nagl] = {dgl1, dgl2};  
            //Labels 
            int l1 = im+Detector::instance()->getPixelXYN()*Detector::instance()->getLayer()[lyr]-1;
            int l2 = im+Detector::instance()->getPixelXYN()*Detector::instance()->getLayer()[lyr]+1000-1;
            int label[nalc] = {l1, l2}; 
            
            //multiple scattering errors (no correlations) (for imodel == 1)
            //add break points multiple scattering later XXX (for imodel == 2)
            //! add 'broken lines' offsets for multiple scattering XXX (for imodel == 3)
           
            float rMeas_mp2 =  generated_line.y_hits[i]; 
            float sigma_mp2 = generated_line.hit_sigmas[i]; 

            //Sanity Plots 
            h_sigma -> Fill(sigma_mp2);
           // h_Xhits-> Fill (); /TODO return from genlin2 as a new container vector<array[4]> ? 
           // h_Yhits-> Fill ();
            h_hits_MP2 -> Fill (rMeas_mp2); 

            
            m.mille(nalc, derlc, nagl, dergl, label, rMeas_mp2, sigma_mp2);

            // For debugging
            if (debugBool){
                debug_mp2 << endl; 
                debug_mp2 << " LC #:              " << nalc << "          LC1 :     " << derlc[0] << "        LC2 :  " << derlc[1] << "       LC3 :    " << derlc[2] << "   LC4 :     " << derlc[3] << endl
                          << " GL #:              " << nagl << "          GL1 :     " << dergl[0] << "       GL2 :  " << dergl[1] << endl
                          << " LB1 :              " << label[0] << "         LB2:      " << label[1] << "            Y Hit: " << rMeas_mp2 << "     Sigma : " << sigma_mp2 << endl;
                //debug_mp2 << endl; 
            }


            if (debugBool) {
                
                            debug<< "nhits: " << generated_line.hit_count << endl
                                 <<  "Hit #: " << i+1 << endl
                                 <<  "ihit: " << generated_line.i_hits[i] << endl
                                 << "x= " << generated_line.x_true[i] << "y= " << generated_line.y_true[i] << endl
                                // << "Local: " << derlc[0] << " " << derlc[1] << " " << derlc[2] << " " << derlc[3] << endl
                                // << "Global: " << dergl[0] << " " << dergl[1] << endl
                                // << "Label: " << label[0] << " " << label[1] << endl
                                  << "xhits = " << generated_line.x_hits[i] << endl
                                 << "yhits Hit= " << rMeas_mp2 << endl
                                 << "Sigma= " << sigma_mp2 << endl << endl;            }


        hitsN++; //count hits
        } // end of hits loop
         
        // XXX additional measurements from MS IF (imodel == 2) THEN

        //IF (imodel >= 3) THEN

        //cout << "Recored passed to bin file" << endl; 
        m.end(); // Write buffer (set of derivatives with same local parameters) to file.
        recordN++; // count records;
       
    } // end of rack count
   

    cout << " " << endl;
    cout << Detector::instance()->getTrackCount() << " tracks generated with " << hitsN << " hits." << endl;
    cout << recordN << " records written." << endl;
    cout << " " << endl;
    cout << "Ready for PEDE alogrithm: ./pede Mp2str.txt" << endl;
    cout << "Sanity Plots: root -l Mptest2.root" << endl; 

    // Close text files
    constraint_file.close();
    steering_file.close();
    debug.close();
    debug_mp2.close();
    //ROOT stuff
    file->Write();
    file->Close(); //good habit!
    return 0; 
} //end of main 
