// TODOs: 
// 0) get simple case working: hitsN not incimentig correctly via genline2()
// 0.5) add logger 
// 1) rename variables to something more sensible, and replace from Line_data. 
// 2) Plot hits in ROOT - sanity plots  [check for loops for <= vs <].
// 3) add x [i=0], generate .bin file from Fortran for simplest case and compare.
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

                              TODO impliment passing arguments to the executable

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
!! - 50 modules per layer (1*2cm)
!! - 10 cm spacing, no B-field
!! - layers 1,4,7,10 have additional +/-5deg stereo modules   TODO simplify this (all detectors the same for now)
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
//           0: 'straight line', ignoring multiple scattering // XXX using this explicitly for now 
//           1: 'straight line', using diagonal of m.s. error matrix
//           2: 'break points'
//           3: 'broken lines', fine
//           4: 'broken lines', coarse (stereo layers combined)
//

!!
**/

#include "Mptest2.h"
#include "Logger.cc" //XXX - need to specicfy debug level and ouput file 

using namespace std; 

//// -- Initialising logger staff -- ///
//add logger methods here
const unsigned int logLevel = 4; // DEBUG
//Logger l; XXX 
//Logger l (logLevel); 




//////----Variable Intialisation-------///////////
int imodel = 0;  //XXX Model type (see above)
int ip = 0;  // verbosity level of genlin2 [0= none, 1=verbose output] XXX 

//arguments for Mille constructor:
const char* outFileName = "Mptest2.bin";
bool asBinary = true; 
bool writeZero = false;
string conFileName = "Mp2con.txt";
string strFileName = "Mp2str.txt";
//output ROOT file
TFile* file = new TFile("Mptest2.root", "recreate");  // recreate = owerwrite if already exisists

///initialsing physics varibles
const int detectorN = 10; //number of detector layers
const int layerN = 14; //number of measurement layers   //XXX why 14 is to do with stereo-angles for 1,4,7,10
const int moduleXN = 10; //number of modules in x direction
const int moduleYN = 5; //number of modules in y direction   //total 50 modules moduleXN * moduleYN = 50

const int modulesTotalN=detectorN*moduleYN*moduleXN; //total number of modules
//  define detector geometry
float arcLength_Plane1= 10.0; // arclength of first plane
float planeDistance= 10.0; // distance between planes //cm / Pede works in cm
float width= 0.02; //thickness/width of plane (X0)
float offset=  0.5;  // offset of stereo modules
float stereoTheta=0.08727;  // stereo angle
float layerSize= 20.0; //size of layers  //cm 
float resolution =0.002;  // <resolution  // 20um = 0.002 cm 

//int hitsN = 0; // number of hits  //XXX passed from genlin2(hit_count)
float scatterError = 0; // multiple scattering error
int layer[layerN];// (detector) layer
//int moduleN[layerN]; // module number //XXX passed from genlin2(i_hits)
float sdevX[modulesTotalN];// shift in x (alignment parameter)
float sdevY[modulesTotalN] ; //shift in y (alignment GLOBAL parameter)
float arcLength[layerN];  // arc length
float resolutionLayer[layerN];   //resolution
float projection[2][layerN]; //projection of measurent direction in (XY)
//float hitsX[layerN];   //position perp. to plane (hit) //XXX passed from genlin2(x_hits)
//float hitsY[layerN];    //measured position in plane (hit)  //XXX passed from genlin2(y_hits)
//float sigma[layerN];    // measurement sigma (hit) //XXX passed from genlin2(hit_sigmas)

// Structure to contain data of a generated line, with the number of hits, their positions, the uncertainty in the positions, and the plane number hit.
struct Line_data {
    int hit_count;  // number of hits
    vector<float> x_hits;
    vector<float> y_hits; 
    vector<float> hit_sigmas; // measurment sigma 
    vector<int> i_hits;   //mdoule number 
};

// Random devices for seeding
random_device uniform_device;
random_device gaus_device;
// Distributions for random numbers
uniform_real_distribution<float> uniform_dist(0.0, 1.0);
normal_distribution<float> gaus_dist(0.0, 1.0);


/////************MAIN***************/////////////
int main(){

    //TODO draw a respectable-looking millipede here 
    cout << "" << endl;
    cout << "$            $"<< endl;
    cout << " $          $ "<< endl;
    cout << "  $        $ "<< endl;
    cout << "   $      $ "<< endl;
    cout << "     $ $ $ "<< endl;
    cout << "" << endl; 

    // Get sequences of seeds for random number generation
    seed_seq uniform_seeds{uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device()}; 
    seed_seq gaus_seeds{gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device()}; 

    // Set up Marsenne Twister random number generators with seeds
    mt19937 uniform_generator(uniform_seeds);
    mt19937 gaus_generator(gaus_seeds);
     
   // Creating .bin, steering, constrating and ROOT files here:
    Mille m (outFileName, asBinary, writeZero);  // call to Mille.cc to create a .bin file
     // Book histograms
    TH1F* h_1 = new TH1F("h_1", "Test",  1000,  -100, 100); // D=double bins, name, title, nBins, Min, Max

    cout << "" << endl;
    cout << "Generating test data for Mp II...";
    cout << "" << endl; 

    // fstreams for str and cons files 
    ofstream constraint_file(conFileName);
    ofstream steering_file(strFileName);
    
    float s=arcLength_Plane1;
    int i_counter = 0;
    float sgn = 1.0;

    for (int layer=1; layer<10; layer++){
        i_counter++;
        layer[i_counter] = layer;  // layer
        arcLength[i_counter] = s;  //arclength
        resolutionLayer[i_counter] = resolution; //resolution
        projection[1][i_counter]=1.0;  // y
        s=s+planeDistance;  // incrimenting distance between detecors 
    }  // end of looping over layers

    //Now misaligning detecors
    float dispym = 0.01; // module displacement in Y .05 mm * N(0,1)

    //so we are only displacing 9/10 detectors? XXX
    for (int i=0; i<detectorN-1; i++){
        for(int k=0; k<=moduleYN-1; k++){
            sdevY[(i*moduleYN+k)+1] = dispym * uniform_dist(uniform_generator);          
        } // end of number of modules in y 
    } // end of layers

    
    //Now writing steering and constraint files
    if(steering_file.is_open()){

        cout << "" << endl;
        cout<< "Writing the steering file" << endl;
        cout << "" << endl;

        steering_file <<  "*            Default test steering file" << endl
        << "Cfiles ! following bin files are C++" << endl   // XXX 
        << "Mp2con.txt   ! constraints text file " << endl
        << "Mptest2.bin   ! binary data file" << endl
        << "Cfiles       ! following bin files are Cfiles" << endl
        << "*outlierrejection 100.0 ! reject if Chi^2/Ndf >" << endl
        << "*outliersuppression 3   ! 3 local_fit iterations" << endl

        << "*hugecut 50.0     !cut factor in iteration 0" << endl
        << "*chisqcut 1.0 1.0 ! cut factor in iterations 1 and 2" << endl
        << "*entries  10 ! lower limit on number of entries/parameter" << endl
        <<  "" << endl 
        <<  "*pairentries 10 ! lower limit on number of parameter pairs"  << endl
        <<    "                ! (not yet!)"            << endl
        << "*printrecord   1  2      ! debug printout for records" << endl
        <<   "" << endl
        <<    "*printrecord  -1 -1      ! debug printout for bad data records" << endl
        <<   "" << endl
        <<   "*outlierdownweighting  2 ! number of internal iterations (> 1)"<< endl
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
    } 

    int moduleXNy = moduleYN; // TODO fix variable name  
    if (constraint_file.is_open()) {
        cout << "" << endl;
        cout << "Writing Constraint File" << endl;
        cout << "" << endl;

        //Evaluation of constraints
       
        int lunt = 9;
        float one = 1.0;
        for (int i = 1; i <= detectorN; i=i+(detectorN-1)){  //XXX 
        constraint_file << "Constraint 0.0" << endl;
            //TODO fix contstrain file output
            for (int k=0; k<=moduleYN-1; k++){
                int labelt=(i*moduleYN+k)+1000-1;
                constraint_file << labelt << " " << fixed << setprecision(7) << one<< endl;
                sdevY[(i-1)*moduleYN+k]=0.0;      // fix center modules at 0.
            } // end of y loop
        } // end of detecors loop 


    } //end of constraints

    //record loop  TODO put this into genlin2
    int ncount = 10000; // = number of generated tracks 
    //int ncount = 2;  XXX for debug // = number of generated tracks 
    int nthits = 0;
    int nrecds=0;

    //Generating particles with energies: 10..100 Gev
    for (int icount=1; icount<=ncount; icount++){
        float p=pow(10.0, 1+uniform_dist(uniform_generator));
        //scatterError=sqrt(width)*0.014/p

        //Generating hits
        Line_data generated_line = genlin2();

        //XXX HACK!!!
        hitsN = 10;

        for (int i=1; i<=hitsN; i++){
            //simple straight line
            int lyr = moduleN[i]/moduleXNy+1;
            int im = moduleN[i]%moduleXNy;
            //const int nalc = 4; // XXX  number of LC paremeters? 
            const int nalc = 2; // XXX  number of LC paremeters? 
            //only hitsX? XXX see line 320 fix this 
            derlc[nalc] = {projection[1][lyr], hitsY[i]*projection[1][lyr]}; ///XXX
            const int nagl = 2;
            dergl[nagl] = {projection[1][lyr], projection[1][lyr]}; //XXX  
            label[nalc] = {im+moduleXNy*layer[lyr], im+moduleXNy*layer[lyr]+1000};  ///XXX
            //add multiple scattering errors later XXX

            if (imodel == 1){
                for (int j=1; j<hitsN; j++){
                    sigma[j] = sqrt(pow(sigma[j],2) + pow(hitsY[j]-hitsY[i],2));  //XXX 
                }
            }

            //add break points multiple scattering later XXX
           
            #if 0
            cout << "derlc1= " << derlc[1] << endl; 
            cout << "derlc2= " << derlc[2] << endl; 
            cout << "dergl1= " << dergl[1] << endl; 
            cout << "dergl2= " << dergl[2] << endl; 
            cout << "label1= " << label[1] << endl; 
            cout << "label2= " << label[2] << endl;
            cout << "sigma= " << sigma[1] << endl; 
            cout << "label2= " << sigma[2] << endl; 
            #endif 

            m.mille(nalc, derlc, nagl, dergl, label, hitsY[i], sigma[i]);
            nthits++; //count hits
        } // end of hits loop

        // XXX additional measurements from MS

        //IF (imodel >= 3) THEN

        cout << "Recored passed to bin file" << endl; 
        m.end(); // Write buffer (set of derivatives with same local parameters) to file.
        nrecds++; // count records;

    } // end of N trials (track count)


    cout << " " << endl;
    cout << ncount << " tracks generated with " << nthits << " hits." << endl;
    cout << nrecds << " records written." << endl;
    cout << " " << endl;
    cout << "Ready for PEDE alogrithm: ./pede Mp2str.txt" << endl; 

    //ROOT stuff
    file->Write();
    file->Close(); //good habit! 
    return 0; 
} //end of main 


///////----------Function Defenition--------------------- ///// 
//Source code for genlin courtesy of John. 
// Function to simulate a linear track through the detector, returning data about detector hits.

Line_data genlin2() {

    // Get sequences of seeds for random number generation
    seed_seq uniform_seeds{uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device(), uniform_device()}; 
    seed_seq gaus_seeds{gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device(), gaus_device()}; 

    // Set up Mersenne Twister random number generators with seeds
    mt19937 uniform_generator(uniform_seeds);
    mt19937 gaus_generator(gaus_seeds);

    // Set up new container for track data, with hit count set to zero
    Line_data line;
    line.hit_count = 0;

    // Track parameters for rand-generated line
    float x_0 = layerSize * (uniform_dist(uniform_generator)-0.5); //uniform vertex
    float y_0 = layerSize * (uniform_dist(uniform_generator)-0.5); //uniform vertex 
    float x_1 = layerSize * (uniform_dist(uniform_generator)-0.5); //uniform exit point: so fitting a line to these two points
    float y_1 = layerSize * (uniform_dist(uniform_generator)-0.5); //uniform exit point: 
    float x_slope=(x_1-x_1)/arcLength[layerN];
    float y_slope=(y_1-y_0)/arcLength[layerN];

    if (ip =! 0){
        cout << "" << endl;
        cout << "Track: " << x_0 <<  y_0 << x_slope << y_slope << endl;
    }
    
    float x = x_0;
    float dx = x_slope;
    float y = y_0;
    float dy = y_slope;
    float sold = 0.0;  // XXX ??? 
    
    for(int i=0; i<=layerN; i++){
        float ds = arcLength[i] - sold;
        sold = arcLength[i];

        //position with parameters 1. hit
        float xs=x_0 + arcLength[i] * x_slope;
        float ys=y_0 + arcLength[i] * y_slope;
        
        //true track position
        x=x+dx*ds;
        y=y+dy*ds;

        //multiple scattering
        dx = dx+ gaus_dist(gaus_generator) * scatterError;
        dy = dy+ gaus_dist(gaus_generator) * scatterError;

        // TODO understand purpose of this part properly 
        float imx=int(x+layerSize*0.5)/layerSize*float(moduleXN);
        if (imx < 0. || imx >= moduleXN) continue;
        float imy=int(y+layerSize*0.5)/layerSize*float(moduleYN);
        if (imy < 0. || imy >= moduleYN) continue;      

        
        int ihit= ((i)*moduleYN+imy)*moduleXN; // XXX i from 0 to 14
        int ioff=((layer[i]-1)*moduleYN+imy)*moduleXN*imx+1;
        //line.i_hits.push_back(ihit);
        line.i_hits.push_back(i); //XXX which one do we need?
        float xl=x-sdevX[ioff];
        float yl=y-sdevY[ioff];
        line.x_hits.push_back(arcLength[i]);
        line.y_hits.push_back((xl-xs)*projection[1][i]+(yl-ys)*projection[2][i]+gaus_dist(gaus_generator)*resolutionLayer[i]);
        line.hit_sigmas.push_back(resolutionLayer[i]);
        line.hit_count++;


        if (ip =! 0){
            cout << "" << endl;
            cout << "Generated Line data: " <<   line.hit_sigmas[i] <<  ine.x_hits[i] << ine.y_hits[i] << ine.hit_sigmas[i] << endl;
        }
    }// end of looping over detector layers
    
    //101 FORMAT(3I3,5F8.4) XXX

    return line; // Return data from simulated track

} // end of genlin2

