// TODOs: 
// 0) get simple case working 
// 1) rename variables to something more sensible 
// 2) Plot hits in ROOT - sanity plots  [check for loops for <= vs <]
// 3) add x [i=0], generate .bin file from Fortran for simplest case and compare
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
!!   Fit 4+2*(\ref mptest2::nmlyr "nmlyr"-2) parameters.
!!   Matrix of corresponding linear equation system is full and solution
!!   is obtained by inversion (time ~ parameters^3).
!! - \c BRLF: Use (fine) broken lines (see \ref ref_sec). Multiple scattering kinks
!!   are described by triplets of offsets at scatterers as track parameters.
!!   Fit 4+2*(\ref mptest2::nmlyr "nmlyr"-2) parameters. Matrix of corresponding
!!   linear equation system has band structure and solution
!!   is obtained by root-free Cholesky decomposition (time ~ parameters).
!! - \c BRLC: Use (coarse) broken lines. Similar to \c BRLF, but with stereo layers
!!   combined into single layer/scatterer. Fit 4+2*(\ref mptest2::nlyr "nlyr"-2) parameters.
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

using namespace std; 

//////----Variable Intialisation-------///////////

//XXX Model type (see above)
int imodel = 0; 

//arguments for Mille constructor:
const char* outFileName = "Mptest2.bin";
bool asBinary = true; 
bool writeZero = false;
string conFileName = "Mp2con.txt";
string strFileName = "Mp2str.txt";

///initialsing physics varibles
const int nlyr = 10; //number of detector layers
const int nmlyr = nlyr; //number of measurement layers   //why 14 is to do with stereo-angles for 1,4,7,10
//const int nmx = 10; //number of modules in x direction
const int nmy = 1; //number of modules in y direction   //total 50 modules nmx * nmy = 50

const int ntot=nlyr*nmy; //total number of modules
//  define detector geometry
float dets= 10.0; // arclength of first plane
float diss= 10.0; // distance between planes //cm / Pede works in cm
float thck= 0.02; //thickness of plane (X0)
//float offs=  0.5;  // offset of stereo modules
//float stereo=0.08727;  // stereo angle
float sizel= 20.0; //size of layers  //cm 
float sigl =0.002;  // <resolution  // 20um = 0.002 cm 

int nhits = 0; // number of hits  //XXX is this passed as argument 
//float the0 = 0; // multiple scattering error
int islyr[nmlyr];// (detector) layer
int ihits[nmlyr]; // module number
//float sdevx[ntot];// shift in x (alignment parameter)
float sdevy[ntot] ; //shift in y (alignment GLOBAL parameter)
float sarc[nmlyr];  // arc length
float ssig[nmlyr];   //resolution
float spro[1][nmlyr]; //projection of measurent direction in (Y) [change 1-> for XY]
//float xhits[nmlyr];   //position perp. to plane (hit)
float yhits[nmlyr];    //measured position in plane (hit)
float sigma[nmlyr];    // measurement sigma (hit)

// Structure to contain data of a generated line, with the number of hits, their positions, the uncertainty in the positions, and the plane number hit.
struct Line_data {
    int hit_count;
    vector<float> x_hits;
    //vector<float> y_hits;
    vector<float> hit_sigmas;
    //vector<float> y_drifts;
    vector<int> i_hits;
};

// Random number generators and distributions, for uniform and gaussian distribution
default_random_engine uniform_generator;
default_random_engine gaus_generator;
uniform_real_distribution<float> uniform_dist(0.0,1.0);
normal_distribution<float> gaus_dist(0.0, 1.0);


///////----------Function Defenition--------------------- ///// 
//Source code for genlin courtesy of John. 
// Function to simulate a linear track through the detector, returning data about detector hits.

Line_data genlin2() {

    // Set up new container for track data, with hit count set to zero
    Line_data line;
    line.hit_count = 0;

    // Track parameters for rand-generated line
    float ynull = sizel * (uniform_dist(uniform_generator)*0.5); //uniform vertex 
    float yexit = sizel * (uniform_dist(uniform_generator)*0.5); //uniform exit point: so fitting a line to these two points
    float yslop=(yexit-ynull)/sarc[nmlyr];

    nhits=0;
    float y = ynull;
    float dy = yslop;
    float sold = 0.0; 
    
    for(int i=0; i<nmlyr; i++){
        float ds = sarc[i] - sold;
        sold = sarc[i];

        //position with parameters 1. hit
        float ys=ynull+sarc[i]*yslop;
        //true track position
        y=y+dy*ds;

        float imy=int(y+sizel*0.5)/sizel*float(nmy);

        if (imy < 0. || imy >= nmy) continue;

        float ihit = ((i-1)*nmy+imy);
        int ioff=((islyr[i]-1)*nmy+imy)+1;
        nhits++;
        ihits[nhits]=ihit;
        float yl=y-sdevy[ioff];
        yhits[nhits]=(yl-ys)*spro[1][i]+gaus_dist(gaus_generator)*ssig[i];
        sigma[nhits]=ssig[i];

    }// end of looping over detector layers
    
    //101 FORMAT(3I3,5F8.4) TODO 

    return line; // Return data from simulated track

} // end of genlin2




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
       
   // Creating .bin, steering, constrating and ROOT files here:
    Mille m (outFileName, asBinary, writeZero);  // call to Mille.cc to create a .bin file
     // Book histograms
    TFile* file = new TFile("Mptest2.root", "recreate");  // recreate = owerwrite if already exisists
    TH1F* h_1 = new TH1F("h_1", "Test",  1000,  -100, 100); // D=double bins, name, title, nBins, Min, Max

    cout << "" << endl;
    cout << "Generating test data for Mp II...";
    cout << "" << endl; 

    // fstreams for str and cons files 
    ofstream constraint_file(conFileName);
    ofstream steering_file(strFileName);
    
    float s=dets;
    int i_counter = 0;
    float sgn = 1.0;

    for (int layer=1; layer<10; layer++){
        i_counter++;
        islyr[i_counter] = layer;  // layer
        sarc[i_counter] = s;  //arclength
        ssig[i_counter] = sigl; //resolution
        spro[1][i_counter]=1.0;  // y
        s=s+diss;  // incrimenting distance between detecors 
    }  // end of looping over layers

    //Now misaligning detecors
    float dispym = 0.01; // module displacement in Y .05 mm * N(0,1)

    //so we are only displacing 9/10 dets.? XXX
    for (int i=0; i<nlyr-1; i++){
        for(int k=0; k<=nmy-1; k++){
            sdevy[(i*nmy+k)+1] = dispym * uniform_dist(uniform_generator);          
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

    int nmxy = nmy; 
    if (constraint_file.is_open()) {
        cout << "" << endl;
        cout << "Writing Constraint File" << endl;
        cout << "" << endl;

        //Evaluation of constraints
       
        int lunt = 9;
        float one = 1.0;
        for (int i = 1; i <= nlyr; i=i+(nlyr-1)){  //XXX 
        constraint_file << "Constraint 0.0" << endl;
            //TODO fix contstrain file output
            for (int k=0; k<=nmy-1; k++){
                int labelt=(i*nmy+k)+1000-1;
                constraint_file << labelt << " " << fixed << setprecision(7) << one<< endl;
                sdevy[(i-1)*nmy+k]=0.0;      // fix center modules at 0.
            } // end of y loop
        } // end of detecors loop 


    } //end of constraints

    //record loop  TODO put this into genlin2
    int ncount = 10000;
    int nthits = 0;
    int nrecds=0;

    //Generating particles with energies: 10..100 Gev
    for (int icount=1; icount<ncount; icount++){
        float p=pow(10.0, 1+uniform_dist(uniform_generator));
        //the0=sqrt(thck)*0.014/p

        //Generating hits
        Line_data generated_line = genlin2();

        for (int i=1; i<nhits; i++){
            //simple straight line
            int lyr = ihits[i]/nmxy+1;
            int im = ihits[i]%nmxy;
            //const int nalc = 4; // XXX  number of LC paremeters? 
            const int nalc = 2; // XXX  number of LC paremeters? 
            //only xhits? XXX see line 320 fix this 
           
            float derlc[nalc] = {spro[1][lyr], yhits[i]*spro[1][lyr]};
            const int nagl = 2;
            float dergl[nagl] = {spro[1][lyr], spro[1][lyr]}; //XXX
            int label[nalc] = {im+nmxy*islyr[lyr], im+nmxy*islyr[lyr]+1000};
            //add multiple scattering errors later XXX

            if (imodel == 1){
                for (int j=1; i<nhits; i++){
                    sigma[j] = sqrt(pow(sigma[j],2) + pow(yhits[j]-yhits[i],2));  //XXX 
                }
            }

            //add break points multiple scattering later XXX
            
            m.mille(2, derlc, 2, dergl, label, yhits[i], sigma[i]);
            nthits++; //count hits
        } // end of hits loop

        // XXX additional measurements from MS

        //IF (imodel >= 3) THEN

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