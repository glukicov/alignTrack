// TODOs: 
//TODO move random intialisation outside of the genlin2() function, oterwise might get the same numbers...
// 0) FIX: Warning: insufficient constraint equations
// A) Generate RND to txt file and use from there
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
//#include "Logger.cc" //XXX - need to specicfy debug level and ouput file 

using namespace std; 

//// -- Initialising logger staff -- /// 
//add logger methods here
const unsigned int logLevel = 4; // DEBUG
//Logger l; XXX 
//Logger l (logLevel); 


//////----Variable Intialisation-------///////////
int imodel = 0;  //XXX Model type (see above)
int ip = 0;  // verbosity level of genlin2 [0= none, 1=verbose output] XXX
int track_count = 10000; // = number of generated tracks (recods) ///XXX

//arguments for Mille constructor:
const char* outFileName = "Mptest2.bin";
bool asBinary = true; 
bool writeZero = false;
string conFileName = "Mp2con.txt";
string strFileName = "Mp2str.txt";
//output ROOT file
TFile* file = new TFile("Mptest2.root", "recreate");  // recreate = owerwrite if already exisists
 // Book histograms
TH1F* h_1 = new TH1F("h_1", "Test",  1000,  -0.002, 0.004); // D=double bins, name, title, nBins, Min, Max
TH2F* h_2 = new TH2F("h_2", "Test",  100,  -2, 99, 100, -0.1 , 0.1); // D=double bins, name, title, nBins, Min, Max

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
float stereoTheta=0.08727;  // stereo angle  // radians (5 deg = 0.087.. rad)  
float layerSize= 20.0; //size of layers  //cm 
float resolution =0.002;  // <resolution  // 20um = 0.002 cm 

float scatterError = 0; // multiple scattering error
int layer[layerN];// (detector) layer
float sdevX[modulesTotalN];// shift in x (alignment parameter)
float sdevY[modulesTotalN] ; //shift in y (alignment GLOBAL parameter)
float arcLength[layerN];  // arc length
float resolutionLayer[layerN];   //resolution
float projection[2][layerN]; //projection of measurent direction in (XY)
//float sigma[layerN];    // measurement sigma (hit) //XXX passed from genlin2(hit_sigmas)

// Structure to contain data of a generated line, with the number of hits, their positions, the uncertainty in the positions, and the plane number hit.
struct Line_data {
    int hit_count;  // number of hits
    vector<float> x_hits;
    vector<float> y_hits; 
    vector<float> hit_sigmas; // measurment sigma 
    vector<int> i_hits;   //mdoule number 
};


///////----------Function Defenition--------------------- ///// 
//Source code for genlin courtesy of John. 
// Function to simulate a linear track through the detector, returning data about detector hits.

Line_data genlin2() {
    
     // TODO fix this! (segmentation fault if taken outstide the fucntion)
    //int seed = 123456;  //pseudorandom 
    //Mersenne Twister generator TR3 root class: 
    //TRandom3* rand_gen = new TRandom3(seed);  // ran->Rndm(); // uniform in [0,1] TRandom::Rndm 


    Line_data line;
    line.hit_count = 0;

    // Track parameters for rand-generated line
    float x_0 = layerSize * (rand_gen->Rndm()-0.5); //uniform vertex
    float y_0 = layerSize * (rand_gen->Rndm()-0.5); //uniform vertex 
    float x_1 = layerSize * (rand_gen->Rndm()-0.5); //uniform exit point: so fitting a line to these two points
    float y_1 = layerSize * (rand_gen->Rndm()-0.5); //uniform exit point: 
    float x_slope=(x_1-x_1)/arcLength[layerN];
    float y_slope=(y_1-y_0)/arcLength[layerN];
     

    if (ip != 0){
        cout << "" << endl;
        cout << "Track: " << "x0= " << x_0 << " y0= " << y_0 << "x_slope = " << x_slope << "y_slope = " << y_slope << endl;
   } 
    
    float x = x_0;
    float dx = x_slope;
    float y = y_0;
    float dy = y_slope;
    float sold = 0.0;  // XXX ??? 
    
    for(int i=0; i<layerN; i++){
        float ds = arcLength[i] - sold;
        sold = arcLength[i];

        //position with parameters 1. hit
        float xs=x_0 + arcLength[i] * x_slope;
        float ys=y_0 + arcLength[i] * y_slope;
        
        //true track position
        x=x+dx*ds;
        y=y+dy*ds;

        //multiple scattering
        dx = dx+ rand_gen-> Gaus(0,1) * scatterError;
        dy = dy+ rand_gen -> Gaus(0,1) * scatterError;

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
        line.y_hits.push_back((xl-xs)*projection[1][i]+(yl-ys)*projection[2][i]+ rand_gen -> Gaus(0,1)*resolutionLayer[i]);
        line.hit_sigmas.push_back(resolutionLayer[i]);
        line.hit_count++;

        h_2 -> Fill(line.x_hits[i], line.y_hits[i]);

        if (ip != 0){
            cout << "" << endl;
            cout << "Generated Line data: " <<  "Sigma= " << line.hit_sigmas[i] << " X hit= "<< line.x_hits[i] << " Y hit =" << line.y_hits[i] << endl;
        }
    }// end of looping over detector layers
    
    return line; // Return data from simulated track

} // end of genlin2




/////************MAIN***************/////////////
int main(){

        
    // Millepede courtesy of John 
   cout << endl;
    cout << "********************************************" << endl;
    cout << "*                 MPTEST 1                 *" << endl;
    cout << "********************************************" << endl;
    cout << endl;
    cout << "    _____________________________  \\  /" << endl;
    cout << "   {_|_|_|_|_|_|_|_|_|_|_|_|_|_|_( ͡° ͜ʖ ͡°) " << endl;
    cout << "    /\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/" << endl;
    cout << endl; 

       
   // Creating .bin, steering, constrating and ROOT files here:
    Mille m (outFileName, asBinary, writeZero);  // call to Mille.cc to create a .bin file
    
    cout << "" << endl;
    cout << "Generating test data for Mp II...";
    cout << "" << endl; 

    //Random number generator
    //int seed = 12345; 
    //TRandom3* rand_gen = new TRandom3(seed); // TODO fix this! (segmentation fault if taken outstide the fucntion)

    // fstreams for str and cons files 
    ofstream constraint_file(conFileName);
    ofstream steering_file(strFileName);
    
    float s=arcLength_Plane1;
    int i_counter = 0;
    float sign = 1.0;

    // Geometry of detecor arrangement 
    for (int layer_i=1; layer_i<=10; layer_i++){
        i_counter++;
        layer[i_counter] = layer_i;  // layer
        arcLength[i_counter] = s;  //arclength
        resolutionLayer[i_counter] = resolution; //resolution
        projection[1][i_counter]=1.0;  // x
        projection[2][i_counter]=0.0;  // y
        //taking care of stereo modules 
        if ((layer_i % 3) == 1){
            i_counter++;
            layer[i_counter] = layer_i;  // layer
            arcLength[i_counter] = s+offset;  //arclength
            resolutionLayer[i_counter] = resolution; //resolution
            projection[1][i_counter]=sqrt(1.0-pow(stereoTheta,2));  // x
            projection[2][i_counter]=stereoTheta*sign;  // y
            sign=-sign;
        }
        s=s+planeDistance;  // incrimenting distance between detecors 
    }  // end of looping over layers

    // XXX: definition of broken lines here in the future

    //Now misaligning detecors
    float dispX = 0.01; // module displacement in Y .05 mm * N(0,1)
    float dispY = 0.01; // module displacement in Y .05 mm * N(0,1)

    //so we are only displacing 9/10 detectors? XXX
    for (int i=0; i<detectorN-1; i++){   /// XXX
        for(int k=0; k<=moduleYN-1; k++){
            for(int l=1; l<=moduleXN; l++){
                sdevX[(i*moduleYN+k)*moduleXN+l] = dispX * rand_gen -> Gaus(0,1);  //XXX where is that used in imodel=0?? 
                sdevY[(i*moduleYN+k)*moduleXN+l] = dispY * rand_gen -> Gaus(0,1);         
            } // // end of number of modules in x
        } // end of number of modules in y 
    } // end of layers

    
    //Now writing steering and constraint files
    if(steering_file.is_open()){

        cout << "" << endl;
        cout<< "Writing the steering file" << endl;
        cout << "" << endl;


        steering_file <<  "*            Default test steering file" << endl
        << "Cfiles ! following bin files are Cfiles" << endl   // XXX 
        << "Mp2con.txt   ! constraints text file " << endl
        << "Mptest2.bin   ! binary data file" << endl
       // << "Cfiles       ! following bin files are Cfiles" << endl
        //<< "*outlierrejection 100.0 ! reject if Chi^2/Ndf >" << endl
        //<< "*outliersuppression 3   ! 3 local_fit iterations" << endl
        << "*hugecut 50.0     !cut factor in iteration 0" << endl
        << "*chisqcut 1.0 1.0 ! cut factor in iterations 1 and 2" << endl
        << "*entries  10 ! lower limit on number of entries/parameter" << endl
        << "" << endl 
        << "*pairentries 10 ! lower limit on number of parameter pairs"  << endl
        << "                ! (not yet!)"            << endl
        << "*printrecord   1  2      ! debug printout for records" << endl
        << "" << endl
        << "*printrecord  -1 -1      ! debug printout for bad data records" << endl
        << "" << endl
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
    } 


    int ncx = (moduleXN+1)/2; 
    int moduleXYN=0;
    int lunt = 9;
    float one = 1.0;

    if (constraint_file.is_open()) {
        cout << "" << endl;
        cout << "Writing Constraint File" << endl;
        cout << "" << endl;

        //Evaluation of constraints
        
        moduleXYN = moduleXN*moduleYN; 
        
        for (int i = 1; i <= detectorN; i=i+(detectorN-1)){  //XXX coorect implimentation of DO i=1,nlyr,nlyr-1
            constraint_file << "Constraint 0.0" << endl;
            for (int k=0; k<=moduleYN-1; k++){
                int labelt=(i*moduleYN+k)*moduleXN+ncx-1;
                constraint_file << labelt << " " << fixed << setprecision(7) << one<< endl;
                sdevX[((i-1)*moduleYN+k)*moduleXN+ncx]=0.0;      // fix center modules at 0.
            } // end of y loop
            constraint_file << "Constraint 0.0" << endl;
            for(int k=0; k<=moduleYN-1; k++){
                int labelt=(i*moduleYN+k)*moduleXN+ncx+1000-1;
                constraint_file << labelt << " " << fixed << setprecision(7) << one<< endl;
                sdevY[((i-1)*moduleYN+k)*moduleXN+ncx]=0.0; // fix center modules at 0.
            } // end of x loop
        } // end of detecors loop 

    } //end of constraints
   

    //Set up counters for hits and records (tracks)
    int hitsN = 0;
    int recordN=0;

    //Generating particles with energies: 10..100 Gev
    // track_count is set manually 
    for (int icount=0; icount<track_count; icount++){
        float p=pow(10.0, 1+rand_gen->Rndm());
        scatterError=sqrt(width)*0.014/p;

        //Generating hits for N=track_count
        
        //TODO fix *** Break *** illegal instruction

        Line_data generated_line = genlin2();

        for (int i=0; i<generated_line.hit_count; i++){
            //simple straight line
            int lyr = generated_line.i_hits[i]/moduleXYN+1;
            int im = generated_line.i_hits[i]%moduleXYN;
            const int nalc = 4; // XXX  number of LC paremeters? 
            const int nagl = 2; //XXX  number of GL paremeters? 
            
            float dlc1=projection[1][lyr];
            float dlc2=projection[2][lyr];
            float dlc3=generated_line.x_hits[i]*projection[1][lyr];
            float dlc4=generated_line.x_hits[i]*projection[2][lyr];  //XXX xhits again? 
            float derlc[nalc] = {dlc1, dlc2, dlc3, dlc4};
            
            float dgl1 = projection[1][lyr];
            float dgl2 = projection[2][lyr];
            float dergl[nagl] = {dgl1, dgl2};  
            
            int l1 = im+moduleXYN*layer[lyr];
            int l2 = im+moduleXYN*layer[lyr]+1000;
            int label[nalc] = {l1, l2};  ///XXX
            
            //multiple scattering errors (no correlations) (for imodel == 1)
            //add break points multiple scattering later XXX (for imodel == 2)
            //! add 'broken lines' offsets for multiple scattering XXX (for imodel == 3)
           
            float rMeas_mp2 =  generated_line.y_hits[i]; 
            float sigma_mp2 = generated_line.hit_sigmas[i]; 

            h_1 -> Fill(sigma_mp2);
            
            m.mille(nalc, derlc, nagl, dergl, label, rMeas_mp2, sigma_mp2);
            hitsN++; //count hits
        } // end of hits loop

        // XXX additional measurements from MS IF (imodel == 2) THEN

        //IF (imodel >= 3) THEN

        //cout << "Recored passed to bin file" << endl; 
        m.end(); // Write buffer (set of derivatives with same local parameters) to file.
        recordN++; // count records;
    } // end of N trials (track count)


    cout << " " << endl;
    cout << track_count << " tracks generated with " << hitsN << " hits." << endl;
    cout << recordN << " records written." << endl;
    cout << " " << endl;
    cout << "Ready for PEDE alogrithm: ./pede Mp2str.txt" << endl;
    cout << "Sanity Plots: root -l Mptest2.root" << endl; 


    //ROOT stuff
    file->Write();
    file->Close(); //good habit!
    delete rand_gen; // Cleaning up 
    return 0; 
} //end of main 
