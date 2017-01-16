/** 
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
!!
**/

//TODO some includes might be redundant
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <random>
#include <cmath> //math class`
#include <TH1D.h> //1D Histo Root class
#include <TFile.h> // data records for ROOT 
#include <TRandom3.h> // rnd generator class for ROOT
#include "Mille.h"
#include "Mille.cc"
//////----Variable Intialisation-------///////////

//arguments for Mille constructor:
const char* outFileName = "Mptest2.bin";
bool asBinary = true; 
bool writeZero = false;

//initilising variables for the data file
int NLC = 0; // # of Local derivatives
int NGL = 0 ;  // # of Global parameters
float rMeas = 0.0;  // this will be a Gaussian random (0,x) for now - later to be filled be simulated residuals
float sigma = 0.0; // for now, for simplicity


///initialsing physics varibles
const int nlyr = 10; //number of detector layers
const int nmlyr = 14; //number of measurement layers
const int nmx = 10; //number of modules in x direction
const int nmy = 10; //number of modules in y direction

int ntot=nlyr*nmx*nmy; //total number of modules
//  define detector geometry
float dets= 10.0; // arclength of first plane
float diss= 10.0; // distance between planes
float thck= 0.02; //thickness of plane (X0)
float offs=  0.5;  // offset of stereo modules
float stereo=0.08727;  // stereo angle
float sizel= 20.0; //size of layers
float sigl =0.002;  // <resolution

int nhits = 0; // number of hits
float the0 = 0; // multiple scattering error
int[nmlyr]islyr;// (detector) layer
int[nmlyr]ihits; // module number
float[ntot]sdevx;// shift in x (alignment parameter)
float[ntot] sdevy; //shift in y (alignment parameter)
float[nmlyr]sarc;  // arc length
float[nmlyr]ssig;   //resolution
float[2][nmlyr]spro; //projection of measurent direction in (XY)
float[nmlyr]xhits;   //position perp. to plane (hit)
float[nmlyr]yhits;    //measured position in plane (hit)
float[nmlyr]sigma;    // measurement sigma (hit)

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

const int imodel = 0;  //straight line 
// Defining function Mptst2 
void Mptst2(imodel){

	float cmbbrl = 0.0;
    float dispxm =0.0;
    float dispym =0.0;
    float dn = 0.0;
    float dp = 0.0;
    float gran = 0.0;
    float one = 0.0;
    float p = 0.0;
    float s = 0.0;
    float sgn = 0.0;
    float sbrl = 0.0;
    float sold = 0.0;
    float uran = 0.0;
    float wbrl = 0.0;
    int i = 0;
    int ibrl = 0;
    int icount = 0;
    int im = 0;
    int ios = 0;
    int ip = 0;
    int j = 0;
    int k = 0;
    int l = 0; 
    int labelt = 0;
    int layer = 0;
    int lb = 0;
    int lbrl = 0;
    int luns = 0;
    int lunt = 0;
    int lyr = 0;
    int nalc  = 0;
    int nbrl = 0;
    int ncount = 0;
    int ncx = 0;
    int nmxy = 0;
    int nrecds = 0;
    int nthits  = 0;


    int size_gl = 1;  // XXX
    int size_lc = 1; // XXX
    float *derlc=new float[size_lc]{mlyr*2+3};  // XXX
    float *dergl=new float[size_gl]{nmlyr*2+3}; // XXX 
    int *label=new int[size_gl]{2};  // XXX
    bool ex1;  // XXX
    bool ex2; // XXX
    bool ex3; //XXX 

   // !     for broken lines: 1=fine, 2=coarse
    // DIMENSION nbrl(2),lbrl(nmlyr,2),sbrl(nmlyr,2),wbrl(nmlyr,2), cmbbrl(2)
    //DATA cmbbrl / 0.0, 1.0 / //! cut for combining layers


}




//Dynamic allocation for arrays (cleaning up in the end!)
//TODO 

//XXX:  --std=c++11 flag is required for pre-initilised arrays
//int *label=new int[size_gl]{1, 2, 3, 4};

//float *derGl=new float[size_gl]{1.0, 2.0, 3.0, 4.0};  


///////////------------Function prototyping----------//////////////////
void cleanUP();

void Mptst2(const int imodel); 


/////************MAIN***************/////////////
int main(){

Mille m (outFileName, asBinary, writeZero);  // call to Mille.cc to create a .bin file

TFile* file = new TFile("Mptest2.root", "recreate");  // recreate = owerwrite if already exisists
// Book histograms
TH1F* h_1 = new TH1F("h_1", "Plane 1",  20,  -0.05, 0.05); // D=doube bins, name, title, nBins, Min, Max


//keep this simple for now... TODO
NLC = size_lc;
NGL = size_gl;
sigma=100;   //TODO





cleanUP(); // TODO: Really (!) clarify how to use dynamic cleanup stuff..

//ROOT stuff
file->Write();
file->Close(); //good habit! 
return 0; 

} //end of main 



////////----------Function Defenition--------------------- ///// 
void cleanUP(){
//TODO cleaning up: read bookmark - see the neat and correct way to do this (e.g pointer array-object?) or if it is really necessary - smart pointers c++11? 
// When done, free memory pointed.
delete [] label; 
delete [] derGl; 
delete [] derLc;  
// Clear a to prevent using invalid memory reference.
derGl = NULL;
derLc = NULL; 
label = NULL;     
}




