/** 
* Translation into C++11 of mptest2.f90 (description of purpose there) 
* pass data as number of local and global parameters, and their derivatives, and residual and their sigma (smearing/accuracy)
* Mille will then pack this in .bin to be used by PEDE routine.
*
*
* Gleb Lukicov 11 Jan 2017 
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
const char* outFileName = "cppMptest2.bin";
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

//TODO continue line 77 













//Dynamic allocation for arrays (cleaning up in the end!)
//TODO 
int size_gl = 1;
int size_lc = 2;
//XXX:  --std=c++11 flag is required for pre-initilised arrays
//int *label=new int[size_gl]{1, 2, 3, 4};
int *label=new int[size_gl]{0};
//float *derGl=new float[size_gl]{1.0, 2.0, 3.0, 4.0};  
float *derGl=new float[size_gl]{0.0};  
float *derLc=new float[size_lc]{0.0, 0.0};   //same for all measurements 

///////////------------Function prototyping----------//////////////////
void cleanUP();

/////************MAIN***************/////////////
int main(){

Mille m (outFileName, asBinary, writeZero);  // call to Mille.cc to create a .bin file

TFile* file = new TFile("cppMptest2.root", "recreate");  // recreate = owerwrite if already exisists
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




