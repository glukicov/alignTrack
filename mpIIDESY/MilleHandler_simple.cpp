/**
* File to interact with Mille: SIMPLE VERISON (1 entry only)
* pass data as number of local and global parameters, and their derivatives, and residual and their sigma (smearing/accuracy)
* Mille will then pack this in .bin to be used by PEDE routine.
*
* clang++ MilleHandler.cpp -o MilleHandler --std=c++11   -flag required (see below)
*
*
* Gleb Lukicov 17 Nov 2016 
**/

#include <fstream>
#include <string.h>
#include <stdlib.h>
#include "Mille.h"
#include "Mille.cc"

//arguments for Mille constructor:
const char* outFileName = "test_simple.bin";
bool asBinary = true; 
bool writeZero = false;

//declaring variables for the data file
int NLC = 0; // # of Local derivatives
int NGL = 0 ;  // # of Global parameters
float rMeas = 0.0;  // this will be a Gaussian random (0,x) for now - later to be filled be simulated residuals
float sigma = 0.0; // for now, for simplicity

//Dynamic allocation for arrays (cleaning up in the end!)
//TODO use of size in practice
int size = 1;
//XXX:  --std=c++11 flag is required for pre-initilised arrays
int *label=new int[size]{1};
float *derGl=new float[size]{0.0};  // 'std::bad_alloc' if 4x4... 
float *derLc=new float[size]{1.0};

///////////Function prototyping

void cleanUP();


int main(){

Mille m (outFileName, asBinary, writeZero);  // call to Mille.cc to create a test.bin file

// TODO add looping over events and planes here to set values

NLC = 2;
NGL = 0;
sigma=1E-10;  //cannot be zero! 
rMeas = 4.0;

// TODO set array value via looping 
//for those which values do change with hits and planes... 

m.mille(NLC, derLc, NGL, derGl, label, rMeas, sigma);   //Add measurement to buffer
label[0]=2;
m.mille(NLC, derLc, NGL, derGl, label, rMeas, sigma);   //Add measurement to buffer
m.end();  // Write buffer (set of derivatives with same local parameters) to file
m.kill(); // Reset buffers, i.e. kill derivatives accumulated for current set

cleanUP(); 

return 0; 

} //end of main 



//////// Function Defenition ///// 
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




