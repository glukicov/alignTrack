/** 
* File to interact with Mille:
* pass data as number of local and global parameters, and their derivatives, and residual and their sigma (smearing/accuracy)
* Mille will then pack this in .bin to be used by PEDE routine.
*
* clang++ MilleHandler.cpp -o MilleHandler --std=c++11   -flag required (see below) - now has a Makefile (gmake)
*
*
* Gleb Lukicov 17 Nov 2016 
**/

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
const char* outFileName = "test.bin";
bool asBinary = true; 
bool writeZero = false;

//initilising variables for the data file
int NLC = 0; // # of Local derivatives
int NGL = 0 ;  // # of Global parameters
float rMeas = 0.0;  // this will be a Gaussian random (0,x) for now - later to be filled be simulated residuals
float sigma = 0.0; // for now, for simplicity

//Dynamic allocation for arrays (cleaning up in the end!)
//TODO use of size in practice
int size_gl = 1;
int size_lc = 2;
//XXX:  --std=c++11 flag is required for pre-initilised arrays
//int *label=new int[size_gl]{1, 2, 3, 4};
int *label=new int[size_gl]{0};
//float *derGl=new float[size_gl]{1.0, 2.0, 3.0, 4.0};  
float *derGl=new float[size_gl]{0.0};  
float *derLc=new float[size_lc]{0.0, 0.0};   //same for all measurements 

//Gausian random generator
std::default_random_engine generator;
std::normal_distribution<float> generate_residual(0.0, 100); //centerred at 0

///////////------------Function prototyping----------//////////////////
void cleanUP();


/////************MAIN***************/////////////
int main(){

Mille m (outFileName, asBinary, writeZero);  // call to Mille.cc to create a test.bin file

TFile* file = new TFile("test.root", "recreate");  // recreate = owerwrite if already exisists
// Book histograms
TH1F* h_1 = new TH1F("h_1", "Plane 1",  20,  -0.05, 0.05); // D=doube bins, name, title, nBins, Min, Max
TH1F* h_2 = new TH1F("h_2", "Plane 2",  20,  -0.05, 0.05); 
TH1F* h_3 = new TH1F("h_3", "Plane 3",  20,  -0.05, 0.05); 
TH1F* h_4 = new TH1F("h_4", "Plane 4",  20,  -0.05, 0.05); 

//keep this simple for now... 
NLC = size_lc;
NGL = size_gl;
sigma=100;   


int nPlanes=4;
int nEvents=1000; 


for (int n=0; n<nEvents; ++n){
//For planes
	for (int i=0; i<nPlanes; ++i) {
		rMeas = generate_residual(generator);   //generate a residual
	
    	//TODO make this effficient! see MWPCPlots_module.cc for example! 
    	if(i==0){h_1->Fill(rMeas);}
    	if(i==1){h_2->Fill(rMeas);}
    	if(i==2){h_3->Fill(rMeas);}
    	if(i==3){h_4->Fill(rMeas);}
		derLc[0] = (i+1)*1E6; 
		derLc[1] = 1; 
		derGl[0] = 1;
		label[0] = i+1;   
		m.mille(NLC, derLc, NGL, derGl, label, rMeas, sigma);   //Add measurement to buffer 
	}//end of planes 
	
	m.end();  // TODO: where to put this w.r.t planes/events. Write buffer (set of derivatives with same local parameters) to file
		
}// end of events


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




