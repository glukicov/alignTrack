#ifndef ALIGNTRACKER
#define ALIGNTRACKER

//This header file declares main function, and imports necessary dependencies for AlignTracker.cpp [description of purpose is there]
#include "Mille.h"  // courtesy of Gero Flucke (DESY) 
#include "Mille.cc" // courtesy of Gero Flucke (DESY) 
#include "AlignTracker_methods.h" // Methods and functions for the main programme (AlignTracker.cpp)
#include "Logger.hh"  // Logger courtesy of Tom Stuttard (UCL) - from gm2trackdaq repository

//XXX Some includes maybe redundant 
#include <iostream>
#include <fstream> 
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <chrono>
#include <string>
#include <vector>
#include <random>
#include <iomanip>
#include <chrono>
#include <typeinfo>
#include <stdlib.h>
#include <stdio.h>
#include <cmath> 
#include <TH1D.h>
#include <TF1.h> 
#include <TH2D.h> //2D Histo Root class
#include <TH3D.h> //3D Histo Root class
#include <TFile.h> // data records for ROOT 
#include <TRandom3.h> // Random generator class for ROOT
#include <TTree.h>
#include <TROOT.h>
#include <TClass.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLine.h>
#include <TApplication.h> //for Displaying Canvas
#include <TDirectory.h>
#include <THStack.h>
#include <TText.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <TRatioPlot.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMarker.h>
#include <TPaveStats.h>


int main(int argc, char* argv[]);  

#endif
