#ifndef ALIGNTRACKER
#define ALIGNTRACKER

//This header file declares main function, and imports necessary dependencies.
#include "Mille.h"  // DESY 
#include "Mille.cc" // DESY
#include "AlignTracker_methods.h" // Methods and fucntions for main programme
#include "Logger.hh"  // Logger courtesy of Tom (from gm2trackdaq repository) 

#include <iostream>
#include <fstream> 
#include <string>
#include <vector>
#include <random>
#include <iomanip>
#include <chrono>
#include <typeinfo>
#include <stdlib.h>
#include <cmath> //math class`
#include <TH1D.h> //1D Histo Root class
#include <TH2D.h> //2D Histo Root class
#include <TH3D.h> //3D Histo Root class
#include <TFile.h> // data records for ROOT 
#include <TRandom3.h> // rnd generator class for ROOT
#include <TTree.h>

int main(int argc, char* argv[]);  

#endif
