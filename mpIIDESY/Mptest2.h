#ifndef MPTEST2
#define MPTEST2

//TODO some includes might be redundant  
#include "Mille.h"
#include "Mille.cc"
#include <iostream>
#include <fstream> 
#include <string>
#include <vector>
#include <random>
#include <iomanip>
#include <stdlib.h>
#include <random>
#include <cmath> //math class`
#include <TH1D.h> //1D Histo Root class
#include <TH2D.h> //1D Histo Root class
#include <TFile.h> // data records for ROOT 
#include <TRandom3.h> // rnd generator class for ROOT
#include <TTree.h>

std::vector<float> genlin2(int&, std::vector<float>&, std::vector<float>&, std::vector<float>&, std::vector<int>&);

int main();

#endif
