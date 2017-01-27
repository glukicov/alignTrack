#ifndef MPTEST1
#define MPTEST1

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <typeinfo>

#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>

#include <Mille.h>
#include <Mille.cc>

std::vector<float> genlin(int&, std::vector<float>&, std::vector<float>&, std::vector<float>&, std::vector<int>&);

int main(int argc, char* argv[]);

#endif
