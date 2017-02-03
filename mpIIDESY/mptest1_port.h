/**
   mptest1_port.h

   Purpose: Port of mptest1.f90 Mille test program. Simulates a plane drift chamber, with variable plane offset and drift velocity. Writes global and local derivatives to a binary file, and writes appropriate steering and constraint files. This header file declares main function, imports necessary dependencies.

   @author John Smeaton
   @version 03/02/2017

 */

#ifndef MPTEST1
#define MPTEST1

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <iomanip>
#include <chrono>
#include <typeinfo>

#include "TFile.h"
#include "TTree.h"

#include "mptest1_port_detector.h"

#include "Mille.h"
#include "Mille.cc"

int main(int argc, char* argv[]);

#endif
