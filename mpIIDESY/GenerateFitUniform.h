#ifndef GENERATEFITUNIFROM
#define GENERATEFITUNIFROM

#include <TF1.h>
#include <TH1.h>
#include <TError.h>
#include <fstream>
#include "random_buffer.h" // courtesy of John Smeaton (UCL)

float generate_uniform(); // using the RandomBuffer class

void set_uniform_file(std::string); // Set filename for uniform random numbers [randomIntGenerator.py - see https://github.com/glukicov/alignTrack]

int main(int argc, char* argv[]);

static float twoR = 2.0;


#endif