/**
   random_buffer.cpp

   Purpose: Reads plaintext files of pre-generated random numbers, with functions to read one random number at a time. This source file contains definitions of various functions in RandomBuffer class.

   @author John Smeaton
   @version 03/02/2017

 */

#include "random_buffer.h"

using namespace std;

// Set up empty pointer for instance of class.
RandomBuffer* RandomBuffer::s_instance = NULL;

/**
   Empty constructor for RandomBuffer class.
 */
RandomBuffer::RandomBuffer() {
}

/**
   Destructor for RandomBuffer class. Closes files of random numbers, to free resources for OS.
 */
RandomBuffer::~RandomBuffer() {
	uniform_file.close();
	gaussian_file.close();
}


/**
   Get pointer to only instance of RandomBuffer class, creating this instance if it doesn't already exist.

   @return Pointer to RandomBuffer instance.
 */
RandomBuffer* RandomBuffer::instance() {

	// Create pointer to class instance if one doesn't exist already, then return that pointer.
	if (s_instance == NULL) s_instance = new RandomBuffer();
	return s_instance;
}


/**
   Opens file of uniformly distributed random numbers, allowing them to be read.

   @param uniform_filename String of filename for the file of random numbers.
 */
void RandomBuffer::open_uniform_file(string uniform_filename) {
	cout << "Opened file of uniform randoms - " << uniform_filename << endl;
	uniform_file.open(uniform_filename.c_str());
}

/**
   Opens file of normally distributed random numbers, allowing them to be read.

   @param gaussian_filename String of filename for the file of random numbers.
 */

void RandomBuffer::open_gaussian_file(string gaussian_filename) {
	cout << "Opened file of gaussian randoms - " << gaussian_filename << endl;
	gaussian_file.open(gaussian_filename.c_str());
}


/**
   Read next line of file of uniformly distributed random numbers, and return number in that line.

   @return Next random number in file.
 */
float RandomBuffer::get_uniform_number() {

	// Test if file open, then read and return random number if it is. Otherwise, throw exception.
	if (uniform_file.is_open()) {
		string rand_num;
		uniform_file >> rand_num;
		return stof(rand_num);
	} else {
		throw "Please open file of uniform random numbers."
	}
}


/**
   Read next line of file of normally distributed random numbers, and return number in that line.

   @return Next random number in file.
 */
float RandomBuffer::get_gaussian_number() {

	// Test if file is open, then read and return random number if it is. Otherwise, throw exception.
	if (gaussian_file.is_open()) {
		string rand_num;
		gaussian_file >> rand_num;
		return stof(rand_num);
	} else {
		throw "Please open file of gaussian random numbers"
	}
}
