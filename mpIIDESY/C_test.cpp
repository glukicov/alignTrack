 //TODO extend to reading random buffer.cpp? 

#include <iostream>
#include <fstream>  
#include <string>
#include <vector>
#include <random>
#include <iomanip>
#include <chrono>
#include <typeinfo>
#include <stdlib.h>
#include <cmath> //math class
#include <TH1D.h> //1D Histo Root class
#include <TH2D.h> //2D Histo Root class
#include <TH3D.h> //3D Histo Root class
#include <TFile.h> // data records for ROOT 
#include <TRandom3.h> // rnd generator class for ROOT
#include <TTree.h>
using namespace std; 

int main(){

	int const test = 7; 
	int const MUL = 1000;

	int const a = 5;
	int const c = 4; 
	
	float const bF  = 1.012003;
	float const dF  = 2.320421;
	double const bD = 1.012003;
	double const dD = 2.320421;

	string debugFileName = "C_test.txt";
	string pFileName = "C_P_test.txt"; 
	ofstream debug(debugFileName);
	ofstream p(pFileName);

	debug<< "Integer test number is " << test << endl;

	debug << std::setprecision(8);
	p<< setprecision(8);

	debug<< "Floats: precsion set to 8 decimal points " << endl;

	int Fcounter = 0;
	float FDiv, FSub, FMul, FAdd;
	float FResult = float(test);
	for (int i=0; i<1000; i++){
		FDiv = FResult/float(a);
		FSub = FDiv - bF; 
		FMul = FSub * float(c);
		FAdd = FMul + dF; 
		FResult = FAdd;
		Fcounter ++; 
	}
	Fcounter=0;
	
	debug<< "After " << Fcounter << " iterations " << " the result is " <<  FResult <<endl;
	p << FResult << endl; 
	
	
	debug << std::setprecision(17);
	p << std::setprecision(17);
	debug<< "Doubles: precsion set to 16 decimal points " << endl;

	int Dcounter = 0;
	double DDiv, DSub, DMul, DAdd;
	double DResult = double(test);
	for (int i=0; i<1000; i++){
		DDiv = DResult/double(a);
		DSub = DDiv - bD; 
		DMul = DSub * double(c);
		DAdd = DMul + dD;
		DResult = DAdd;  
		Dcounter ++; 
	}
	Dcounter=0;
	
	debug<< "After " << Dcounter << " iterations " << " the result is " <<  DResult <<endl;

	debug<< "" << endl;
	p << DResult<<endl;   


	debug<< "Integer test number is " << test << " multiplication factor is " <<  MUL <<endl;

	debug << std::setprecision(8);
	p << std::setprecision(8);

	debug<< "Floats: precsion set to 8 decimal points " << endl;

	Fcounter = 0;
	FDiv, FSub, FMul, FAdd = 0;
	FResult = float(test)*MUL;
	for (int i=0; i<1000; i++){
		FDiv = FResult/float(a);
		FSub = FDiv - bF*MUL; 
		FMul = FSub * float(c);
		FAdd = FMul + dF*MUL; 
		FResult = FAdd;
		Fcounter ++; 
	}

	FResult = FResult/MUL;
	
	debug<< "After " << Fcounter << " iterations " << " the result is " <<  FResult <<endl;
	p << FResult<<endl; 
	
	
	debug << std::setprecision(17);
	p << std::setprecision(17);
	debug<< "Doubles: precsion set to 16 decimal points " << endl;

	Dcounter = 0;
	DDiv, DSub, DMul, DAdd = 0;
	DResult = double(test)*MUL;
	for (int i=0; i<1000; i++){
		DDiv = DResult/double(a);
		DSub = DDiv - bD*MUL; 
		DMul = DSub * double(c);
		DAdd = DMul + dD*MUL;
		DResult = DAdd;  
		Dcounter ++; 
	}
	
	DResult = DResult/MUL;

	debug<< "After " << Dcounter << " iterations " << " the result is " <<  DResult <<endl;
	p<<DResult<<endl; 
	



	debug.close();
	return 0;
}