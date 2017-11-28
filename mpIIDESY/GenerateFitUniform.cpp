#include "GenerateFitUniform.h"

using namespace std;

//Function Definition 

float generate_uniform() {
	float uniform = (( RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()));
	return uniform;
}

void set_uniform_file(string uniform_filename) {
	RandomBuffer::instance()->open_uniform_file(uniform_filename);
}


int main(int argc, char* argv[]) {
	//gErrorIgnoreLevel = kWarning; // Display ROOT Warning and above messages [i.e. suppress info]

	float nEvents = stoi(argv[1]);

	//File to continuously write Chi2/ndf
	ofstream Chi2Uni("Chi2Uni.txt", ios_base::app);

	//Open existing uniform rand file
	set_uniform_file("uniform_ran.txt");

	//Book histogram
	int hBinNumber = 80; float binEdge = 0.02; 
	TH1F *h_uniform = new TH1F("h_uniform", "h_uniform",  hBinNumber,  -binEdge, binEdge);
	h_uniform->SetDirectory(0); // to decouple it from the open file directory

	//Fill
	float functionEdge=0.015;
	for (int i_count = 0; i_count < nEvents; i_count++) {
		float signXSlope;
		if (generate_uniform() >= 0.5) {
			signXSlope = 1.0;
		}
		else {
			signXSlope = -1.0;
		}
		float xSlope = (generate_uniform() * signXSlope) * functionEdge;
		h_uniform->Fill(xSlope);

	} // end of fill

	//Set canvas and stat parameters
	//TCanvas *cUnifrom = new TCanvas("cUnifrom", "cUnifrom", 700, 700);

	//Set function to fit
	TF1* lineF = new TF1("lineF", "pol 0", -functionEdge, functionEdge);
	//h_uniform->Draw("E1"); //Set errors on all bins
	h_uniform->Fit("lineF", "Q");

	float Chi2 = lineF->GetChisquare();
	float ndf = lineF->GetNDF();
	float Chi2_ndf = Chi2/ndf;

	Chi2Uni << Chi2_ndf << endl;
	//cUnifrom->Print("cUnifrom.root");

	return 0; //Goodbye! 
}
