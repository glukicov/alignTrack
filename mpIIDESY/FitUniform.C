{
	// Open ROOT file of interest
	TFile fFile("TrackerAlignment.root");

	//Open histogram of interest [by coping it as an object]
	TH1F *h_uniform = (TH1F*)fFile.Get("TrackerAlignment/Tracks/h_pValue");

	h_uniform->SetDirectory(0); // to decouple it from the open file directory

	fFile.Close();

	//Set canvas and stat parameters
	TCanvas *cUnifrom = new TCanvas("cUnifrom", "cUnifrom", 700, 700);

	//Get parameters from the histogram
	float hBinMin = h_uniform -> FindFirstBinAbove(0, 1);
	float hBinMax = h_uniform -> FindLastBinAbove(0, 1);
	cout << "hBinMin= " << hBinMin << " hBinMax " << hBinMax << endl;
	float hBinNumber = hBinMax - hBinMin; // number of non-zero bins
	cout << " hBinNumber= " << hBinNumber;

	//Calculate mean value of all bins

	float hsum = 0.0;
	for (int i_bin=hBinMin; i_bin<=hBinMax; i_bin++){
		hsum += h_uniform->GetBinContent(i_bin);
	}

	float hMean = hsum/hBinNumber;
	cout << " hMean= " << hMean;

	TAxis *xaxis = h_uniform->GetXaxis();
	float minF = xaxis->GetBinLowEdge(hBinMin);
	float maxF = xaxis->GetBinUpEdge(hBinMax);

	cout << " minF= " << minF << " maxF= " << maxF << endl;

	//Set function to fit
	TF1* lineF = new TF1("lineF", "pol 0", minF, maxF);
	//lineF->SetParameters(0.0, hMean);
	//Fit function
	h_uniform->Draw("E1"); //Set errors on all bins
	h_uniform->Fit("lineF");

	//Save canvas as .png file
	//cUnifrom->Write();
	cUnifrom->Print("cUnifrom.png");
	
	// delete lineF;
	// delete h_uniform;
	// delete cUnifrom;
}