void Mom() {

    gStyle->SetOptStat(0);

    TCanvas *can = new TCanvas(" ", " ", 1600, 700);
    can->Divide(2);

    //open files
    TFile *file0 = TFile::Open("/Users/gleb/software/alignTrack/mpIIDESY/Curve/Truth/gm2tracker_ana.root");
    TFile *file1 = TFile::Open("/Users/gleb/software/alignTrack/mpIIDESY/Curve/Plus/gm2tracker_ana.root");
    TFile *file2 = TFile::Open("/Users/gleb/software/alignTrack/mpIIDESY/Curve/Minus/gm2tracker_ana.root");

    //get histos
    // TH1F* pTruth = (TH1F*)file0->Get("TrackSummaryS12/FitResults/P");
    // TH1F* pPlus = (TH1F*)file1->Get("TrackSummaryS12/FitResults/P");
    // TH1F* pMinus = (TH1F*)file2->Get("TrackSummaryS12/FitResults/P");

    TH1F* pTruth = (TH1F*)file0->Get("TrackSummaryS12/FitResults/pValues");
    TH1F* pPlus = (TH1F*)file1->Get("TrackSummaryS12/FitResults/pValues");
    TH1F* pMinus = (TH1F*)file2->Get("TrackSummaryS12/FitResults/pValues");

    //rebin
    pTruth->Rebin(5);
    pPlus->Rebin(5);
    pMinus->Rebin(5);

    // clone histors and devide
    TH1F* hTruth_pMinus = (TH1F*)pTruth->Clone();
    TH1F* hTruth_pPlus = (TH1F*)pTruth->Clone();

    // divide pmax/pmin
    hTruth_pPlus->Divide(pPlus);
    hTruth_pMinus->Divide(pMinus);

    // loop over
    std::string labels[] = {"Truth/Plus", "Truth/Minus"};
    double rangeMin[] = { 0.5, 0.5};
    double rangeMax[] = { 1.5, 1.5};
    TH1F* histos[] = {hTruth_pPlus, hTruth_pMinus};

    for (int i = 0; i < 2; i++) {
        can->cd(i+1);
        histos[i]->SetMarkerStyle(8);
        histos[i]->SetMarkerColor(8);
        histos[i]->GetYaxis()->SetRangeUser(rangeMin[i], rangeMax[i]);
        histos[i]->GetYaxis()->CenterTitle();
        double binWidth = histos[i]->GetBinWidth(1);
        std::cout << "binWidth=" << binWidth << "\n";
        labels[i] = labels[i] + " /" +std::to_string(binWidth) + " ";
        histos[i]->GetYaxis()->SetTitle(labels[i].c_str());
        histos[i]->GetYaxis()->SetTitleOffset(1.2);
        histos[i]->GetXaxis()->CenterTitle();
        histos[i]->Draw("P");
    }

    can->Draw();
    can->Print("Mom.png");
}