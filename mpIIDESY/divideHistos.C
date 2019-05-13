// opens and divides histos by truth/nominal

void divideHistos() {

    gStyle->SetOptStat(0); // no stats

    //iterations
    // states are the folder names also (!)
    vector<std::string> states{"Truth", "a=+1e-6", "a=-1e-6", "a=+0.5e-6", "a=-0.5e-6",  "a=+1e-7", "a=-1e-7"};
    // vector<std::string> states{"Truth", "a=+1e-6", "a=-1e-6"};
    int stateN = states.size();
    int splitN = 2; // Plus and Minus
    std::string labels[] = {"#frac{Positive}{Truth}", "#frac{Negative}{Truth}"}; // split labels on y-axis
    int curveN = 3; // each split has N/2-truth curvatures
    int colors[] = {2, 8, 9}; // red, green, blue
    int iterable_ids[] = {1, 3, 5}; // truth =0, "+"=1,3,5 "-"=2,4,6
    // int iterable_ids[] = {1}; // truth =0, "+"=1,3,5 "-"=2,4,6

    //Plotting constants
    int rebinFactor = 4;
    int mrkStyle = 8;

    //main dir
    std::string mainDir = "/Users/gleb/software/alignTrack/mpIIDESY/Curve/";
    std::string fileName = "/gm2tracker_ana.root";
    std::stringstream fullPath; // iterable stringstr. var.

    //define plots of interest
    std::string plotPath = "TrackSummaryS12/FitResults/";
    vector<std::string> plotName{"P", "pValues"};
    vector<std::string> plotUnits{"MeV", " "};
    int plotN = plotName.size();

    //keep all TFile and histos open
    TFile* fileArray[stateN];
    vector<TH1F*> truthArray;
    vector<TLegend*> legendArray;

    //open all files and push into array
    for (int i_state = 0; i_state < stateN; i_state++) {
        fullPath.str(""); fullPath << mainDir << states[i_state] << fileName;
        fileArray[i_state] = TFile::Open(fullPath.str().c_str());
    }

    // new canvas
    TCanvas *can = new TCanvas(" ", " ", 1200, 700);
    can->Divide(plotN, splitN);
    int i_pad = 1; // to cd into a pad

    // user range
    double rangeMin[] = {0.8, 0.3, 0.8, 0.8};
    double rangeMax[] = {2.3, 1.6, 1.4, 1.4};
    int i_total = 0; // range iterator

    // loop over plots and states 
    for (int i_plot = 0; i_plot < plotN; i_plot++) {
        //get truth histos
        fullPath.str(""); fullPath << plotPath << plotName[i_plot];
        std::cout << "fullPath: " << fullPath.str() << "\n";
        TH1F* pTruth = (TH1F*)fileArray[0]->Get(fullPath.str().c_str());
        pTruth->Rebin(rebinFactor);
        truthArray.push_back(pTruth);
        for (int i_split = 0; i_split < splitN; i_split++) {
            can->cd(i_pad); // move onto next pad 
            vector<TH1F*> splitArray; // temp storage 
            vector<TH1F*> dHistArray;
            TLegend* legend =  new TLegend(0.65, 0.88, 0.15, 0.60);
            legendArray.push_back(legend);
            for (int i_id = 0; i_id < curveN; i_id++) {
                int fileID = i_split + iterable_ids[i_id];
                std::cout << "fileID: " << fileID << "\n";
                splitArray.push_back( (TH1F*)fileArray[fileID]->Get(fullPath.str().c_str()) );
                TH1F* dHist = (TH1F*)splitArray[i_id]->Clone();
                dHist->Rebin(rebinFactor);
                dHistArray.push_back(dHist);
                dHistArray[i_id]->Divide(truthArray[i_plot]);
                dHistArray[i_id]->SetMarkerStyle(mrkStyle);
                dHistArray[i_id]->SetMarkerColor(colors[i_id]);
                dHistArray[i_id]->SetLineColor(colors[i_id]);
                dHistArray[i_id]->GetYaxis()->SetRangeUser(rangeMin[i_total], rangeMax[i_total]);
                dHistArray[i_id]->GetYaxis()->SetTitleSize(.048);
                dHistArray[i_id]->GetYaxis()->CenterTitle();
                dHistArray[i_id]->GetXaxis()->SetTitleSize(.048);
                dHistArray[i_id]->GetXaxis()->CenterTitle();
                double binWidth = dHistArray[i_id]->GetBinWidth(1);
                stringstream yTitle;  yTitle << labels[i_split] << fixed << setprecision(2) <<  "/ " << binWidth << " " << plotUnits[i_plot];
                dHistArray[i_id]->GetYaxis()->SetTitle(yTitle.str().c_str());
                legend->AddEntry(dHistArray[i_id], states[i_split + iterable_ids[i_id]].c_str() , "P"); 
                if(i_plot == 0) dHistArray[i_id]->GetXaxis()->SetRangeUser(500, 3060);
                if (i_id == 0) dHistArray[i_id]->Draw("P");
                if (i_id > 0) dHistArray[i_id]->Draw("P same");
            } // per on of  +/-
            i_pad++; // 1->4
            i_total++; // 0->3
            //draw legend once per pad
            legend->SetTextSize(.048);
            legend->Draw("same");
        } //per split +/-

    } // per plot

    can->Draw();
    can->Print("Mom.png");
}