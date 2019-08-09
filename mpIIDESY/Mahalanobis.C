void Mahalanobis() {


//Define constants
    const int stationN = 2;
    const int moduleN = 8;
    std::string labels[2] = {"S12", "S18"};
    const int colors[2] = {2, 4}; // red, blue
    const int markerStyle = 8; // dot
    int roundTo = 3;
    gStyle->SetOptFit(0); gStyle->SetOptStat(0); gROOT->ForceStyle(); // no default stats box

//Define survey meas.
    double Z_s12[moduleN-1] = {0.0, 138.4261909, 275.7030146, 412.8463143, 684.9748242, 820.6609492, 956.1864672};
    double Z_s18[moduleN] = {0.0, 138.4261909, 275.7030146, 412.8463143, 551.9417241, 684.9748242, 820.6609492, 956.1864672};
    double errX_s12[moduleN-1] = { 0 };  // all zeroes
    double errX_s18[moduleN] = { 0 };  // all zeroes
    double* errX[stationN] = {errX_s12, errX_s18}; 
    double measuredR_s12[moduleN-1] = {1.72, 2.08, 2.15, 2.29, 2.56, 2.73, 2.92};
    double measuredR_s18[moduleN] = {3.05, 3.43, 3.81, 4.06, 4.46, 4.75, 5.05, 5.38};
    double* measuredR[stationN] = {measuredR_s12, measuredR_s18};
    double* Z[stationN] = {Z_s12, Z_s18};
    double erR = 0.2; // 200 um / 0.2 mmm
    double errR_s12[moduleN-1]; std::fill_n(errR_s12, 7, erR); // fill with constant value
    double errR_s18[moduleN]; std::fill_n(errR_s18, moduleN, erR); // fill with constant value
    double* errR[stationN] = {errR_s12, errR_s18}; 

    // Fit function is'
    double pivotPoint = (Z_s18[3] + Z_s18[4])/2; 
    // double pivotPoint = 0.0; 
    std::cout << "pivotPoint= "<< pivotPoint << "\n";
    // ax^2 + bx +c or bx +c
    std::string curve = "[0]*(x-"+to_string(pivotPoint)+")*(x-"+to_string(pivotPoint)+") + [1]*(x-"+to_string(pivotPoint)+") + [2]"; // bx + c
    // std::string curve = "[0]*x*x + [1]*x + [2]"; // bx + c
    const int parN = 3; // curvature, slope, intercept (a, b, c)

    std::cout << "RShift from Curve: ";
    for (int i_module=0; i_module < moduleN; i_module++){
        // std::cout << 1.0*pow(10, -6)*(Z[i_module]*Z[i_module]) << ", ";
        std::cout << 0.2*pow(10, -6)*(Z_s18[i_module]-pivotPoint)*(Z_s18[i_module]-pivotPoint) << " ";
    }
    std::cout<<"\n";

//Make new canvas and legend with plotting range
    TCanvas* canvas_survey = new TCanvas("canvas_survey", "", 900, 800);
    TLegend* legend = new TLegend(0.1, 0.89, 0.45, 0.4);
    double rmin(2.3), rmax(6.0), zmin(-500.0), zmax(980.0);

    std::vector<TGraphErrors*> tge_vector; // keep all graphs in scope
    std::vector<TFitResultPtr> fit_vector; // keep all fits results in scope

// Loop over the two stations and do the initial fit
    for (int i_station = 0; i_station < stationN; i_station++) {

        int totalModules = moduleN;
        if (i_station == 0 ) totalModules=totalModules-1;

        //Plot data -
        TGraphErrors* tge = new TGraphErrors(totalModules, Z[i_station], measuredR[i_station], errX[i_station], errR[i_station]);
        tge_vector.push_back(tge);
        tge->SetTitle("; Z-position of modules [mm]; #Delta R (error = #pm 0.2) [mm]");
        tge->GetXaxis()->CenterTitle(); tge->GetYaxis()->CenterTitle();
        tge->GetYaxis()->SetRangeUser(rmin, rmax); tge->GetXaxis()->SetRangeUser(zmin, zmax);
        tge->SetLineColor(colors[i_station]); tge->SetMarkerColor(colors[i_station]); tge->SetMarkerStyle(markerStyle);
        if (i_station == 0) tge->Draw("AP");
        if (i_station != 0) tge->Draw("P same");

        //Fit data
        std::stringstream function_name; function_name << "tf_" << i_station; //unique function name
        TF1* tf = new TF1(function_name.str().c_str(), curve.c_str(), Z[i_station][0], Z[i_station][totalModules - 1]) ; // line fit in Z range
        TFitResultPtr fit_result = tge->Fit(function_name.str().c_str(), "QRS"); // Quite Chi2 fit over Range of the function, with pointer return (S)
        fit_vector.push_back(fit_result);
        tge->GetFunction(function_name.str().c_str())->SetLineColor(colors[i_station]);
        // double slope = tge->GetFunction(function_name.str().c_str())->GetParameter(0);
        // double slope_error = tge->GetFunction(function_name.str().c_str())->GetParError(0);
        // double angle_deg = atan(slope) * (180 / TMath::Pi());
        // double angle_deg_error = slope_error * (180 / TMath::Pi());
        double a = tge->GetFunction(function_name.str().c_str())->GetParameter(0);
        double a_error = tge->GetFunction(function_name.str().c_str())->GetParError(0);
        double b = tge->GetFunction(function_name.str().c_str())->GetParameter(1);
        double b_error = tge->GetFunction(function_name.str().c_str())->GetParError(1);
        double c = tge->GetFunction(function_name.str().c_str())->GetParameter(2);
        double c_error = tge->GetFunction(function_name.str().c_str())->GetParError(2);

        //Display fit results
        legend->SetBorderSize(0);
        // std::stringstream legend_string; legend_string << std::setprecision(roundTo) << labels[i_station] << " #theta= " << angle_deg << " #pm " << angle_deg_error;
        std::stringstream legend_string; legend_string << std::setprecision(roundTo) << labels[i_station] << " a= " << a << " #pm " << a_error << " b= " << b << " #pm " << b_error << " c= " << c << " #pm" << c_error;
        std::cout << legend_string.str() << "\n";
        legend->AddEntry(tge, legend_string.str().c_str(), "l");
        gStyle->SetLegendTextSize(0.028);

    } // per station loop

    // legend->Draw("SAME");
    canvas_survey->Draw();
    canvas_survey->Print("Survey.png");

// Mahalanobis fits
    int col[] = {3, 7, 6, 7, 5, 2, 4, 1, 12, 3, 7, 6, 7, 5, 2, 4, 1, 12, 3, 7, 6, 7, 5, 2, 4, 1, 12, 3, 7, 6, 7, 5, 2, 4, 1, 12, 3, 7, 6, 7, 5, 2, 4, 1, 12};

// Pick 1 sigma points per station
    for  (int i_station = 0; i_station < stationN; i_station++) {

        int totalModules = moduleN;
        if (i_station == 0 ) totalModules=totalModules-1;

        //Clear previous canvas
        canvas_survey->Clear();
        legend->Clear();

        //get parameter values
        TVectorD fitVals(parN);
        std::stringstream function_name; function_name << "tf_" << i_station; //unique function name
        for (int i_par = 0; i_par < parN; i_par++) fitVals[i_par] = tge_vector[i_station]->GetFunction(function_name.str().c_str())->GetParameter(i_par);

        // Covariance matrix
        TMatrixD covMatrix = fit_vector[i_station]->GetCovarianceMatrix();
        TDecompChol decompCholCov(covMatrix);
        decompCholCov.Decompose();
        TMatrixD matrixCovI = decompCholCov.GetU();
        TMatrixD matrixCov = TMatrixD(parN, parN);
        matrixCov.Transpose(matrixCovI);

        double CL =  0.317311; // 1 - 68% (one sigma)
        double r2 = 0; // r^2 (see https://upload.wikimedia.org/wikipedia/commons/a/a2/Cumulative_function_n_dimensional_Gaussians_12.2013.pdf)

        //compute the Mahalanobis distance for number for par
        while (TMath::Prob(r2, parN) > CL) {
            r2 += 0.00001;
        }
        double r = sqrt(r2); // This is the Mahalanobis distance threshold under which CL % of points fall below
        std::cout << "Mahalanobis distance is " << r << " in " << labels[i_station] << "\n";

        //draw the nominal line
        tge_vector[i_station]->Draw("AP");

        int i_state = 0; // this is a counter, expect 3^parN states: 3 base state (+1, 0, -1) with ^parN permutations
        std::vector<double> a, b, c; // containers to store slopes and intercepts

        //loop over the 1-sigma changes
        for (int i_curve = -1; i_curve < 2; i_curve++) {
        for (int i_slope = -1; i_slope < 2; i_slope++) {
            for (int i_intercept = -1; i_intercept < 2; i_intercept++) {
                TVectorD stateVector(parN);
                stateVector[0] = i_curve; stateVector[1] = i_slope; stateVector[2] = i_intercept;

                double scale = stateVector.Norm2Sqr() > 0 ? r / sqrt(stateVector.Norm2Sqr()) : 1;
                for (int i_par = 0; i_par < parN; i_par++) stateVector[i_par] *= scale;

                // z = r*u*M + mean
                // calculate the new state parameters
                TVectorD newParameters = matrixCov * stateVector + fitVals;
                for (int i_par = 0; i_par < parN; i_par++) std::cout << "newParameter= " << newParameters[i_par] << " ";
                std::cout << endl;

                //set new state parameters for the fit
                TF1* stateFit = new TF1(Form("State_%d_%d_%d", i_curve, i_slope, i_intercept), curve.c_str(), Z[i_station][0]-10, Z[i_station][totalModules - 1]+10);
                for (int i_par = 0; i_par < parN; i_par++) stateFit->SetParameter(i_par, newParameters[i_par]);

                // bold dashed lines for state fits
                stateFit->SetLineWidth(3); stateFit->SetLineStyle(2); stateFit->SetLineColor(col[i_state]);
                stateFit->Draw("same");

                //get new slopes and intercepts
                double a_state = stateFit->GetParameter(0);
                double b_state = stateFit->GetParameter(1);
                double c_state = stateFit->GetParameter(2);
                // double b_state_deg = atan(b_state) * ( 180 / TMath::Pi() );
                 a.push_back(a_state); b.push_back(b_state); c.push_back(c_state);

                //print keeping the +/- with showpos
                std::stringstream legend_string; legend_string << std::showpos << std::setprecision(roundTo) << i_curve << i_slope << i_intercept << " " << labels[i_station] << " a:" << a_state << " b: " << b_state << " c: " << c_state;
                legend->AddEntry(stateFit, legend_string.str().c_str(), "l");

                i_state++;

            } // intercept
        } // slope
    } // curve
        //save per station
        // legend->Draw("same");
        canvas_survey->Draw();
        std::stringstream plotName; plotName << "Mach_line_" << labels[i_station] << ".png";
        canvas_survey->Print(plotName.str().c_str());

        //Clear previous canvas
        canvas_survey->Clear();
        legend->Clear();

        // Summary plots
        TGraph* summary_ab = new TGraph(i_state, a.data(), b.data());
        summary_ab->SetTitle(Form("%s; a; b", labels[i_station].c_str()));
        summary_ab->GetXaxis()->SetTitle("a");
        summary_ab->GetYaxis()->SetTitle("b");
        summary_ab->GetXaxis()->CenterTitle();
        summary_ab->GetYaxis()->CenterTitle();
        summary_ab->SetMarkerColor(colors[i_station]);
        summary_ab->SetMarkerStyle(markerStyle);
        summary_ab->SetMarkerSize(2);
        summary_ab->Draw("AP");
        plotName.str(""); plotName << "Mach_point_ab_" << labels[i_station] << ".png";
        canvas_survey->Print(plotName.str().c_str());
        canvas_survey->Clear();

        TGraph* summary_ac = new TGraph(i_state, a.data(), c.data());
        summary_ac->SetTitle(Form("%s; a; c", labels[i_station].c_str()));
        summary_ac->GetXaxis()->CenterTitle();
        summary_ac->GetYaxis()->CenterTitle();
        summary_ac->SetMarkerColor(colors[i_station]);
        summary_ac->SetMarkerStyle(markerStyle);
        summary_ac->SetMarkerSize(2);
        summary_ac->Draw("AP");
        plotName.str(""); plotName << "Mach_point_ac_" << labels[i_station] << ".png";
        canvas_survey->Print(plotName.str().c_str());
        canvas_survey->Clear();

        TGraph* summary_bc = new TGraph(i_state, b.data(), c.data());
        summary_bc->SetTitle(Form("%s; b; c", labels[i_station].c_str()));
        summary_bc->GetXaxis()->CenterTitle();
        summary_bc->GetYaxis()->CenterTitle();
        summary_bc->SetMarkerColor(colors[i_station]);
        summary_bc->SetMarkerStyle(markerStyle);
        summary_bc->SetMarkerSize(2);
        summary_bc->Draw("AP");
        plotName.str(""); plotName << "Mach_point_bc_" << labels[i_station] << ".png";
        canvas_survey->Print(plotName.str().c_str());
        canvas_survey->Clear();

    } // Mah per station

} // main