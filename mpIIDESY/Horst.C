//numbers taken from Tracker-locations-v1.xslx
void Horst() {

  double distanceFromM00[8] = {0.0, 138.4261909, 275.7030146, 412.8463143, 551.9417241, 684.9748242, 820.6609492, 956.1864672};
  double errX[8] = {0., 0., 0., 0., 0., 0., 0., 0.};

  double measuredR_s12[8] = {1.72, 2.08, 2.15, 2.29, 2.92, 2.56, 2.73, 2.92};
  double measuredR_s18[8] = {3.05, 3.43, 3.81, 4.06, 4.46, 4.75, 5.05, 5.38};

  double measuredY_s12[8] = {0.573, 0.504, 0.553, 0.615, 0.605, 0.506, 0.501, 0.663};
  double measuredY_s18[8] = {0.679, 0.545, 0.507, 0.361, 0.350, 0.423, 0.208, 0.291};

  double eR = 0.2;
  double errR[8] = {eR, eR, eR, eR, eR, eR, eR, eR};

  double eY = 0.2;
  double errY[8] = {eY, eY, eY, eY, eY, eY, eY, eY};

  TGraphErrors* gR_12 = new TGraphErrors(8, distanceFromM00, measuredR_s12, errX, errR );
  TGraphErrors* gR_18 = new TGraphErrors(8, distanceFromM00, measuredR_s18, errX, errR );
  TGraphErrors* gY_12 = new TGraphErrors(8, distanceFromM00, measuredY_s12, errX, errY );
  TGraphErrors* gY_18 = new TGraphErrors(8, distanceFromM00, measuredY_s18, errX, errY );

  TCanvas* c1 = new TCanvas("c1", "", 800, 600);
  double rmin = 0.0;
  double rmax = 6.0;
  double ymin = -0.5;
  double ymax = 1.0;
  double xmin = -350.;
  double xmax = 1050.;

  gR_12->SetTitle("Radial difference from sim (e = #pm 0.2mm); distance around ring [mm]; #Delta R [mm]");
  gR_12->GetXaxis()->CenterTitle();
  gR_12->GetYaxis()->CenterTitle();
  gR_12->GetYaxis()->SetRangeUser(rmin, rmax);
  gR_12->GetXaxis()->SetRangeUser(xmin, xmax);
  gR_12->SetLineColor(2);
  gR_12->SetMarkerColor(2);
  gR_12->SetMarkerStyle(8);

  gR_18->SetTitle("Radial difference from sim (e = #pm 0.2mm); distance around ring [mm]; #Delta R [mm]");
  gR_18->GetXaxis()->CenterTitle();
  gR_18->GetYaxis()->CenterTitle();
  gR_18->GetYaxis()->SetRangeUser(rmin, rmax);
  gR_18->GetXaxis()->SetRangeUser(xmin, xmax);
  gR_18->SetMarkerStyle(8);
  gR_18->SetLineColor(4);
  gR_18->SetMarkerColor(4);

  TF1* fR = new TF1("fR", "[0]*x + [1]", 0.0, distanceFromM00[7]);
  TFitResultPtr fitResult_r12 = gR_12->Fit("fR", "RS");
  TFitResultPtr fitResult_r18 = gR_18->Fit("fR", "RS");

  gR_12->GetFunction("fR")->SetLineColor(2);
  gR_18->GetFunction("fR")->SetLineColor(4);

  gR_12->Draw("AP");
  gR_18->Draw("P SAME");

  double t12 = gR_12->GetFunction("fR")->GetParameter(0);
  double e12 = gR_12->GetFunction("fR")->GetParError(0);

  double t12_d = atan(gR_12->GetFunction("fR")->GetParameter(0)) * (180 / TMath::Pi());
  double e12_d = (e12 / t12) * t12_d;

  double t18 = gR_18->GetFunction("fR")->GetParameter(0);
  double e18 = gR_18->GetFunction("fR")->GetParError(0);

  double t18_d = atan(gR_18->GetFunction("fR")->GetParameter(0)) * (180 / TMath::Pi());
  double e18_d = (e18 / t18) * t18_d;

  stringstream ss;
  ss << std::setprecision(3) << t12_d << " +/- " << e12_d;
  string theta12 = ss.str();

  ss.str("");
  ss << std::setprecision(3) << t18_d << " +/- " << e18_d;
  string theta18 = ss.str();

  TLegend* leg = new TLegend(0.4, 0.11, 0.89, 0.25);
  leg->SetBorderSize(0);
  leg->AddEntry(gR_12, Form("station 12, #theta = %s", theta12.c_str()), "l");
  leg->AddEntry(gR_18, Form("station 18, #theta = %s", theta18.c_str()), "l");
  leg->Draw("SAME");

  c1->SaveAs("SystematicPlots/radialFit.png");
  c1->Clear();

  ////NOW DO VERTICAL
  gY_12->SetTitle("Vertical difference from sim (e = #pm 0.2mm); distance around ring [mm]; #Delta y [mm]");
  gY_12->GetXaxis()->CenterTitle();
  gY_12->GetYaxis()->CenterTitle();
  gY_12->GetYaxis()->SetRangeUser(ymin, ymax);
  gY_12->GetXaxis()->SetRangeUser(xmin, xmax);
  gY_12->SetLineColor(2);
  gY_12->SetMarkerColor(2);
  gY_12->SetMarkerStyle(8);

  gY_18->SetTitle("Vertical difference from sim (e = #pm 0.2mm); distance around ring [mm]; #Delta y [mm]");
  gY_18->GetXaxis()->CenterTitle();
  gY_18->GetYaxis()->CenterTitle();
  gY_18->GetYaxis()->SetRangeUser(ymin, ymax);
  gY_18->GetXaxis()->SetRangeUser(xmin, xmax);
  gY_18->SetMarkerStyle(8);
  gY_18->SetLineColor(4);
  gY_18->SetMarkerColor(4);

  TF1* fY = new TF1("fY", "[0]*x + [1]", 0.0, distanceFromM00[7]);
  TFitResultPtr fitResult_y12 = gY_12->Fit("fY", "RS");
  TFitResultPtr fitResult_y18 = gY_18->Fit("fY", "RS");

  gY_12->GetFunction("fY")->SetLineColor(2);
  gY_18->GetFunction("fY")->SetLineColor(4);

  gY_12->Draw("AP");
  gY_18->Draw("P SAME");

  t12 = gY_12->GetFunction("fY")->GetParameter(0);
  e12 = gY_12->GetFunction("fY")->GetParError(0);

  t12_d = atan(gY_12->GetFunction("fY")->GetParameter(0)) * (180 / TMath::Pi());
  e12_d = (e12 / t12) * t12_d;

  t18 = gY_18->GetFunction("fY")->GetParameter(0);
  e18 = gY_18->GetFunction("fY")->GetParError(0);

  t18_d = atan(gY_18->GetFunction("fY")->GetParameter(0)) * (180 / TMath::Pi());
  e18_d = (e18 / t18) * t18_d;

  ss.str("");
  ss << std::setprecision(3) << t12_d << " +/- " << e12_d;
  theta12 = ss.str();

  ss.str("");
  ss << std::setprecision(3) << t18_d << " +/- " << e18_d;
  theta18 = ss.str();

  TLegend* leg2 = new TLegend(0.4, 0.11, 0.89, 0.25);
  leg2->SetBorderSize(0);
  leg2->AddEntry(gY_12, Form("station 12, #theta = %s", theta12.c_str()), "l");
  leg2->AddEntry(gY_18, Form("station 18, #theta = %s", theta18.c_str()), "l");
  leg2->Draw("SAME");

  c1->SaveAs("SystematicPlots/verticalFit.png");
  c1->Clear();

  ////Mahalanobis
  TGraph* graph[4] = {gR_12, gR_18, gY_12, gY_18};
  TFitResultPtr fitResult[4] = {fitResult_r12, fitResult_r18, fitResult_y12, fitResult_y18};
  string outname[4] = {"radialBands_s12", "radialBands_s18", "verticalBands_s12", "verticalBands_s18"};
  string legname[4];
  int col[8] = {3, kOrange + 7, 6, 7, 5, 2, 4, 1};
  legname[0] = "station 12, #theta_{r} = ";
  legname[1] = "station 18, #theta_{r} = ";
  legname[2] = "station 12, #theta_{y} = ";
  legname[3] = "station 18, #theta_{y} = ";

  for (int iFit(0); iFit < 4; iFit++) {

    c1->Clear();
    // Get correlation matrix
    TMatrixD corrMatrix = fitResult[iFit]->GetCorrelationMatrix();

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Pick 1 sigma points!
    /////////////////////////////////////////////////////////////////////////////////////////////////

    int nPars = 2;
    string gName = (iFit < 2) ? "fR" : "fY";

    // Mean values
    TVectorD meanVals(nPars);
    for (int n = 0; n < nPars; n++) {
      meanVals[n] = graph[iFit]->GetFunction(gName.c_str())->GetParameter(n);
    }

    // Covariance matrix
    TMatrixD covMatrix = fitResult[iFit]->GetCovarianceMatrix();
    TDecompChol decompCholCov(covMatrix);
    decompCholCov.Decompose();
    TMatrixD matrixCovI = decompCholCov.GetU();
    TMatrixD matrixCov(nPars, nPars);
    matrixCov.Transpose(matrixCovI);

    // Calculate Mahalanobis distance for one sigma CL
    double CL = 0.317311; // 1 - 68% (one sigma)
    //double CL = 0.05; // 1 - 68% (one sigma)
    double r2 = 0; // r^2 (see https://upload.wikimedia.org/wikipedia/commons/a/a2/Cumulative_function_n_dimensional_Gaussians_12.2013.pdf)
    while (TMath::Prob(r2, nPars) > CL) {
      r2 += 0.00001;
    }
    double r = sqrt(r2); // This is the Mahalanobis distance threshold under which CL % of points fall below
    cout << "Mahalanobis distance = " << r << endl;

    graph[iFit]->Draw("AP");

    TLegend* legTmp = new TLegend(0.3, 0.11, 0.87, 0.35);
    legTmp->SetBorderSize(0);
    int iLine = 0;

    vector<double> p0, p1;

    for (int i = -1; i < 2; i++) {
      for (int j = -1; j < 2; j++) {
        TVectorD u(nPars);
        u[0] = i;
        u[1] = j;

        double scale = u.Norm2Sqr() > 0 ? r / sqrt(u.Norm2Sqr()) : 1;
        for (int par = 0; par < nPars; par++) u[par] *= scale;

        // z = r*u*M + mean
        TVectorD z = matrixCov * u + meanVals;

        for (int par = 0; par < nPars; par++) cout << z[par] << " ";
        cout << endl;

        TF1* cmCombo = new TF1(Form("cmCombo_%d_%d", i, j), "[1] + [0]*x", 0.0, distanceFromM00[7]);
        for (int par = 0; par < nPars; par++) cmCombo->SetParameter(par, z[par]);

        cmCombo->SetLineWidth(3);
        cmCombo->SetLineStyle(2);
        cmCombo->SetLineColor(col[iLine]);
        cmCombo->Draw("SAME");

        double mtmp = cmCombo->GetParameter(0);
        double ctmp = cmCombo->GetParameter(1);

        double mtmp_d = atan(mtmp) * (180 / TMath::Pi());
        p0.push_back(mtmp_d);
        p1.push_back(ctmp);

        ss.str("");
        ss << std::showpos << i << j << " " << std::setprecision(3) << legname[iFit] << " " << mtmp_d << ", intercept = " << ctmp;
        legTmp->AddEntry(cmCombo, Form("%s", ss.str().c_str() ), "l");

        iLine++;
      }
    }
    legTmp->Draw("SAME");

    stringstream ss1;
    ss1 << "SystematicPlots/" << outname[iFit] << ".png";
    c1->SaveAs(ss1.str().c_str());

    TGraph* g = new TGraph(9, p0.data(), p1.data());
    g->SetTitle(Form("%s; #theta [degrees]; intercept", outname[iFit].c_str()));
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    g->SetMarkerColor(2);
    g->SetMarkerStyle(8);
    g->Draw("AP");
    stringstream ss2;
    ss2 << "SystematicPlots/" << "scan_" << outname[iFit] << ".png";
    c1->SaveAs(ss2.str().c_str());

  }
}
