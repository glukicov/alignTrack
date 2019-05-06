// Author: Joe
// Modified by Gleb  

#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include <iostream>
#include <vector>
#include <string>
using std::cout;

double momDiff(double curvature, double mom, bool Plot = false);

int main(int argc, char** argv){

  cout << "argc: " << argc << "\n";
  
  if (argc < 3){
    cout << "arguments: " << argc << "\n";
    cout << "usage: ./tmp.exe <curvature> <mom>\n"; 
    return 1;
  }

  double curvature = atof(argv[1]);
  double mom = atof(argv[2]);

  momDiff(curvature, mom, true);

  TCanvas* c2 = new TCanvas("c2","", 800, 600);
  std::string opt = "AP";
  int icol = 2;
  std::vector<double> curvatures = {1.0E-6, 0.5E-6, 1.0E-7};
  //std::vector<TGraph*> graphs;

  TLegend* leg = new TLegend(0.11, 0.68, 0.49, 0.89);
  leg->SetBorderSize(0);

  for (auto& curve: curvatures){

    std::vector<double> pTrue = {};
    std::vector<double> pDiff = {};
    double testMom = 1.0;
    int n = 0;
    while (testMom < 3.2){
      pTrue.push_back(testMom);
      pDiff.push_back( momDiff(curve, testMom, false) * 1000.0);
      testMom += 0.1; 
      n++;
    }

    int j = 0;
    for (auto& p: pDiff ) {
      cout << "pTrue: " << pTrue.at(j) << " p: " << p << "\n";
      j++;
    }
    TGraph* g = new TGraph(n, pTrue.data(), pDiff.data());

    c2->cd();
    //c2->GetYaxis()->SetNdivisions(555);
    // g->SetGrid(1, 1); g->Update();
    g->SetTitle("");
    g->GetXaxis()->SetTitle("True P [GeV]");
    g->GetYaxis()->SetTitle("#Delta P  [MeV]");
    g->GetYaxis()->CenterTitle();
    g->GetXaxis()->CenterTitle();
    g->SetMarkerStyle(8);
    g->SetMarkerColor(icol);

    g->GetYaxis()->SetRangeUser(2E-1, 200.0);
    g->Draw(opt.c_str());
    leg->AddEntry(g, Form("Curvature: %.2f E-6", curve/1.0E-6), "P");
    opt = "SAME P";
    icol++;
  }

  leg->Draw("SAME");
  c2->SetGrid();
  c2->SaveAs("curvatureEffect.png");
  c2->SetLogy();
  c2->SaveAs("curvatureEffect_log.png");

  return 0;
}

double momDiff(double curvature, double mom, bool Plot){

  TCanvas* c1 = new TCanvas(Form("c1%.6f",mom/curvature),"", 800, 600);

  //for momentum 1GeV = 5.36E-19 kgm/s
  double c = 2.99792458 * 1E8; //m/s
  double q = 1.60217 * 1E-19; 
  double GeV2SI = q * 1E9 / c;
  double p_SI = mom *  GeV2SI;
  double B = 1.45147; //Tesla
  double q_SI = 1.602E-19;
  double radius = p_SI / (q_SI * B);

  radius *= 1000.0; //to get in mm

  // Gleb's parameterisation
  TF1* f1 = new TF1("f1", "[0]*(x-[1])^2 + [2]", 0, 1000);
  f1->SetParameter(0, curvature);
  f1->SetParameter(1, 956.0/2.0);
  f1->SetParameter(2, 0.0);

  f1->SetTitle("; module position [mm]; change in R [mm]");

  f1->GetYaxis()->SetRangeUser(-1.0, 1.0);
  f1->Draw();
  if (Plot)c1->SaveAs("tmp.png");
  
  double modDiff = 956.0 / 8.0;
  //plot points on a cirlce
  double xmod[8] = {1.0, 2*modDiff, 3*modDiff, 4*modDiff, 5*modDiff, 6*modDiff, 7*modDiff, 8*modDiff}; 
  double xpos[32] = {0.0};
  double ypos[32] = {0.0};
  double ypos2[32] = {0.0};//, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double ydiff[32] = {0.0};
  double xerr[32] = {0.0};
  double yerr[32] = {0.0};

  for (int i(0); i < 8; i++){
    xpos[i*4+0] = xmod[i];
    xpos[i*4+1] = xmod[i]+6.0;
    xpos[i*4+2] = xmod[i]+12.0;
    xpos[i*4+3] = xmod[i]+18.0;
  }

  for (int i(0); i < 32; i++){
    ypos[i] = sqrt(radius*radius - xpos[i]*xpos[i]);
    ypos2[i] = ypos[i] + f1->Eval(xpos[i]); 
    ydiff[i] = f1->Eval(xpos[i]); 
    xerr[i] = 0.001;
    yerr[i] = 0.1;
  }
  
  TF1* circ = new TF1("circ", "sqrt([0]^2 - (x - [1])^2) + [2]", 0.5, 960);
  circ->SetParameter(0, radius);
  circ->SetParameter(1, 0.0);
  circ->SetParameter(2, 0.0);

  TGraphErrors* g = new TGraphErrors(32, xpos, ypos, xerr, yerr);
  g->SetTitle(";module x position [mm]; raduis of centre of module [mm]");
  g->SetMarkerStyle(8);
  g->SetMarkerColor(2);
  g->Fit("circ", "RQME");  
  g->Draw("AP");

  TF1* circ2 = new TF1("circ2", "sqrt([0]^2 - (x - [1])^2) + [2]", 0.5, 960);
  circ2->SetParameter(0, radius);
  circ2->SetParameter(1, 0.0);
  circ2->SetParameter(2, 0.0);

  TGraphErrors* g2 = new TGraphErrors(32, xpos, ypos2, xerr, yerr);
  g2->SetMarkerStyle(8);
  g2->SetMarkerColor(3);
  g2->Fit("circ2", "RQME");
  g2->GetFunction("circ2")->SetLineColor(3);
  g2->Draw("P SAME");

  double pFit1 = g->GetFunction("circ")->GetParameter(0)*(q_SI * B) / (1000.0 * GeV2SI);
  double pFit2 = g2->GetFunction("circ2")->GetParameter(0)*(q_SI * B) / (1000.0 * GeV2SI);
  if (Plot){
    cout << "r1: " << g->GetFunction("circ")->GetParameter(0) << " mm, gives momentum: " << pFit1 << " GeV\n";
    cout << "r2: " << g2->GetFunction("circ2")->GetParameter(0) << " mm, gives momentum: " << pFit2 << " GeV\n";
  }

  TLegend* leg = new TLegend(0.5, 0.7, 0.89, 0.89);
  leg->SetBorderSize(0);
  leg->AddEntry(g->GetFunction("circ"), Form("no misalignment, p = %.3f GeV", pFit1), "l");
  leg->AddEntry(g2->GetFunction("circ2"), Form("misalignment, p = %.3f GeV", pFit2), "l");
  leg->Draw("SAME");
  if (Plot) {
    c1->SaveAs("circle.png");
    cout << "when momentum = " << mom << " GeV, radius of curvature = " << radius << " mm \n";
    cout << "fit gives mom: " << pFit1 << " GeV, fit with added curvature of " << curvature << " gives mom: " << pFit2 << " GeV, diff: " << pFit2 - pFit1 << " GeV \n"; 
  }

  //TGraph* gdiff = new TGraph(8, xpos, ydiff);
  //gdiff->SetMarkerStyle(8);
  //gdiff->SetMarkerColor(3);
  //
  ////an ellipse where yc = a
  //TF1* ellipse = new TF1("ellipse", "-1 * sqrt([0]^2 * ( 1 -  ( (x - [1])^2) / [2]^2 ) ) + [0]", 0, 900);
  //ellipse->SetParameter(0, 5.0); //a
  //ellipse->SetParameter(1, 500.0); //xc
  //ellipse->SetParameter(2, 500.0); //b 
  //
  //gdiff->Fit("ellipse");
  //gdiff->Draw("AP");
  //c1->SaveAs("elipseDiff.png");
  ////ellipse->Draw();
  //
  //TF1* ellipse2 = new TF1("ellipse2", "1 * sqrt([0]^2 * ( 1 -  ( (x - [1])^2) / [2]^2 ) ) + [3]", 0, 900);
  //ellipse2->SetParameter(0, 7100.0); //a
  //ellipse2->SetParameter(1, 0.0); //xc
  //ellipse2->SetParameter(2, 7100.0); //b 
  //ellipse2->SetParameter(3, 0.0); //yc 
  //g2->Fit("ellipse2");  
  //g2->Draw("AP");
  //ellipse2->Draw("SAME");
  ////  ellipse2->Draw("SAME");
  //
  //c1->SaveAs("elipseFit.png");

  return pFit2 - mom;
  //return pFit2 - pFit1;
}
