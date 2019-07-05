#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "TList.h"
#include "TObject.h"
#include "TColor.h"
#include "TPaveStats.h"
#include "TKey.h"
#include "TROOT.h"
#include "TVirtualFFT.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <sstream>
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TClass.h"
#include "TLegendEntry.h"
#include "TGaxis.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TDirectoryFile.h"
//#include "boost/format.hpp"
#include <exception>
#include <cstdlib>
#include <iomanip>

using std::string;
using std::stringstream;
using std::vector;
using std::map;

TH2F* moduloPlot(TH2F* hist, double period, double xmin, double xmax) {

  int nbinsx = hist->GetXaxis()->GetNbins();
  int nbinsy = hist->GetYaxis()->GetNbins();
  double maxy = hist->GetYaxis()->GetBinUpEdge(nbinsy);
  double miny = hist->GetYaxis()->GetBinLowEdge(1);
  
  double binw = hist->GetXaxis()->GetBinWidth(1);
  double t_tot = nbinsx * binw;
  double frac = t_tot / period;
  int n_osc_tot = int(frac / binw);
  int nbins_mod = int(period / binw);
  
  int n_osc(1);
  int newBin = 1;
  
  TH2F* h2mod = new TH2F("h2mod","",nbins_mod,0,nbins_mod*binw,nbinsy,miny,maxy);

  for (int ibinx(1); ibinx < nbinsx; ibinx++) {
    for (int ibiny(1); ibiny < nbinsy; ibiny++) {
      int binc = hist->GetBinContent(ibinx,ibiny);
      double x = hist->GetXaxis()->GetBinCenter(ibinx);
      double y = hist->GetYaxis()->GetBinCenter(ibiny);

      if (x < xmin || x > xmax) continue;

      int current_binc = h2mod->GetBinContent(newBin,ibiny);
      h2mod->SetBinContent(newBin,ibiny,current_binc+binc);
    }
    if (ibinx % (nbins_mod * n_osc) == 0) {
      newBin = 0;
      n_osc++;
    }
    newBin++;
  }
  return h2mod;
}

//TGraph must have equally spaced points
TH1F* GraphToHist(TGraphErrors* g){
  
  double n = g->GetN();
  //if (n < 2) //exception;  
  double* xArray = g->GetX();
  double binSize = xArray[1] - xArray[0];
  TH1F* h = new TH1F(g->GetName(),"",n, xArray[0] - binSize/2.0, xArray[0]+ n*binSize);

  for (int i(0); i < n; i++){
    double x(0.0), y(0.0);
    g->GetPoint(i, x, y);
    double yErr = g->GetErrorY(i);
    
    double binx = h->FindBin(x);
    h->SetBinContent(binx, y);
    h->SetBinError(binx, yErr);
  }
  return h;
}

TH1F* moduloPlot(TH1* hist, double period, double xmin, double xmax) {
  
  int nbins = hist->GetXaxis()->GetNbins();
  double binw = hist->GetXaxis()->GetBinWidth(1);
  double t_tot = nbins * binw;
  double frac = t_tot / period;
  int n_osc_tot = int(frac / binw);
  int nbins_mod = int(period / binw);
  
  TH1F* hmod = new TH1F("hmod","",nbins_mod,0,period);
  int n_osc(1);
  int newBin = 1;
  
  for (int ibin(1); ibin < nbins; ibin++) {
    int binc = hist->GetBinContent(ibin);
    double x = hist->GetXaxis()->GetBinCenter(ibin);

    if (x < xmin || x > xmax) {
      continue;
    }
    
    int current_binc = hmod->GetBinContent(newBin);
    hmod->SetBinContent(newBin,current_binc+binc);
    //cout << "ibin: " << ibin << " x: " << x << " has content: " << binc << "\n";

    if (ibin % (nbins_mod * n_osc) == 0) {
      newBin = 0;
      n_osc++;
    }
    newBin++;    
  }
  return hmod;
}

void radial(){

  gStyle->SetOptStat(0);
  string outFolderName = "Radial_60hr";
  
  //60 hr with track quality
  TFile* f = new TFile("trackRecoPlots_15921_15969.root");
  double fieldIndex = 0.108; 

  //uncomment for 9 day
  //TFile* f  = new TFile("/gm2/data/g2be/Production/Plots/BLAH.root");
  //fieldIndex = 0.1215;
  //outFolderName = "VerticalCBO_9day"

  // Run2
  //outFolderName = "VerticalCBO_Run2";
  //TFile* f  = new TFile("/gm2/data/g2be/Production/Plots/v9_17_01/Run2/trackRecoPlots_23745_23756.root");

  string app = ".png";

  TDirectoryFile* extrap    = (TDirectoryFile*)f->Get("Extrapolation");
  TDirectoryFile* vertices  = (TDirectoryFile*)extrap->Get("vertices");
  TDirectoryFile* station12 = (TDirectoryFile*)vertices->Get("station12");
  TDirectoryFile* station18 = (TDirectoryFile*)vertices->Get("station18");

  map<string, TH1*> replot;
  map<string, TGraph*> replotG;
  double average = 0.0;
  double averageV = 0.0;

  for (auto station : {station12, station18}){
    //  for (auto station : {station18}){

    string stationName = station->GetName();
    cout << "looping over " << stationName << "\n";
    
    TDirectoryFile* pVal = (TDirectoryFile*)station;
    
    TCanvas* c1 = new TCanvas("c1","c1", 2000, 2000);
    c1->SetCanvasSize(1200, 850);
    c1->SetLeftMargin(0.13);
    c1->SetBottomMargin(0.13);
    c1->SetFixedAspectRatio();
    
    //BeamSpot
    TH2F* beam = (TH2F*)pVal->Get("h_vertexPosSpread");
    beam->GetXaxis()->SetRangeUser(-45.0, 45.0);
    beam->GetYaxis()->SetRangeUser(-45.0, 45.0);
    beam->Draw("COLZ");
    beam->GetXaxis()->SetTitle("Radial Beam Position [mm]");
    beam->GetYaxis()->SetTitle("Vertical Beam Position [mm]");
    beam->GetYaxis()->SetTitleOffset(1.0);
    beam->GetXaxis()->SetTitleOffset(1.0);
    beam->GetYaxis()->CenterTitle();
    beam->GetXaxis()->CenterTitle();
    beam->SetTitle("");
    //c1->SaveAs(Form("%s/beam_%s%s",outFolderName.c_str(), stationName.c_str(), app.c_str()));
    c1->Print("Spot.png");

    // //rootManager_->Add(subDir, new TH2F("h_radialPos_vs_time_fine"  ,"Radial pos vs mean time;Track time [us];Radial pos [mm]",20000,0,200,140,-70,70));
    // TH2F* h1 = (TH2F*)pVal->Get("h_radialPos_vs_time");
    // h1->GetXaxis()->SetRangeUser(0.0, 600.0);
    // h1->GetYaxis()->SetRangeUser(-60, 60.0);
    // h1->Draw("COLZ");
    // c1->SaveAs(Form("%s/h1_%s%s",outFolderName.c_str(), stationName.c_str(), app.c_str()));

    // h1->SetTitle(";Track time [us]; Average Radial Position [mm]");
    // h1->GetXaxis()->SetRangeUser(30.0, 100.0);
    // h1->GetYaxis()->SetRangeUser(-60, 60.0);
    // h1->Draw("COLZ");
    // c1->SaveAs(Form("%s/h1_30-100_%s%s",outFolderName.c_str(), stationName.c_str(), app.c_str()));

    // TH2F* h2 = (TH2F*)pVal->Get("h_verticalPos_vs_time_fine");
    // h2->RebinX(10);
    // h2->GetXaxis()->SetRangeUser(30.0, 50.0);
    // h2->GetYaxis()->SetRangeUser(-60, 60.0);
    // h2->SetTitle(";Track time [us]; Average Vertical Position [mm]");
    // h2->Draw("COLZ");
    // c1->SaveAs(Form("%s/hy_30-100_%s%s",outFolderName.c_str(), stationName.c_str(), app.c_str()));

    ////uncomment for shift
    //double cboPeriod = 3.8;
    //double shift = (stationName == "station12")? 1 * cboPeriod / 2 : 3 * cboPeriod /4;
    //cout << "SHIFT: " << shift << " for " << stationName << "\n";
    //TH2F* h1_shifted = new TH2F("shifted", ";adjusted time [us]; radial position [mm]", h1->GetXaxis()->GetNbins(), 0., 200., h1->GetYaxis()->GetNbins(), -70, 70);
    //
    //for (int ix(1); ix < h1->GetXaxis()->GetNbins(); ix++){
    //  double t1 = h1->GetXaxis()->GetBinCenter(ix);
    //  double t2 = t1 - shift;
    //  if (t2 > 200.0) continue;
    //  int bx = h1_shifted->GetXaxis()->FindBin(t2);
    //  //cout << bx << " new time for shifted: " << t2 << " original time: " << t1 << " from bin " << ix << "\n";
    //  for (int iy(1); iy < h1->GetYaxis()->GetNbins(); iy++){
    //	double y = h1->GetYaxis()->GetBinCenter(iy);
    //	int by = h1_shifted->GetYaxis()->FindBin(y);
    //	h1_shifted->SetBinContent(bx,by,h1->GetBinContent(ix,iy));
    //  }
    //}
    
 //    //TProfile it and look for 
 //    TProfile* tp = (TProfile*)h1->ProfileX("_pfx");
 //    tp->Draw();
 //    c1->SaveAs(Form("%s/tp_%s%s",outFolderName.c_str(), stationName.c_str(), app.c_str()));

 //    //TProfile with time cut
 //    int binLow = h1->GetXaxis()->FindBin(30.0);
 //    int binHigh = h1->GetXaxis()->FindBin(600.0);
 //    tp = (TProfile*)h1->ProfileX("_pfx", binLow, binHigh);
 //    tp->Draw();
 //    c1->SaveAs(Form("%s/tp_after30_%s%s",outFolderName.c_str(), stationName.c_str(), app.c_str()));
    
 //    TH1F* proj = (TH1F*)h1->ProjectionY("_pfy");
 //    proj->GetXaxis()->SetRangeUser(-60.0, 60.0);
 //    proj->Draw();
 //    c1->SaveAs(Form("%s/proj_%s%s",outFolderName.c_str(), stationName.c_str(), app.c_str()));

 //    proj = (TH1F*)h1->ProjectionY("_pfy", binLow, binHigh);
 //    proj->GetXaxis()->SetRangeUser(-60.0, 60.0);
 //    proj->Draw();
 //    c1->SaveAs(Form("%s/proj_after30_%s%s",outFolderName.c_str(), stationName.c_str(), app.c_str()));


 //    //fix for all
 //    tp->GetYaxis()->SetRangeUser(-25, 35);

 //    //make 0.2us steps
 //    vector<double> windows;
 //    for (int i(0); i < -1; i++) {
 //      double t0 = i * 0.2;
 //      windows.push_back(t0);
 //    }

 //    int iplot = 0;
 //    for (auto& ti: windows){
 //      double tf = ti + 5.0;
 //      tp->GetXaxis()->SetRangeUser(ti, tf);
 //      tp->Draw();
 //      c1->SaveAs(Form("%s/tpZoomed_%.4i_%.4f-%.4f_%s%s",outFolderName.c_str(), iplot, ti, tf, stationName.c_str(), app.c_str()));

 //      if (replot[Form("tpZoomed_%.4f", ti)]){
	// replot[Form("tpZoomed_%.4f", ti)]->SetLineColor(4);
	// replot[Form("tpZoomed_%.4f", ti)]->Draw();
	// tp->SetLineColor(2);
	// tp->Draw("SAME");
	
	// TLegend* legend = new TLegend(0.65,0.79,0.85,0.89);
	// legend->SetBorderSize(0);
	// legend->AddEntry(replot[Form("tpZoomed_%.4f", ti)],"station 12","l");
	// legend->AddEntry(tp,"station 18","l");
	// legend->Draw("SAME");
	
	// //c1->SaveAs(Form("%s/tpZoomed_%.4i_%.4f-%.4f_both%s",outFolderName.c_str(), iplot, ti, tf, app.c_str()));
	// c1->SaveAs(Form("%s/tpZoomed_%.4i_both%s",outFolderName.c_str(), iplot, app.c_str()));
 //      }
 //      else replot[Form("tpZoomed_%.4f", ti)] = (TH1D*)tp->Clone();
 //      iplot++;
 //    }


  // }
  }
}
