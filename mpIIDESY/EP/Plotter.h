//Simple code to make plots from the Europa OmegaA ntuples
// gavin.hesketh@ucl.ac.uk

#ifndef Plotter_h
#define Plotter_h

#include "Reader.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

class Plotter {

public:
  Plotter(TString file, TString output_name="");
  ~Plotter();

  void Run();
  void SetOutputName(TString name) { output_name = name;}

  // private:
  //  bool NextTrEvent(){return tr->NextEvent();}
  bool NextallmuonsEvent(){return am->NextEvent();}
  bool NextClEvent(){return cr->NextEvent();}
  bool NextTrEvent(){return tr->NextEvent();}
  bool NextClTrClEvent(){return ctc->NextEvent();}
  bool NextClTrTrEvent(){return ctt->NextEvent();}
  void InitTrees(TString input_file);
  void InitHistos();
  void plot1D(TString name, int nbins, float xlow, float xhigh, TString xtitle="x axis", TString ytitle="", int col=2);  
  TH1* plot1DMaker(TString name, int nbins, float xlow, float xhigh, TString xtitle="x axis", TString ytitle="", int col=2);  
  TH2* plot2DMaker(TString name, int nbinsx, float xlow, float xhigh, int nbinsy, float ylow, float yhigh, TString xtitle="x axis", TString ytitle="", int col=2);  
  void Fill1D(TString name, double val, double weight=1);
  void plot2D(TString name, int xnbins, float xlow, float xhigh, int ynbins, float ylow, float yhigh, TString xtitle="x axis", TString ytitle="", int col=2);
  void Fill2D(TString name, double xval, double yval, double weight=1);

  void plot1DTr(TString name, int nbins, float xlow, float xhigh, TString xtitle="x axis", TString ytitle="", int col=2);  
  void Fill1DTr(TString name, double val, double weight=1);
  void plot2DTr(TString name, int xnbins, float xlow, float xhigh, int ynbins, float ylow, float yhigh, TString xtitle="x axis", TString ytitle="", int col=2);
  void Fill2DTr(TString name, double xval, double yval, double weight=1);

  double SetCaloY(int CaloNum, double Y_raw);

  
  std::vector<TH1*> plots1D;
  std::vector<TH2*> plots2D;
  TString output_name;
  double y_13[6], y_19[6];
  double height_13[6], height_19[6];

  void Extrapolate(double &x, double &y, double distance=64.9);
  double ExtrapolateToExit(const double x, const double y, bool xonly=true);

  bool FiducialXtal(const double x, const double y);
  bool FiducialMain(const double x, const double y);

  //could add any number of readers here:
  clusterTrackerTrackReader *ctt;
  clusterTrackerClusterReader *ctc;
  clusterReader *cr;
  trackerReader *tr;
  allmuonsReader *am;
};

#endif

//==============================================================================

#ifdef Plotter_C

Plotter::Plotter(TString input_file, TString output_file) : output_name(output_file) {

  InitTrees(input_file);
  
  InitHistos();

  //build the output name if not provided:
  if(output_name.Sizeof()<2) {
    //assumes files in the format eg /scratch/gm2/Omega/DATASET/gm2rootTrees_ana_RUN.root
    //dataset can be eg 60hr
    //   cout << "input name: " << input_file << "\n";
    //separate file name and path
    size_t found = input_file.Last('/');
    TString path = input_file(0,found);
    TString file = input_file(found+1, input_file.Sizeof()-(found+2) );

    //get dataset type from path
    found = path.Last('/');
    TString dataset = path(found+1, path.Sizeof()-(found+2) );
    // cout << "from " << path << " extracted dataset: " << dataset << "\n";
    
    //get run number from path
    found = file.Last('_');
    TString run = file(found+1,5);
    // cout << "from " << file << " extracted run num: " << run << "\n";
    
    output_name="plots_"+dataset+"_"+run+".root";
  }


  height_13[0] = 25.2;
  height_13[1] = 25.25;
  height_13[2] = 25.18;
  height_13[3] = 25.14;
  height_13[4] = 25.19;
  height_13[5] = 24.75;
  
  y_13[0] = -50.07 - 0.05 - height_13[0];
  y_13[1] = -24.77 - 0.05 - height_13[1];
  y_13[2] =  0.46 - 0.05 - height_13[2];
  y_13[3] = 25.65 - 0.05 - height_13[3];
  y_13[4] =  50.89 - 0.05 - height_13[4];
  y_13[5] =  75.69 - 0.05 - height_13[5];
  
  height_19[0] = 25.24;
  height_19[1] = 25.19;
  height_19[2] = 25.07;
  height_19[3] = 25.22;
  height_19[4] = 25.02;
  height_19[5] = 24.99;
  y_19[0] = -49.35 - 0.05 - height_19[0];
  y_19[1] = -24.11 - 0.05 - height_19[1];
  y_19[2] = 1.01 - 0.05 -  height_19[2];
  y_19[3] = 26.28 - 0.05 - height_19[3];
  y_19[4] = 51.53 - 0.05 - height_19[4];
  y_19[5] = 76.39 - 0.05 - height_19[5];
  

  
};

//----------------------------------

Plotter::~Plotter() {
  
  TFile * output_file = new TFile(output_name, "RECREATE");
  output_file->cd();
  
  for(int i=0; i<plots1D.size(); i++) {
    if(plots1D[i]->GetEntries()==0) cout <<"1D Plot not used: "<<plots1D[i]->GetName()<<endl;
    else plots1D[i]->Write();
  }
  for(int i=0; i<plots2D.size(); i++) {
    if(plots2D[i]->GetEntries()==0) cout <<"2D Plot not used: "<<plots2D[i]->GetName()<<endl;
    else plots2D[i]->Write();
  }
  
  output_file->Close();
  cout<<"Plots written to "<<output_file->GetName()<<endl;
  delete output_file;
  
  for(int i=0; i<plots1D.size(); i++) delete plots1D[i];
  for(int i=0; i<plots2D.size(); i++) delete plots2D[i];
  plots1D.clear();
  plots2D.clear();

}

//----------------------------------

TH1 * Plotter::plot1DMaker(TString name, int nbins, float xlow, float xhigh, TString xtitle, TString ytitle, int col){
  TH1 * p = new TH1D(name, name, nbins, xlow, xhigh);
  p->Sumw2();
  p->SetDirectory(0);
  p->SetXTitle(xtitle);
  p->SetYTitle(ytitle);
  p->GetXaxis()->SetLabelSize(0.04);
  p->GetYaxis()->SetLabelSize(0.04);
  p->GetXaxis()->SetTitleOffset(1.2);
  p->GetYaxis()->SetTitleOffset(1.2);
  p->GetXaxis()->SetLabelFont(42);
  p->GetYaxis()->SetLabelFont(42);
  p->GetXaxis()->SetTitleFont(42);
  p->GetYaxis()->SetTitleFont(42);
  p->SetTitle(0);
  p->SetLineWidth(2);
  p->SetLineColor(col);
  p->SetMarkerColor(col);
  p->SetMarkerSize(1.5);
  return p;
}

void Plotter::plot1D(TString name, int nbins, float xlow, float xhigh, TString xtitle, TString ytitle, int col){
  plots1D.push_back(plot1DMaker(name, nbins, xlow, xhigh, xtitle, ytitle, col));
}

void Plotter::plot1DTr(TString name, int nbins, float xlow, float xhigh, TString xtitle, TString ytitle, int col){
  plot1D("St12_"+name, nbins, xlow, xhigh, xtitle, ytitle, col);
  plot1D("St18_"+name, nbins, xlow, xhigh, xtitle, ytitle, col);
}


TH2 * Plotter::plot2DMaker(TString name, int nbinsx, float xlow, float xhigh, int nbinsy, float ylow, float yhigh, TString xtitle, TString ytitle, int col){
  TH2 * p = new TH2D(name, name, nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
  p->Sumw2();
  p->SetDirectory(0);
  p->SetXTitle(xtitle);
  p->SetYTitle(ytitle);
  p->GetXaxis()->SetLabelSize(0.04);
  p->GetYaxis()->SetLabelSize(0.04);
  p->GetXaxis()->SetTitleOffset(1.2);
  p->GetYaxis()->SetTitleOffset(1.2);
  p->GetXaxis()->SetLabelFont(42);
  p->GetYaxis()->SetLabelFont(42);
  p->GetXaxis()->SetTitleFont(42);
  p->GetYaxis()->SetTitleFont(42);
  p->SetTitle(0);
  p->SetLineWidth(2);
  p->SetLineColor(col);
  p->SetMarkerColor(col);
  p->SetMarkerSize(1.5);
  return p;
}

void Plotter::plot2D(TString name, int xnbins, float xlow, float xhigh, int ynbins, float ylow, float yhigh, TString xtitle, TString ytitle, int col){
  plots2D.push_back(plot2DMaker(name, xnbins, xlow, xhigh, ynbins, ylow, yhigh, xtitle, ytitle, col));
}

void Plotter::plot2DTr(TString name, int xnbins, float xlow, float xhigh, int ynbins, float ylow, float yhigh, TString xtitle, TString ytitle, int col){
  plot2D("St12_"+name, xnbins, xlow, xhigh, ynbins, ylow, yhigh, xtitle, ytitle, col);
  plot2D("St18_"+name, xnbins, xlow, xhigh, ynbins, ylow, yhigh, xtitle, ytitle, col);
}  

//----------------------------------

void Plotter::Fill1D(TString name, double val, double weight) {
  for(int i=0; i<plots1D.size(); i++){
    if(name.CompareTo(plots1D[i]->GetName())) continue;
    plots1D[i]->Fill(val, weight);
    return;
  }
  cout<<"Plot1D "<<name<<" not found!"<<endl;
}

void Plotter::Fill1DTr(TString name, double val, double weight) {
  Fill1D( "St"+std::to_string(ctt->station)+"_"+name, val, weight);
}

void Plotter::Fill2D(TString name, double xval, double yval, double weight) {
  for(int i=0; i<plots2D.size(); i++){
    if(name.CompareTo(plots2D[i]->GetName())) continue;
    plots2D[i]->Fill(xval, yval, weight);
    return;
  }
  cout<<"Plot2D "<<name<<" not found!"<<endl;
}

void Plotter::Fill2DTr(TString name, double xval, double yval, double weight) {
  Fill2D( "St"+std::to_string(ctt->station)+"_"+name, xval, yval, weight);
}



//=========================================================


double Plotter::SetCaloY(int CaloNum, double Y_raw) {
  double y_int=0;
  double y_frac=modf(Y_raw, &y_int);

  if( CaloNum==13 ) return  y_13[(int)y_int] + y_frac * height_13[(int)y_int];
  return  y_19[(int)y_int] + y_frac * height_19[(int)y_int];

}


//=========================================================


void Plotter::Extrapolate(double &x, double &y, double dist) {
  double trans_dist = dist * TMath::Sin(ctt->trackMomentumTheta);

  x += trans_dist * TMath::Cos(ctt->trackMomentumPhi);
  y += trans_dist * TMath::Sin(ctt->trackMomentumPhi);
  return;
}

double Plotter::ExtrapolateToExit(const double x, const double y, bool xonly) {

  //extrapolate until x= +- 112, or y=+-75

  // x = gradient * z
  const double gradx = TMath::Sin(ctt->trackMomentumTheta) * TMath::Cos(ctt->trackMomentumPhi);
  const double grady = TMath::Sin(ctt->trackMomentumTheta) * TMath::Sin(ctt->trackMomentumPhi);

  double zx = 0, zy=0;
  if(gradx>0) { //extrapolate to x = +112
    zx = (112.5-x) / gradx;
  }
  else zx = (-112.5-x) / gradx;

  if(grady>0) { //extrapolate to y = +75
    zy = (75-y) / grady;
  }
  else zy = (-75-y) / grady;

  double z = zx;
  if(!xonly)  z = zx<zy ? zx : zy;
  if (z>140) return 70;
  else return z-70;

}

//=========================================================


bool Plotter::FiducialXtal(const double x, const double y) {

  //crystal is 25 x 25
  //fiducial centre is 0.707 * 25

  const double fc = 0.707 * 25; //the fiducial centre
  const double step = 25;
  
  double fx = fabs(x) ; //start at low edge of fiducial volume
  double fy = fabs(y);
  bool fidx=false, fidy=false;

  if( fx < fc/2) fidx=true;
  fx-=step-fc/2;
  while (fx >0){
    if( fx < fc) fidx=true;
    fx-=step;
  }

  fy-= 0.5*(step-fc);
  while (fy >0){
    if( fy < fc) fidy=true;
    fy-=step;
  }

  
  
  return fidx && fidy;
}
// Different definition of fiducial, take the central crystals including edges
// Take inner region of 5x4 crystals
bool Plotter::FiducialMain(const double x, const double y) {

  bool fid = false;
  bool fidy = false;
  bool fidx = false;
 
  //  const double fx = x;
  // const double fy = y;
  // const double xCentre = (9*25)/2;
  //const double yCentre = (6*25)/2;
  const double xLo = -62.5;//-80;
  const double xHi = 12.5;
  const double yLo = -25;
  const double yHi = 25;
  // const double yBound = 50;

  if( yLo < y && y < yHi) fidy = true;
  if( xLo < x && x < xHi) fidx = true;

  if(fidy && fidx) {   
    fid = true;
  }
   
  return fid;
}




//=========================================================


// Can call the Plotter with a list of files to read.
// Defaults to a single hard-coded file
// Will make a separate output file for each input file.

int main(int argc, char *argv[]){

  vector<TString> fileNames;
  
  if (argc < 2){
    cout<<"Using hard-coded file list"<<endl;
    fileNames.push_back("../../60hr_v9_17_01/gm2rootTrees_ana_15922.root");
    //fileNames.push_back("/scratch/folders/junk/gm2rootTrees_ana.root");
  }
  
  for (int i = 1; i < argc; ++i){
    fileNames.push_back(argv[i]);
  }
  cout << "Reading a total of " << fileNames.size() << " files \n";
  for (int i = 0; i < fileNames.size(); ++i){
    //write out individual file for each input file
    Plotter *p = new Plotter(fileNames[i]);
    p->Run();
    delete p;
  }
  
  return 1;
}
  


#endif
