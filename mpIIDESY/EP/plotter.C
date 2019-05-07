//Simple code to make plots from the Europa OmegaA ntuples
// gavin.hesketh@ucl.ac.uk
// modified by: samuel.grant.18@ucl.ac.uk

//to use a particular branch, make sure it is uncommented in Reader.C
//branch variables are listed in Reader.h

#define Plotter_C
#include "Plotter.h"
#include "TMath.h"
#include "TFile.h"


int CaloNum(int caloX, int caloY) {
  return caloX + 9 * caloY;
}

void Plotter::InitTrees(TString input_file) {
  //initialise the trees you want to read
  //then enable the relevant branches in Reader.C
  // ctt = new clusterTrackerTrackReader(input_file, "clusterTracker");
  // cr = new clusterTrackerClusterReader(input_file);
  // cl = new clusterReader(input_file);
  //  tr = new trackerReader(input_file);
  am = new allmuonsReader(input_file);
}

void Plotter::InitHistos() {


  plot2D("Ep_vs_t_early", 50, 0, 210, 200, 0, 2, "In Fill Time [us]", "E/p");
  plot2D("Ep_vs_t_all", 500, 0, 810, 2000, 0, 4, "In Fill Time [us]", "E/p");

  plot2D("E_vs_p_all", 500, 0, 3100, 500, 0, 3100, "P [MeV]", "E [MeV]");

  plot2D("E_vs_p_60us", 500, 0, 3100, 500, 0, 3100, "P [MeV]", "E [MeV]");
  
  plot2D("E_vs_p_60us_S12", 105, 1000, 3100, 105, 1000, 3100, "Track P/ 20 MeV", "Calorimeter E/ 20 MeV");
  plot2D("E_vs_p_60us_S18", 105, 1000, 3100, 105, 1000, 3100, "Track P/ 20 MeV", "Calorimeter E/ 20 MeV");

   plot2D("E_vs_p_120us", 500, 0, 3100, 500, 0, 3100, "P [MeV]", "E [MeV]");
  
  plot2D("E_vs_p_150us", 500, 0, 3100, 500, 0, 3100, "P [MeV]", "E [MeV]");

  plot2D("E_vs_p_positrons", 500, 0, 3100, 500, 0, 3100, "P [MeV]", "E [MeV]");




}

//=========================================================

//loop over the entries in the tree, making plots:

void Plotter::Run() {


  while ( NextallmuonsEvent() ) {

    //loop over the matches in this event:
    for (int i = 0; i < am->nmatches; i++) {

      const double t = am->cluTime[i]*1e-3; // ns -> us 

      const double p = sqrt(am->trkMomX[i] * am->trkMomX[i] + am->trkMomY[i] * am->trkMomY[i] + am->trkMomZ[i] * am->trkMomZ[i]); // MeV

      const double E = am->cluEne[i]; // MeV

      const double Ep = E / p;

      const int station =am->trkStationNum[i];
  
      //apply track quality to all
      if (am->trkPassTrackQuality[i] == false) continue; // 0/false=fail track quality 

      Fill2D("Ep_vs_t_early", t, Ep); // us 

      Fill2D("Ep_vs_t_all", t, Ep); // us 

      Fill2D("E_vs_p_all", p, E);

      //late in fill 
      if (t > 60) Fill2D("E_vs_p_60us", p, E);
      if (t > 120) Fill2D("E_vs_p_120us", p, E);
      if (t > 150) Fill2D("E_vs_p_150us", p, E);

      //per station and late in fill 
      if (t > 60 and station==12) Fill2D("E_vs_p_60us_S12", p, E);
      if (t > 60 and station==18) Fill2D("E_vs_p_60us_S18", p, E);

      // Positrons: -0.3 < Log(E/p) < 0.2
      if ( -0.3 < log(Ep) < 0.2 ) Fill2D("E_vs_p_positrons", p, E);


    }

  }


  return;

}
