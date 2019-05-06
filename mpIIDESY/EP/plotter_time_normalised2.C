//Simple code to make plots from the Europa OmegaA ntuples
// gavin.hesketh@ucl.ac.uk
// modified by: samuel.grant.18@ucl.ac.uk

//to use a particular branch, make sure it is uncommented in Reader.C
//branch variables are listed in Reader.h


//Just going with one region



#define Plotter_C
#include "Plotter.h"
#include "TMath.h"
#include "TFile.h"


int CaloNum(int caloX, int caloY) {
  return caloX+9*caloY;
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

  const double ymin = 0.5;
  const double ymax = 1.5;
 
  for (int stn = 13; stn < 20 ; stn = stn + 6) {
    for (int cut = 0; cut < 4; cut++) {
      for (int brd = 1; brd < 3; brd++) {
	  plot2D("St"+std::to_string(stn)+"_Ep_vs_t_early_"+std::to_string(cut)+"_"+std::to_string(brd),50,0,4200*50,200,ymin,ymax,"In Fill Time [ns]", "E/p");
	  plot2D("St"+std::to_string(stn)+"_xy_calo_"+std::to_string(cut)+"_"+std::to_string(brd),500,-150,150,500,-120,120,"Track position at calo face X [mm]", "Track position at calo face Y [mm]");
	  plot1D("St"+std::to_string(stn)+"_efrac_"+std::to_string(cut)+"_"+std::to_string(brd),101,-0.005,1.005,"Cluster Energy in Crystal / Cluster Energy","N");  

      }
    }
  }

}

//=========================================================

//loop over the entries in the tree, making plots:

void Plotter::Run() {

 TFile *input1 = TFile::Open("../makePlots2/fitted_Ep_xtal.root");  
     
  //loop over the clusterTracker/tracker tree:

 while( NextallmuonsEvent() ) {

    //loop over the matches in this event:
    for(int i=0; i<am->nmatches; i++) {
      //  cout<<i<<endl;
      // require quality cut pass
      if(am->trkPassTrackQuality[i] == false) continue;
      
      double p = sqrt(am->trkMomX[i]*am->trkMomX[i] + am->trkMomY[i]*am->trkMomY[i] + am->trkMomZ[i]*am->trkMomZ[i]);
      
      const double logEop = log(am->EovP[i]);
      const double dt = am->Tdiff[i];
 
      const int caloSt = am->cluCaloNum[i];
      if(caloSt==13){
	if (dt<-8 || dt>3 ) continue;
      }
      else {
	if(dt<-9 || dt>1) continue;        
      }
      
      //calo range from 0->8 in x (needs flipping), 0->5 in y
      //tracker range from -120 -> 100 in x, -50 to 50 in y
      // x scaling ~9/220=0.041; y scaling ~ 6/100 = 0.06
      
      const double caloX_raw = am->cluX[i];
      const double caloY_raw = am->cluY[i];
      const int xtal = CaloNum(caloX_raw, caloY_raw);
      //Convert to mm from calo number I believe
      const double caloX = 112.5 - 25*(caloX_raw);
      const double caloY = SetCaloY(caloSt, caloY_raw);
      // Tracker decay vertices
      double trX = am->vX[i];
      double trY = am->vY[i];

      const double dX = caloX - trX;
      const double dY = caloY - trY;
      const double dR = sqrt(dX*dX + dY*dY);
      
      if(dR>30) continue;

      double t = (am -> decayTime[i]);
      if (t < 4200) continue;

      double E = am->cluEne[i];
      double Ep = E/p;

      //========================== SNIP!
  
      bool region[3] = {false};
      // Positrons 
      if(logEop>-0.3 && logEop<0.2 ) region[0]=true;
      // High Flux - Energy - All Tracsk
      if (1200 < E && E < 2400) region[1]=true;
      // High Flux - Energy - Positrons
      if(region[0]==true && region[1]==true) region[2]=true;   
      // Select high flux positrons
      if(!region[2]) continue;

      ////////////////////////////////////////////
      //Energy Fraction Cut

      const double efrac = am->efracmaxclu[i];
      //  if(efrac < 0.605) continue;
      
      ////////////////////////////////////////////
      int shortLifeXtal[22] = {0,9,10,11,14,15,18,19,20,23,24,27,30,31,34,35,36,39,40,43,44,45};
      // cout<<shortLifeXtal[0]<<endl;
      int brd;
      for (int j = 0 ; j < 22 ; j++ ) {
	//	cout<<shortLifeXtal[i]<<" "<<xtal<<endl;
	if (xtal == shortLifeXtal[j]) {
	  brd = 1;
	  break;
	  
	}
	else {
	  brd = 2;
	}
      }

      /////////////////////////
      //Normalisation
      std::string h1 = "St"+std::to_string(caloSt)+"_fit_Ep_vs_xtal_0";
       
      TH1D *scale1 = (TH1D*)input1->Get(h1.c_str());
      if(scale1 == 0) continue;
     
      double sf1 = scale1->GetBinContent(xtal+1);
       
      double Ep1 = Ep * (1 / sf1 );
      
      Fill2D("St"+std::to_string(caloSt)+"_Ep_vs_t_early_0_"+std::to_string(brd),t,Ep1);
      Fill2D("St"+std::to_string(caloSt)+"_xy_calo_0_"+std::to_string(brd),trX,trY);
      Fill1D("St"+std::to_string(caloSt)+"_efrac_0_"+std::to_string(brd),efrac);
	
      if (efrac > .605) {
	Fill2D("St"+std::to_string(caloSt)+"_Ep_vs_t_early_1_"+std::to_string(brd),t,Ep1);
	Fill2D("St"+std::to_string(caloSt)+"_xy_calo_1_"+std::to_string(brd),trX,trY);
	Fill1D("St"+std::to_string(caloSt)+"_efrac_1_"+std::to_string(brd),efrac);
      }

      if (efrac > .755) {
	Fill2D("St"+std::to_string(caloSt)+"_Ep_vs_t_early_2_"+std::to_string(brd),t,Ep1);
	Fill2D("St"+std::to_string(caloSt)+"_xy_calo_2_"+std::to_string(brd),caloX,caloY);
	Fill1D("St"+std::to_string(caloSt)+"_efrac_2_"+std::to_string(brd),efrac);
      }

      if (efrac > .995) {
	
	Fill2D("St"+std::to_string(caloSt)+"_Ep_vs_t_early_3_"+std::to_string(brd),t,Ep1);
	Fill2D("St"+std::to_string(caloSt)+"_xy_calo_3_"+std::to_string(brd),caloX,caloY);
	Fill1D("St"+std::to_string(caloSt)+"_efrac_3_"+std::to_string(brd),efrac);
      }

     
    }
    
 }
  
 input1 -> Close();
  
 return;
  
}
