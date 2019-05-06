//Simple code to read the Europa OmegaA ntuples
// gavin.hesketh@ucl.ac.uk
//
// You shouldn't have to edit this, except to declare new branches
//    or a reader for a new kind of tree.


#ifndef Reader_h
#define Reader_h

#include "TString.h"
#include <iostream>
#include <vector>
#include "TMatrixD.h"

class TBranch;
class TTree;


//====================================================================
//=== Generic tree reader
//====================================================================
//this has some general functions to open the file, load next event, etc
class Reader {

 public: 
  Reader(TString filename, TString folder, TString treename);
  ~Reader(){};
  bool NextEvent();
  long Entries() {return nentries_;}

 protected:
  template<class VAR>
    void LoadBranch(TString name, VAR &var, TBranch *&branch);
  template<class VAR>
    void LoadBranchVector(TString name, VAR &var, TBranch *&branch);
  virtual void Init(){};
  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  int           fCurrent; //!current Tree number in a TChain
  long nentries_;
  long jentry_;
};



//====================================================================
//=== clusterTracker Track reader
//====================================================================
//contains the branches specific to the clusterTracker track tree

class clusterTrackerTrackReader : public Reader {

 public :
 clusterTrackerTrackReader(TString filename, TString treename="clusterTracker"):
  Reader(filename, treename, "tracker"){
    Init();
    std::cout<<"==================================="<<std::endl<<std::endl;
    
  }


   // Declaration of leaf types
   Int_t           runNum;
   Int_t           subRunNum;
   Int_t           eventNum;
   Int_t           bunchNum;
   Int_t           clusterFillNum;
   Int_t           clusterCaloNum;
   Int_t           clusterIslandNum;
   Double_t        clusterX;
   Double_t        clusterY;
   Int_t           clusterHits;
   Double_t        clusterTime;
   Double_t        clusterE;
   Double_t        trackMomentum;
   Double_t        trackMomentumX;
   Double_t        trackMomentumY;
   Double_t        trackMomentumZ;
   Double_t        trackMomentumPhi;
   Double_t        trackMomentumTheta;
   Double_t        caloVertexPosX;
   Double_t        caloVertexPosY;
   Double_t        caloVertexPosZ;
   Double_t        caloVertexPosPhi;
   Double_t        caloVertexPosTheta;
   Double_t        decayVertexPosX;
   Double_t        decayVertexPosY;
   Double_t        decayVertexPosZ;
   Double_t        decayVertexMomX;
   Double_t        decayVertexMomY;
   Double_t        decayVertexMomZ;
   Double_t        clusterEoverP;
   Double_t        trackTimeDiff;
   Double_t        decayTime;
   Bool_t          hitVolume;
   Double_t        trackPValue;
   Int_t           station;

   
 private:

    // List of branches
   TBranch        *b_runNum;   //!
   TBranch        *b_subRunNum;   //!
   TBranch        *b_eventNum;   //!
   TBranch        *b_bunchNum;   //!
   TBranch        *b_clusterFillNum;   //!
   TBranch        *b_clusterCaloNum;   //!
   TBranch        *b_clusterIslandNum;   //!
   TBranch        *b_clusterX;   //!
   TBranch        *b_clusterY;   //!
   TBranch        *b_clusterHits;   //!
   TBranch        *b_clusterTime;   //!
   TBranch        *b_clusterE;   //!
   TBranch        *b_trackMomentum;   //!
   TBranch        *b_trackMomentumX;   //!
   TBranch        *b_trackMomentumY;   //!
   TBranch        *b_trackMomentumZ;   //!
   TBranch        *b_trackMomentumPhi;   //!
   TBranch        *b_trackMomentumTheta;   //!
   TBranch        *b_caloVertexPosX;   //!
   TBranch        *b_caloVertexPosY;   //!
   TBranch        *b_caloVertexPosZ;   //!
   TBranch        *b_caloVertexPosPhi;   //!
   TBranch        *b_caloVertexPosTheta;   //!
   TBranch        *b_decayVertexPosX;   //!
   TBranch        *b_decayVertexPosY;   //!
   TBranch        *b_decayVertexPosZ;   //!
   TBranch        *b_decayVertexMomX;   //!
   TBranch        *b_decayVertexMomY;   //!
   TBranch        *b_decayVertexMomZ;   //!
   TBranch        *b_clusterEoverP;   //!
   TBranch        *b_trackTimeDiff;   //!
   TBranch        *b_decayTime;   //!
   TBranch        *b_hitVolume;   //!
   TBranch        *b_trackPValue;   //!
   TBranch        *b_station;   //!
   
   void Init();

   //  TMatrixD w2t(3,3);
   //  TMatrixD w2c(3,3);

};


//====================================================================
//=== clusterTracker cluster reader
//====================================================================


class clusterTrackerClusterReader : public Reader {

 public :
 clusterTrackerClusterReader(TString filename):
  Reader(filename, "clusterTracker", "clusters"){
    Init();
    std::cout<<"==================================="<<std::endl<<std::endl;
      
  }

  // Declaration of leaf types
   Int_t           runNum;
   Int_t           subRunNum;
   Int_t           eventNum;
   Int_t           bunchNum;
   Int_t           islandNum1;
   Int_t           islandNum2;
   Int_t           islandNum3;
   Int_t           calo1;
   Int_t           calo2;
   Int_t           calo3;
   Int_t           xtal1;
   Int_t           xtal2;
   Int_t           xtal3;
   Double_t        x1;
   Double_t        x2;
   Double_t        x3;
   Double_t        y1;
   Double_t        y2;
   Double_t        y3;
   Double_t        time1;
   Double_t        energy2;
   Double_t        time2;
   Double_t        energy1;
   Double_t        energy3;
   Double_t        time3;
   Double_t        energy1_xtal;
   Double_t        Tcoin21;
   Double_t        deltaE21;
   Double_t        Tcoin31;
   Double_t        deltaE31;
   Double_t        Tcoin32;
   Double_t        deltaE32;
   Int_t           dcoin;
   Int_t           tcoin;
   Int_t           nHit1;
   Int_t           nHit2;
   Int_t           nHit3;
   Double_t        dE_low;
   Double_t        dE_high;
   Double_t        dT_low;
   Double_t        dT_high;
   Double_t        T_spl;
  
  
 private:


   // List of branches
   TBranch        *b_runNum;   //!
   TBranch        *b_subRunNum;   //!
   TBranch        *b_eventNum;   //!
   TBranch        *b_bunchNum;   //!
   TBranch        *b_islandNum1;   //!
   TBranch        *b_islandNum2;   //!
   TBranch        *b_islandNum3;   //!
   TBranch        *b_calo1;   //!
   TBranch        *b_calo2;   //!
   TBranch        *b_calo3;   //!
   TBranch        *b_xtal1;   //!
   TBranch        *b_xtal2;   //!
   TBranch        *b_xtal3;   //!
   TBranch        *b_x1;   //!
   TBranch        *b_x2;   //!
   TBranch        *b_x3;   //!
   TBranch        *b_y1;   //!
   TBranch        *b_y2;   //!
   TBranch        *b_y3;   //!
   TBranch        *b_time1;   //!
   TBranch        *b_time2;   //!
   TBranch        *b_time3;   //!
   TBranch        *b_energy1_xtal;   //!
   TBranch        *b_energy1;   //!
   TBranch        *b_energy2;   //!
   TBranch        *b_energy3;   //!
   TBranch        *b_deltaE21;   //!
   TBranch        *b_deltaE31;   //!
   TBranch        *b_deltaE32;   //!
   TBranch        *b_Tcoin21;   //!
   TBranch        *b_Tcoin31;   //!
   TBranch        *b_Tcoin32;   //!
   TBranch        *b_dcoin;   //!
   TBranch        *b_tcoin;   //!
   TBranch        *b_nHit1;   //!
   TBranch        *b_nHit2;   //!
   TBranch        *b_nHit3;   //!
   TBranch        *b_dE_low;   //!
   TBranch        *b_dE_high;   //!
   TBranch        *b_dT_low;   //!
   TBranch        *b_dT_high;   //!
   TBranch        *b_T_spl;   //!

   void Init();

};


//====================================================================
//=== cluster (clusterTree/clusters) reader
//====================================================================
//contains the branches specific to the clusterTree/cluster tree

class clusterReader : public Reader {

 public :
 clusterReader(TString filename):
  Reader(filename, "clusterTree", "clusters"){
    Init();
    std::cout<<"==================================="<<std::endl<<std::endl;
  }

  // Declaration of leaf types
   Double_t        energy;
   Double_t        time;
   Double_t        x;
   Double_t        y;
   std::vector<double>  *xtalNums;
   std::vector<double>  *xtalEnergies;
   std::vector<double>  *xtalTimes;
   UInt_t          nHit;
   UInt_t          caloNum;
   UInt_t          xtalNum;
   UInt_t          islandNum;
   UInt_t          eventNum;
   UInt_t          bunchNum;
   UInt_t          midasSerialNum;
   UInt_t          subRunNum;
   UInt_t          runNum;

 private:

   // List of branches
   TBranch        *b_energy;   //!
   TBranch        *b_time;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_xtalNums;   //!
   TBranch        *b_xtalEnergies;   //!
   TBranch        *b_xtalTimes;   //!
   TBranch        *b_nHit;   //!
   TBranch        *b_caloNum;   //!
   TBranch        *b_xtalNum;   //!
   TBranch        *b_islandNum;   //!
   TBranch        *b_eventNum;   //!
   TBranch        *b_bunchNum;   //!
   TBranch        *b_midasSerialNum;   //!
   TBranch        *b_subRunNum;   //!
   TBranch        *b_runNum;   //
   
   void Init();

};


//====================================================================
//=== tracker (trackerNTup/tracker) reader
//====================================================================
//contains the branches specific to the trackerNtup/tracker tree

class trackerReader : public Reader {

 public :
 trackerReader(TString filename):
  Reader(filename, "trackerNTup", "tracker"){
    Init();
    std::cout<<"==================================="<<std::endl<<std::endl;
  }
  // Declaration of leaf types
   Int_t           runNum;
   Int_t           subRunNum;
   Int_t           eventNum;
   Float_t         trackMomentum;
   Float_t         trackMomentumX;
   Float_t         trackMomentumY;
   Float_t         trackMomentumZ;
   Float_t         trackMomentumUnc;
   Float_t         decayVertexPosX;
   Float_t         decayVertexPosY;
   Float_t         decayVertexPosZ;
   Float_t         decayVertexMomX;
   Float_t         decayVertexMomY;
   Float_t         decayVertexMomZ;
   Float_t         decayVertexUncR;
   Float_t         decayVertexUncY;
   Float_t         decayVertexUncPR;
   Float_t         decayVertexUncPY;
   Float_t         caloVertexPosX;
   Float_t         caloVertexPosY;
   Float_t         caloVertexPosZ;
   Float_t         caloVertexMomX;
   Float_t         caloVertexMomY;
   Float_t         caloVertexMomZ;
   Float_t         caloVertexUncX;
   Float_t         caloVertexUncY;
   Float_t         caloVertexUncPX;
   Float_t         caloVertexUncPY;
   Float_t         trackT0;
   Float_t         time;
   Float_t         decayTime;
   Bool_t          hitVolume;
   Float_t         trackPValue;
   Int_t           station;
   Int_t           nHits;
   Int_t           nUHits;
   Int_t           nVHits;
   Float_t         missedLayersFrac;
   Float_t         minDriftTime;
   Float_t         maxDriftTime;
   Float_t         maxResidual;
   Float_t         extrapolatedDistance;

 private:
   // List of branches
   TBranch        *b_runNum;   //!
   TBranch        *b_subRunNum;   //!
   TBranch        *b_eventNum;   //!
   TBranch        *b_trackMomentum;   //!
   TBranch        *b_trackMomentumX;   //!
   TBranch        *b_trackMomentumY;   //!
   TBranch        *b_trackMomentumZ;   //!
   TBranch        *b_trackMomentumUnc;   //!
   TBranch        *b_decayVertexPosX;   //!
   TBranch        *b_decayVertexPosY;   //!
   TBranch        *b_decayVertexPosZ;   //!
   TBranch        *b_decayVertexMomX;   //!
   TBranch        *b_decayVertexMomY;   //!
   TBranch        *b_decayVertexMomZ;   //!
   TBranch        *b_decayVertexUncR;   //!
   TBranch        *b_decayVertexUncY;   //!
   TBranch        *b_decayVertexUncPR;   //!
   TBranch        *b_decayVertexUncPY;   //!
   TBranch        *b_caloVertexPosX;   //!
   TBranch        *b_caloVertexPosY;   //!
   TBranch        *b_caloVertexPosZ;   //!
   TBranch        *b_caloVertexMomX;   //!
   TBranch        *b_caloVertexMomY;   //!
   TBranch        *b_caloVertexMomZ;   //!
   TBranch        *b_caloVertexUncX;   //!
   TBranch        *b_caloVertexUncY;   //!
   TBranch        *b_caloVertexUncPX;   //!
   TBranch        *b_caloVertexUncPY;   //!
   TBranch        *b_trackT0;   //!
   TBranch        *b_time;   //!
   TBranch        *b_decayTime;   //!
   TBranch        *b_hitVolume;   //!
   TBranch        *b_trackPValue;   //!
   TBranch        *b_station;   //!
   TBranch        *b_nHits;   //!
   TBranch        *b_nUHits;   //!
   TBranch        *b_nVHits;   //!
   TBranch        *b_missedLayersFrac;   //!
   TBranch        *b_minDriftTime;   //!
   TBranch        *b_maxDriftTime;   //!
   TBranch        *b_maxResidual;   //!
   TBranch        *b_extrapolatedDistance;   //!

   
   void Init();

};







//====================================================================
//=== tracker (trackerNTup/tracker) reader
//====================================================================
//contains the branches specific to the trackerNtup/tracker tree

class allmuonsReader : public Reader {

 public :
 allmuonsReader(TString filename):
  Reader(filename, "allmuons", "tree"){
    Init();
    std::cout<<"==================================="<<std::endl<<std::endl;
  }

   // Declaration of leaf types
   UInt_t          midasSerialNum;
   UInt_t          runNum;
   UInt_t          subRunNum;
   UInt_t          eventNum;
   UInt_t          bunchNum;
   UInt_t          ncoin;
   UInt_t          firstclu[108];   //[ncoin]
   UInt_t          nclu_in_coin[108];   //[ncoin]
   UInt_t          ncluster;
   UInt_t          calonum[257];   //[ncluster]
   UInt_t          nhits[257];   //[ncluster]
   Float_t         eneclu[257];   //[ncluster]
   Float_t         efracmaxclu[257];   //[ncluster]
   Float_t         tclu[257];   //[ncluster]
   Float_t         xclu[257];   //[ncluster]
   Float_t         yclu[257];   //[ncluster]
   UInt_t          addrclu[257];   //[ncluster]
   UInt_t          nmatches;
   UInt_t          nvertices;
   UInt_t          cluCaloNum[110];   //[nmatches]
   Float_t         cluX[110];   //[nmatches]
   Float_t         cluY[110];   //[nmatches]
   UInt_t          cluNhit[110];   //[nmatches]
   Float_t         cluTime[110];   //[nmatches]
   Float_t         cluEne[110];   //[nmatches]
   Float_t         EovP[110];   //[nmatches]
   Float_t         Tdiff[110];   //[nmatches]
   UInt_t          trkStationNum[110];   //[nmatches]
   Float_t         trkMomX[110];   //[nmatches]
   Float_t         trkMomY[110];   //[nmatches]
   Float_t         trkMomZ[110];   //[nmatches]
   Float_t         trkPvalue[110];   //[nmatches]
   Int_t           trkNHits;
   Int_t           trkNUHits;
   Int_t           trkNVHits;
   Float_t         trkT0[110];   //[nmatches]
   Bool_t          trkPassCandidateQuality[110];   //[nmatches]
   Bool_t          trkPassTrackQuality[110];   //[nmatches]
   Bool_t          trkPassVertexQuality[110];   //[nmatches]
   Float_t         vX[110];   //[nmatches]
   Float_t         vY[110];   //[nmatches]
   Float_t         vZ[110];   //[nmatches]
   Float_t         vPX[110];   //[nmatches]
   Float_t         vPY[110];   //[nmatches]
   Float_t         vPZ[110];   //[nmatches]
   Float_t         vP[110];   //[nmatches]
   Float_t         decayvX[110];   //[nmatches]
   Float_t         decayvY[110];   //[nmatches]
   Float_t         decayvZ[110];   //[nmatches]
   Float_t         decayMomX[110];   //[nmatches]
   Float_t         decayMomY[110];   //[nmatches]
   Float_t         decayMomZ[110];   //[nmatches]
   Float_t         decayTime[110];   //[nmatches]
   Bool_t          decayHitVolume[110];   //[nmatches]
   UInt_t          nMatchedClusters[221];   //[nvertices]

 private:
   // List of branches
   TBranch        *b_midasSerialNum;   //!
   TBranch        *b_runNum;   //!
   TBranch        *b_subRunNum;   //!
   TBranch        *b_eventNum;   //!
   TBranch        *b_bunchNum;   //!
   TBranch        *b_ncoin;   //!
   TBranch        *b_firstclu;   //!
   TBranch        *b_nclu_in_coin;   //!
   TBranch        *b_ncluster;   //!
   TBranch        *b_calonum;   //!
   TBranch        *b_nhits;   //!
   TBranch        *b_eneclu;   //!
   TBranch        *b_efracmaxclu;   //!
   TBranch        *b_tclu;   //!
   TBranch        *b_xclu;   //!
   TBranch        *b_yclu;   //!
   TBranch        *b_addrclu;   //!
   TBranch        *b_nmatches;   //!
   TBranch        *b_nvertices;   //!
   TBranch        *b_cluCaloNum;   //!
   TBranch        *b_cluX;   //!
   TBranch        *b_cluY;   //!
   TBranch        *b_cluNhit;   //!
   TBranch        *b_cluTime;   //!
   TBranch        *b_cluEne;   //!
   TBranch        *b_EovP;   //!
   TBranch        *b_Tdiff;   //!
   TBranch        *b_trkStationNum;   //!
   TBranch        *b_trkMomX;   //!
   TBranch        *b_trkMomY;   //!
   TBranch        *b_trkMomZ;   //!
   TBranch        *b_trkPvalue;   //!
   TBranch        *b_trkNHits;   //!
   TBranch        *b_trkNUHits;   //!
   TBranch        *b_trkNVHits;   //!
   TBranch        *b_trkT0;   //!
   TBranch        *b_trkPassCandidateQuality;   //!
   TBranch        *b_trkPassTrackQuality;   //!
   TBranch        *b_trkPassVertexQuality;   //!
   TBranch        *b_vX;   //!
   TBranch        *b_vY;   //!
   TBranch        *b_vZ;   //!
   TBranch        *b_vPX;   //!
   TBranch        *b_vPY;   //!
   TBranch        *b_vPZ;   //!
   TBranch        *b_vP;   //!
   TBranch        *b_decayvX;   //!
   TBranch        *b_decayvY;   //!
   TBranch        *b_decayvZ;   //!
   TBranch        *b_decayMomX;   //!
   TBranch        *b_decayMomY;   //!
   TBranch        *b_decayMomZ;   //!
   TBranch        *b_decayTime;   //!
   TBranch        *b_decayHitVolume;   //!
   TBranch        *b_nMatchedClusters;   //!

      void Init();

};

#endif
