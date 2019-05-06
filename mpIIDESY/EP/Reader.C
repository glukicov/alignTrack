//Simple code to read the Europa OmegaA ntuples
// gavin.hesketh@ucl.ac.uk
//
// Edit this code simply to un-comment any branches you want to use
// Can uncomment them all, but it will make the code slower.

#include "Reader.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//================================================
// trackReader
//================================================

void clusterTrackerTrackReader::Init() {

  //comment out lines to decativate branches you don't need.

  //LoadBranch("runNum", runNum, b_runNum);
  //  LoadBranch("subRunNum", subRunNum, b_subRunNum);
  LoadBranch("eventNum", eventNum, b_eventNum);
  //  LoadBranch("bunchNum", bunchNum, b_bunchNum);
  //  LoadBranch("clusterFillNum", clusterFillNum, b_clusterFillNum);
 LoadBranch("clusterCaloNum", clusterCaloNum, b_clusterCaloNum);
  //  LoadBranch("clusterIslandNum", clusterIslandNum, b_clusterIslandNum);
  LoadBranch("clusterX", clusterX, b_clusterX);
  LoadBranch("clusterY", clusterY, b_clusterY);
  LoadBranch("clusterHits", clusterHits, b_clusterHits);
  LoadBranch("clusterTime", clusterTime, b_clusterTime);
  LoadBranch("clusterE", clusterE, b_clusterE);  
  LoadBranch("trackMomentum", trackMomentum, b_trackMomentum); 
  // LoadBranch("trackMomentumX", trackMomentumX, b_trackMomentumX);
  // LoadBranch("trackMomentumY", trackMomentumY, b_trackMomentumY);
  //LoadBranch("trackMomentumZ", trackMomentumZ, b_trackMomentumZ);
  LoadBranch("trackMomentumPhi", trackMomentumPhi, b_trackMomentumPhi);
  LoadBranch("trackMomentumTheta", trackMomentumTheta, b_trackMomentumTheta); 
  LoadBranch("caloVertexPosX", caloVertexPosX, b_caloVertexPosX);
  LoadBranch("caloVertexPosY", caloVertexPosY, b_caloVertexPosY);
  //  LoadBranch("caloVertexPosZ", caloVertexPosZ, b_caloVertexPosZ);
  // LoadBranch("caloVertexPosPhi", caloVertexPosPhi, b_caloVertexPosPhi);
  //LoadBranch("caloVertexPosTheta", caloVertexPosTheta, b_caloVertexPosTheta);
  //LoadBranch("decayVertexPosX", decayVertexPosX, b_decayVertexPosX);
  LoadBranch("decayVertexPosY", decayVertexPosY, b_decayVertexPosY);
  //LoadBranch("decayVertexPosZ", decayVertexPosZ, b_decayVertexPosZ);
  //LoadBranch("decayVertexMomX", decayVertexMomX, b_decayVertexMomX);
  //LoadBranch("decayVertexMomY", decayVertexMomY, b_decayVertexMomY);
  //LoadBranch("decayVertexMomZ", decayVertexMomZ, b_decayVertexMomZ);
  LoadBranch("decayTime", decayTime, b_decayTime);
  // LoadBranch("hitVolume", hitVolume, b_hitVolume);
  LoadBranch("trackPValue", trackPValue, b_trackPValue);
  LoadBranch("station", station, b_station);
  LoadBranch("clusterEoverP", clusterEoverP, b_clusterEoverP);
  LoadBranch("trackTimeDiff", trackTimeDiff, b_trackTimeDiff);
  
}



//================================================
// clusterReader (in clusterTracker folder)
//================================================

void clusterTrackerClusterReader::Init() {
  //comment out lines to decativate branches you don't need
  //LoadBranch("runNum", runNum, b_runNum);
  //LoadBranch("subRunNum", subRunNum, b_subRunNum);
  //LoadBranch("eventNum", eventNum, b_eventNum);
  //LoadBranch("bunchNum", bunchNum, b_bunchNum);
  //LoadBranch("islandNum1", islandNum1, b_islandNum1);
  //LoadBranch("islandNum2", islandNum2, b_islandNum2);
  //LoadBranch("islandNum3", islandNum3, b_islandNum3);
  //LoadBranch("calo1", calo1, b_calo1);
  //LoadBranch("calo2", calo2, b_calo2);
  //LoadBranch("calo3", calo2, b_calo3);
  //LoadBranch("xtal1", xtal1, b_xtal1);
  //LoadBranch("xtal2", xtal2, b_xtal2);
  //LoadBranch("xtal3", xtal2, b_xtal3);
  //LoadBranch("x1", x1, b_x1);
  //LoadBranch("x2", x2, b_x1);
  //LoadBranch("x3", x2, b_x3);
  //LoadBranch("y1", y1, b_y1);
  //LoadBranch("y2", y2, b_y2);
  //LoadBranch("y3", y2, b_y3);
  //LoadBranch("time1", time1, b_time1);
  //LoadBranch("time2", time2, b_time2);
  //LoadBranch("time3", time3, b_time3);
  //LoadBranch("energy1", energy1, b_energy1);
  //LoadBranch("energy1_xtal", energy1_xtal, b_energy1_xtal);
  //LoadBranch("energy2", energy2, b_energy2);
  //LoadBranch("energy3", energy3, b_energy3);
  //LoadBranch("Tcoin21", Tcoin21, b_Tcoin21);
  //LoadBranch("Tcoin31", Tcoin31, b_Tcoin31);
  //LoadBranch("Tcoin32", Tcoin32, b_Tcoin32);
  //LoadBranch("deltaE21", deltaE21, b_deltaE21);
  //LoadBranch("deltaE31", deltaE31, b_deltaE31);
  //LoadBranch("deltaE32", deltaE32, b_deltaE32);
  //LoadBranch("dcoin", dcoin, b_dcoin);
  //LoadBranch("tcoin", tcoin, b_tcoin);
  //LoadBranch("nHit1", nHit1, b_nHit1);
  //LoadBranch("nHit2", nHit2, b_nHit2);
  //LoadBranch("nHit3", nHit3, b_nHit3);
  //LoadBranch("dE_low", dE_low, b_dE_low);
  //LoadBranch("dE_high", dE_high, b_dE_high);
  //LoadBranch("dT_low", dT_low, b_dT_low);
  //LoadBranch("dT_high", dT_high, b_dT_high);
  //LoadBranch("T_spl", T_spl, b_T_spl);
}





//================================================
// cluster Reader (in clusterTree folder)
//================================================

void clusterReader::Init() {
  //comment out lines to decativate branches you don't need

  LoadBranch("energy", energy, b_energy);
  LoadBranch("time", time, b_time);
  LoadBranch("x", x, b_x);
  LoadBranch("y", y, b_y);
  //  LoadBranch("xtalNums", xtalNums, b_xtalNums);  //NOT IN NEW NTUPLE!
  //  LoadBranch("xtalEnergies", xtalEnergies, b_xtalEnergies);  //NOT IN NEW NTUPLE!
  //  LoadBranch("xtalTimes", xtalTimes, b_xtalTimes);  //NOT IN NEW NTUPLE!
  LoadBranch("nHit", nHit, b_nHit);
  LoadBranch("caloNum", caloNum, b_caloNum);
  //LoadBranch("xtalNum", xtalNum, b_xtalNum);
  LoadBranch("islandNum", islandNum, b_islandNum);
  LoadBranch("eventNum", eventNum, b_eventNum);
  //LoadBranch("bunchNum", bunchNum, b_bunchNum);
  //LoadBranch("midasSerialNum", midasSerialNum, b_midasSerialNum);
  LoadBranch("subRunNum", subRunNum, b_subRunNum);
  LoadBranch("runNum", runNum, b_runNum);
}




//================================================
// cluster Reader (in clusterTree folder)
//================================================

void trackerReader::Init() {
  // Declaration of leaf types


   LoadBranch("runNum", runNum, b_runNum);
   LoadBranch("subRunNum", subRunNum, b_subRunNum);
   LoadBranch("eventNum", eventNum, b_eventNum);
   LoadBranch("trackMomentum", trackMomentum, b_trackMomentum);
   LoadBranch("trackMomentumX", trackMomentumX, b_trackMomentumX);
   LoadBranch("trackMomentumY", trackMomentumY, b_trackMomentumY);
   LoadBranch("trackMomentumZ", trackMomentumZ, b_trackMomentumZ);
   LoadBranch("trackMomentumUnc", trackMomentumUnc, b_trackMomentumUnc);
   LoadBranch("decayVertexPosX", decayVertexPosX, b_decayVertexPosX);
   LoadBranch("decayVertexPosY", decayVertexPosY, b_decayVertexPosY);
   LoadBranch("decayVertexPosZ", decayVertexPosZ, b_decayVertexPosZ);
   LoadBranch("decayVertexMomX", decayVertexMomX, b_decayVertexMomX);
   LoadBranch("decayVertexMomY", decayVertexMomY, b_decayVertexMomY);
   LoadBranch("decayVertexMomZ", decayVertexMomZ, b_decayVertexMomZ);
   LoadBranch("decayVertexUncR", decayVertexUncR, b_decayVertexUncR);
   LoadBranch("decayVertexUncY", decayVertexUncY, b_decayVertexUncY);
   LoadBranch("decayVertexUncPR", decayVertexUncPR, b_decayVertexUncPR);
   LoadBranch("decayVertexUncPY", decayVertexUncPY, b_decayVertexUncPY);
   LoadBranch("caloVertexPosX", caloVertexPosX, b_caloVertexPosX);
   LoadBranch("caloVertexPosY", caloVertexPosY, b_caloVertexPosY);
   LoadBranch("caloVertexPosZ", caloVertexPosZ, b_caloVertexPosZ);
   LoadBranch("caloVertexMomX", caloVertexMomX, b_caloVertexMomX);
   LoadBranch("caloVertexMomY", caloVertexMomY, b_caloVertexMomY);
   LoadBranch("caloVertexMomZ", caloVertexMomZ, b_caloVertexMomZ);
   LoadBranch("caloVertexUncX", caloVertexUncX, b_caloVertexUncX);
   LoadBranch("caloVertexUncY", caloVertexUncY, b_caloVertexUncY);
//    LoadBranch("caloVertexMomX", caloVertexMomX, b_caloVertexUncPX);
//    LoadBranch("caloVertexMomY", caloVertexMomY, b_caloVertexUncPY);
   LoadBranch("trackT0", trackT0, b_trackT0);
   LoadBranch("time", time, b_time);
   LoadBranch("decayTime", decayTime, b_decayTime);
   LoadBranch("hitVolume", hitVolume, b_hitVolume);
   LoadBranch("trackPValue", trackPValue, b_trackPValue);
   LoadBranch("station", station, b_station);
   LoadBranch("nHits", nHits, b_nHits);
   LoadBranch("nUHits", nUHits, b_nUHits);
   LoadBranch("nVHits", nVHits, b_nVHits);
   LoadBranch("missedLayersFrac", missedLayersFrac, b_missedLayersFrac);
   LoadBranch("minDriftTime", minDriftTime, b_minDriftTime);
   LoadBranch("maxDriftTime", maxDriftTime, b_maxDriftTime);
   LoadBranch("maxResidual", maxResidual, b_maxResidual);
   LoadBranch("extrapolatedDistance", extrapolatedDistance, b_extrapolatedDistance);
   
  
   
}


void allmuonsReader::Init() {


  
  // LoadBranch("midasSerialNum", midasSerialNum, b_midasSerialNum);
  // LoadBranch("runNum", runNum, b_runNum);
  //  LoadBranch("subRunNum", subRunNum, b_subRunNum);
  //  LoadBranch("eventNum", eventNum, b_eventNum);
  // LoadBranch("bunchNum", bunchNum, b_bunchNum);
   // number of coincidences
  // LoadBranch("trkPassCandidateQuality",trkPassCandidateQuality,b_trkPassCandidateQuality);
   LoadBranchVector("trkPassTrackQuality",trkPassTrackQuality,b_trkPassTrackQuality);
  // LoadBranch("trkPassVertexQuality",trkPassVertexQuality,b_trkPassVertexQuality);
  LoadBranch("ncoin", ncoin, b_ncoin);
   LoadBranchVector("firstclu", firstclu, b_firstclu);
   LoadBranchVector("nclu_in_coin", nclu_in_coin, b_nclu_in_coin);

   LoadBranch("ncluster", ncluster, b_ncluster);
   LoadBranchVector("calonum", calonum, b_calonum);
   LoadBranchVector("nhits", nhits, b_nhits);
   LoadBranchVector("eneclu", eneclu, b_eneclu);
   LoadBranchVector("efracmaxclu", efracmaxclu, b_efracmaxclu);
   LoadBranchVector("tclu", tclu, b_tclu);
   LoadBranchVector("xclu", xclu, b_xclu);
   LoadBranchVector("yclu", yclu, b_yclu);
   //   LoadBranchVector("addrclu", addrclu, b_addrclu);

   LoadBranch("nmatches", nmatches, b_nmatches);
   LoadBranchVector("cluCaloNum", cluCaloNum, b_cluCaloNum);
   LoadBranchVector("cluX", cluX, b_cluX);
   LoadBranchVector("cluY", cluY, b_cluY);
   LoadBranchVector("cluNhit", cluNhit, b_cluNhit);  // number of crystals per cluster
   LoadBranchVector("cluTime", cluTime, b_cluTime);
   LoadBranchVector("cluEne", cluEne, b_cluEne);
   LoadBranchVector("EovP", EovP, b_EovP);
   LoadBranchVector("Tdiff", Tdiff, b_Tdiff);
   LoadBranchVector("trkStationNum", trkStationNum, b_trkStationNum);
   LoadBranchVector("trkMomX", trkMomX, b_trkMomX);
   LoadBranchVector("trkMomY", trkMomY, b_trkMomY);
   LoadBranchVector("trkMomZ", trkMomZ, b_trkMomZ);
   LoadBranchVector("trkPvalue", trkPvalue, b_trkPvalue);
   //   LoadBranchVector("trkNHits", trkNHits, b_trkNHits);
   //  LoadBranchVector("trkNUHits", trkNUHits, b_trkNUHits);
   //  LoadBranchVector("trkNVHits", trkNVHits, b_trkNVHits);
   // LoadBranchVector("trkT0", trkT0, b_trkT0);
   // LoadBranchVector("trkPassCandidateQuality", trkPassCandidateQuality, b_trkPassCandidateQuality);
   //LoadBranchVector("trkPassTrackQuality", trkPassTrackQuality, b_trkPassTrackQuality);
   //LoadBranchVector("trkPassVertexQuality", trkPassVertexQuality, b_trkPassVertexQuality);
   LoadBranchVector("vX", vX, b_vX);
   LoadBranchVector("vY", vY, b_vY);
   LoadBranchVector("vZ", vZ, b_vZ);
   //  LoadBranchVector("vPX", vPX, b_vPX);
   // LoadBranchVector("vPY", vPY, b_vPY);
   //LoadBranchVector("vPZ", vPZ, b_vPZ);
   // LoadBranchVector("vP", vP, b_vP);
   LoadBranchVector("decayvX", decayvX, b_decayvX);
   LoadBranchVector("decayvY", decayvY, b_decayvY);
   LoadBranchVector("decayvZ", decayvZ, b_decayvZ);
   LoadBranchVector("decayMomX", decayMomX, b_decayMomX);
   LoadBranchVector("decayMomY", decayMomY, b_decayMomY);
   LoadBranchVector("decayMomZ", decayMomZ, b_decayMomZ);
   LoadBranchVector("decayTime", decayTime, b_decayTime);
   LoadBranchVector("decayHitVolume", decayHitVolume, b_decayHitVolume);

   // LoadBranch("nvertices", nvertices, b_nvertices);
   // LoadBranchVector("nMatchedClusters", nMatchedClusters, b_nMatchedClusters);



   
}



//================================================
//General Reader:
//================================================

Reader::Reader(TString filename, TString folder, TString treename) :
  fChain(0), fCurrent(-1), jentry_(0) {  
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
  if (!f || !f->IsOpen()) {
    f = new TFile(filename);
  }
  TDirectory * dir = (TDirectory*)f->Get(filename+":/"+folder);
  dir->GetObject(treename,fChain);
  if (!fChain) std::cout<<"Tree Not Found "<<filename<<" "<<filename+":/"+folder+"/"+treename<<std::endl;
  fChain->SetMakeClass(1);
  fChain->SetBranchStatus("*",0);  // disable all branches by default
  
  nentries_ = fChain->GetEntriesFast();
  
  std::cout<<std::endl<<"==================================="<<std::endl;
  std::cout<<"Loaded tree "<<folder<<"/"<<fChain->GetName()<<" from "<<filename<<" with "<<nentries_<<" entries"<<std::endl;
}




template<class VAR>
void Reader::LoadBranch(TString name, VAR &var, TBranch *&branch) { 
    std::cout<<"Activating branch "<<name<<std::endl;
    var=0;
    fChain->SetBranchStatus(name,1);
    fChain->SetBranchAddress(name, &var, &branch);
}
template<class VAR>
void Reader::LoadBranchVector(TString name, VAR &var, TBranch *&branch) { 
    std::cout<<"Activating branch "<<name<<std::endl;
    // var=0;
    fChain->SetBranchStatus(name,1);
    fChain->SetBranchAddress(name, var, &branch);
}

bool Reader::NextEvent() {
  //  std::cout<<"Loading entry: "<<jentry_<<std::endl;
  Long64_t centry = fChain->LoadTree(jentry_);
  if (centry < 0) return false;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
  }
  fChain->GetEntry(jentry_);   
  jentry_++;
 
  return true;
}

