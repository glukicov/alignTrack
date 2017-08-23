//Straight line tracker - takes time islands and forms StraightLineTrackArtRecords
//James Mott (19th September 2016)

// Include needed ART headers
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

//Geant4
#include "Geant4/G4SystemOfUnits.hh"

//art records
#include "gm2dataproducts/strawtracker/StrawTimeIslandArtRecord.hh"
#include "gm2dataproducts/strawtracker/StraightLineTrackArtRecord.hh"

//Straw Geometry
#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"

//Coord systems
#include "gm2geom/coordSystems/CoordSystemsStore_service.hh"
#include "gm2geom/coordSystems/CoordSystem3Vector.hh"  

//Utils
#include "gm2geom/common/Gm2Constants_service.hh"
#include "gm2util/common/dataModuleDefs.hh"
#include "gm2util/common/RootManager.hh"
#include "gm2tracker/utils/StrawGeomUtils.hh"

//C++ includes
#include <iostream>
#include <string>
#include <vector>
#include <math.h> 

namespace gm2strawtracker {

  //
  // Class declaration
  //
  class StraightLineTracker : public art::EDProducer {

  public:

    explicit StraightLineTracker(fhicl::ParameterSet const& pset);

    //Override desired art::EDProducer functions
    void produce(art::Event& event ) override;
    void beginRun(art::Run& run) override;

  private:

    //TODO: SWAP THESE FOR OTHER FUNCTIONS (MEMBERS OF HIT & TRACKS, I GUESS)
    //Forward declare these structs for the next two functions - we'll fill them further down
    struct DriftCircle;
    struct Tangent;
    ////////////////////////////////////////////////////////////////////////////

    //Function to return tangents to given drift circles
    vector<Tangent> calculateTangents(vector<DriftCircle>& circles);

    //Function to return track from intersection of planes described by tangents
    StraightLineTrackArtRecord getTrackFromTangents(Tangent* uTangent, Tangent* vTangent);

    //Producer labels
    std::string islandModuleLabel_;
    std::string islandInstanceName_;
    std::string tracksInstanceName_;
    std::string allTracksInstanceName_;

    //Flags to decide whether to use truth information for LR ambiguity & DCA
    bool useTruthLR_;
    bool useTruthDCA_;

    //Straw geometry
    gm2geom::StrawTrackerGeometry geom_;

    //Helper tools
    StrawGeomUtils geomUtils_;
    gm2geom::CoordSystemsStoreData cs_;

  }; //End of class StraightLineTracker


  //
  // Class implementation
  //

  StraightLineTracker::StraightLineTracker(fhicl::ParameterSet const& pset)
    : islandModuleLabel_( pset.get<std::string>("islandModuleLabel") )
    , islandInstanceName_( pset.get<std::string>("islandInstanceName") )
    , tracksInstanceName_( pset.get<std::string>("tracksInstanceName") )
    , allTracksInstanceName_( pset.get<std::string>("allTracksInstanceName") )
    , useTruthLR_( pset.get<bool>("useTruthLR", false) )
    , useTruthDCA_( pset.get<bool>("useTruthDCA", false) )
    , geom_()
    , geomUtils_()
    , cs_()
  { 
    produces<StraightLineTrackArtRecordCollection> (tracksInstanceName_);
    produces<StraightLineTrackArtRecordCollection> (allTracksInstanceName_);
  }

  // Drift circle struct definition - need these to draw tangents between hits
  struct StraightLineTracker::DriftCircle {

    double uv; // U or V position
    double z;  // Z position
    double r;  // Radius
      
    DriftCircle() : uv(), z(), r(){}
    
    DriftCircle(const double uv, const double z, const double r)
      : uv(uv)
      , z(z)
      , r(r)
    {
    }

    bool operator<(const DriftCircle& circle) {
      return z < circle.z;
    }
  };
  
  // Tangent defintion - for tangents between circles
  struct StraightLineTracker::Tangent {

    double m; // Gradient
    double c; // Intercept
    bool hit0_left; // LR ambiguity for hit 0
    bool hit1_left;
    
    Tangent() : m(), c(), hit0_left(), hit1_left() {}
    
    Tangent(const double m, const double c, const bool hit0_left, const bool hit1_left)
      : m(m)
      , c(c)
      , hit0_left(hit0_left)
      , hit1_left(hit1_left)
    {
    }
  }; // Tangent
  

  void StraightLineTracker::produce(art::Event& event) {
    

    //Create the output collections ready
    std::unique_ptr<StraightLineTrackArtRecordCollection> straightLineTracks(new StraightLineTrackArtRecordCollection);
    std::unique_ptr<StraightLineTrackArtRecordCollection> allStraightLineTracks(new StraightLineTrackArtRecordCollection);

    // Check we're not trying to use truth on real data
    if (event.isRealData() and (useTruthLR_ or useTruthDCA_)) {
      throw cet::exception("StraightLineTracker") << "Module configued to access simulation information, but event is real data.\n";
      return;
    }

    //Get islands
    art::Handle<gm2strawtracker::StrawTimeIslandArtRecordCollection> islandDataHandle;
    bool foundTimeIslands = event.getByLabel(islandModuleLabel_,islandInstanceName_,islandDataHandle);
    if( ! foundTimeIslands ) {
      throw cet::exception("StraightLineTracker") << "No islands in this event (\"" << islandModuleLabel_ << "\":\"" << islandInstanceName_ << "\")\n";
      return;
    }
    gm2strawtracker::StrawTimeIslandArtRecordCollection islands = *islandDataHandle; //Resolve handle to get collection

    // Loop over all islands in event - use old-school for loop as I want the i_island when making art::Ptr    
    for(unsigned int i_island = 0; i_island < islands.size(); i_island++){
      
      // Get this island
      auto island = islands.at(i_island);

      // Loop over digits, fill some histograms and separate out drift circles into U and V
      vector< DriftCircle > uCircles, vCircles;

      // Loop over digits in island
      for(auto digit : island.strawDigits) {

	// In these coordinates, cosmics are travelling in Z, X runs across straws and Y is vertical within module
	gm2geom::CoordSystem3Vector digitWireCentre = digit->wireID.getCentreInStation(cs_);

	// Split based on whether u or v view
	if(digit->wireID.getView() == gm2strawtracker::u_view){

	  // Convert to U coordinate from X,Y values
	  double wirePosU = geomUtils_.getUfromXY(geom_, digitWireCentre.x(),digitWireCentre.y());

	  // Get DCA (either truth or data)
	  double dca = useTruthDCA_ ? digit->strawMCDigit.dca : digit->dca;

	  // Pair holding centre and radius of circle - add sign if we want to use truth LR data
	  if(useTruthLR_){
	    gm2geom::CoordSystem3Vector dcaTruePos = digit->strawMCDigit.position.transform(cs_, WireID::getStationCoordSysName(0)); // Station Num
	    double dcaTruePosU = geomUtils_.getUfromXY(geom_, dcaTruePos.x(),dcaTruePos.y());
	    uCircles.push_back(DriftCircle(wirePosU, digitWireCentre.z(), (dcaTruePosU-wirePosU > 0) ? dca : -dca));
	  } else {
	    uCircles.push_back(DriftCircle(wirePosU, digitWireCentre.z(), dca));
	  }

	} else if(digit->wireID.getView() == gm2strawtracker::v_view){

	  // Convert to V coordinate from X,Y values
	  double wirePosV = geomUtils_.getVfromXY(geom_, digitWireCentre.x(),digitWireCentre.y());
	  
	  // Get DCA (either truth or data)
	  double dca = useTruthDCA_ ? digit->strawMCDigit.dca : digit->dca;

	  // Pair holding centre and radius of circle
	  if(useTruthLR_){
	    gm2geom::CoordSystem3Vector dcaTruePos = digit->strawMCDigit.position.transform(cs_, WireID::getStationCoordSysName(0)); // Station Num
	    double dcaTruePosV = geomUtils_.getVfromXY(geom_, dcaTruePos.x(),dcaTruePos.y());
	    vCircles.push_back(DriftCircle(wirePosV, digitWireCentre.z(), (dcaTruePosV-wirePosV > 0) ? dca : -dca));
	  } else {
	    vCircles.push_back(DriftCircle(wirePosV, digitWireCentre.z(), dca));
	  }

	} // if(view)

      } // it_digit

      //
      // Get true track position and momentum
      //
      gm2geom::CoordSystem3Vector trueTrackVertex;
      gm2geom::CoordSystem3Vector trueTrackMomentum;
      if(!event.isRealData()){

	// Loop over digits and set true track information from hit in lowest straw layer we encounter
	bool firstDigit = true;
	WireID lowestWireID;
	for ( auto digit : island.strawDigits ) {
	  
	  // Get coordinate system for station
	  std::string stationCoordSysName = WireID::getStationCoordSysName(0); // Station Num

	  // If first hit, then store as trueTrack without comparing to anything
	  if(firstDigit) {
	    trueTrackVertex = digit->strawMCDigit.strawMCHits.front()->position.transform(cs_, stationCoordSysName);
	    trueTrackMomentum = digit->strawMCDigit.strawMCHits.front()->momentum.transform(cs_, stationCoordSysName, true); // true -> Momentum-like transform
	    lowestWireID = digit->wireID;
	    firstDigit = false;
    
          // If wire is earlier, then store this one instead
	  } else if (digit->wireID < lowestWireID){
	    trueTrackVertex = digit->strawMCDigit.strawMCHits.front()->position.transform(cs_, stationCoordSysName);
	    trueTrackMomentum = digit->strawMCDigit.strawMCHits.front()->momentum.transform(cs_, stationCoordSysName, true); // true -> Momentum-like transform
	    lowestWireID = digit->wireID;
	  }

	} // for digit

      } // if (event.isRealData())

      //
      // Calculate tangents and see what they overlap with
      //

      // We know there's only 2 U and 2 V hits right now so panic if this is wrong
      if(uCircles.size() != 2 || vCircles.size() != 2){
	throw cet::exception("StraightLineTracker") << "Number of u hits (" << uCircles.size() << ") or number v hits (" << vCircles.size() << ") does not equal expected value (2)\n";
      }

      // Sort drift circle vectors by z value (as it's easier to record LR side)
      std::sort(uCircles.begin(),uCircles.end());
      std::sort(vCircles.begin(),vCircles.end());

      // Get tangents to circles using helper function
      vector<Tangent> uTangents = calculateTangents(uCircles);
      vector<Tangent> vTangents = calculateTangents(vCircles);

      // If we're using truth info for LR ambiguity, then remove all the false combinations from the vectors
      if(useTruthLR_){
	
	// Loop through uTangents and store index of correct combination
	int LRindex = -1;
	for (unsigned int i = 0; i < uTangents.size(); i++){
	  bool hit0_left = (uCircles.at(0).r > 0);
	  bool hit1_left = (uCircles.at(1).r > 0);
	  if( (hit0_left and uTangents.at(i).hit0_left) or (!hit0_left and !uTangents.at(i).hit0_left)){
	    if ( (hit1_left and uTangents.at(i).hit1_left) or (!hit1_left  and !uTangents.at(i).hit1_left) ){
	      LRindex = i;
	    } 
	  }
	}

	// Remove all other tangents - first beginning to element, then after element (which is now at start) to end
	uTangents.erase(uTangents.begin(), uTangents.begin() + LRindex);
	uTangents.erase(uTangents.begin()+1, uTangents.end());

	// Loop through vTangents and store index of correct combination
	for (unsigned int i = 0; i < vTangents.size(); i++){
	  bool hit0_left = (vCircles.at(0).r > 0);
	  bool hit1_left = (vCircles.at(1).r > 0);
	  if( (hit0_left and vTangents.at(i).hit0_left) or (!hit0_left and !vTangents.at(i).hit0_left)){
	    if ( (hit1_left and vTangents.at(i).hit1_left) or (!hit1_left  and !vTangents.at(i).hit1_left) ){
	      LRindex = i;
	    } 
	  }
	}

	// Remove all other tangents - first beginning to element, then after element (which is now at start) to end
	vTangents.erase(vTangents.begin(), vTangents.begin() + LRindex);
	vTangents.erase(vTangents.begin()+1, vTangents.end());

      } // if(useTruthLR_)

      // Vector that we'll populate with tracks that pass cuts
      vector<StraightLineTrackArtRecord> selectedTracks;

      // Loop over all different combinations of tangents
      for(unsigned int i_uTang = 0; i_uTang < uTangents.size(); i_uTang++){
        for(unsigned int i_vTang = 0; i_vTang < vTangents.size(); i_vTang++){

	  // Get this tangent combination
	  Tangent* uTangent = &uTangents.at(i_uTang);
	  Tangent* vTangent = &vTangents.at(i_vTang);

	  // Get track formed by intersection of two planes defined by tangents (one track per combination)
	  StraightLineTrackArtRecord fitTrack = getTrackFromTangents(uTangent, vTangent);

	  // Put island in fitTrack
	  art::Ptr<gm2strawtracker::StrawTimeIslandArtRecord> islandPtr(islandDataHandle, i_island);
	  fitTrack.island = islandPtr;

	  // Convert true track to be defined at z = 0 and with pz = 1 (as this is what fit track has)
	  fitTrack.trueMom = trueTrackMomentum / trueTrackMomentum.z();  
	  fitTrack.truePos = trueTrackVertex - trueTrackVertex.z()*fitTrack.trueMom;

	  // Put into all tracks collection
	  allStraightLineTracks->push_back(fitTrack);

	  // Check for tracks reconstructed outside of tracker region

	  // Some variables we'll fill and flag for cut
	  double xPos, yPos;
	  bool trackInStraw = true;

	  // Cut based on 2 sigma from edge 41.5 - four different straw layers
	  geomUtils_.getTrackXY(fitTrack, xPos, yPos, -12.7); if(fabs (yPos) > 51.9) trackInStraw = false;
	  geomUtils_.getTrackXY(fitTrack, xPos, yPos, -7.55); if(fabs (yPos) > 51.9) trackInStraw = false;
	  geomUtils_.getTrackXY(fitTrack, xPos, yPos, 7.5);   if(fabs (yPos) > 51.9) trackInStraw = false;
	  geomUtils_.getTrackXY(fitTrack, xPos, yPos, 12.69); if(fabs (yPos) > 51.9) trackInStraw = false;

	  // If track was in straws for all modules, then we'll store it in selected tracks
	  if (trackInStraw) straightLineTracks->push_back(fitTrack);
	  
	} //i_vTang
      } //i_uTang

    } //island
    
    // Put the output collections in the event
    event.put(std::move(straightLineTracks),tracksInstanceName_);
    event.put(std::move(allStraightLineTracks),allTracksInstanceName_);

  }//analyze

  
  void StraightLineTracker::beginRun(art::Run & run) {

    //Get coord systems
    cs_ = artg4::dataFromRunOrService<gm2geom::CoordSystemsStoreData, gm2geom::CoordSystemsStore>
          ( run, dataModuleDefs::coordSysModuleLabel(),dataModuleDefs::coordSysInstanceLabel() );
    if( cs_.size() == 0 ) {
      mf::LogWarning("StraightLineTracker") << "This run does not contain any data associated with the coordinate system\n";
    }

  }//beginRun


  vector<StraightLineTracker::Tangent> StraightLineTracker::calculateTangents(vector<DriftCircle>& circles) {

    vector<Tangent> calcTangents;

    // Get two circles
    DriftCircle* circle0 = &circles.at(0);
    DriftCircle* circle1 = &circles.at(1);
    
    // Parameters describing positions of circles - put in terms of u here, but equally applicable for v
    double du = circle1->uv - circle0->uv;
    double dz = circle1->z - circle0->z;
    double r0 = circle0->r;
    double r1 = circle1->r;
    double du_dz = du/dz;
    double centreDist = sqrt( du*du + dz*dz );

    // Have four tangents to 2 circles
    int nTangents = 4;

    // Unless we've passed just the wire values when there's only one
    if(r0 == 0 and r1 == 0) nTangents = 1;

    // Calculate each one by one in this four loop (shameless stolen from Wikipedia article)
    for (int iTang = 0; iTang < nTangents; iTang++){

      // These two parameters define which tangent we're taking
      int sgn;
      double radDiff;
      switch (iTang) {
      case 0:
	sgn = 1;
	radDiff = fabs(r1) - fabs(r0);
	break;
      case 1:
	sgn = -1;
	radDiff = fabs(r1) - fabs(r0);
	break;
      case 2:
	sgn = -1;
	radDiff = -fabs(r1) - fabs(r0);
	break;
      case 3:
	sgn = 1;
	radDiff = -fabs(r1) - fabs(r0);
	break;
      }

      // This sqrt shows up a lot - just redefined it to neaten things up
      double sqrtRad = sqrt(pow(centreDist/radDiff,2)-1);
      
      // Calculate gradient and intercept of line (U = m*Z + c)
      double m = (1 + sgn*du_dz*sqrtRad) / (sgn*sqrtRad - du_dz);
      double c = circle0->uv - m*circle0->z - fabs(r0)/(radDiff/(centreDist*centreDist) * (du-sgn*dz*sqrtRad));
      
      // Fix special case where radii are equal (and looking at external tangents)
      if(radDiff == 0){
	m = du_dz;
	c = circle0->uv - m*circle0->z - fabs(r0)/(sgn*dz/centreDist);
      }

      // Calculate whether we passed on the left or right side of each wire
      bool hit0_left = false;
      bool hit1_left = false;
      if (m*circle0->z + c > circle0->uv) hit0_left = true;
      if (m*circle1->z + c > circle1->uv) hit1_left = true;

      // Push back this gradient to vector we're going to return
      calcTangents.push_back(Tangent(m, c, hit0_left, hit1_left));
    }

    return calcTangents;
  }

  
  StraightLineTrackArtRecord StraightLineTracker::getTrackFromTangents(Tangent* uTangent, Tangent* vTangent) {

    // Each tangent defines a plane and intersection of these two places gives track that satisfies both tangents

    // Track that we'll return
    StraightLineTrackArtRecord track;

    // Coordinate system we'll use
    std::string stationCoordSysName = WireID::getStationCoordSysName(0); // Station Num

    // Get point on track by considering z = 0 (where x0, y0 is and u- and v-tangents have c values)
    double x0, y0;
    geomUtils_.getXYfromUV(geom_, x0, y0, uTangent->c, vTangent->c);
    gm2geom::CoordSystem3Vector trackPos(x0, y0, 0, stationCoordSysName);
    track.pos = trackPos;

    // Direction is given by considering vector defined by gradients in u and v and translating into x,y (assume pz = 1)
    double px, py;
    geomUtils_.getXYfromUV(geom_, px, py, uTangent->m, vTangent->m);
    gm2geom::CoordSystem3Vector trackMom(px, py, 1, stationCoordSysName);
    track.mom = trackMom;

    return track;
  }

} // End of namespace gm2strawtracker

//
// Extras
//

// These are some necessary boilerplate for the ROOT persistency system
using gm2strawtracker::StraightLineTracker;
DEFINE_ART_MODULE(StraightLineTracker)
