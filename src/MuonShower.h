/**
 *  Class: MuonShower
 *
 *  Description: Object class storing EM shower parameters
 *  in the muon system. 
 *  Parameters stored are:
 *  * the transverse hit-cluster size
 *  * the number of hits not used by the recsegments (uncorrelated hits)
 *
 *  Authors :
 *  A. Svyatkovskiy, Purdue University
 *
 **/

// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include <iostream>

namespace edm {class Event;}
namespace reco {class TransientTrack;}

//
// class declaration
//

class MuonShower {

   public:

	  typedef TransientTrackingRecHit::ConstRecHitContainer ConstRecHitContainer;
          typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;

   public:

	  MuonShower() {};
	  MuonShower(std::vector<double>&, std::vector<int>&);
          ~MuonShower();

         std::vector<double> getStationShowerSizes() const {return theStationShowerSizes;}
         std::vector<int> getStationUncorrelatedHits() const {return theStationUncorrelatedHits;}
	 void setStationShowerSizes(std::vector<double>& allShowerSizes) {theStationShowerSizes = allShowerSizes;}
	 void setStationUncorrelatedHits(std::vector<int>& allUncorrelatedHits) {theStationUncorrelatedHits = allUncorrelatedHits;}

   private:

         std::vector<double> theStationShowerSizes;
         std::vector<int> theStationUncorrelatedHits;

};
