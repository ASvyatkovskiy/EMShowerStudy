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

#include "RecoMuon/GlobalTrackingTools/interface/MuonShower.h"

using namespace std;
using namespace edm;


//
//Constructor
//
MuonShower::MuonShower(std::vector<double>& stationShowerSizes, std::vector<int>& stationUncorrelatedHits) {

	 setStationShowerSizes(stationShowerSizes);
	 setStationUncorrelatedHits(stationUncorrelatedHits);

}

//
//Destructor
//
MuonShower::~MuonShower() {}
