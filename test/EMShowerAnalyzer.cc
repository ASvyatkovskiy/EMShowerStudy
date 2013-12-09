// -*- C++ -*-
//
// Package:    EMShowerAnalyzer
// Class:      EMShowerAnalyzer
// 
/**\class EMShowerAnalyzer EMShowerAnalyzer.cc EMShowers/EMShowerAnalyzer/src/EMShowerAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alexey Svyatkovskiy
//         Created:  Thu Aug  5 22:58:40 EDT 2010
// $Id: EMShowerAnalyzer.cc,v 1.3 2010/09/13 15:48:03 asvyatko Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "RecoMuon/GlobalTrackingTools/interface/MuonShowerFinder.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

using namespace edm;
using namespace std;

//
// class declaration
//

class EMShowerAnalyzer : public edm::EDAnalyzer {
   public:
      explicit EMShowerAnalyzer(const edm::ParameterSet&);
      ~EMShowerAnalyzer();

      typedef TransientTrackingRecHit::ConstRecHitContainer ConstRecHitContainer;

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
   MuonServiceProxy * theService;
   MuonShowerFinder* theShowerFinder;
   Handle<reco::TrackCollection> muons;
   InputTag trackLabel_;

   bool debug_;

   TH2D* h_mom_ShowerSize[4];
   TH2D* h_uncorr_ShowerSize[4];
   TH2D* h_uncorr_ShowerSize_prob5[4];
   TH2D* h_uncorr_ShowerSize_noprob5[4];

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EMShowerAnalyzer::EMShowerAnalyzer(const edm::ParameterSet& iConfig)

{
   // the service parameters
   ParameterSet serviceParameters = iConfig.getParameter<ParameterSet>("ServiceParameters");
   theService = new MuonServiceProxy(serviceParameters);

   ParameterSet showerParameters = iConfig.getParameter<ParameterSet>("EMShowerParameters");
   theShowerFinder = new MuonShowerFinder(theService,showerParameters);

   //now do what ever initialization is needed
   trackLabel_ = iConfig.getParameter<edm::InputTag>("TrackLabel");
   debug_ = iConfig.getParameter<bool>("debug");

}


EMShowerAnalyzer::~EMShowerAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

   if (theService) delete theService;
   if (theShowerFinder) delete theShowerFinder;

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
EMShowerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

      theService->update(iSetup); 
      theShowerFinder->setEvent(iEvent);
      theShowerFinder->setServices(theService->eventSetup());

      iEvent.getByLabel(trackLabel_, muons);
      if( !muons.failedToGet() || muons->size() != 0) {
       for(reco::TrackCollection::const_iterator imuon = muons->begin(); imuon != muons->end(); ++imuon ) {

          vector<int> uncorrelatedHits(4);
          vector<double> showerSizes(4);
          vector<bool> hasShowers(4);

          showerSizes = theShowerFinder->stationShowerSizes(*imuon);
          uncorrelatedHits = theShowerFinder->stationUncorrelatedHits(*imuon);
          hasShowers = theShowerFinder->hasStationShower(uncorrelatedHits, showerSizes);

/*          cout << (theShowerFinder->hasStationShower(uncorrelatedHits, showerSizes)).at(0) << " "
                              << (theShowerFinder->hasStationShower(uncorrelatedHits, showerSizes)).at(1) << " " 
                              << (theShowerFinder->hasStationShower(uncorrelatedHits, showerSizes)).at(2) << " " 
                              << (theShowerFinder->hasStationShower(uncorrelatedHits, showerSizes)).at(3) << endl;
*/

          //fill the histograms
          for (int i = 0; i < 4; i++) {

          cout << "size " << theShowerFinder->stationShowerSizes(*imuon).at(i) << endl;
          cout << "hits " << theShowerFinder->stationUncorrelatedHits(*imuon).at(i) << endl;


             h_mom_ShowerSize[i]->Fill(imuon->p(), showerSizes.at(i));
             h_uncorr_ShowerSize[i]->Fill(uncorrelatedHits.at(i), showerSizes.at(i));

             if (hasShowers.at(i)) {
                h_uncorr_ShowerSize_prob5[i]->Fill(uncorrelatedHits.at(i), showerSizes.at(i)); 
              } else {
                h_uncorr_ShowerSize_noprob5[i]->Fill(uncorrelatedHits.at(i), showerSizes.at(i));
             }
         }
      }
   }



}


// ------------ method called once each job just before starting event loop  ------------
void 
EMShowerAnalyzer::beginJob()
{
   edm::Service<TFileService> fs;

   std::ostringstream pprint;
   string histoName;

  for (int i = 0; i != 4; i++) {

    pprint.str("");
    pprint<<i;
    string iname = pprint.str();


    histoName = "h_mom_ShowerSize"+iname;
    h_mom_ShowerSize[i] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(), 100, 0, 1000, 1000, 0, 500);  

    histoName = "h_uncorr_ShowerSize"+iname;
    h_uncorr_ShowerSize[i] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(), 200, 0, 200, 1000, 0, 500);

    histoName = "h_uncorr_ShowerSize_prob5"+iname;
    h_uncorr_ShowerSize_prob5[i] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),200,0,200,1000,0,500);

    histoName = "h_uncorr_ShowerSize_noprob5"+iname;
    h_uncorr_ShowerSize_noprob5[i] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),200,0,200,1000,0,500);
   }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EMShowerAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EMShowerAnalyzer);
