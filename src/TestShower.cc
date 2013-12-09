// -*- C++ -*-
//
// Package:    TestShower
// Class:      TestShower
// 
/**\class TestShower TestShower.cc TestShower/TestShower/src/TestShower.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alexey Svyatkovskiy
//         Created:  Mon Nov 22 11:25:17 EST 2010
// $Id: TestShower.cc,v 1.1 2011/03/31 23:58:48 asvyatko Exp $
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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/MuonReco/interface/MuonShower.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/PluginManager/interface/PluginManager.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

#include <iostream>
#include <string>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

//
// class declaration
//

using namespace std;
using namespace edm;

class TestShower : public edm::EDAnalyzer {
   public:
      explicit TestShower(const edm::ParameterSet&);
      ~TestShower();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      //station hits used by a track
      int StationHitsUsedByTrack(const reco::Track&, int);

      // ----------member data ---------------------------
      //histograms for shower related vars by station
      TH1D* Muon_Station_showerT[4];
      TH1D* Muon_Station_deltaR[4];
      TH1D* Muon_Station_AllHits[4];
      TH1D* Muon_Station_CorrHits[4];
      TH1D* Muon_Station_UncorrHits[4];
      TH1D* Muon_Station_HitsNotUsedByTrack[4];

      //as a function of inner p, pT by station
      TH2D* Muon_Station_showerT_pT[4];
      TH2D* Muon_Station_deltaR_pT[4];
      TH2D* Muon_Station_AllHits_pT[4];
      TH2D* Muon_Station_CorrHits_pT[4];
      TH2D* Muon_Station_UncorrHits_pT[4];
      TH2D* Muon_Station_HitsNotUsedByTrack_pT[4];

      TH2D* Muon_Station_showerT_p[4];
      TH2D* Muon_Station_deltaR_p[4];
      TH2D* Muon_Station_AllHits_p[4];
      TH2D* Muon_Station_CorrHits_p[4];
      TH2D* Muon_Station_UncorrHits_p[4];
      TH2D* Muon_Station_HitsNotUsedByTrack_p[4];
    
      //as a function of eta by station
      TH2D* Muon_Station_showerT_eta[4];
      TH2D* Muon_Station_deltaR_eta[4];
      TH2D* Muon_Station_AllHits_eta[4];
      TH2D* Muon_Station_CorrHits_eta[4];
      TH2D* Muon_Station_UncorrHits_eta[4];
      TH2D* Muon_Station_HitsNotUsedByTrack_eta[4];

      //simple monitoring histograms
      TH1D* Muon_inner_pT;
      TH1D* Muon_inner_p;
      TH1D* Muon_eta;

      TH1D* Muon_inner_pT_preselected;
      TH1D* Muon_inner_p_preselected;
      TH1D* Muon_eta_preselected;

      //histograms for shower rates
      //by station , in principle the same histo, do we need it that way?
      TH1D* Muon_Station_inner_pT_prob5[4];
      TH1D* Muon_Station_inner_p_prob5[4];
      TH1D* Muon_Station_eta_prob5[4];

      std::string theRootFileName;
      edm::InputTag theMuonLabel;
      edm::InputTag theTrackLabel;
      edm::InputTag inputMuonShowerInformationValueMap_;
      
};

//
// constructors and destructor
//
TestShower::TestShower(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  theRootFileName                   = iConfig.getUntrackedParameter<string>("out", "output.root");

  theMuonLabel                      = iConfig.getUntrackedParameter<edm::InputTag>("Muon", edm::InputTag("selectedPatMuonsTriggerMatch"));
  theTrackLabel                     = iConfig.getUntrackedParameter<edm::InputTag>("Track", edm::InputTag("generalTracks"));

  inputMuonShowerInformationValueMap_ = iConfig.getUntrackedParameter<edm::InputTag>("inputMuonShowerInformationValueMap");

}


TestShower::~TestShower()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TestShower::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
   cout << "Analyzing event: " << iEvent.id() << endl;

   using namespace edm;

   // call pat objects
   edm::Handle< reco::MuonCollection > muonHandle;
   iEvent.getByLabel("muons", muonHandle);

   edm::Handle<reco::VertexCollection> pvHandle;
   iEvent.getByLabel("offlinePrimaryVertices",pvHandle);

   edm::Handle<edm::ValueMap<reco::MuonShower> > muonShowerInformationValueMapH_;
   iEvent.getByLabel(inputMuonShowerInformationValueMap_.label(), muonShowerInformationValueMapH_);

//
//Determine the leadig muon momentum (for collision events)
//
   double leadingMu = 0.;
   for(reco::MuonCollection::const_iterator imuon = muonHandle->begin(); imuon != muonHandle->end(); ++imuon ) {

       reco::TrackRef tk = imuon->innerTrack();
       if (!tk.isNonnull()) continue; 
       if (tk->p() > leadingMu) leadingMu = tk->p();
   }

//
//Fill
//
   int imucount1 = -1;
   for(reco::MuonCollection::const_iterator iMuon = muonHandle->begin(); iMuon != muonHandle->end(); ++iMuon) {

       imucount1++;
       reco::MuonRef muonRef(muonHandle, imucount1);
       reco::MuonShower muonShowerInformation = (*muonShowerInformationValueMapH_)[muonRef];

       //shower preselection
       reco::TrackRef tk = iMuon->innerTrack();
       reco::TrackRef gm = iMuon->globalTrack();
       reco::TrackRef ot = iMuon->outerTrack();

       if (gm.isNonnull() && tk.isNonnull() && ot.isNonnull()) {

       ///FIXME remove this afterwards 
       if (fabs(tk->eta()) > 0.9) continue;
       if (tk->p() != leadingMu || tk->p() < 150) continue;

       //monitor
       Muon_inner_pT->Fill(tk->pt());
       Muon_inner_p->Fill(tk->p());
       Muon_eta->Fill(tk->eta());


       //offline selection for collision muons         
       if (pvHandle->empty() || pvHandle->front().isFake()) continue;
       const reco::Vertex &vtx = pvHandle->front();
       if (fabs(tk->dxy(vtx.position())) > 0.2) continue;
       int pixelHits = tk->hitPattern().numberOfValidPixelHits();
       int trackerHits = tk->hitPattern().numberOfValidTrackerHits();
       int muonHits = gm->hitPattern().numberOfValidMuonHits(); 
       double tknormalizedChi2 = tk->normalizedChi2();
       double glbnormalizedChi2 = gm->normalizedChi2();
       if (iMuon->numberOfMatches() < 2 || glbnormalizedChi2 > 10 || tknormalizedChi2 > 5 ||  muonHits < 1 || pixelHits < 1 || trackerHits < 11 || (iMuon->isolationR03().sumPt + iMuon->isolationR03().hadEt + iMuon->isolationR03().emEt) > 2.5) continue;

       //monitor
       Muon_inner_pT_preselected->Fill(tk->pt());
       Muon_inner_p_preselected->Fill(tk->p());
       Muon_eta_preselected->Fill(tk->eta());

       //by station
       for (int istat = 0; istat < 4; istat++) {

          if (istat != 3 && (((muonShowerInformation.stationShowerSizeT).at(istat) > 20 && ((muonShowerInformation.nStationHits).at(istat) - 2*StationHitsUsedByTrack(*ot, istat+1)) > 40) || ((muonShowerInformation.stationShowerSizeT).at(istat) > 15 && ((muonShowerInformation.nStationHits).at(istat) - 2*StationHitsUsedByTrack(*ot, istat+1)) > 60))) {
             cout << "++++++++++++++STAT: " << istat << " ++++++++++++++++++" << endl;
             cout << "Muon p/pT/eta " << tk->p() << " " << tk->pt() << " " << tk->eta() << endl;
             cout << "SizeT " << (muonShowerInformation.stationShowerSizeT).at(istat) << " Station Hits used by track " << 2*StationHitsUsedByTrack(*ot, istat+1) << " out of Station Hits " << (muonShowerInformation.nStationHits).at(istat) << endl;  

             Muon_Station_inner_pT_prob5[istat]->Fill(tk->pt());
             Muon_Station_inner_p_prob5[istat]->Fill(tk->p());
             Muon_Station_eta_prob5[istat]->Fill(tk->eta());
          } else if (istat == 3 && ((muonShowerInformation.stationShowerSizeT).at(istat) > 20 && ((muonShowerInformation.nStationHits).at(istat) - 2*StationHitsUsedByTrack(*ot, istat+1)) > 25)) {
             cout << "++++++++++++++STAT: " << istat << " ++++++++++++++++++" << endl;
             cout << "Muon pT/eta " << tk->pt() << " " << tk->eta() << endl;
             cout << "SizeT " << (muonShowerInformation.stationShowerSizeT).at(istat) << " Station Hits used by track " << 2*StationHitsUsedByTrack(*ot, istat+1) << " out of Station Hits " << (muonShowerInformation.nStationHits).at(istat) << endl;
             Muon_Station_inner_pT_prob5[istat]->Fill(tk->pt());
             Muon_Station_inner_p_prob5[istat]->Fill(tk->p());
             Muon_Station_eta_prob5[istat]->Fill(tk->eta());
}

          Muon_Station_showerT[istat]->Fill((muonShowerInformation.stationShowerSizeT).at(istat));
          Muon_Station_deltaR[istat]->Fill((muonShowerInformation.stationShowerDeltaR).at(istat));
          Muon_Station_AllHits[istat]->Fill((muonShowerInformation.nStationHits).at(istat));
          Muon_Station_CorrHits[istat]->Fill((muonShowerInformation.nStationCorrelatedHits).at(istat));
          Muon_Station_UncorrHits[istat]->Fill((muonShowerInformation.nStationHits).at(istat) - (muonShowerInformation.nStationCorrelatedHits).at(istat));
          Muon_Station_HitsNotUsedByTrack[istat]->Fill((muonShowerInformation.nStationHits).at(istat) - 2*StationHitsUsedByTrack(*ot, istat+1));
          //vs pT
          Muon_Station_showerT_pT[istat]->Fill(tk->pt(),(muonShowerInformation.stationShowerSizeT).at(istat));
          Muon_Station_deltaR_pT[istat]->Fill(tk->pt(),(muonShowerInformation.stationShowerDeltaR).at(istat));
          Muon_Station_AllHits_pT[istat]->Fill(tk->pt(),(muonShowerInformation.nStationHits).at(istat));
          Muon_Station_CorrHits_pT[istat]->Fill(tk->pt(),(muonShowerInformation.nStationCorrelatedHits).at(istat));
          Muon_Station_UncorrHits_pT[istat]->Fill(tk->pt(),(muonShowerInformation.nStationHits).at(istat) - (muonShowerInformation.nStationCorrelatedHits).at(istat));
          Muon_Station_HitsNotUsedByTrack_pT[istat]->Fill(tk->pt(),(muonShowerInformation.nStationHits).at(istat) - 2*StationHitsUsedByTrack(*ot, istat+1));
          //vs p
          Muon_Station_showerT_p[istat]->Fill(tk->p(),(muonShowerInformation.stationShowerSizeT).at(istat));
          Muon_Station_deltaR_p[istat]->Fill(tk->p(),(muonShowerInformation.stationShowerDeltaR).at(istat));
          Muon_Station_AllHits_p[istat]->Fill(tk->p(),(muonShowerInformation.nStationHits).at(istat));
          Muon_Station_CorrHits_p[istat]->Fill(tk->p(),(muonShowerInformation.nStationCorrelatedHits).at(istat));
          Muon_Station_UncorrHits_p[istat]->Fill(tk->p(),(muonShowerInformation.nStationHits).at(istat) - (muonShowerInformation.nStationCorrelatedHits).at(istat));
          Muon_Station_HitsNotUsedByTrack_p[istat]->Fill(tk->p(),(muonShowerInformation.nStationHits).at(istat) - 2*StationHitsUsedByTrack(*ot, istat+1));
          //vs eta
          Muon_Station_showerT_eta[istat]->Fill(tk->eta(),(muonShowerInformation.stationShowerSizeT).at(istat));
          Muon_Station_deltaR_eta[istat]->Fill(tk->eta(),(muonShowerInformation.stationShowerDeltaR).at(istat));
          Muon_Station_AllHits_eta[istat]->Fill(tk->eta(),(muonShowerInformation.nStationHits).at(istat));
          Muon_Station_CorrHits_eta[istat]->Fill(tk->eta(),(muonShowerInformation.nStationCorrelatedHits).at(istat));
          Muon_Station_UncorrHits_eta[istat]->Fill(tk->eta(),(muonShowerInformation.nStationHits).at(istat) - (muonShowerInformation.nStationCorrelatedHits).at(istat));
          Muon_Station_HitsNotUsedByTrack_eta[istat]->Fill(tk->eta(),(muonShowerInformation.nStationHits).at(istat) - 2*StationHitsUsedByTrack(*ot, istat+1));
        }//end by station
     }
  }
}

//
//Station hits usec by a track
//
int TestShower::StationHitsUsedByTrack(const reco::Track& muon, int station) {

    int hitsInThisStation = 0;
    //start numbering from 0
      for( trackingRecHit_iterator ihit = muon.recHitsBegin(); ihit != muon.recHitsEnd(); ihit++ ) {
          if( (*ihit)->isValid() ) {
              DetId id((*ihit)->geographicalId());
              if (id.subdetId() == MuonSubdetId::DT) {
                DTChamberId chamberId(id.rawId());
                if (chamberId.station() == station) hitsInThisStation++;
              } else if (id.subdetId() == MuonSubdetId::CSC) {
                CSCDetId did(id.rawId());
                if (did.station() == station) hitsInThisStation++;
              } else {
                //do nothing, RPC
              }
          }
     }
   
    return hitsInThisStation;
}

// ------------ method called once each job just before starting event loop  ------------
void 
TestShower::beginJob()
{
    edm::Service<TFileService> fs;
    std::ostringstream pprint;
    string histoName;

    histoName = "h_Muon_inner_pT";
    Muon_inner_pT = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),1000,0,1000);
    histoName = "h_Muon_inner_p";
    Muon_inner_p = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),1000,0,1000);
    histoName = "h_Muon_eta";
    Muon_eta = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),60,-3.0,3.0);

    histoName = "h_Muon_inner_pT_preselected";
    Muon_inner_pT_preselected= fs->make<TH1D>(histoName.c_str(),histoName.c_str(),1000,0,1000);
    histoName = "h_Muon_inner_p_preselected";
    Muon_inner_p_preselected = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),1000,0,1000);
    histoName = "h_Muon_eta_preselected";
    Muon_eta_preselected = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),60,-3.0,3.0);

    for (int istat = 0; istat < 4; istat++) {

       pprint.str("");
       pprint<<istat;
       string iname = pprint.str();

       histoName = "h_Muon_Station_inner_pT_prob5"+iname;
       Muon_Station_inner_pT_prob5[istat] = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),1000,0,1000);
       histoName = "h_Muon_Station_inner_p_prob5"+iname;
       Muon_Station_inner_p_prob5[istat] = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),1000,0,1000);
       histoName = "h_Muon_Station_eta_prob5"+iname;
       Muon_Station_eta_prob5[istat] = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),60,-3.0,3.0);

       histoName = "h_Muon_Station_showerT"+iname;
       Muon_Station_showerT[istat] = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),1000,0,500);
       histoName = "h_Muon_Station_deltaR"+iname;
       Muon_Station_deltaR[istat] = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),5000,0,0.5);
       histoName = "h_Muon_Station_AllHits"+iname;
       Muon_Station_AllHits[istat] = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),400,0,400);
       histoName = "h_Muon_Station_CorrHits"+iname;
       Muon_Station_CorrHits[istat] = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),400,0,400);
       histoName = "h_Muon_Station_UncorrHits"+iname;
       Muon_Station_UncorrHits[istat] = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),400,0,400);
       histoName = "h_Muon_Station_HitsNotUsedByTrack"+iname;
       Muon_Station_HitsNotUsedByTrack[istat] = fs->make<TH1D>(histoName.c_str(),histoName.c_str(),400,0,400);

       histoName = "h_Muon_Station_showerT_pT"+iname;
       Muon_Station_showerT_pT[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),1000,0,1000,1000,0,500);
       histoName = "h_Muon_Station_deltaR_pT"+iname;
       Muon_Station_deltaR_pT[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),1000,0,1000,5000,0,0.5);
       histoName = "h_Muon_Station_AllHits_pT"+iname;
       Muon_Station_AllHits_pT[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),1000,0,1000,400,0,400);
       histoName = "h_Muon_Station_CorrHits_pT"+iname;
       Muon_Station_CorrHits_pT[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),1000,0,1000,400,0,400);
       histoName = "h_Muon_Station_UncorrHits_pT"+iname;
       Muon_Station_UncorrHits_pT[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),1000,0,1000,400,0,400);
       histoName = "h_Muon_Station_HitsNotUsedByTrack_pT"+iname;
       Muon_Station_HitsNotUsedByTrack_pT[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),1000,0,1000,400,0,400);

       histoName = "h_Muon_Station_showerT_p"+iname;
       Muon_Station_showerT_p[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),1000,0,1000,1000,0,500);
       histoName = "h_Muon_Station_deltaR_p"+iname;
       Muon_Station_deltaR_p[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),1000,0,1000,5000,0,0.5);
       histoName = "h_Muon_Station_AllHits_p"+iname;
       Muon_Station_AllHits_p[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),1000,0,1000,400,0,400);
       histoName = "h_Muon_Station_CorrHits_p"+iname;
       Muon_Station_CorrHits_p[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),1000,0,1000,400,0,400);
       histoName = "h_Muon_Station_UncorrHits_p"+iname;
       Muon_Station_UncorrHits_p[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),1000,0,1000,400,0,400);
       histoName = "h_Muon_Station_HitsNotUsedByTrack_p"+iname;
       Muon_Station_HitsNotUsedByTrack_p[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),1000,0,1000,400,0,400);

       histoName = "h_Muon_Station_showerT_eta"+iname;
       Muon_Station_showerT_eta[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),60,-3,3,1000,0,500);
       histoName = "h_Muon_Station_deltaR_eta"+iname;
       Muon_Station_deltaR_eta[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),60,-3,3,5000,0,0.5);
       histoName = "h_Muon_Station_AllHits_eta"+iname;
       Muon_Station_AllHits_eta[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),60,-3,3,400,0,400);
       histoName = "h_Muon_Station_CorrHits_eta"+iname;
       Muon_Station_CorrHits_eta[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),60,-3,3,400,0,400);
       histoName = "h_Muon_Station_UncorrHits_eta"+iname;
       Muon_Station_UncorrHits_eta[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),60,-3,3,400,0,400);
       histoName = "h_Muon_Station_HitsNotUsedByTrack_eta"+iname;
       Muon_Station_HitsNotUsedByTrack_eta[istat] = fs->make<TH2D>(histoName.c_str(),histoName.c_str(),60,-3,3,400,0,400);

    }//end by station

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TestShower::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TestShower);
