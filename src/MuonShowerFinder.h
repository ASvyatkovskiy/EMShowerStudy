/**
 *  \class: MuonShowerFinder
 *
 *  Description: class for the identification of EM showers in the muon system
 *
 *  $Date: 2011/01/10 01:39:24 $
 *  $Revision: 1.9 $
 *
 *  \author : A. Svyatkovskiy, Purdue University
 *
 **/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"

namespace edm {class ParameterSet; class Event; class EventSetup;}
namespace reco {class TransientTrack;}

class MuonServiceProxy;
class Trajectory;
class Cylinder;
class BoundDisk;
class BarrelDetLayer;
class ForwardDetLayer;
class TransientTrackingRecHitBuilder;
class GeometricSearchTracker;
class GlobalTrackingGeometry;
class MuonDetLayerGeometry;

//
// class declaration
//

class MuonShowerFinder {

  public:

    typedef TransientTrackingRecHit::ConstRecHitContainer ConstRecHitContainer;
    typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;

  public:

    MuonShowerFinder() {};
    MuonShowerFinder(const MuonServiceProxy*, const edm::ParameterSet&);
    ~MuonShowerFinder();

    /// pass the Event to the algorithm at each event
    virtual void setEvent(const edm::Event&);

    /// set the services needed by the TrackTransformer
    void setServices(const edm::EventSetup&);

    MuonRecHitContainer findPhiCluster(MuonRecHitContainer&, const GlobalPoint&) const;
    MuonRecHitContainer findThetaCluster(MuonRecHitContainer&, const GlobalPoint&) const;
    MuonRecHitContainer recHits4D(const GeomDet*,edm::Handle<DTRecSegment4DCollection>, edm::Handle<CSCSegmentCollection>) const;
    int numberOfCorrelatedHits(const MuonRecHitContainer&) const;
    std::vector<double> stationShowerSizes(const reco::Track&) const;
    std::vector<double> stationShowerRSizes(const reco::Track&) const;
    std::vector<int> stationUncorrelatedHits(const reco::Track&) const;
    std::vector<const GeomDet*> getCompatibleDets(const reco::Track&) const;
    std::vector<bool> hasStationShower(std::vector<int>& uncorrelatedHits, std::vector<double>& showerSizes) const;


   protected:

     const MuonServiceProxy* getService() const { return theService; }

   private:

     const MuonServiceProxy* theService;

//     GlobalPoint getClusteringRefpoint(MuonRecHitContainer&, const reco::Track&) const;
     GlobalPoint crossingPoint(const GlobalPoint&, const GlobalPoint&, const BarrelDetLayer* ) const;
     GlobalPoint crossingPoint(const GlobalPoint&, const GlobalPoint&, const Cylinder& ) const;
     GlobalPoint crossingPoint(const GlobalPoint&, const GlobalPoint&, const ForwardDetLayer* ) const;
     GlobalPoint crossingPoint(const GlobalPoint&, const GlobalPoint&, const BoundDisk& ) const;
     std::vector<const GeomDet*> dtPositionToDets(const GlobalPoint&) const;
     std::vector<const GeomDet*> cscPositionToDets(const GlobalPoint&) const;
     std::vector<MuonRecHitContainer> fillHitsByStation(const reco::Track&) const;
     MuonRecHitContainer findPerpCluster(MuonRecHitContainer& muonRecHits) const;

     struct LessMag {
       LessMag(const GlobalPoint& point) : thePoint(point) {}
       bool operator()(const GlobalPoint& lhs,
		       const GlobalPoint& rhs) const{ 
	    return (lhs - thePoint).mag() < (rhs -thePoint).mag();
	}
        bool operator()(const MuonTransientTrackingRecHit::MuonRecHitPointer& lhs,
                       const MuonTransientTrackingRecHit::MuonRecHitPointer& rhs) const{
           return (lhs->globalPosition() - thePoint).mag() < (rhs->globalPosition() -thePoint).mag();
        }
	GlobalPoint thePoint;
      };

      struct LessDPhi {
        LessDPhi(const GlobalPoint& point) : thePoint(point) {}
        bool operator()(const MuonTransientTrackingRecHit::MuonRecHitPointer& lhs,
                       const MuonTransientTrackingRecHit::MuonRecHitPointer& rhs) const{
           return deltaPhi(lhs->globalPosition().phi(), thePoint.phi()) < deltaPhi(rhs->globalPosition().phi(), thePoint.phi());
        }

        GlobalPoint thePoint;
      };

      struct AbsLessDPhi {
        AbsLessDPhi(const GlobalPoint& point) : thePoint(point) {}
        bool operator()(const MuonTransientTrackingRecHit::MuonRecHitPointer& lhs,
                       const MuonTransientTrackingRecHit::MuonRecHitPointer& rhs) const{
           return ( fabs(deltaPhi(lhs->globalPosition().phi(), thePoint.phi())) < fabs(deltaPhi(rhs->globalPosition().phi(), thePoint.phi())) );
        }

        GlobalPoint thePoint;

      };

      struct AbsLessDTheta {
        AbsLessDTheta(const GlobalPoint& point) : thePoint(point) {}
        bool operator()(const MuonTransientTrackingRecHit::MuonRecHitPointer& lhs,
                       const MuonTransientTrackingRecHit::MuonRecHitPointer& rhs) const{
           return ( fabs(lhs->globalPosition().phi() - thePoint.phi()) < fabs(rhs->globalPosition().phi() - thePoint.phi()) );
        }

        GlobalPoint thePoint;

      };

      //sort by phi
      struct LessPhi {
        LessPhi() : thePoint(0,0,0) {}
        bool operator()(const MuonTransientTrackingRecHit::MuonRecHitPointer& lhs,
                       const MuonTransientTrackingRecHit::MuonRecHitPointer& rhs) const{
           return (lhs->globalPosition().phi() < rhs->globalPosition().phi());
        }

        GlobalPoint thePoint;
      };


      // sort by perp
      struct LessPerp {
        LessPerp() : thePoint(0,0,0) {}
        bool operator()(const MuonTransientTrackingRecHit::MuonRecHitPointer& lhs,
                       const MuonTransientTrackingRecHit::MuonRecHitPointer& rhs) const{
           return (lhs->globalPosition().perp() < rhs->globalPosition().perp());
        }

        GlobalPoint thePoint;
      };
 
      //sort 
      struct LessAbsMag {
        LessAbsMag() : thePoint(0,0,0) {}
        bool operator()(const MuonTransientTrackingRecHit::MuonRecHitPointer& lhs,
                       const MuonTransientTrackingRecHit::MuonRecHitPointer& rhs) const{
           return (lhs->globalPosition().mag() < rhs->globalPosition().mag());
        }

        GlobalPoint thePoint;
      };

      std::string category_;
      double sizeThreshold1_, sizeThreshold2_;
      int hitThreshold1_, hitThreshold2_;

      unsigned long long theCacheId_TRH;
      unsigned long long theCacheId_MT;

      std::string theTrackerRecHitBuilderName;
      edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;

      std::string theMuonRecHitBuilderName;
      edm::ESHandle<TransientTrackingRecHitBuilder> theMuonRecHitBuilder;

      edm::InputTag theDTRecHitLabel;
      edm::InputTag theCSCSegmentsLabel;
      edm::InputTag theCSCRecHitLabel;
      edm::InputTag theDT4DRecSegmentLabel;
      edm::Handle<DTRecHitCollection> theDTRecHits;
      edm::Handle<CSCSegmentCollection> theCSCSegments;
      edm::Handle<DTRecSegment4DCollection> theDT4DRecSegments;
      edm::Handle<CSCRecHit2DCollection> theCSCRecHits;

      // geometry
      edm::ESHandle<GeometricSearchTracker> theTracker;
      edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
      edm::ESHandle<MagneticField> theField;
      edm::ESHandle<MuonDetLayerGeometry> theMuonGeometry;

};
