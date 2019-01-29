#ifndef RecoHGCal_HGCTracking_HGCTkTrajectoryBuilder_h
#define RecoHGCal_HGCTracking_HGCTkTrajectoryBuilder_h

/// Class that does the trajectory building starting from a seed

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"
#include "TrackingTools/PatternTools/interface/TempTrajectory.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"
#include "TrackingTools/TrajectoryCleaning/interface/TrajectoryCleaner.h"
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"

#include "RecoHGCal/HGCTracking/interface/HGCTrackingRecHit.h"
#include "RecoHGCal/HGCTracking/interface/HGCTrackingClusteringRecHit.h"
#include "RecoHGCal/HGCTracking/interface/HGCTracker.h"
#include "RecoHGCal/HGCTracking/interface/HGCTrackingBasicCPE.h"
#include "RecoHGCal/HGCTracking/interface/TrajectoryCleanerBySharedEndpoints.h"
#include "RecoHGCal/HGCTracking/interface/HGCTrackingData.h"
#include "RecoHGCal/HGCTracking/interface/hgcdebug.h"

class HGCTkTrajectoryBuilder {
    public:
        HGCTkTrajectoryBuilder(const edm::ParameterSet& ps, edm::ConsumesCollector && c );

        // no copy, no move, no assign (they would be all bad ideas)
        HGCTkTrajectoryBuilder(const HGCTkTrajectoryBuilder &other) = delete;
        HGCTkTrajectoryBuilder(HGCTkTrajectoryBuilder && other) = delete;
        HGCTkTrajectoryBuilder& operator=(const HGCTkTrajectoryBuilder &other) = delete;
        HGCTkTrajectoryBuilder& operator=(HGCTkTrajectoryBuilder && other) = delete;

        /// to be called at the beginning of the event
        void init(const edm::Event &event, const edm::EventSetup &iSetup) ;
        /// to be called at the end of the event
        void done() ;

        /// Build a trajectory from a track reference
        unsigned int trajectories(const reco::TrackRef &track, std::vector<Trajectory> &out, PropagationDirection direction=alongMomentum) const ; 
        unsigned int trajectories(const FreeTrajectoryState &start, std::vector<Trajectory> &out, PropagationDirection direction=alongMomentum) const ; 

        /// Ambiguity resolution
        void cleanTrajectories(std::vector<Trajectory> &trajectories) const ;

        /// Backwards fit
        Trajectory bwrefit(const Trajectory &traj, float scaleErrors = 100.) const ;

        //// ---- for debugging ----
        // mc truth association map 
        void setTruth(const CaloTruthRevMap *truth) { truthMap_ = truth; }
        // we keep the template private, and define public
        void printTraj(const Trajectory &traj) const { printTraj_(traj); }
        void printTraj(const TempTrajectory &traj) const { printTraj_(traj); }
        std::vector<std::pair<const CaloParticle *, float>> truthMatch(const Trajectory &traj) const { return truthMatch_(traj); }
        std::vector<std::pair<const CaloParticle *, float>> truthMatch(const TempTrajectory &traj) const { return truthMatch_(traj); }

    protected:
        enum PatternRecoAlgo { SingleHitAlgo, ClusterizingAlgo, MixedAlgo };

        // --- Inputs ---
        const edm::EDGetTokenT<HGCRecHitCollection> srcEE_, srcFH_;  //, srcBH_;
        const edm::EDGetTokenT<reco::CaloClusterCollection> srcClusters_;

        // --- Configuration ---
        PatternRecoAlgo algo_;
        const std::string propName_, propNameOppo_;
        const std::string estimatorName_;
        const std::string updatorName_;
        const std::string trajectoryCleanerName_;

        const double lostHitRescaleErrorFactor_;
        const bool fastCleaner_, endpointCleaner_;
        const unsigned int theMaxCand;
        const double theFoundHitBonus, theLostHitPenalty;
        const unsigned int theMaxStartingEmptyLayers, theLayersBeforeCleaning;
        const bool bestHitOnly_;
        const double deltaChiSquareForHits_, minChi2ForInvalidHit_;
        const double clusterRadius_;
        const bool lostHitsOnBH_;

        // --- Event Setup stuff (or similar) ---
        uint32_t geomCacheId_;
        std::unique_ptr<HGCTracker> hgcTracker_;
        edm::ESHandle<CaloGeometry> caloGeom_;
        edm::ESHandle<GlobalTrackingGeometry> trkGeom_;
        edm::ESHandle<Propagator> prop_, propOppo_;
        edm::ESHandle<MagneticField> bfield_;
        edm::ESHandle<Chi2MeasurementEstimatorBase> estimator_;
        edm::ESHandle<TrajectoryStateUpdator> updator_;
        edm::ESHandle<TrajectoryCleaner> trajectoryCleaner_;
        std::unique_ptr<TrajectoryFilter> trajFilter_;
        std::unique_ptr<HGCTrackingBasicCPE> cpe_;

        // --- Event Data ---
        // note that handles have to be kept, since they're needed to build refs
        edm::Handle<HGCRecHitCollection> srcEE, srcFH; //, srcBH;
        edm::Handle<reco::CaloClusterCollection> srcClusters;
        // data re-arranged for tracking
        std::unique_ptr<HGCTrackingData> data_;
        // for debug
        const CaloTruthRevMap *truthMap_;

        /// --- Pattern reco ---
        template<class Start>
        std::vector<TempTrajectory> advanceOneLayer(const Start &start, const TempTrajectory &traj, const HGCDiskGeomDet *disk, bool bestHitOnly) const ;

        /// --- Utilities ---
        void trim(std::vector<TempTrajectory> &tempTrajectories) const {
            tempTrajectories.erase(std::remove_if( tempTrajectories.begin(),tempTrajectories.end(), std::not1(std::mem_fun_ref(&TempTrajectory::isValid))), tempTrajectories.end());
        }
        void trim(std::vector<Trajectory> &trajectories) const {
            trajectories.erase(std::remove_if( trajectories.begin(),trajectories.end(), std::not1(std::mem_fun_ref(&Trajectory::isValid))), trajectories.end());
        }

        /// --- Debug ---
        template<typename Traj> void printTraj_(const Traj &t) const ;
        template<typename Traj> std::vector<std::pair<const CaloParticle *, float>> truthMatch_(const Traj &t) const ;
};

#endif
