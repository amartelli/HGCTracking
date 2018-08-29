#ifndef RecoParticleFlow_HGCTracking_interface_HGCTrackingDiskData_h
#define RecoParticleFlow_HGCTracking_interface_HGCTrackingDiskData_h

/// Class that holds the rechits on one HGC Layer, and does the lookup
//  Almost like a MeasurementDetWithData in the tracker
//
//  note 1: the class is typedeffed on the rechit collection, but is not a template
//  note 2: the class is cheap to move, but expensive to copy


#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingRecHit.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingClusteringRecHit.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTracker.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingBasicCPE.h"

#include "TrackingTools/DetLayers/interface/MeasurementEstimator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"

#include "RecoParticleFlow/HGCTracking/interface/hgcdebug.h"


class HGCTrackingDiskData {
    public:
        typedef HGCRecHitCollection TColl;
        typedef TColl::value_type value_type;
        typedef TColl::const_iterator const_iterator;
        typedef edm::Ref<TColl> ref_type;
        typedef edm::RefVector<TColl> refvector_type;
        typedef HGCTrackingRecHit<ref_type> rechit_type;
        typedef HGCTrackingClusteringRecHit rechitcluster_type;

        HGCTrackingDiskData() {}
        HGCTrackingDiskData(const edm::Handle<TColl> &data, int subdet, int zside, int layer, const HGCTrackingBasicCPE *cpe) ;

        void addClusters(const edm::Handle<reco::CaloClusterCollection> &data, int subdet, int zside, int layer) ;

        unsigned int size() const { return index_.size(); }
        unsigned int nclusters() const { return clusterIndex_.size(); }

        std::vector<TrajectoryMeasurement> measurements(TrajectoryStateOnSurface &tsos, const MeasurementEstimator &mest, bool cutSoN=false)  const ;
        std::vector<TrajectoryMeasurement> clusterizedMeasurements(TrajectoryStateOnSurface &tsos, const MeasurementEstimator &mest, float rCut) const ;
        std::vector<TrajectoryMeasurement> clusterMeasurements(TrajectoryStateOnSurface &tsos, const MeasurementEstimator &mest) const ;

        // for debugging
        void setTruth(const CaloTruthRevMap *truth) const { truthMap_ = truth; }
        void setTruthBis(const std::map<uint32_t, float > *truthX, const std::map<uint32_t, float > *truthY, const std::map<uint32_t, float > *truthZ) const {
	  truthMapX_ = truthX; truthMapY_ = truthY; truthMapZ_ = truthZ;}

    private:
        const edm::Handle<TColl> *alldata_;
        std::vector<std::pair<HGCTrackingBasicCPE::PositionHint,const_iterator>> index_;

        const edm::Handle<reco::CaloClusterCollection> *clusterData_;
        std::vector<std::pair<HGCTrackingBasicCPE::PositionHint,reco::CaloClusterCollection::const_iterator>> clusterIndex_;

        const HGCTrackingBasicCPE *cpe_;

        // for debugging
        mutable const CaloTruthRevMap *truthMap_;
        mutable const std::map<uint32_t, float > *truthMapX_;
        mutable const std::map<uint32_t, float > *truthMapY_;
        mutable const std::map<uint32_t, float > *truthMapZ_;

	void buildIndex_() {
            // could do some sorting here
        }
    
};

#endif
