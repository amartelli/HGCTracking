#ifndef RecoHGCal_HGCTracking_interface_HGCTrackingDiskData_h
#define RecoHGCal_HGCTracking_interface_HGCTrackingDiskData_h

/// Class that holds the rechits on one HGC Layer, and does the lookup
//  Almost like a MeasurementDetWithData in the tracker
//
//  note 1: the class is typedeffed on the rechit collection, but is not a template
//  note 2: the class is cheap to move, but expensive to copy


#include "RecoHGCal/HGCTracking/interface/HGCTrackingRecHit.h"
#include "RecoHGCal/HGCTracking/interface/HGCTrackingClusteringRecHit.h"
#include "RecoHGCal/HGCTracking/interface/HGCTracker.h"
#include "RecoHGCal/HGCTracking/interface/HGCTrackingBasicCPE.h"

#include "TrackingTools/DetLayers/interface/MeasurementEstimator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"

#include "RecoHGCal/HGCTracking/interface/hgcdebug.h"


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

        std::vector<TrajectoryMeasurement> measurements(TrajectoryStateOnSurface &tsos, const MeasurementEstimator &mest)  const ;
        std::vector<TrajectoryMeasurement> clusterizedMeasurements(TrajectoryStateOnSurface &tsos, const MeasurementEstimator &mest, float rCut) const ;
        std::vector<TrajectoryMeasurement> clusterMeasurements(TrajectoryStateOnSurface &tsos, const MeasurementEstimator &mest) const ;

        // for debugging
        void setTruth(const CaloTruthRevMap *truth) const { truthMap_ = truth; }

    private:
        const edm::Handle<TColl> *alldata_;
        std::vector<std::pair<HGCTrackingBasicCPE::PositionHint,const_iterator>> index_;

        const edm::Handle<reco::CaloClusterCollection> *clusterData_;
        std::vector<std::pair<HGCTrackingBasicCPE::PositionHint,reco::CaloClusterCollection::const_iterator>> clusterIndex_;

        const HGCTrackingBasicCPE *cpe_;

        // for debugging
        mutable const CaloTruthRevMap *truthMap_;

        void buildIndex_() {
            // could do some sorting here
        }
    
};

#endif
