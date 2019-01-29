#include "RecoHGCal/HGCTracking/interface/HGCTrackingRecHit.h"
#include "RecoHGCal/HGCTracking/interface/HGCTrackingClusteringRecHit.h"

template<>
bool HGCTrackingRecHit<HGCRecHitRef>::sharesInput( const TrackingRecHit* other, SharedInputType what) const 
{
    if (typeid(*other) == typeid(HGCTrackingRecHit<HGCRecHitRef>)) {
        return objRef() == (static_cast<const HGCTrackingRecHit<HGCRecHitRef> *>(other))->objRef();
    } else if (typeid(*other) == typeid(HGCTrackingClusteringRecHit)) {
        if (what == TrackingRecHit::some) {
            const auto &refv = (static_cast<const HGCTrackingClusteringRecHit *>(other))->objRefs();
            return std::find(refv.begin(), refv.end(), objRef()) != refv.end();
        } else {
            return objRef() == (static_cast<const HGCTrackingClusteringRecHit *>(other))->objRef();
        }
    } else if (typeid(*other) == typeid(HGCTrackingRecHit<reco::CaloClusterPtr>)) {
        return other->sharesInput( this, what ); // see below
    } else {
        return false;
    }
}

template<>
bool HGCTrackingRecHit<reco::CaloClusterPtr>::sharesInput( const TrackingRecHit* other, SharedInputType what) const 
{
    if (typeid(*other) == typeid(HGCTrackingRecHit<reco::CaloClusterPtr>)) {
        return objRef() == (static_cast<const HGCTrackingRecHit<reco::CaloClusterPtr> *>(other))->objRef();
    } else if (typeid(*other) == typeid(HGCTrackingRecHit<HGCRecHitRef>)) {
        DetId otherId = (static_cast<const HGCTrackingRecHit<HGCRecHitRef> *>(other))->objRef()->id();
        for (const auto & pair : objref_->hitsAndFractions()) {
            if (pair.second == 0) continue;
            if (pair.first == otherId) return true;
        }
        return false;
    } else if (typeid(*other) == typeid(HGCTrackingClusteringRecHit)) {
        const auto &refv = (static_cast<const HGCTrackingClusteringRecHit *>(other))->objRefs();
        for (HGCRecHitRef ref : refv)  {
            DetId otherId = ref->id();
            for (const auto & pair : objref_->hitsAndFractions()) {
                if (pair.second == 0) continue;
                if (pair.first == otherId) return true;
            }
        }
        return false;
    } else {
        return false;
    }
}


