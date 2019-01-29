#include "RecoHGCal/HGCTracking/interface/HGCTrackingRecHit.h"
#include "RecoHGCal/HGCTracking/interface/HGCTrackingClusteringRecHit.h"

bool HGCTrackingClusteringRecHit::sharesInput( const TrackingRecHit* other, SharedInputType what) const 
{
    if (typeid(*other) == typeid(HGCTrackingClusteringRecHit)) {
        const auto & refv = (static_cast<const HGCTrackingClusteringRecHit *>(other))->objRefs();
        if (refv.id() != objrefs_.id()) return false; // unlikely
        const auto &keys1 = objrefs_.refVector().keys(), &keys2 = refv.refVector().keys();
        if (keys1.front() == keys2.front()) return true; // same seed
        if (what == TrackingRecHit::some) {
            return std::find_first_of(keys1.begin(), keys1.end(), keys2.begin(), keys2.end()) != keys1.end();
        } else { // require each to contain the seed of the other
            return (std::find(keys1.begin(), keys1.end(), keys2.front()) != keys1.end()) &&
                (std::find(keys2.begin(), keys2.end(), keys1.front()) != keys2.end());
        }
    } else if (typeid(*other) == typeid(HGCTrackingRecHit<ObjRef>)) {
        const auto & ref = (static_cast<const HGCTrackingRecHit<ObjRef> *>(other))->objRef();
        if (what == TrackingRecHit::some) {
            return std::find(objrefs_.begin(), objrefs_.end(), ref) != objrefs_.end();
        } else {
            return objRef() == ref; 
        }
    } else if (typeid(*other) == typeid(HGCTrackingRecHit<reco::CaloClusterPtr>)) {
        return other->sharesInput(this, what); // logic is implemeted there
    } else {
        return false;
    }
}

