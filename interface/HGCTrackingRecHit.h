#ifndef RecoHGCal_HGCTracking_HGCTrackingRecHit_h
#define RecoHGCal_HGCTracking_HGCTrackingRecHit_h

/// Basic template class for a RecHit wrapping a Ref to an object

#include <cassert>
#include "DataFormats/TrackingRecHit/interface/RecHit2DLocalPos.h"

template<typename ObjRef>
class HGCTrackingRecHit : public RecHit2DLocalPos {
    public:
        typedef ObjRef ref_type;
        typedef typename ObjRef::value_type obj_type;

        HGCTrackingRecHit() {}
        HGCTrackingRecHit(DetId id, const ObjRef &objref, const LocalPoint &pos, const LocalError &err):
            RecHit2DLocalPos(id),
            objref_(objref),
            pos_(pos), err_(err) {}

        virtual HGCTrackingRecHit<ObjRef> * clone() const override { return new HGCTrackingRecHit<ObjRef>(*this); }

        virtual LocalPoint localPosition() const override { return pos_; }
        virtual LocalError localPositionError() const override { return err_; }

        const ObjRef & objRef() const { return objref_; }

        float energy() const { return objref_->energy(); }

        virtual bool sharesInput( const TrackingRecHit* other, SharedInputType what) const override { assert(false); }
    protected:
        ObjRef objref_;
        LocalPoint pos_;
        LocalError err_;
};


// Instantiations and specializations for HGCRecHitRef and reco::CaloClusterPtr
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
typedef HGCTrackingRecHit<HGCRecHitRef> HGCTrackingRecHitFromHit;
typedef HGCTrackingRecHit<reco::CaloClusterPtr> HGCTrackingRecHitFromCluster;

template<>
bool HGCTrackingRecHit<HGCRecHitRef>::sharesInput( const TrackingRecHit* other, SharedInputType what) const ;

template<>
bool HGCTrackingRecHit<reco::CaloClusterPtr>::sharesInput( const TrackingRecHit* other, SharedInputType what) const ;

#endif
