#ifndef RecoHGCal_HGCTracking_HGCTrackingClusteringRecHit_h
#define RecoHGCal_HGCTracking_HGCTrackingClusteringRecHit_h

/// Class for a TrackingRecHit made from an on-the-fly clustering of HGCRecHits

#include "DataFormats/TrackingRecHit/interface/RecHit2DLocalPos.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

class HGCTrackingClusteringRecHit : public RecHit2DLocalPos {
    public:
        typedef HGCRecHitRef ObjRef;
        typedef edm::RefVector<ObjRef::product_type, ObjRef::value_type, ObjRef::finder_type> ObjRefVector;

        HGCTrackingClusteringRecHit() {}
        HGCTrackingClusteringRecHit(DetId id, const ObjRefVector &objrefs, float etot, const LocalPoint &pos, const LocalError &err):
            RecHit2DLocalPos(id),
            objrefs_(objrefs), etot_(etot),
            pos_(pos), err_(err) {}

        virtual HGCTrackingClusteringRecHit * clone() const override { return new HGCTrackingClusteringRecHit(*this); }

        virtual LocalPoint localPosition() const override { return pos_; }
        virtual LocalError localPositionError() const override { return err_; }

        ObjRef objRef() const { return objrefs_[0]; }
        const ObjRefVector & objRefs() const { return objrefs_; }

        float energy() const { return etot_; }

        virtual bool sharesInput( const TrackingRecHit* other, SharedInputType what) const override ;

    protected:
        ObjRefVector objrefs_;
        float etot_;
        LocalPoint pos_;
        LocalError err_;
};

#endif

