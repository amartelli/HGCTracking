#ifndef RecoHGCal_HGCTracking_TrajectorySeedFromTrack_h
#define RecoHGCal_HGCTracking_TrajectorySeedFromTrack_h

/// Basic class  for a TrajectorySeed made from a reco::Track.
//  FIXME: not HGC specific, so should go in DataFormats/TrajectorySeed or similar

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"

class TrajectorySeedFromTrack : public TrajectorySeed  {
    public:
        TrajectorySeedFromTrack(const reco::TrackRef &tkRef, PTrajectoryStateOnDet const &ptsos, PropagationDirection dir) :
             TrajectorySeed(ptsos,TrajectorySeed::recHitContainer(),dir),
             trkRef_(tkRef) {}
        const reco::TrackRef & track() const { return trkRef_; }
    protected:
        reco::TrackRef trkRef_;
};


#endif
