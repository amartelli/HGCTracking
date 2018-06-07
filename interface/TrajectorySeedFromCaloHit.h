#ifndef RecoParticleFlow_HGCTracking_TrajectorySeedFromCaloHit_h
#define RecoParticleFlow_HGCTracking_TrajectorySeedFromCaloHit_h

/// Basic class  for a TrajectorySeed made from a reco::Track.
//  FIXME: not HGC specific, so should go in DataFormats/TrajectorySeed or similar

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingRecHit.h"


class TrajectorySeedFromCaloHit : public TrajectorySeed  {
 public:
 TrajectorySeedFromCaloHit(HGCTrackingRecHitFromCluster seedCl, int seedL, PTrajectoryStateOnDet const &ptsos, PropagationDirection dir):
  sclRef_(seedCl.objRef()), seedLayer_(seedL)
  {
    TrajectorySeed::recHitContainer vecHit; vecHit.push_back(seedCl); 
    //TrajectorySeed(ptsos, TrajectorySeed::recHitContainer(), dir);
    TrajectorySeed(ptsos, vecHit, dir);
  }
  
  const edm::Ptr<reco::CaloCluster> & clusterRef() const { return sclRef_; }
  
  int getSeedL() const {return seedLayer_;}

 protected:
  const edm::Ptr<reco::CaloCluster> sclRef_;
  int seedLayer_;
};

  
#endif
