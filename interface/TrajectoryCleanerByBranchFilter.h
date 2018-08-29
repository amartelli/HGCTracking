#ifndef RecoParticleFlow_HGCTracking_TrajectoryCleanerByBranchFilter_h
#define RecoParticleFlow_HGCTracking_TrajectoryCleanerByBranchFilter_h

/// Trajectory cleaner that merges trajectories based on sharing of measurements

#include <vector>
#include <set>
#include <map>

#include "TrackingTools/PatternTools/interface/Trajectory.h"
//#include "TrackingTools/PatternTools/interface/TempTrajectory.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingRecHit.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingClusteringRecHit.h"


class TrajectoryCleanerByBranchFilter {
 public:

 TrajectoryCleanerByBranchFilter() :
  theFoundHitBonus(0.), theLostHitPenalty(0.) {}
  
 TrajectoryCleanerByBranchFilter(float foundHitBonus, float lostHitPenalty) :
  theFoundHitBonus(foundHitBonus), theLostHitPenalty(lostHitPenalty) {}
  
  std::map<unsigned int, std::vector<TrajectoryMeasurement> > mergeTraj(const std::vector<Trajectory> &trajs);
  
 private:
  float theFoundHitBonus, theLostHitPenalty;
  
};

#endif
