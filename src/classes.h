#include "RecoHGCal/HGCTracking/interface/HGCTrackingRecHit.h"
#include "RecoHGCal/HGCTracking/interface/HGCTrackingClusteringRecHit.h"
#include "RecoHGCal/HGCTracking/interface/TrajectorySeedFromTrack.h"

namespace RecoHGCal_HGCTracking {
    struct dictionary {
        HGCTrackingRecHitFromHit dummy_hit;
        HGCTrackingRecHitFromCluster dummy_hit_cluster;
        HGCTrackingClusteringRecHit dummy_clustering_hit; 
    };
}
