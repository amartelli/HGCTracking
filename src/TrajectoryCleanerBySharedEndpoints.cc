#include "RecoHGCal/HGCTracking/interface/TrajectoryCleanerBySharedEndpoints.h"

void TrajectoryCleanerBySharedEndpoints::clean(std::vector<Trajectory> &trajs) const 
{ 
    clean_(trajs); 
}

void TrajectoryCleanerBySharedEndpoints::clean(std::vector<TempTrajectory> &trajs) const 
{ 
    clean_(trajs); 
}

template<typename Traj>
void TrajectoryCleanerBySharedEndpoints::clean_(std::vector<Traj> &trajs) const 
{
    // stupid O(n^2) algo
    for (unsigned int i = 0, nt = trajs.size(); i < nt; ++i) {
        Traj &t1 = trajs[i]; if (!t1.isValid() || t1.foundHits() < 2) continue;
        //const auto & first1 = t1.firstMeasurement();
        auto last1 = t1.measurements().rbegin();
        while (!last1->recHit()->isValid()) ++last1;
        for (unsigned int j = i+1; j < nt; ++j) {
            Traj &t2 = trajs[j]; if (!t2.isValid() || t2.foundHits() < 2) continue;
            //const auto & first2 = t2.firstMeasurement();
            auto last2 = t2.measurements().rbegin();
            while (!last2->recHit()->isValid()) ++last2;
            if (last2->recHit()->sharesInput(&*last1->recHit(), TrackingRecHit::all)) {
                //printf("\t\t\tfound two trajectories sharing last hit %12d\n", last1->recHit()->geographicalId()());
                //++last1; ++last2;
                //while (!last1->recHit()->isValid()) ++last1;
                //while (!last2->recHit()->isValid()) ++last2;
                //if (last2->recHit()->sharesInput(&*last1->recHit(), TrackingRecHit::all)) {
                //    printf("\t\t\t\tthey also share next-to-last hit %12d --> will clean them\n", last1->recHit()->geographicalId()());
                if (cmp(t1,t2)) t2.invalidate(); else t1.invalidate();
                //} else if (first2.recHit()->sharesInput(&*first1.recHit(), TrackingRecHit::all)) {
                //    printf("\t\t\t\tthey also share the first hit %12d --> will clean them\n", first1.recHit()->geographicalId()());
                //    if (cmp(t1,t2)) t2.invalidate(); else t1.invalidate();
                //}
            }
        }
    }  
}


