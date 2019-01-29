#ifndef RecoHGCal_HGCTracking_HGCTracker_h
#define RecoHGCal_HGCTracking_HGCTracker_h

/// Class holding all the layers of the HGC Tracker
/// With some more effort it could inherit from TrackingGeometry etc
// - holds (and owns) all the disks
// - implements navigation and lookup by id

#include "RecoHGCal/HGCTracking/interface/HGCDiskGeomDet.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

class HGCTracker {
    public:
        HGCTracker(const CaloGeometry *geom) ;
        ~HGCTracker() { 
            for (HGCDiskGeomDet * disk : disksPos_) delete disk;
            for (HGCDiskGeomDet * disk : disksNeg_) delete disk;
        }

        // disable copy constructor and assignment
        HGCTracker(const HGCTracker &other) = delete;
        HGCTracker & operator=(const HGCTracker &other) = delete;

        unsigned int numberOfDisks(int zside) const { return (zside > 0 ? disksPos_ : disksNeg_).size(); }
        const HGCDiskGeomDet * disk(int zside, int disk) const { return (zside > 0 ? disksPos_ : disksNeg_).at(disk); }
        const HGCDiskGeomDet * firstDisk(int zside) const { return (zside > 0 ? disksPos_ : disksNeg_).front(); }
        const HGCDiskGeomDet * lastDisk(int zside) const { return (zside > 0 ? disksPos_ : disksNeg_).back(); }
        const HGCDiskGeomDet * firstDisk(int zside, PropagationDirection direction) const { return direction == alongMomentum ? firstDisk(zside) : lastDisk(zside); }
        const HGCDiskGeomDet * lastDisk(int zside, PropagationDirection direction) const { return direction == alongMomentum ? lastDisk(zside) : firstDisk(zside); }

        // FIXME: interface to be extended for outside-in 
        const HGCDiskGeomDet * nextDisk(const HGCDiskGeomDet *from, PropagationDirection direction=alongMomentum) const ;

        const HGCDiskGeomDet * idToDet(DetId id) const ;

    protected:
        const CaloGeometry *geom_;
        std::vector<HGCDiskGeomDet *> disksPos_, disksNeg_;

        void makeDisks(int subdet, int disks) ;
        void addDisk(HGCDiskGeomDet *disk) { 
            (disk->zside() > 0 ? disksPos_ : disksNeg_).push_back(disk);
        }
};

#endif
