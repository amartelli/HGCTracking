#ifndef RecoHGCal_HGCTracking_HGCDiskGeomDet_h
#define RecoHGCal_HGCTracking_HGCDiskGeomDet_h

/// Class corresponding to one layer of HGC

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"

class HGCDiskGeomDet : public GeomDet {
    public:
        HGCDiskGeomDet(int subdet, int zside, int layer, float z, float rmin, float rmax, float radlen, float xi) ;

        int subdet() const { return subdet_; }
        int zside() const { return zside_; }
        int layer() const { return layer_; }

        bool operator<(const HGCDiskGeomDet &other) const { 
            return (subdet_ == other.subdet_ ? layer_ < other.layer_ : subdet_ < other.subdet_);
        }

    protected:
        const int subdet_, zside_, layer_;
};

#endif
