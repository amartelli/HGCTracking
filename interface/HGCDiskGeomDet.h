#ifndef RecoParticleFlow_HGCTracking_HGCDiskGeomDet_h
#define RecoParticleFlow_HGCTracking_HGCDiskGeomDet_h

/// Class corresponding to one layer of HGC

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

class HGCDiskGeomDet : public GeomDet {
    public:
        HGCDiskGeomDet(int subdet, int zside, int layer, float z, float rmin, float rmax, float radlen, float xi) ;

        int subdet() const { return subdet_; }
        int zside() const { return zside_; }
        int layer() const { return layer_; }

	int progressiveLayer() const {
	  int offset = (subdet_ > 3) ? 28 : 0;
	  offset += (subdet_ > 4) ? 12 : 0;
	  return layer_ + offset;
	}

        bool operator<(const HGCDiskGeomDet &other) const { 
            return (subdet_ == other.subdet_ ? layer_ < other.layer_ : subdet_ < other.subdet_);
        }

    protected:
        const int subdet_, zside_, layer_;
};

#endif
