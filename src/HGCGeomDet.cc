#include "RecoHGCal/HGCTracking/interface/HGCDiskGeomDet.h"
#include "DataFormats/GeometrySurface/interface/MediumProperties.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

HGCDiskGeomDet::HGCDiskGeomDet(int subdet, int zside, int layer, float z, float rmin, float rmax, float radlen, float xi) :
            GeomDet( Disk::build(Disk::PositionType(0,0,z), Disk::RotationType(), SimpleDiskBounds(rmin, rmax, -20, 20)).get() ),
            subdet_(subdet), zside_(zside), layer_(layer) 
{
    // if (subdet == 3 || subdet == 4) {
    //     setDetId(HGCalDetId(ForwardSubdetector(subdet),zside,layer,0,0,0));
    // } else {
    //     setDetId(HcalDetId(HcalSubdetector::HcalEndcap, zside>0 ? +30 : -30, 0, layer));
    // }
    if(subdet == 8 || subdet == 9) setDetId(HGCSiliconDetId(DetId::Detector(subdet), z, 0, layer, 0, 0, 0, 0));
    if (radlen > 0) {
        (const_cast<Plane &>(surface())).setMediumProperties(MediumProperties(radlen,xi));
    }
}

