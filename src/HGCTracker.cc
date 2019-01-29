#include "RecoHGCal/HGCTracking/interface/HGCTracker.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "FWCore/Utilities/interface/Exception.h"


HGCTracker::HGCTracker(const CaloGeometry *geom) :
    geom_(geom)
{
    //if (hgctracking::g_debuglevel > 0) std::cout << "Making the HGCTracker" << std::endl;
    makeDisks(8, 28);
    makeDisks(9, 24);
    //makeDisks(10, 12);

    auto ptrSort = [](const HGCDiskGeomDet *a, const HGCDiskGeomDet *b) -> bool { return (*a) < (*b); };
    std::sort(disksPos_.begin(), disksPos_.end(), ptrSort);
    std::sort(disksNeg_.begin(), disksNeg_.end(), ptrSort);
}

void HGCTracker::makeDisks(int subdet, int disks)
{
    // const CaloSubdetectorGeometry *subGeom = subdet < 5 ? geom_->getSubdetectorGeometry(DetId::Forward, subdet) :
    //                                                       geom_->getSubdetectorGeometry(DetId::Hcal, 2);

    const CaloSubdetectorGeometry *subGeom = geom_->getSubdetectorGeometry(DetId::Detector(subdet), ForwardSubdetector::ForwardEmpty);

    std::vector<float>  rmax(disks, 0), rmin(disks, 9e9);
    std::vector<double> zsumPos(disks), zsumNeg(disks);
    std::vector<int>    countPos(disks), countNeg(disks);
    const std::vector<DetId> & ids = subGeom->getValidDetIds();
    //if (hgctracking::g_debuglevel > 0) std::cout << "on subdet " << subdet << " I got a total of " << ids.size() << " det ids " << std::endl;
    for (auto & i : ids) {
        const GlobalPoint & pos = geom_->getPosition(i); 
        float z = pos.z();
        float rho = pos.perp();
        int side = z > 0 ? +1 : -1;
        // int layer = (subdet < 5 ? HGCalDetId(i).layer()-1 : HcalDetId(i).depth()-1 ); 
        // if (subdet == 5) {
        //     if (std::abs(z) < 400 || std::abs(z) > 600 || (layer == 4 && rho > 240)) continue; // bad, not BH
        // }

	int layer = std::numeric_limits<unsigned int>::max();
	if (i.det() == DetId::HGCalEE)    layer = HGCSiliconDetId(i).layer() - 1;
	else if(i.det() == DetId::HGCalHSi) layer = HGCSiliconDetId(i).layer() - 1;
	else if (i.det() == DetId::HGCalHSc)  layer = HGCScintillatorDetId(i).layer() - 1;

        (side > 0 ? zsumPos : zsumNeg)[layer] += z;
        (side > 0 ? countPos : countNeg)[layer]++;
        if (rho > rmax[layer]) rmax[layer] = rho;
        if (rho < rmin[layer]) rmin[layer] = rho;
    }
    for (int i = 0; i < disks; ++i) {
        float radlen=-1, xi=-1; // see DataFormats/GeometrySurface/interface/MediumProperties.h
        switch(subdet) {
           case 8: 
                // Phase II TDR page 98 for radlength in units of X0 (ignoring the silicon); 
                //  xi = radlen * X0 * 0.307075 * Z/A * 1/2 = radlen * 6.76 * 0.307075 * 0.402 * 0.5 = radlen * 0.42 // for W
                if (i < 10) { radlen = 0.65; }
                else if (i < 20) { radlen = 0.88; }
                else { radlen = 1.26; }
                xi = radlen * 0.42e-3; // e-3 because the 0.307075 was in MeV
                break;
            case 9:
                // radlen = 3.5 * lambda / X0 / 12 = 3.5 * 137 / 13 / 12 = 3.1 // for Cu
                // xi = radlen * 13 * 0.307075 * 0.456 * 0.5 = radlen * 0.9 // for Cu
                radlen = 3.1; xi = radlen * 0.9e-3; // e-3 because the 0.307075 was in MeV
                break;
            case 10: 
                radlen = 5; xi = radlen * 0.9e-3; // just scale FH with the number of lambdas 3.5 -> 5
                // who knows?
                break;
        }
        if (countPos[i]) {
            //printf("Positive disk %2d at z = %+7.2f   %6.1f <= rho <= %6.1f\n", i+1, zsumPos[i]/countPos[i], rmin[i], rmax[i]);
            addDisk(new HGCDiskGeomDet(subdet, +1, i+1, zsumPos[i]/countPos[i], rmin[i], rmax[i], radlen, xi));
        }
        if (countNeg[i]) {
            //printf("Negative disk %2d at z = %+7.2f   %6.1f <= rho <= %6.1f\n", i+1, zsumNeg[i]/countPos[i], rmin[i], rmax[i]);
            addDisk(new HGCDiskGeomDet(subdet, -1, i+1, zsumNeg[i]/countNeg[i], rmin[i], rmax[i], radlen, xi));
        }
    }
}


const HGCDiskGeomDet * HGCTracker::nextDisk(const HGCDiskGeomDet *from, PropagationDirection direction) const 
{
    const std::vector<HGCDiskGeomDet *> & vec = (from->zside() > 0 ? disksPos_ : disksNeg_);
    auto it = std::find(vec.begin(), vec.end(), from);
    if (it == vec.end()) throw cms::Exception("LogicError", "nextDisk called with invalid starting disk");
    if (direction == alongMomentum) {
        if (*it == vec.back()) return nullptr;
        return *(++it);
    } else {
        if (it == vec.begin()) return nullptr;
        return *(--it);
    }
}

const HGCDiskGeomDet * HGCTracker::idToDet(DetId id) const 
{
    // FIXME: can be made less stupid
    for (const HGCDiskGeomDet * disk : disksPos_) {
        if (disk->geographicalId() == id) return disk;
    }
    for (const HGCDiskGeomDet * disk : disksNeg_) {
        if (disk->geographicalId() == id) return disk;
    }
    return nullptr;
}



