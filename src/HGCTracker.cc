#include "RecoHGCal/HGCTracking/interface/HGCTracker.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "FWCore/Utilities/interface/Exception.h"


HGCTracker::HGCTracker(const CaloGeometry *geom) :
    geom_(geom)
{
  computeAbsorbers();

    //if (hgctracking::g_debuglevel > 0) std::cout << "Making the HGCTracker" << std::endl;
    makeDisks(8, 28);
    makeDisks(9, 24);
    //makeDisks(10, 24);

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
	  if (i%2 == 0) {
	    radlen = 0.748 * xi_["Pb"] + 0.068 * xi_["Fe"] + 0.014 * xi_["Cu"];
	    xi = radlen / (0.748 + 0.068 + 0.014) * 1.e-3;
	  }
	  else{
	    radlen = 0.648 * xi_["WCu"] + 0.417 * xi_["Cu"];
	    xi = radlen / (0.648 + 0.417) * 1.e-3;
	  }
	  break;
	case 9:
	  if (i == 0){
	    radlen = 0.374 * xi_["Pb"] + (0.007+0.07) * xi_["Cu"] + (0.034+2.277) * xi_["Fe"];
	    xi = radlen / (0.374 + 0.007+0.07 + 0.034+2.277) * 1.e-3;
	  }
	  else if(i < 12){
	    radlen = 0.2315 * xi_["WCu"] + 0.487 * xi_["Cu"] + 1.992 * xi_["Fe"];
	    xi = radlen / (0.2315 + 0.487 + 1.992) * 1.e-3;
	  }
	  else{
	    radlen = 0.2315 * xi_["WCu"] + 0.487 * xi_["Cu"] + 3.870 * xi_["Fe"];
	    xi = radlen / (0.2315 + 0.487 + 3.870) * 1.e-3;
	  }
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


void HGCTracker::computeAbsorbers()
{
  std::map<std::string, float> X0_;
  std::map<std::string, float> lambda_;
  std::map<std::string, float> ZoA_;

  X0_["Fe"] = 13.84;
  X0_["Pb"] = 6.37;
  X0_["Cu"] = 12.86;
  X0_["W"] = 6.76;
  X0_["WCu"] = combineX0(0.75, X0_["W"], 0.25, X0_["Cu"]);
  //X0_["WCu"] = combinedEdX(0.75, X0_["W"], 0.25, X0_["Cu"]);

  lambda_["Fe"] = 132.1;
  lambda_["Pb"] = 199.6;
  lambda_["Cu"] = 137.3;
  lambda_["W"] = 191.9;
  lambda_["WCu"] = combineX0(0.75, lambda_["W"], 0.25, lambda_["Cu"]);
  //lambda_["WCu"] = combinedEdX(0.75, lambda_["W"], 0.25, lambda_["Cu"]);

  ZoA_["Fe"] = 0.466;
  ZoA_["Pb"] = 0.396;
  ZoA_["Cu"] = 0.456;
  ZoA_["W"] = 0.403;
  ZoA_["WCu"] = combinedEdX(0.75, ZoA_["W"], 0.25, ZoA_["Cu"]);

  // see DataFormats/GeometrySurface/interface/MediumProperties.h
  for(auto ij : X0_){
    xi_[ij.first] = X0_[ij.first] * 0.307075 * ZoA_[ij.first] * 0.5;
    std::cout << " " << ij.first << " xi = " << xi_[ij.first] << std::endl;
  }


  //from TDR pag 18
  //radlen = number of X0
  //EE odd layers: 0.748 Pb + 0.068 Fe
  //EE even layers: 0.648 WCu + 0.417 Cu

  //FH (Had first 12) L28: 0.374 Pb + 0.034 Fe + 0.007 Cu + 2.277 Fe + 0.07 Cu
  //FH (Had first 12) L>29: 0.2315 WCu + 1.992 Fe + 0.487 Cu

  //BH (Had last 12): 0.2315 WCu + 3.87 Fe + 0.487 Cu

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



