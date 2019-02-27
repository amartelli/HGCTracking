#include "RecoHGCal/HGCTracking/interface/HGCTrackingBasicCPE.h"

HGCTrackingBasicCPE::HGCTrackingBasicCPE(const CaloGeometry* geom, float fixedError, float fixedErrorBH)  : 
    geom_(geom), 
    geomEE_(static_cast<const HGCalGeometry*>(geom_->getSubdetectorGeometry(DetId::HGCalEE, ForwardSubdetector::ForwardEmpty))),
    geomFH_(static_cast<const HGCalGeometry*>(geom_->getSubdetectorGeometry(DetId::HGCalHSi, ForwardSubdetector::ForwardEmpty))),
    geomBH_(static_cast<const HGCalGeometry*>(geom_->getSubdetectorGeometry(DetId::HGCalHSc, ForwardSubdetector::ForwardEmpty))),
    fixedError_(fixedError), fixedError2_(fixedError*fixedError),
    fixedErrorBH_(fixedErrorBH), fixedErrorBH2_(fixedErrorBH*fixedErrorBH) 
{
}

