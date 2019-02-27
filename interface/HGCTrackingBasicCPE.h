#ifndef RecoHGCal_HGCTracking_HGCTrackingBasicCPE_h
#define RecoHGCal_HGCTracking_HGCTrackingBasicCPE_h

/// Class that estimates positions and uncertainties for a hit or cluster
//  Similar interface to the ClusterParameterEstimators of the tracking
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "RecoHGCal/HGCTracking/interface/HGCTrackingRecHit.h"
#include "RecoHGCal/HGCTracking/interface/HGCTrackingClusteringRecHit.h"

class HGCTrackingBasicCPE {
    public:
        struct PositionHint {
            PositionHint(float x_, float y_, float size_) : x(x_), y(y_), size(size_) {}
            const float x, y, size;
        };

        // For a circle, E<x^2> = 1/4 r^2 = A / (4 pi); HGcal sensors have 1cm area --> sqrt(E<x^2>) = 0.28
        HGCTrackingBasicCPE(const CaloGeometry* geom, float fixedError=0.28*1.2, float fixedErrorBH=3) ;

       std::pair<LocalPoint,LocalError> localParameters(const HGCRecHit &obj, const Surface &surf) const {
	 float loce2 = ( (obj.id().det() == 8 || getLayer(obj) < 13) ? fixedError2_ : fixedErrorBH2_);
            return std::make_pair(
                        surf.toLocal(getPosition(obj)),
                        LocalError(loce2,0,loce2)
                   );
        }
        std::pair<LocalPoint,LocalError> localParameters(const reco::CaloCluster & obj, const Surface &surf) const {
	  float loce2 = ((obj.hitsAndFractions().front().first.det() == 8 ||
			  getLayer(obj.hitsAndFractions().front().first) < 13) ? fixedError2_ : fixedErrorBH2_);
            return std::make_pair(
                        surf.toLocal(GlobalPoint(obj.position().X(), obj.position().Y(), obj.position().Z())),
                        LocalError(loce2 * obj.size(),0,loce2 * obj.size())
                   );
        }
        PositionHint hint(const HGCTrackingRecHitFromHit & hit) const { return hint(*hit.objRef()); }
        PositionHint hint(const HGCTrackingRecHitFromCluster & hit) const { return hint(*hit.objRef()); }
        PositionHint hint(const HGCRecHit & obj) const {
	  float loce = ((obj.id().det() == 8 || getLayer(obj) < 13) ? fixedError_ : fixedErrorBH_);
            const GlobalPoint & gp = getPosition(obj);
            return PositionHint(gp.x(), gp.y(), loce);
        }
        PositionHint hint(const reco::CaloCluster & obj) const {
	  float loce = ((obj.hitsAndFractions().front().first.det() == 8 ||
			 getLayer(obj.hitsAndFractions().front().first) < 13) ? fixedError_ : fixedErrorBH_);
            const math::XYZPoint &gp = obj.position();
            return PositionHint(gp.X(), gp.Y(), loce * sqrt(obj.size()));
        }
        GlobalPoint getPosition(const HGCRecHit &obj) const { return getPosition(obj.id()); }
        GlobalPoint getPosition(const DetId id) const {
	  if(id.det() == 8) return geomEE_->getPosition(id);
	  else if (id.det() == 9) return geomFH_->getPosition(id);
	  else if(id.det() == 10 ) return geomBH_->getPosition(id);
	  else throw cms::Exception("BadDetId") << "ERROR, unknown detid " << int(id.det()) << ":" << id.subdetId() << "\n" ;	  
        }

	int getLayer(const HGCRecHit &obj) const {return getLayer(obj.id()); }
	int getLayer(const DetId id) const {
	  if(id.det() == 8) return HGCSiliconDetId(id).layer();
	  else if (id.det() == 9) return HGCSiliconDetId(id).layer();
	  else if(id.det() == 10 ) return HGCScintillatorDetId(id).layer();
	  else throw cms::Exception("BadDetId") << "ERROR, unknown detid " << int(id.det()) << ":" << id.subdetId() << "\n" ;
	}

    private:
        const CaloGeometry *geom_;
        const HGCalGeometry *geomEE_, *geomFH_, *geomBH_;

        const float fixedError_, fixedError2_;
        const float fixedErrorBH_, fixedErrorBH2_;
};

#endif
