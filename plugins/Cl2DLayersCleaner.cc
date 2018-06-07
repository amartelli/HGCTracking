/** \class Cl2DLayersCleaner
 *   produce Candidates for each RecHit **/

#include <cassert>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <cmath>

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingRecHit.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingClusteringRecHit.h"
#include "RecoParticleFlow/HGCTracking/interface/TrajectorySeedFromTrack.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCTkTrajectoryBuilder.h"
#include "RecoParticleFlow/HGCTracking/interface/hgcdebug.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"


class Cl2DLayersCleaner : public edm::stream::EDProducer<> {
    public:
        explicit Cl2DLayersCleaner(const edm::ParameterSet& ps);
        ~Cl2DLayersCleaner() {}

        virtual void produce(edm::Event& evt, const edm::EventSetup& es) override;

    private:
        const edm::EDGetTokenT<reco::CaloClusterCollection> srcClusters_;
  
        hgcal::RecHitTools rhtools_;

};


Cl2DLayersCleaner::Cl2DLayersCleaner(const edm::ParameterSet& ps) :
    srcClusters_(consumes<reco::CaloClusterCollection>(ps.getParameter<edm::InputTag>("srcClusters")))
{
  produces<std::vector<reco::CaloCluster> > ();
}


void
Cl2DLayersCleaner::produce(edm::Event& evt, const edm::EventSetup& es)
{
  rhtools_.getEventSetup(es);

    // Get 2d clusters
    edm::Handle<reco::CaloClusterCollection> sClusters;
    evt.getByToken(srcClusters_, sClusters);

    std::unique_ptr<std::vector<reco::CaloCluster> > outCl(new std::vector<reco::CaloCluster>());

    for (const reco::CaloCluster &sCl : *sClusters) {
      reco::CaloCluster new2D;
      //      float energyNew2D = 0.;
      /*
      float x = 0;
      float y = 0;
      float z = 0;
      float totW = 0;
      */
      for (auto & p : sCl.hitsAndFractions()) {
	if (p.second != 0) {
	  new2D.addHitAndFraction(p.first, p.second);
	  //  energyNew2D += p.second;
	  /*
	  const GlobalPoint position( std::move( rhtools_.getPosition( p.first ) ) );
	  float localE = p.second;
	  x += position.x() * localE;
	  y += position.y() * localE;
	  z += position.z() * localE;
	  totW += localE;
	  
	  float dist = sqrt( pow(position.x() - sCl.position().x(), 2) + pow(position.y() - sCl.position().y(), 2) );
	  unsigned int layer = rhtools_.getLayerWithOffset(p.first);
	  if(dist > 2.) std::cout << " xH = " << position.x() << " yH = " << position.y() 
				  << " scL x = " << sCl.position().x() << " scL y = " << sCl.position().y() << " scL E = " << sCl.energy() 
				  << " hit E = " << p.second 
				  << " scL size = " << new2D.size() << " layer = " << layer << std::endl;
	  */
	}
      }
      if (new2D.size() == 0) continue; 

      new2D.setEnergy(sCl.energy());                                                                                                   
      new2D.setCorrectedEnergy(sCl.energy());                                                                                                   

      new2D.setPosition(sCl.position());
      outCl->push_back(new2D);                                                                                                   
    }//loop over 2Dcl

    evt.put(std::move(outCl));
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( Cl2DLayersCleaner );
