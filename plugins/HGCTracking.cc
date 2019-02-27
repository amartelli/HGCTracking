/** \class HGCTracking
 *   produce Candidates for each RecHit **/

#include <cassert>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <cmath>

#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "RecoHGCal/HGCTracking/interface/HGCTrackingRecHit.h"
#include "RecoHGCal/HGCTracking/interface/HGCTrackingClusteringRecHit.h"
#include "RecoHGCal/HGCTracking/interface/TrajectorySeedFromTrack.h"
#include "RecoHGCal/HGCTracking/interface/HGCTkTrajectoryBuilder.h"
#include "RecoHGCal/HGCTracking/interface/hgcdebug.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"


class HGCTracking : public edm::stream::EDProducer<> {
    public:
        explicit HGCTracking(const edm::ParameterSet& ps);
        ~HGCTracking() {}
        virtual void produce(edm::Event& evt, const edm::EventSetup& es) override;

    private:
        /// workhorse 
        HGCTkTrajectoryBuilder builder_;

        /// for seeding
        const edm::EDGetTokenT<reco::TrackCollection> srcTk_;
        const StringCutObjectSelector<reco::Track> cutTk_;

        /// do the backwards fit
        const bool doBackwardsRefit_;

        /// for the debug output
        const edm::EDGetTokenT<HGCRecHitCollection> srcEE_, srcFH_;  //, srcBH_;

        /// for debugging
        const bool hasMCTruth_;
        const edm::EDGetTokenT<std::vector<CaloParticle>> srcTruth_;
        CaloTruthRevMap revTruthMap_;

        void declareOutput(const std::string &label); 
        void makeOutput(const std::vector<Trajectory> &trajs, edm::Event &out, const std::string &label); 
        void writeAllHitsByLayer(const HGCRecHitCollection &hits, edm::Event &out, const std::string &label, int layers);
};


HGCTracking::HGCTracking(const edm::ParameterSet& ps) :
    builder_(ps, consumesCollector()),
    srcTk_(consumes<reco::TrackCollection>(ps.getParameter<edm::InputTag>("srcTk"))),
    cutTk_(ps.getParameter<std::string>("cutTk")),
    doBackwardsRefit_(ps.getParameter<bool>("doBackwardsRefit")),
    srcEE_(consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("srcEE"))),
    srcFH_(consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("srcFH"))),
    //srcBH_(consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("srcBH"))),
    hasMCTruth_(ps.existsAs<edm::InputTag>("srcTruth")),
    srcTruth_(hasMCTruth_ ? consumes<std::vector<CaloParticle>>(ps.getParameter<edm::InputTag>("srcTruth")) : edm::EDGetTokenT<std::vector<CaloParticle>>())
{
    hgctracking::g_debuglevel = ps.getUntrackedParameter<uint32_t>("debugLevel",0);

    declareOutput("");
    //if (doBackwardsRefit_) declareOutput("bw");
    if (hgctracking::g_debuglevel > 0) {
        produces<std::vector<reco::CaloCluster>>("EE");
        produces<std::vector<reco::CaloCluster>>("FH");
        //produces<std::vector<reco::CaloCluster>>("BH");
    }
}

void
HGCTracking::produce(edm::Event& evt, const edm::EventSetup& es)
{
    builder_.init(evt, es);

    edm::Handle<std::vector<CaloParticle>> truth;
    if (hasMCTruth_) {
        evt.getByToken(srcTruth_, truth);
        for (const CaloParticle &p : *truth) {
            if (hgctracking::g_debuglevel > 0) {
                if (p.eventId().event() == 0 && p.eventId().bunchCrossing() == 0) {
                    printf("Considering truth particle of pdgId %+6d eid %d/%d pt %7.1f eta %+5.2f phi %+5.2f simclusters %4d \n", 
                            p.pdgId(), p.eventId().event(), p.eventId().bunchCrossing(), p.pt(), p.eta(), p.phi(), int(p.simClusters().size()));
                }
            }
            for (const auto & scr : p.simClusters()) {
                //printf("    simcluster pt %7.1f eta %+5.2f phi %+5.2f rechits %4d simhits %d \n",
                //    scr->pt(), scr->eta(), scr->phi(), scr->numberOfRecHits(), scr->numberOfSimHits());
                for (const auto & pair : scr->hits_and_fractions()) {
		    revTruthMap_.emplace(pair.first, std::make_pair(&p, pair.second));
                }
            }
        }
        builder_.setTruth(&revTruthMap_);
    }
    

    if (hgctracking::g_debuglevel > 0) {
      edm::Handle<HGCRecHitCollection> srcEE, srcFH; //, srcBH;
        evt.getByToken(srcEE_, srcEE); 
        evt.getByToken(srcFH_, srcFH); 
        //evt.getByToken(srcBH_, srcBH); 
        writeAllHitsByLayer(*srcEE, evt, "EE", 28);
        writeAllHitsByLayer(*srcFH, evt, "FH", 24);
        //writeAllHitsByLayer(*srcBH, evt, "BH", 24);
    }

    // Get tracks
    edm::Handle<reco::TrackCollection> tracks;
    evt.getByToken(srcTk_, tracks);

    // Loop on tracks, make one or more trajectories per track
    std::vector<Trajectory> finaltrajectories;
    unsigned int itrack = 0;
    for (const reco::Track &tk : *tracks) { ++itrack;
        if (!cutTk_(tk)) continue;
        if (hgctracking::g_debuglevel > 1) {
            printf("\n\nConsidering track pt %7.1f eta %+5.2f phi %+5.2f valid hits %d outer lost hits %d highPurity %1d\n", tk.pt(), tk.eta(), tk.phi(), tk.hitPattern().numberOfValidHits(), tk.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS), tk.quality(reco::Track::highPurity));
        }
        builder_.trajectories( reco::TrackRef(tracks,itrack-1), finaltrajectories, alongMomentum );
    }
    if (hgctracking::g_debuglevel > 0) {
        printf("\n\nA total of %lu trajectories found in the event\n",finaltrajectories.size());
    }

    // Global cleaning (should not do much except in case of collimated tracks)
    unsigned int size_pre = finaltrajectories.size();
    builder_.cleanTrajectories(finaltrajectories);
    if (hgctracking::g_debuglevel > 0 && size_pre != finaltrajectories.size()) {
        printf("A total of %lu trajectories after multiple cleaner\n", finaltrajectories.size());
    }

    // Debug Printout 
    if (hgctracking::g_debuglevel > 0)  {
        std::set<const reco::Track *> done;
        for (const Trajectory &t : finaltrajectories) {
            printf("- Trajectory     "); builder_.printTraj(t);
            const reco::Track & tk  = * (dynamic_cast<const TrajectorySeedFromTrack &>(t.seed())).track();
            printf("\t from track pt %7.1f eta %+5.2f phi %+5.2f valid hits %2d outer lost hits %d highPurity %1d\n", tk.pt(), tk.eta(), tk.phi(), tk.hitPattern().numberOfValidHits(), tk.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS), tk.quality(reco::Track::highPurity));
            done.insert(&tk);
        }
        printf("Seed tracks not reconstructed:\n");
        for (const reco::Track &tk : *tracks) {
            if (!cutTk_(tk)) continue;
            if (done.count(&tk)) continue; // reconstructed
            printf("- Track pt %7.1f eta %+5.2f phi %+5.2f valid hits %2d outer lost hits %d highPurity %1d\n", tk.pt(), tk.eta(), tk.phi(), tk.hitPattern().numberOfValidHits(), tk.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS), tk.quality(reco::Track::highPurity));
        }
        printf("done.\n");
    }

    // Output (direct fit)
    makeOutput(finaltrajectories, evt, "");

    // Make backwards fit
    if (doBackwardsRefit_) {
        // Run
        if (hgctracking::g_debuglevel > 0) printf("Now going backwards\n");
        std::vector<Trajectory> bwfits;
        for (const Trajectory &t : finaltrajectories) {
            Trajectory && tbw = builder_.bwrefit(t);
            if (tbw.isValid()) bwfits.push_back(tbw);
        }
        finaltrajectories.swap(bwfits);

        if (hgctracking::g_debuglevel > 0)  {
            for (const Trajectory &t : finaltrajectories) {
                printf("- TrajectoryBW "); builder_.printTraj(t);
            }
        }

        // Output (direct fit)
        //makeOutput(finaltrajectories, evt, "");
    }

    builder_.done();
    revTruthMap_.clear();
}

void HGCTracking::declareOutput(const std::string &label) 
{
    produces<TrackingRecHitCollection>(label);
    produces<std::vector<reco::TrackExtra>>(label);
    produces<std::vector<reco::Track>>(label);
    produces<std::vector<reco::Track>>(label+"Seed");
    produces<std::vector<reco::CaloCluster>>(label);
}

void HGCTracking::makeOutput(const std::vector<Trajectory> &trajs, edm::Event &out, const std::string &label) 
{
    std::unique_ptr<TrackingRecHitCollection> outHits(new TrackingRecHitCollection());
    std::unique_ptr<std::vector<reco::TrackExtra>> outTkEx(new std::vector<reco::TrackExtra>());
    std::unique_ptr<std::vector<reco::Track>> outTk(new std::vector<reco::Track>());
    std::unique_ptr<std::vector<reco::CaloCluster>> outCl(new std::vector<reco::CaloCluster>());
    std::unique_ptr<std::vector<reco::Track>> outTkSeed(new std::vector<reco::Track>());
    auto rTrackExtras = out.getRefBeforePut<reco::TrackExtraCollection>(label);
    auto rHits = out.getRefBeforePut<TrackingRecHitCollection>(label);

    for (const Trajectory &t : trajs) {
        reco::CaloCluster hits;
        bool first = true;
        float energy = 0;
        GlobalPoint pos; GlobalVector mom; int q = 1;
        GlobalPoint outpos; GlobalVector outmom; 
        unsigned int ifirst = outHits->size(), nhits = 0;
        for (const auto & tm : t.measurements()) {
            outHits->push_back(tm.recHit()->clone()); nhits++;
            if (!tm.recHit()->isValid()) continue;
            if (first) {
                pos = tm.updatedState().globalPosition();
                mom = tm.updatedState().globalMomentum();
                q = tm.updatedState().charge();
                hits.setPosition(math::XYZPoint(pos.x(), pos.y(), pos.z()));
                first = false;
            } else if (tm.updatedState().isValid()) {
                outpos = tm.updatedState().globalPosition();
                outmom = tm.updatedState().globalMomentum();
            }
            if (typeid(*tm.recHit()) == typeid(HGCTrackingRecHitFromCluster)) {
                const reco::CaloCluster &cl = *(dynamic_cast<const HGCTrackingRecHitFromCluster&>(*tm.recHit())).objRef();
                for (auto & p : cl.hitsAndFractions()) {
                    if (p.second > 0) hits.addHitAndFraction(p.first, p.second);
                }
                energy += cl.energy();
            } else if (typeid(*tm.recHit()) == typeid(HGCTrackingClusteringRecHit)) {
                for (auto & hr : (dynamic_cast<const HGCTrackingClusteringRecHit&>(*tm.recHit())).objRefs()) {
                    hits.addHitAndFraction(hr->id(), 1.0);
                }
                energy += ((dynamic_cast<const HGCTrackingClusteringRecHit&>(*tm.recHit())).energy());
            } else {
                energy += ((dynamic_cast<const HGCTrackingRecHitFromHit&>(*tm.recHit())).energy());
                hits.addHitAndFraction(tm.recHit()->geographicalId(), 1.0);
            }
        }
        hits.setEnergy(energy);
        hits.setCorrectedEnergy(energy);
        hits.setPosition(math::XYZPoint(0.5*(pos.x()+outpos.x()), 0.5*(pos.y()+outpos.y()), 0.5*(pos.z()+outpos.z())));
        outCl->push_back(hits);

        reco::TrackExtra extra(reco::Track::Point(outpos.x(), outpos.y(), outpos.z()), reco::Track::Vector(outmom.x(), outmom.y(), outmom.z()), true,  
                               reco::Track::Point(pos.x(), pos.y(), pos.z()), reco::Track::Vector(mom.x(), mom.y(), mom.z()), true,
                               reco::Track::CovarianceMatrix(), 0,  reco::Track::CovarianceMatrix(), 0, t.direction());
        extra.setHits(rHits, ifirst, nhits);
        outTkEx->push_back(extra);

        reco::Track trk(t.chiSquared(), t.foundHits()*2, hits.position(),
                        reco::Track::Vector(mom.x(), mom.y(), mom.z()), q, reco::Track::CovarianceMatrix());
        for (const auto & tm : t.measurements()) {
	    // int subdet = 5, layer = 0;
            // if (tm.recHit()->geographicalId().det() == DetId::Forward) {
            //     subdet = tm.recHit()->geographicalId().subdetId();
            //     layer = HGCalDetId(tm.recHit()->geographicalId()).layer();
            // } else {
            //     layer = HcalDetId(tm.recHit()->geographicalId()).depth();
            // }
	    int subdet = tm.recHit()->geographicalId().det();
	    int layer = 0;

            // convert subdet to something representable in a hitpattern
            switch (subdet) {
	    case 8: subdet = PixelSubdetector::PixelEndcap; layer = HGCSiliconDetId(tm.recHit()->geographicalId()).layer(); break;
	    case 9: subdet = StripSubdetector::TID; layer = HGCSiliconDetId(tm.recHit()->geographicalId()).layer(); break;
	    case 10: subdet = StripSubdetector::TEC; layer = HGCScintillatorDetId(tm.recHit()->geographicalId()).layer(); break;
            }
            trk.appendTrackerHitPattern(subdet, layer/2, layer%2, tm.recHit()->type());
        }
        if (hasMCTruth_) {
            auto mresult = builder_.truthMatch(t);
            trk.setNLoops(mresult.empty() ? 0 : ceil(mresult.front().second));
            if (!mresult.empty() && mresult.front().first->eventId().event()==0 && mresult.front().first->eventId().bunchCrossing()==0) {
                trk.setQuality(reco::Track::confirmed);
            }
        }
        outTk->push_back(trk);

        const reco::Track & seedtk  = * (dynamic_cast<const TrajectorySeedFromTrack &>(t.seed())).track();
        outTkSeed->push_back(seedtk);
    }
    for (unsigned int i = 0, n = outTk->size(); i < n; ++i) {
        (*outTk)[i].setExtra(reco::TrackExtraRef(rTrackExtras, i));
    }
    out.put(std::move(outHits), label);
    out.put(std::move(outTkEx), label);
    out.put(std::move(outTk), label);
    out.put(std::move(outCl), label);
    out.put(std::move(outTkSeed), label+"Seed");
}

void HGCTracking::writeAllHitsByLayer(const HGCRecHitCollection &hits, edm::Event &out, const std::string &label, int layers) 
{
    std::unique_ptr<std::vector<reco::CaloCluster>> outCl(new std::vector<reco::CaloCluster>(2*layers));

    for (const HGCRecHit &hit : hits) {
        int zside, layer;

	if (hit.id().det() == DetId::HGCalEE){
	  layer = HGCSiliconDetId(hit.id()).layer();
	  zside = HGCSiliconDetId(hit.id()).zside();
	} else if(hit.id().det() == DetId::HGCalHSi) {
	  layer = HGCSiliconDetId(hit.id()).layer();
	  zside = HGCSiliconDetId(hit.id()).zside();
	} else if (hit.id().det() == DetId::HGCalHSc) {
	  layer = HGCScintillatorDetId(hit.id()).layer();
	  zside = HGCScintillatorDetId(hit.id()).zside();
	}
        (*outCl)[(zside>0)+2*(layer-1)].addHitAndFraction(hit.id(), 1.0);
    }
    for (int ilayer = 0; ilayer < layers; ++ilayer) {
       (*outCl)[2*ilayer  ].setEnergy(ilayer+1);
       (*outCl)[2*ilayer+1].setEnergy(ilayer+1.5);
       math::RhoEtaPhiVectorF pointP(1, +0.1*(ilayer+1), 0);
       math::RhoEtaPhiVectorF pointN(1, -0.1*(ilayer+1), 0);
       (*outCl)[2*ilayer  ].setPosition(math::XYZPoint(pointN.x(),pointN.y(),pointN.z()));
       (*outCl)[2*ilayer+1].setPosition(math::XYZPoint(pointP.x(),pointP.y(),pointP.z()));
    }

    out.put(std::move(outCl), label);
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( HGCTracking );
