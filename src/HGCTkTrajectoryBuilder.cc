#include "RecoHGCal/HGCTracking/interface/HGCTkTrajectoryBuilder.h"
#include "RecoHGCal/HGCTracking/interface/TrajectorySeedFromTrack.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/PatternTools/interface/TempTrajectory.h"
#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"
#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryCleaning/interface/FastTrajectoryCleaner.h"
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilterFactory.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"


HGCTkTrajectoryBuilder::HGCTkTrajectoryBuilder(const edm::ParameterSet& ps, edm::ConsumesCollector && c ) :
    srcEE_(c.consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("srcEE"))),
    srcFH_(c.consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("srcFH"))),
    //srcBH_(c.consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("srcBH"))),
    srcClusters_(c.consumes<reco::CaloClusterCollection>(ps.getParameter<edm::InputTag>("srcClusters"))),
    propName_(ps.getParameter<std::string>("propagator")),
    propNameOppo_(ps.getParameter<std::string>("propagatorOpposite")),
    estimatorName_(ps.getParameter<std::string>("estimator")),
    updatorName_(ps.getParameter<std::string>("updator")),
    trajectoryCleanerName_(ps.getParameter<std::string>("trajectoryCleaner")),
    thrSoN_(ps.getParameter<double>("singleHitSoNthreshold")),
    lostHitRescaleErrorFactor_(ps.getParameter<double>("lostHitRescaleErrorFactor")),
    fastCleaner_(ps.getParameter<bool>("fastCleaner")), 
    endpointCleaner_(ps.getParameter<bool>("endpointCleaner")), 
    theMaxCand(ps.getParameter<uint32_t>("maxCand")),
    theFoundHitBonus(ps.getParameter<double>("foundHitBonus")),
    theLostHitPenalty(ps.getParameter<double>("lostHitPenalty")),
    theMaxStartingEmptyLayers(ps.getParameter<uint32_t>("maxStartingEmptyLayers")),
    theLayersBeforeCleaning(ps.getParameter<uint32_t>("layersBeforeCleaning")),
    bestHitOnly_(ps.getParameter<bool>("bestHitOnly")), 
    deltaChiSquareForHits_(ps.getParameter<double>("deltaChiSquareForHits")),
    minChi2ForInvalidHit_(ps.getParameter<double>("minChi2ForInvalidHit")),
    clusterRadius_(ps.getParameter<double>("clusterRadius")),
    lostHitsOnBH_(ps.getParameter<bool>("lostHitsOnBH")),
    geomCacheId_(0),
    truthMap_(nullptr)
{
    std::string palgo = ps.getParameter<std::string>("patternRecoAlgo");
    if (palgo == "singleHit") {
        algo_ = SingleHitAlgo;
    } else if (palgo == "hitsAndClusters") {
        algo_ = MixedAlgo;
    } else if (palgo == "clusterizing") {
        algo_ = ClusterizingAlgo;
    } else throw cms::Exception("Configuration") << "Unknown algo: knowns are 'singleHit', 'clusterizing'\n";

    auto pset = ps.getParameter<edm::ParameterSet>("trajectoryFilter");
    trajFilter_.reset( TrajectoryFilterFactory::get()->create(pset.getParameter<std::string>("ComponentType"), pset, c) );
}


void
HGCTkTrajectoryBuilder::init(const edm::Event& evt, const edm::EventSetup& es)
{
    //Get Calo Geometry
    if (es.get<CaloGeometryRecord>().cacheIdentifier() != geomCacheId_) {
        es.get<CaloGeometryRecord>().get(caloGeom_);
        geomCacheId_ = es.get<CaloGeometryRecord>().cacheIdentifier();
        hgcTracker_.reset(new HGCTracker(caloGeom_.product()));
        cpe_.reset(new HGCTrackingBasicCPE(&*caloGeom_)); // FIXME better
    } 

    es.get<GlobalTrackingGeometryRecord>().get(trkGeom_);
    es.get<IdealMagneticFieldRecord>().get(bfield_);
    es.get<TrackingComponentsRecord>().get(propName_, prop_);
    if (!propNameOppo_.empty()) es.get<TrackingComponentsRecord>().get(propNameOppo_, propOppo_);
    es.get<TrackingComponentsRecord>().get(estimatorName_,estimator_);
    es.get<TrackingComponentsRecord>().get(updatorName_,updator_);

    if (!trajectoryCleanerName_.empty()) es.get<TrajectoryCleaner::Record>().get(trajectoryCleanerName_, trajectoryCleaner_);
   
    trajFilter_->setEvent(evt, es); 

    data_.reset(new HGCTrackingData(*hgcTracker_, &*cpe_));

    evt.getByToken(srcEE_, srcEE); data_->addData(srcEE, 8, thrSoN_);
    evt.getByToken(srcFH_, srcFH); data_->addData(srcFH, 9, thrSoN_);
    //evt.getByToken(srcBH_, srcBH); data_->addData(srcBH, 10);

    evt.getByToken(srcClusters_, srcClusters); 
    data_->addClusters(srcClusters);
}


void 
HGCTkTrajectoryBuilder::done()
{
    data_.reset();
}


unsigned int 
HGCTkTrajectoryBuilder::trajectories(const reco::TrackRef &tk, std::vector<Trajectory> &out, PropagationDirection direction) const 
{
    assert(direction == alongMomentum); // the rest is not implemented yet

    FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(*tk, bfield_.product());
    fts.rescaleError(1.0 + lostHitRescaleErrorFactor_ * tk->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS));

    std::vector<Trajectory> myOut;
    trajectories(fts, myOut, direction);

    if (!myOut.empty()) {
        auto ostate = trajectoryStateTransform::outerStateOnSurface(*tk, *trkGeom_, &*bfield_);
        auto pstate = trajectoryStateTransform::persistentState(ostate, DetId(tk->outerDetId()));
        boost::shared_ptr<TrajectorySeed> seed(new TrajectorySeedFromTrack(tk, pstate, alongMomentum));
        for (Trajectory &t : myOut) { 
            t.setSharedSeed(seed);
            out.push_back(std::move(t));
        }
    }

    return myOut.size();
}


unsigned int 
HGCTkTrajectoryBuilder::trajectories(const FreeTrajectoryState &fts, std::vector<Trajectory> &out, PropagationDirection direction) const 
{
    auto trajCandLess = [&](TempTrajectory const & a, TempTrajectory const & b) {
        return  (a.chiSquared() + a.lostHits()*theLostHitPenalty - a.foundHits()*theFoundHitBonus)  <
            (b.chiSquared() + b.lostHits()*theLostHitPenalty - b.foundHits()*theFoundHitBonus);
    };

    int zside = fts.momentum().eta() > 0 ? +1 : -1;
    const HGCDiskGeomDet *disk = hgcTracker_->firstDisk(zside,direction);
    std::vector<TempTrajectory> myfinaltrajectories;
    std::vector<TempTrajectory> trajectories = advanceOneLayer(fts, TempTrajectory(direction, 0), disk, false);
    unsigned int depth = 2;
    for (disk = hgcTracker_->nextDisk(disk,direction); disk != nullptr; disk = hgcTracker_->nextDisk(disk,direction), ++depth) {
        if (trajectories.empty()) continue;
        if (hgctracking::g_debuglevel > 1) {
            printf("   New destination: disk subdet %d, zside %+1d, layer %2d, z = %+8.2f\n", disk->subdet(), disk->zside(), disk->layer(), disk->toGlobal(LocalPoint(0,0,0)).z());
        }
        //printf("   Starting candidates: %lu\n", trajectories.size());
        std::vector<TempTrajectory> newCands;
        int icand = 0;
        for (TempTrajectory & cand : trajectories) {
            if (hgctracking::g_debuglevel > 1) {
                printf("    Processing candidate %2d/%lu with ", ++icand, trajectories.size()); printTraj(cand); //found hits %3d, lost hits %3d, chi2 %8.1f\n", cand.foundHits(), cand.lostHits(), cand.chiSquared());
            }
            TrajectoryStateOnSurface start = cand.lastMeasurement().updatedState();
            if (!start.isValid()) start = cand.lastMeasurement().predictedState();
            bool bestHitOnly = bestHitOnly_ && (depth > theLayersBeforeCleaning);
            std::vector<TempTrajectory> hisTrajs = advanceOneLayer(start, cand, disk, bestHitOnly);
            if (hisTrajs.empty()) {
                if (hgctracking::g_debuglevel > 1) printf("     --> stops here\n");
                if (trajFilter_->qualityFilter(cand)) {
                    if (hgctracking::g_debuglevel > 1) printf("          --> passes filter, being retained for the moment \n");
                    myfinaltrajectories.push_back(cand);
                }
            } else {
                //printf("---- produced %lu trajectories \n", hisTrajs.size());
                for ( TempTrajectory & t : hisTrajs ) { 
                    if (t.foundHits() == 0 && depth > theMaxStartingEmptyLayers) continue;
                    if (trajFilter_->toBeContinued(t)) {
                        newCands.push_back(t);
                    } else {
                        //printf("------    --> one is not to be continued.\n");
                        if (trajFilter_->qualityFilter(t)) {
                            //printf("------         --> but it is to be saved.\n");
                            myfinaltrajectories.push_back(std::move(t));
                        }
                    }
                }
            }
        }
        //printf("   A total of %lu trajectories after this step\n", newCands.size());
        if (endpointCleaner_ && depth > theLayersBeforeCleaning) {
            unsigned int oldsize = newCands.size();
            TrajectoryCleanerBySharedEndpoints endpointCleaner(theFoundHitBonus, theLostHitPenalty);
            endpointCleaner.clean(newCands); trim(newCands);
            if (hgctracking::g_debuglevel > 1) if (oldsize != newCands.size()) printf("    Reduced from %u to %lu trajectories after TrajectoryCleanerBySharedEndpoints\n", oldsize, newCands.size());
        }
        if (newCands.size() > theMaxCand && depth > theLayersBeforeCleaning) {
            std::sort(newCands.begin(), newCands.end(), trajCandLess);
            newCands.resize(theMaxCand);
            //printf("    Reduced to %lu trajectories after sorting and trimming\n", newCands.size());
        }
        trajectories.swap(newCands);
    }
    if (!trajectories.empty()) {
        if (hgctracking::g_debuglevel > 1) printf("A total of %lu trajectories reached the end of the tracker from this track\n", trajectories.size());
        for (TempTrajectory & t : trajectories) {
            if (t.foundHits() > 0 && trajFilter_->qualityFilter(t)) {
                myfinaltrajectories.push_back(std::move(t));
            }
        }
    }
    if (hgctracking::g_debuglevel > 1) printf("A total of %lu trajectories found from this track\n", myfinaltrajectories.size());
    if (fastCleaner_) {
        FastTrajectoryCleaner cleaner(theFoundHitBonus/2,theLostHitPenalty); // the factor 1/2 is because the FastTrajectoryCleaner uses the ndof instead of the number of found hits
        cleaner.clean(myfinaltrajectories); trim(myfinaltrajectories);
        if (hgctracking::g_debuglevel > 1) printf("A total of %lu trajectories after FastTrajectoryCleaner\n", myfinaltrajectories.size());
    }

    std::vector<Trajectory> toPromote;
    for (TempTrajectory &t : myfinaltrajectories) {
        while (!t.lastMeasurement().recHit()->isValid()) t.pop();
        if (hgctracking::g_debuglevel > 1) { printf("- TempTrajectory "); printTraj(t); }
        toPromote.push_back(t.toTrajectory());
    }
    if (toPromote.size() > 1 && !trajectoryCleanerName_.empty()) {
        trajectoryCleaner_->clean(toPromote); trim(toPromote);
        if (hgctracking::g_debuglevel > 1) printf("A total of %lu trajectories after multiple cleaner\n", toPromote.size());
    }
    for (Trajectory &t : toPromote) { 
        if (hgctracking::g_debuglevel > 1) { printf("- Trajectory "); printTraj(t); }
        out.push_back(std::move(t));
    }

    return toPromote.size();
}

template<class Start>
std::vector<TempTrajectory> 
HGCTkTrajectoryBuilder::advanceOneLayer(const Start &start, const TempTrajectory &traj, const HGCDiskGeomDet *disk, bool bestHitOnly) const  
{
    std::vector<TempTrajectory> ret;

    const Propagator &prop = (traj.direction() == alongMomentum ? *prop_ : *propOppo_);
    // propagate to the plane of the layer
    TrajectoryStateOnSurface tsos = prop.propagate(start, disk->surface());
    if (!tsos.isValid()) { 
        if (hgctracking::g_debuglevel > 0)  {
            printf("        Destination disk subdet %d, zside %+1d, layer %2d, z = %+8.2f\n", disk->subdet(), disk->zside(), disk->layer(), disk->toGlobal(LocalPoint(0,0,0)).z());
            printf("         --> propagation failed.\n"); 
        }
        return ret; 
    }

    // check if inside the bounds of the layer
    if (!disk->surface().bounds().inside(tsos.localPosition())) {
        if (hgctracking::g_debuglevel > 0)  {
            GlobalPoint gp = tsos.globalPosition();
            printf("        Prop point: global eta %+5.2f phi %+5.2f  x = %+8.2f y = %+8.2f z = %+8.2f rho = %8.2f\n", gp.eta(), float(gp.phi()), gp.x(), gp.y(), gp.z(), gp.perp());
            //LocalPoint lp = tsos.localPosition();
            //printf("            local                       x = %+8.2f y = %+8.2f z = %+8.2f rho = %8.2f\n", lp.x(), lp.y(), lp.z(), lp.perp());
            printf("         --> outside the bounds.\n"); 
        }
        return ret; 
    }

    // get the data on the layer
    const auto & diskData = data_->diskData(disk);

    // for debug
    if (truthMap_) diskData.setTruth(truthMap_); 

    //printf("        Looking for hits on disk subdet %d, zside %+1d, layer %2d: %u total hits\n", disk->subdet(), disk->zside(), disk->layer(), diskData.size());
    std::vector<TrajectoryMeasurement> meas;
    switch(algo_) {
        case SingleHitAlgo:    meas = diskData.measurements(tsos, *estimator_); break;
        case ClusterizingAlgo: meas = diskData.clusterizedMeasurements(tsos, *estimator_, clusterRadius_); break;
        case MixedAlgo:    
                meas = diskData.clusterMeasurements(tsos, *estimator_); 
                std::vector<TrajectoryMeasurement> meas2 = diskData.measurements(tsos, *estimator_); 
                if (meas.empty()) meas2.swap(meas);
            break;
    };

    // sort hits from better to worse
    std::sort(meas.begin(), meas.end(), TrajMeasLessEstim());
    if (hgctracking::g_debuglevel > 1)  printf("        Compatible hits: %lu\n", meas.size());

    // for each, make a new trajectory candidate
    for (const TrajectoryMeasurement &tm : meas) {
        if (deltaChiSquareForHits_ > 0) {
            if (meas.size() > 1  && !ret.empty() && tm.estimate() > meas.front().estimate() + deltaChiSquareForHits_) {
                if (hgctracking::g_debuglevel > 3) printf("        stop after the first %lu hits, since this chi2 of %.1f is too bad wrt the best one of %.1f\n", ret.size(), tm.estimate(), meas.front().estimate());
                break;
            }
        }
        TrajectoryStateOnSurface updated = updator_->update(tm.forwardPredictedState(), *tm.recHit());
        if (!updated.isValid()) { 
            if (hgctracking::g_debuglevel > 0)  {
                std::cout << "          Hit with chi2 = " << tm.estimate() << std::endl;
                std::cout << "              track state     " << tm.forwardPredictedState().localPosition() << std::endl;
                std::cout << "              rechit position " << tm.recHit()->localPosition() << std::endl;
                std::cout << "               --> failed update state" << std::endl;
            }
            continue;
        }

        // Add a valid hit
        ret.push_back(traj.foundHits() ? traj : TempTrajectory(traj.direction(),0)); // don't start with a lost hit
        ret.back().push(TrajectoryMeasurement(tm.forwardPredictedState(),
                                              updated,
                                              tm.recHit(),
                                              tm.estimate()), 
                        tm.estimate());
    
        // fast return
        if (bestHitOnly) return ret;
    }

    // Possibly add an invalid hit, for the hypothesis that the track didn't leave a valid hit
    if (minChi2ForInvalidHit_ > 0) {
        if (meas.size() > 0  && !ret.empty() && meas.front().estimate() < minChi2ForInvalidHit_) {
            if (hgctracking::g_debuglevel > 3) printf("        will not add the invalid hit after %lu valid hits, as the best valid hit has chi2 of %.1f\n", ret.size(), meas.front().estimate());
            return ret;
        }
    }

    //if subdet == 9 and layer > 12 is in last 12 layers = ex BH
    auto missing = (disk->subdet() == 8 || disk->layer() < 13 || lostHitsOnBH_) ? TrackingRecHit::missing : TrackingRecHit::inactive;
    ret.push_back(traj.foundHits() ? traj : TempTrajectory(traj.direction(),0)); // either just one lost hit, or a trajectory not starting on a lost hit
    ret.back().push(TrajectoryMeasurement(tsos, std::make_shared<InvalidTrackingRecHit>(*disk, missing)));
    return ret;
}

void
HGCTkTrajectoryBuilder::cleanTrajectories(std::vector<Trajectory> &trajectories) const {
    if (!trajectoryCleanerName_.empty()) {
        trajectoryCleaner_->clean(trajectories);
        trim(trajectories);
    }
}

Trajectory
HGCTkTrajectoryBuilder::bwrefit(const Trajectory &traj, float scaleErrors) const {
    Trajectory ret(oppositeToMomentum);

    const Trajectory::DataContainer & tms = traj.measurements();

    ret.reserve(tms.size());

    TrajectoryStateOnSurface tsos = tms.back().updatedState();
    tsos.rescaleError(scaleErrors);

    const Propagator &propOppo = (traj.direction() == alongMomentum ? *propOppo_ : *prop_);

    for (int i = tms.size()-1; i >= 0; --i) {
        const HGCDiskGeomDet * det = hgcTracker_->idToDet(tms[i].recHit()->geographicalId());
        if (det == 0) {
            if (hgctracking::g_debuglevel > 0)  {
                printf(" ---> failure in finding det for step %d on det %d, subdet %d\n",i,tms[i].recHit()->geographicalId().det(),tms[i].recHit()->geographicalId().subdetId());
            }
            ret.invalidate();
            return ret;
        }
        TrajectoryStateOnSurface prop = propOppo.propagate(tsos, det->surface());
        if (!prop.isValid()) {
            if (hgctracking::g_debuglevel > 0)  {
                printf(" ---> failure in propagation for step %d\n",i);
            }
            ret.invalidate();
            return ret;
        }
        if (tms[i].recHit()->isValid()) {
           auto pair = estimator_->estimate( prop, *tms[i].recHit() );
           if (hgctracking::g_debuglevel > 1)  {
               if (i % 7 == 0) {
               printf("   for step %2d pt %6.1f +- %5.1f   q/p %+7.4f +- %6.4f   x = %+7.2f +- %4.2f  y = %+7.2f +- %4.2f   dxdz = %+5.3f +- %4.3f dydz = %+5.3f +- %4.3f    adding hit %+7.2f %+7.2f  results in chi2 = %6.1f (old: %6.1f)\n", i, 
                    prop.globalMomentum().perp(), prop.globalMomentum().perp() * prop.globalMomentum().mag() * sqrt(prop.localError().matrix()(0,0)),
                    prop.globalParameters().signedInverseMomentum(), sqrt(prop.localError().matrix()(0,0)),
                    prop.localPosition().x(), sqrt(prop.localError().matrix()(3,3)), prop.localPosition().y(),  sqrt(prop.localError().matrix()(4,4)),
                    prop.localParameters().dxdz(), sqrt(prop.localError().matrix()(1,1)), prop.localParameters().dydz(), sqrt(prop.localError().matrix()(2,2)),
                    tms[i].recHit()->localPosition().x(), tms[i].recHit()->localPosition().y(),
                    pair.second, tms[i].estimate());
               }
           }
           TrajectoryStateOnSurface updated = updator_->update( prop, *tms[i].recHit() );
           if (!updated.isValid()) {
                if (hgctracking::g_debuglevel > 0) printf(" ---> fail in backwards update for step %d\n",i);
                ret.invalidate(); 
                return ret;
           } 
           //printf("       updated pt %6.1f +- %5.1f   q/p %+7.4f +- %6.4f   x = %+7.2f +- %4.2f  y = %+7.2f +- %4.2f   dxdz = %+5.3f +- %4.3f dydz = %+5.3f +- %4.3f\n",  
           //     updated.globalMomentum().perp(), updated.globalMomentum().perp() * updated.globalMomentum().mag() * sqrt(updated.localError().matrix()(0,0)),
           //     updated.globalParameters().signedInverseMomentum(), sqrt(updated.localError().matrix()(0,0)),
           //     updated.localPosition().x(), sqrt(updated.localError().matrix()(3,3)), updated.localPosition().y(),  sqrt(updated.localError().matrix()(4,4)),
           //     updated.localParameters().dxdz(), sqrt(updated.localError().matrix()(1,1)), updated.localParameters().dydz(), sqrt(updated.localError().matrix()(2,2)));
           ret.push( TrajectoryMeasurement(prop, updated, tms[i].recHit(), pair.second) ); 
           tsos = updated;
        } else {
           ret.push( TrajectoryMeasurement(prop, tms[i].recHit()) ); 
           tsos = prop;
        }
    }
    //printf("Reversed trajectory: "); printTraj(ret);
    return ret;
}


template<typename Traj> 
void 
HGCTkTrajectoryBuilder::printTraj_(const Traj &t) const 
{
    TrajectoryStateOnSurface finalTSOS = t.lastMeasurement().updatedState();
    if (!finalTSOS.isValid()) finalTSOS = t.lastMeasurement().forwardPredictedState();
    GlobalPoint pos = finalTSOS.globalPosition(); GlobalVector mom = finalTSOS.globalMomentum();
    LocalTrajectoryError lte = finalTSOS.localError();
    //printf("found hits %3d, lost hits %3d, chi2 %8.1f ndf %3d outermost state pt = %6.1f +- %5.1f   eta %+5.2f phi %+5.2f x = %+7.2f +- %4.2f  y = %+7.2f +- %4.2f   z = %+7.2f",
    printf("found hits %3d, lost hits %3d, chi2 %8.1f ndf %3d outermost pt = %6.1f x = %+7.2f y = %+7.2f z = %+7.2f last-id %12d first-id %12d ",
            t.foundHits(), t.lostHits(), t.chiSquared(), t.foundHits()*2-5, 
            mom.perp(), //mom.perp() * mom.mag() * sqrt(lte.matrix()(0,0)),
            //pos.eta(), float(pos.phi()), 
            pos.x(), /*sqrt(lte.matrix()(3,3)),*/ pos.y(), /*sqrt(lte.matrix()(4,4)),*/
            pos.z(),
            t.lastMeasurement().recHit()->isValid() ? t.lastMeasurement().recHit()->geographicalId()() : 0,
            t.firstMeasurement().recHit()->isValid() ? t.firstMeasurement().recHit()->geographicalId()() : 0
          );
    // median energy over the last 4 hits
    if (t.foundHits() >= 3) {
        float etotal = 0;
        std::vector<float> energies;
        const auto &meas = t.measurements();
        for (int i = meas.size()-1; i >= 0; --i) {
            if (meas[i].recHit()->isValid()) {
                float energy = -1;
                if (typeid(*meas[i].recHit()) == typeid(HGCTrackingRecHitFromHit)) {
                    energy = ((dynamic_cast<const HGCTrackingRecHitFromHit&>(*meas[i].recHit())).energy());
                } else if (typeid(*meas[i].recHit()) == typeid(HGCTrackingRecHitFromCluster)) {
                    energy = ((dynamic_cast<const HGCTrackingRecHitFromCluster&>(*meas[i].recHit())).energy());
                } else if (typeid(*meas[i].recHit()) == typeid(HGCTrackingClusteringRecHit)) {
                    energy = ((dynamic_cast<const HGCTrackingClusteringRecHit&>(*meas[i].recHit())).energy());
                } else {
                    throw cms::Exception("Invalid valid hit", typeid(*meas[i].recHit()).name());
                }
                if (energies.size() < 4) energies.push_back(energy);
                etotal += energy;
            }
        }
        if (energies.size()) {
            std::sort(energies.begin(),energies.end());
            float median = (energies.size() % 2 == 1) ? energies[energies.size()/2] : 0.5*(energies[energies.size()/2]+energies[energies.size()/2-1]);
            //printf("   energy median %7.3f (%lu)", median, energies.size());
            printf("   energy %6.3f et_tot %6.2f  ", median, etotal * pos.perp() / pos.mag());
        }
    }
    if (truthMap_) {
        for (const auto & pair : truthMatch(t)) {
            printf("  %.1f hits from %s pdgId %+d eid %d/%+d pt %.1f eta %+5.2f phi %+5.2f    ", pair.second, 
                    (pair.first->eventId().event()==0&&pair.first->eventId().bunchCrossing()==0 ? "SIGNAL" : "pileup"),
                    pair.first->pdgId(), pair.first->eventId().event(), pair.first->eventId().bunchCrossing(), pair.first->pt(), pair.first->eta(), pair.first->phi());
        }
    }
    printf("\n");
}

template<typename Traj>
std::vector<std::pair<const CaloParticle *, float>> 
HGCTkTrajectoryBuilder::truthMatch_(const Traj &t) const 
{
    std::vector<std::pair<const CaloParticle *, float>> ret;
    if (truthMap_) { 
        std::map<const CaloParticle *,float> scores;       
        std::vector<const CaloParticle *> keys;
        std::vector<DetId> ids;
        for (const auto & tm : t.measurements()) {
            if (!tm.recHit()->isValid()) continue;
            ids.clear();
            if (typeid(*tm.recHit()) == typeid(HGCTrackingRecHitFromCluster)) {
                for (auto & p : (dynamic_cast<const HGCTrackingRecHitFromCluster&>(*tm.recHit())).objRef()->hitsAndFractions()) {
                    if (p.second > 0) ids.push_back(p.first);
                }
            } else {
                ids.push_back(tm.recHit()->geographicalId());
            }
            for (DetId detid: ids) {
                auto range = truthMap_->equal_range(detid.rawId());
                for (; range.first != range.second; ++range.first) {
                    const auto &pair = *range.first;
                    if (std::find(keys.begin(), keys.end(), pair.second.first) == keys.end()) {
                        keys.push_back(pair.second.first);
                        scores[pair.second.first] += 1;
                    }
                    //scores[pair.second.first] += pair.second.second > 0.0 ? 1 : pair.second.second;
                }   
            }
            keys.clear();
        }
        keys.clear();
        for (auto & pair : scores) keys.push_back(pair.first);
        std::sort(keys.begin(), keys.end(), [&scores](const CaloParticle *a, const CaloParticle *b) -> bool { return scores[b] < scores[a]; });
        for (auto & key : keys) {
            ret.emplace_back(key, scores[key]);
        }
    }
    return ret;
}


