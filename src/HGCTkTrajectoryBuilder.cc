#include "RecoParticleFlow/HGCTracking/interface/HGCTkTrajectoryBuilder.h"
#include "RecoParticleFlow/HGCTracking/interface/TrajectorySeedFromTrack.h"
#include "RecoParticleFlow/HGCTracking/interface/TrajectorySeedFromCaloHit.h"


#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
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
#include "RecoParticleFlow/HGCTracking/interface/TempCaloTrajectory.h"

#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"
#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryCleaning/interface/FastTrajectoryCleaner.h"
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilter.h"
#include "TrackingTools/TrajectoryFiltering/interface/TrajectoryFilterFactory.h"
#include "TrackingTools/TrajectoryState/interface/BasicTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

#include "DataFormats/TrackReco/interface/TrajectoryStopReasons.h"

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"

HGCTkTrajectoryBuilder::HGCTkTrajectoryBuilder(const edm::ParameterSet& ps, edm::ConsumesCollector && c ) :
    srcEE_(c.consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("srcEE"))),
    srcFH_(c.consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("srcFH"))),
    srcBH_(c.consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("srcBH"))),
    srcClusters_(c.consumes<reco::CaloClusterCollection>(ps.getParameter<edm::InputTag>("srcClusters"))),
    propName_(ps.getParameter<std::string>("propagator")),
    propNameOppo_(ps.getParameter<std::string>("propagatorOpposite")),
    estimatorName_(ps.getParameter<std::string>("estimator")),
    updatorName_(ps.getParameter<std::string>("updator")),
    trajectoryCleanerName_(ps.getParameter<std::string>("trajectoryCleaner")),
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
    } else if (palgo == "singleCluster") {
        algo_ = SingleClusterAlgo;
    } else throw cms::Exception("Configuration") << "Unknown algo: knowns are 'singleHit', 'clusterizing'\n";

    auto pset = ps.getParameter<edm::ParameterSet>("trajectoryFilter");
    trajFilter_.reset( TrajectoryFilterFactory::get()->create(pset.getParameter<std::string>("ComponentType"), pset, c) );

    printDEBUG = true;


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

    evt.getByToken(srcEE_, srcEE); data_->addData(srcEE, 3);
    evt.getByToken(srcFH_, srcFH); data_->addData(srcFH, 4);
    evt.getByToken(srcBH_, srcBH); data_->addData(srcBH, 5);

    evt.getByToken(srcClusters_, srcClusters); 
    data_->addClusters(srcClusters);

    evoState = new HGCClusterBuilder();
}


void 
HGCTkTrajectoryBuilder::done()
{
    data_.reset();
    delete evoState;
}


unsigned int 
HGCTkTrajectoryBuilder::trajectories(const reco::TrackRef &tk, std::vector<Trajectory> &out, PropagationDirection direction) const 
{
    assert(direction == alongMomentum); // the rest is not implemented yet

    FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(*tk, bfield_.product());

    fts.rescaleError(1.0 + lostHitRescaleErrorFactor_ * tk->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS));


    auto ostate = trajectoryStateTransform::outerStateOnSurface(*tk, *trkGeom_, &*bfield_);
    auto pstate = trajectoryStateTransform::persistentState(ostate, DetId(tk->outerDetId()));
    boost::shared_ptr<TrajectorySeed> seed(new TrajectorySeedFromTrack(tk, pstate, alongMomentum));

    std::vector<Trajectory> myOut;
    trajectories(*seed, fts, myOut, direction);

    if (!myOut.empty()) {
        for (Trajectory &t : myOut) { 
            t.setSharedSeed(seed);
            out.push_back(std::move(t));
        }
    }

    return myOut.size();
}


unsigned int 
HGCTkTrajectoryBuilder::trajectories(const reco::CaloClusterPtr& sCl, std::vector<Trajectory> &out, PropagationDirection direction, int seedLayer) const 
{
     assert(direction == alongMomentum); // the rest is not implemented yet

    float withinCellError = 0.28*1.2;

    float Xpos = sCl->position().x();
    float Ypos = sCl->position().y();
    float Zpos = sCl->position().z();
    int zside = (Zpos > 0) ? +1 : -1;

    const LocalPoint posCl(sCl->position().x(), sCl->position().y(), sCl->position().z());
    const LocalError errCl(withinCellError, 0., withinCellError);
    float norm = sqrt(Xpos*Xpos+Ypos*Ypos+Zpos*Zpos);
    const LocalVector pCl(sCl->x()/norm, sCl->y()/norm, sCl->z()/norm);

    if(printDEBUG) std::cout << " injected momentum = " << pCl << std::endl;

    //charge 0 works but the pT in the update of the trajectory explodes => multiple scattering becomes negligible and error too small
    // tried with charge 1 problems of pT estimate is ok 
    // (best would be to force the pT to be high enough...straight track) => to-be-checked

    //LocalTrajectoryParameters sClTraj(posCl, pCl, TrackCharge(0));
    //LocalTrajectoryError sClTrajErr(withinCellError, withinCellError, std::abs(withinCellError/Ypos), std::abs(withinCellError/Xpos), 1./sCl->energy());

    LocalTrajectoryParameters sClTraj(posCl, pCl, TrackCharge(1));
    LocalTrajectoryError sClTrajErr(withinCellError, withinCellError, std::abs(withinCellError/Ypos), std::abs(withinCellError/Xpos), 0.5);
    if(printDEBUG) std::cout << " obtained  momentum = " << sClTraj.momentum() << std::endl;


    TrajectoryStateOnSurface tsosCl(sClTraj, sClTrajErr, (hgcTracker_->disk(zside, seedLayer))->surface(), bfield_.product(), SurfaceSideDefinition::atCenterOfSurface);

    if(printDEBUG)
      std::cout << " >>> start HGCTkTrajectoryBuilder::trajectories tsosCl.globalMomentum() = " << tsosCl.globalMomentum() 
		<< " pt " << tsosCl.globalMomentum().perp()  
		<< " tsosCl.globalPosition() = " << tsosCl.globalPosition() 
		<< " trajectoryerror . localerror = " << tsosCl.localError().positionError() << std::endl;
    
    FreeTrajectoryState fts = *(tsosCl.freeState());
    

    if(printDEBUG) std::cout << " >>> HGCTkTrajectoryBuilder::trajectories (sCl) seedL = " << seedLayer+1 << " energy = " << sCl->energy() << std::endl;

    auto pstate = trajectoryStateTransform::persistentState(tsosCl, DetId(sCl->seed()));
    HGCTrackingRecHitFromCluster seedCl(DetId(sCl->seed()), sCl, posCl, errCl);
    boost::shared_ptr<TrajectorySeed> seed(new TrajectorySeedFromCaloHit(seedCl, seedLayer+1, pstate, alongMomentum));
    
    std::vector<Trajectory> myOut;
    trajectories(*seed, fts, myOut, direction, seedLayer);
   
    if (!myOut.empty()) {
      for (Trajectory &t : myOut) { 
	t.setSharedSeed(seed);
	out.push_back(std::move(t));
      }
    }

    return myOut.size();
}


unsigned int 
HGCTkTrajectoryBuilder::trajectories(const TrajectorySeed& trjSeed, const FreeTrajectoryState &fts, std::vector<Trajectory> &out, PropagationDirection direction, int seedLayer) const 
{
    auto trajCandLess = [&](TempCaloTrajectory const & a, TempCaloTrajectory const & b) {
      return  (a.chiSquared() + a.lostHits()*theLostHitPenalty - a.foundHits()*theFoundHitBonus)  <
      (b.chiSquared() + b.lostHits()*theLostHitPenalty - b.foundHits()*theFoundHitBonus);
    };

    int zside = fts.momentum().eta() > 0 ? +1 : -1;

    float startEnergy = 0;
    int startSize = 0;
    int seedL = 0;
    bool trkSeed = 1;
    if (typeid(trjSeed) == typeid(TrajectorySeedFromCaloHit)){
      const reco::CaloCluster& cl = *dynamic_cast<const TrajectorySeedFromCaloHit &>(trjSeed).clusterRef();
      startEnergy = cl.energy();
      startSize = cl.size();
      seedL = seedLayer+1;
      trkSeed = 0;
    }
    evoState->initState(startEnergy, startSize, seedL, (trkSeed == 0));

    //std::cout << " >>> start from seed startEnergy = " << startEnergy << " startSize = " << startSize << " seedL = " << seedL << std::endl;

    //RA FIX check that first layer is considered and not double counted
    const HGCDiskGeomDet *disk = (seedLayer != -1) ? hgcTracker_->disk(zside, seedLayer) : hgcTracker_->firstDisk(zside,direction);
    if(printDEBUG)   std::cout << " >>> start with disk subD + layer = " << disk->subdet() << " " << disk->layer() << std::endl;
    std::vector<TempCaloTrajectory> myfinaltrajectories;
    std::vector<TempCaloTrajectory> trajectories = advanceOneLayer(startEnergy, startSize, seedL, fts, TempCaloTrajectory(direction, 0), disk, false);


    //depth = 2 from original code => seem reasonable | not major arguments
    if(printDEBUG)    std::cout << " start iterative loops "<< std::endl;
    unsigned int depth = 2;
    for (disk = hgcTracker_->nextDisk(disk,direction); disk != nullptr; disk = hgcTracker_->nextDisk(disk,direction), ++depth) {
      if(printDEBUG)
	{
	  std::cout << " >>> main loop in loop disk subD + layer = " << disk->subdet() << " " << disk->layer() << std::endl;
	  std::cout << " depth = " << depth << std::endl;
	}
        if (trajectories.empty()) continue;
        if (hgctracking::g_debuglevel > 1) {
            printf("   New destination in trajectories: disk subdet %d, zside %+1d, layer %2d, z = %+8.2f\n", disk->subdet(), disk->zside(), disk->layer(), disk->toGlobal(LocalPoint(0,0,0)).z());
        }
	if(printDEBUG)         printf("   Starting candidates: %lu\n", trajectories.size());
        std::vector<TempCaloTrajectory> newCands;
        int icand = 0;
        for (TempCaloTrajectory & cand : trajectories) {
            if (hgctracking::g_debuglevel > 1) {
	      printf("    Processing candidate %2d/%lu with ", ++icand, trajectories.size()); 
            }
            TrajectoryStateOnSurface start = cand.getTempTraj()->lastMeasurement().updatedState();
            if (!start.isValid()) start = cand.getTempTraj()->lastMeasurement().predictedState();
            bool bestHitOnly = bestHitOnly_ && (depth > theLayersBeforeCleaning);
	    if(printDEBUG)   std::cout << " >>> again advanceOnelayer " << std::endl;

	    TrajectoryMeasurement localTM = cand.getTempTraj()->lastMeasurement();	    

	    if (typeid(*localTM.recHit()) == typeid(HGCTrackingRecHitFromCluster)){
	      startEnergy = (dynamic_cast<const HGCTrackingRecHitFromCluster &>(*localTM.recHit())).objRef()->energy();
	      startSize = (dynamic_cast<const HGCTrackingRecHitFromCluster&>(*localTM.recHit())).objRef()->size();
	      seedL = disk->progressiveLayer()-1;
	    }
	    else if (typeid(*localTM.recHit()) == typeid(HGCTrackingClusteringRecHit)) {
	      // do something
	    }
	    else if (typeid(*localTM.recHit()) == typeid(HGCTrackingRecHitFromHit)){
	      startEnergy = (dynamic_cast<const HGCTrackingRecHitFromHit&>(*localTM.recHit())).energy();
	      startSize = 1;
	      seedL = disk->progressiveLayer()-1;
	    }
	    else{
	      //invalid hit
	      startEnergy = -1.;
	      startSize = -1;
	      seedL = disk->progressiveLayer()-1;
	      //std::cout << "HGCTkTrajectoryBuilder.cc Invalid valid hit" << typeid(*localTM.recHit()).name() << std::endl;
	      //throw cms::Exception("Invalid valid hit", typeid(*localTM.recHit()).name());
	    }
	    
            std::vector<TempCaloTrajectory> hisTrajs = advanceOneLayer(startEnergy, startSize, seedL, start, cand, disk, bestHitOnly);
	    if(printDEBUG)	    std::cout << " 2nd loop size trajectories = " << hisTrajs.size() << std::endl;
            if (hisTrajs.empty()) {
                if (hgctracking::g_debuglevel > 1) printf("     --> stops here\n");
                if (trajFilter_->qualityFilter(*(cand.getTempTraj()))) {
                    if (hgctracking::g_debuglevel > 1) printf("          --> passes filter, being retained for the moment \n");
                    myfinaltrajectories.push_back(cand);
                }
            } else {
                //printf("---- produced %lu trajectories \n", hisTrajs.size());
                for ( TempCaloTrajectory & tc : hisTrajs ) { 
		  auto t = *(tc.getTempTraj());
		  //RA FIX: check theMaxStartingEmptyLayers
                    if (t.foundHits() == 0 && depth > theMaxStartingEmptyLayers) continue;
                    if (trajFilter_->toBeContinued(t)) {
		      std::cout << " trajs to be continued " << std::endl;
		      newCands.push_back(tc);
                    } else {
		      double stopstop = static_cast<double>(t.stopReason());
		      std::cout << " stop reason = " << stopstop << std::endl;
		      std::cout << " trajs fails quality filters " << std::endl;
                        //printf("------    --> one is not to be continued.\n");
                        if (trajFilter_->qualityFilter(t)) {
			  std::cout << " trajs to be saved " << std::endl;
			  //printf("------         --> but it is to be saved.\n");
			  myfinaltrajectories.push_back(std::move(tc));
                        }
                    }
                }
            }
        }//advance oneLayer wrt first meas found
        //printf("   A total of %lu trajectories after this step\n", newCands.size());

	if(printDEBUG)	std::cout << " end 2nd loop in main one new candidate size = " << newCands.size() << " trajectories size = " << trajectories.size()<< std::endl;
        if (endpointCleaner_ && depth > theLayersBeforeCleaning) {
	    //RA replaced cleaning with branch duplicator...
	    // works also for TempCaloTrajectory (track-calo state)
	    unsigned int oldsize = newCands.size();
            unsigned int oldFsize = myfinaltrajectories.size();

	    TrajectoryCleanerBySharedEndpoints endpointCleaner(theFoundHitBonus, theLostHitPenalty);
	    endpointCleaner.clean(newCands, myfinaltrajectories); 
	    
	    trim(newCands);
            if (hgctracking::g_debuglevel > 0) if (oldsize != newCands.size()) {
		printf("    Reduced from %u to %lu trajectories after TrajectoryCleanerBySharedEndpoints\n", oldsize, newCands.size());
		printf("    Increased from %u to %lu fianlTrajectories after TrajectoryCleanerBySharedEndpoints\n", oldFsize, myfinaltrajectories.size());
	      }
	}
        if (newCands.size() > theMaxCand && depth > theLayersBeforeCleaning) {
            std::sort(newCands.begin(), newCands.end(), trajCandLess);
            newCands.resize(theMaxCand);
            if(printDEBUG) printf("    Reduced to %lu trajectories after sorting and trimming\n", newCands.size());
        }
	trajectories.swap(newCands);
    }//itarative loop for tracking

    if(printDEBUG)    std::cout << " end main loop new trajectories size = " << trajectories.size() << std::endl;
    if (!trajectories.empty()) {
        if (hgctracking::g_debuglevel > 1) printf("A total of %lu trajectories reached the end of the tracker from this track\n", trajectories.size());
        for (TempCaloTrajectory & tc : trajectories) {
	  auto t = *(tc.getTempTraj());
            if (t.foundHits() > 0 && trajFilter_->qualityFilter(t)) {
                myfinaltrajectories.push_back(std::move(tc));
            }
        }
    }
    if (hgctracking::g_debuglevel > 1) printf("A total of %lu trajectories found from this track\n", myfinaltrajectories.size());

    //currently OFF
    if (fastCleaner_) {
        FastTrajectoryCleaner cleaner(theFoundHitBonus/2,theLostHitPenalty); // the factor 1/2 is because the FastTrajectoryCleaner uses the ndof instead of the number of found hits

	for(auto ij : myfinaltrajectories) std::cout << " myfinaltrajectories ij.isValid = " << ij.isValid() << std::endl;
	std::vector<TempTrajectory> tmpTraj = getTempTrajVec(myfinaltrajectories);
	cleaner.clean(tmpTraj); 
	std::vector<TempCaloTrajectory> tmpCaloTraj = getTempCaloTrajVec(myfinaltrajectories, tmpTraj);
	myfinaltrajectories.swap(tmpCaloTraj);
	for(auto ij : myfinaltrajectories) std::cout << " myfinaltrajectories ij.isValid = " << ij.isValid() << std::endl;
	
	trim(myfinaltrajectories);
        if (hgctracking::g_debuglevel > 1) printf("A total of %lu trajectories after FastTrajectoryCleaner\n", myfinaltrajectories.size());
    }

    std::vector<Trajectory> toPromote;
    if(printDEBUG)    std::cout << " >>> in HGCTkTrajectoryBuilder.cc toPromote " << std::endl;
    for (TempCaloTrajectory &tc : myfinaltrajectories) {
      auto t = *(tc.getTempTraj());
      while (!t.lastMeasurement().recHit()->isValid()) { 
	if(printDEBUG) std::cout << " pop hit " << std::endl; 
	t.pop(); }
      if(printDEBUG) std::cout << " >>> post pop: foundHits =  " << t.foundHits() << " lost H = " << t.lostHits() << std::endl;
        if (hgctracking::g_debuglevel > 1) { printf("- TempTrajectory "); printTraj(t); }
        toPromote.push_back(t.toTrajectory());
    }
    /*
    //RA FIX clenaing by shared hits -- removed
    if (toPromote.size() > 1 && !trajectoryCleanerName_.empty()) {
      //RA removed: for the moment replace with tree branch filter
      //        trajectoryCleaner_->clean(toPromote); trim(toPromote);
        if (hgctracking::g_debuglevel > 1) printf("A total of %lu trajectories after multiple cleaner\n", toPromote.size());
    }
    */
    for (Trajectory &t : toPromote) { 
        if (hgctracking::g_debuglevel > 1) { printf("- Trajectory "); printTraj(t); }
        out.push_back(std::move(t));
    }

    return toPromote.size();
}


template<class Start>
std::vector<TempCaloTrajectory> 
HGCTkTrajectoryBuilder::advanceOneLayer(float startE, int startS, int startSeed, const Start &start, const TempCaloTrajectory &traj, const HGCDiskGeomDet *disk, bool bestHitOnly) const  
{
  if(printDEBUG)  std::cout << " >>> advanceOneLayer " << std::endl;

  auto sortBestChi2 = [&](std::pair<int, float> const & a, std::pair<int, float> const & b) {
    return  (a.second  < b.second);
    };

    std::vector<TempCaloTrajectory> ret;


    const Propagator &prop = (traj.direction() == alongMomentum ? *prop_ : *propOppo_);
    // propagate to the plane of the layer
    TrajectoryStateOnSurface tsos = prop.propagate(start, disk->surface());
    if (!tsos.isValid()) { 
        if (hgctracking::g_debuglevel > 0)  {
            printf("         Destination disk subdet %d, zside %+1d, layer %2d, z = %+8.2f\n", disk->subdet(), disk->zside(), disk->layer(), disk->toGlobal(LocalPoint(0,0,0)).z());
            printf("         --> propagation failed.\n"); 
        }
        return ret; 
    }

    //for low pt tracks can be useful to increase errors at the beginning of the propagation
    //tsos.rescaleError(1.*startS); with some criteria...
    //to-be-checked


    if(printDEBUG)    printf("  in AdvanceOneLayer      Destination disk subdet %d, zside %+1d, layer %2d\n", disk->subdet(), disk->zside(), disk->layer());

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
    if (truthMap_) {
      diskData.setTruth(truthMap_); 
      diskData.setTruthBis(truthMapX_, truthMapY_, truthMapZ_); 
    }


    std::cout << " >>> traj.measurements().size() = " << traj.measurements().size() << " traj.foundHits() = " << traj.foundHits() 
	      << " evoState->is2DClusterSeed() = " << evoState->is2DClusterSeed() << std::endl;


    //printf("        Looking for hits on disk subdet %d, zside %+1d, layer %2d: %u total hits\n", disk->subdet(), disk->zside(), disk->layer(), diskData.size());
    //    bool foundHits = false;
    std::vector<TrajectoryMeasurement> measA;
    switch(algo_) {
        case SingleHitAlgo:       measA = diskData.measurements(tsos, *estimator_); break;
        case ClusterizingAlgo:    measA = diskData.clusterizedMeasurements(tsos, *estimator_, clusterRadius_); break;
        case SingleClusterAlgo:   measA = diskData.clusterMeasurements(tsos, *estimator_); break;
        case MixedAlgo:    
                measA = diskData.clusterMeasurements(tsos, *estimator_); 
		//RA seems better to look for recHits only if with 2DCl it's empty... no?
                // std::vector<TrajectoryMeasurement> meas2 = diskData.measurements(tsos, *estimator_); 
                // if (measA.empty()) { meas2.swap(measA); }

		//criteria to cut on lowSoN: is2DClusterSeed
		//assume MIP will be track seeded... 
		//maybe add criteria on #2DCl already found...other?
		if (measA.empty()) measA = diskData.measurements(tsos, *estimator_, (evoState->is2DClusterSeed())); 
		break;
    };


    //RACHECKCODEFINOAQUI
    float sumELayer = 0.;
    int sumSLayer = 0;
    int currentLayer = disk->progressiveLayer();

    if(printDEBUG) std::cout << " >>>>> start seedL = " << startSeed << " E = " << startE << " S " << startS  
			     << " found hits = " << measA.size() << " currentLayer = " << currentLayer << std::endl;
    if(printDEBUG)  {
      std::cout << " evoState->getTotalStateEnergyL = " << evoState->getTotalStateEnergyL(currentLayer) 
		<< " evoState->getSeedLayer = " << evoState->getSeedLayer() 
		<< " evoState->getSeedEnergy = " << evoState->getSeedEnergy() << std::endl; 
    }


    //RA option for cluster evolution
      
    //expState: expected  <energy size> at that layer
    //no hits => (0,0)
    //<3hits => previous layer
    //>3hits => interpolation among last 3
    //only use evoState to access the seed state
    //this allows to check every trajectory being built wrt the hits found one-by-one
    std::pair<float, float> expState = traj.getPredState(evoState, currentLayer);
    std::cout << " expState traj E = " << expState.first << " S = " << expState.second << std::endl;

    //RA clone measurement vector with corrected Chi2
    std::vector<std::pair<float, float> > measEnergyChi2;
    std::vector<std::pair<int, float> > measSizeChi2;    
    std::vector<std::pair<int, float> > measIndexChi2;
    std::vector<TrajectoryMeasurement> meas;
    int indexMeas = 0;
    for (const TrajectoryMeasurement &tm : measA) {

      //current refers to the current measurement
      if(printDEBUG)      std::cout << "disk progressive layer " << currentLayer  <<  std::endl;
      float currentEnergy = 0.;
      int currentSize = 0;
      if (typeid(*tm.recHit()) == typeid(HGCTrackingRecHitFromCluster)) {
	currentEnergy = (dynamic_cast<const HGCTrackingRecHitFromCluster&>(*tm.recHit())).objRef()->energy();
        currentSize = (dynamic_cast<const HGCTrackingRecHitFromCluster&>(*tm.recHit())).objRef()->size();
      }
      else if (typeid(*tm.recHit()) == typeid(HGCTrackingClusteringRecHit)) {
      }
      else if (typeid(*tm.recHit()) == typeid(HGCTrackingRecHitFromHit)){
	currentEnergy = ((dynamic_cast<const HGCTrackingRecHitFromHit&>(*tm.recHit())).energy());
	currentSize = 1;
      }
      else{
	currentEnergy = -1.;
	currentSize = -1.;
      }
      if(printDEBUG){    
	std::cout << " layer " << currentLayer << " E = " << currentEnergy << " S = " << currentSize <<  " startE = " << startE << std::endl;
      }

      //RA FIX heck how to deal with corrected Chi2
      // nb anche se current == previous

      float energyChi2 = (expState.first <= 0) ? 0 : std::abs(expState.first - currentEnergy)/expState.first;
      float sizeChi2 = (expState.second <= 0) ? 0 : std::abs(expState.second - currentSize)/expState.second;
      float updatedChi2 = 0.;
      sumELayer += currentEnergy;
      sumSLayer += currentSize;

      if(currentEnergy == -1. || startE == -1.){
	updatedChi2 = (tm.estimate());	
      }
      else{
	updatedChi2 = sqrt(energyChi2*energyChi2 + sizeChi2*sizeChi2) + tm.estimate()/2.;
	//also tested => not optimized
	//updatedChi2 = (tm.estimate() + energyChi2*energyChi2 + sizeChi2*sizeChi2)/3.;

	if(printDEBUG){  
	  std::cout << " energyChi2 = " << energyChi2 << " sizeChi2 " << sizeChi2 << " original = " << tm.estimate() << std::endl;
	  std::cout << " sumELayer = " << sumELayer << " sumSLayer = " << sumSLayer << std::endl;
	}
      }
      if(printDEBUG)
	std::cout << " original = " << tm.estimate() << "updateChi2 = " << updatedChi2 << std::endl;


      //RA option for cluster evolution
      //for evoState: gloablChi2 and evolveState to consider all measurements from a seed
      //start from tem traj being built and check compatibility with all hits alltogether - as clustering
      //std::pair<float, float> gloabalChi2 = evoState->evaluateESChi2(sumELayer, sumSLayer, currentLayer, startSeed, startS);
      //possibility to stop including hits when the globalChi2 is minimised?
      //not really used for the moment
      //std::cout << " gloabalChi2 E = " << gloabalChi2.first << "  S = " << gloabalChi2.second << std::endl;
      //if(gloabalChi2 > 2) break;
      
      meas.push_back(TrajectoryMeasurement(tm.forwardPredictedState(), tm.backwardPredictedState(), tm.updatedState(), tm.recHitP(), updatedChi2, tm.layer()));
      measEnergyChi2.push_back(std::pair<float, float>(currentEnergy, updatedChi2));
      measSizeChi2.push_back(std::pair<int, float>(currentSize, updatedChi2));
      measIndexChi2.push_back(std::pair<int, float>(indexMeas, updatedChi2));
      ++indexMeas;
    }
    //see above comment
    evoState->evolveState(sumELayer, sumSLayer, currentLayer); 

    if(printDEBUG)    std::cout << "now sort " << std::endl;

    // sort hits from better to worse
    //std::sort(meas.begin(), meas.end(), TrajMeasLessEstim());
    std::sort(measIndexChi2.begin(), measIndexChi2.end(), sortBestChi2);
    std::vector<std::pair<float, float> > measEnergyChi2_sort;
    std::vector<std::pair<int, float> > measSizeChi2_sort;
    measA.clear();
    for(auto ij:measIndexChi2){
      measEnergyChi2_sort.push_back(measEnergyChi2.at(ij.first));
      measSizeChi2_sort.push_back(measSizeChi2.at(ij.first));
      measA.push_back(meas.at(ij.first));    
    }
    meas.swap(measA);
    if (hgctracking::g_debuglevel > 1)  printf("        Compatible hits: %lu\n", meas.size());

    if(printDEBUG)    std::cout << "now add " << std::endl;
    // for each, make a new trajectory candidate
    int countMeas = 0;
    bool returnBoolTrue = false;
    for (const TrajectoryMeasurement &tm : meas) {

      if (deltaChiSquareForHits_ > 0) {
	if (meas.size() > 1  && !ret.empty() && tm.estimate() > meas.front().estimate() + deltaChiSquareForHits_) {
	  if (hgctracking::g_debuglevel > 3) printf("        stop after the first %lu hits, since this chi2 of %.1f is too bad wrt the best one of %.1f\n", ret.size(), tm.estimate(), meas.front().estimate());
	  break;
	}
      }


      TrajectoryStateOnSurface updated = updator_->update(tm.forwardPredictedState(), *tm.recHit());
      //RA => increase search window on next layer
      //seem solved with charge != 0 => need pT to handle multiple scattering
      //updated.rescaleError((deltaSize >= 0) ? (1+0.1*sqrt(deltaSize)) : 1.);
      //updated.rescaleError((deltaSize >= 0) ? (1+0.1*(deltaSize)/2.) : 1.);

      std::cout << " pt post update " << updated.globalMomentum().perp() << std::endl;

      if (!updated.isValid()) { 
	if (hgctracking::g_debuglevel > 0)  {
	  std::cout << "          Hit with chi2 = " << tm.estimate() << std::endl;
	  std::cout << "              track state     " << tm.forwardPredictedState().localPosition() << std::endl;
	  std::cout << "              rechit position " << tm.recHit()->localPosition() << std::endl;
	  std::cout << "               --> failed update state" << std::endl;
	}
            continue;
      }

      //NB in case of propagation of several candidates traj gets duplicated
      // Add a valid hit

      ret.push_back(traj.foundHits() ? traj : TempCaloTrajectory(traj.direction(),0)); // don't start with a lost hit
      ret.back().addMeasurement(TrajectoryMeasurement(tm.forwardPredictedState(),
						       updated,
						       tm.recHit(),
						      tm.estimate()), 
				tm.estimate());

      ret.back().updateState(measEnergyChi2_sort.at(countMeas).first, measSizeChi2_sort.at(countMeas).first, currentLayer);
      ++countMeas;      

      // fast return
      if (bestHitOnly) return ret;
      else returnBoolTrue = true;
      if(printDEBUG) std::cout << " *** build parallel branches " << std::endl;
    }//loop over meas
    if(returnBoolTrue) return ret;

    //RA missing in original => check 
    //add return if !bestHitOnly


    if(printDEBUG)    std::cout << " valuta se trj ok " << std::endl;
    // Possibly add an invalid hit, for the hypothesis that the track didn't leave a valid hit
    if (minChi2ForInvalidHit_ > 0) {
        if (meas.size() > 0  && !ret.empty() && meas.front().estimate() < minChi2ForInvalidHit_) {
            if (hgctracking::g_debuglevel > 3) printf("        will not add the invalid hit after %lu valid hits, as the best valid hit has chi2 of %.1f\n", ret.size(), meas.front().estimate());
            return ret;
        }
    }

    // now add lost hits
    if(printDEBUG)    std::cout << " adding lost hit " << std::endl;
    if(printDEBUG) std::cout << " >>pre ****  ret.size() = " << ret.size() << " last nFound = " << traj.foundHits() 
			     << " last lost = " <<  traj.lostHits() << std::endl;
    auto missing = (disk->subdet() != 5 || lostHitsOnBH_) ? TrackingRecHit::missing : TrackingRecHit::inactive;
    ret.push_back(traj.foundHits() ? traj : TempCaloTrajectory(traj.direction(),0)); // either just one lost hit, or a trajectory not starting on a lost hit
    ret.back().addMeasurement(TrajectoryMeasurement(tsos, std::make_shared<InvalidTrackingRecHit>(*disk, missing)));
    if(printDEBUG) std::cout << " >>post **** ret.size() = " << ret.size() << " last nFound = " << ret.back().foundHits() 
			     << " last lost = " <<  ret.back().lostHits() << std::endl;

    if(printDEBUG) std::cout << " advanceOneLayer found n trajectories = " << ret.size() << std::endl;
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
  if(printDEBUG)  std::cout << " >>> HGCTkTrajectoryBuilder::printTraj_ " << std::endl;
    TrajectoryStateOnSurface finalTSOS = t.lastMeasurement().updatedState();

    std::cout << ">>> fino qui ok  post t.lastMeasurement().updatedState(); " << std::endl;

    if (!finalTSOS.isValid()) finalTSOS = t.lastMeasurement().forwardPredictedState();
    GlobalPoint pos = finalTSOS.globalPosition(); GlobalVector mom = finalTSOS.globalMomentum();
    LocalTrajectoryError lte = finalTSOS.localError();
    //printf("found hits %3d, lost hits %3d, chi2 %8.1f ndf %3d outermost state pt = %6.1f +- %5.1f   eta %+5.2f phi %+5.2f x = %+7.2f +- %4.2f  y = %+7.2f +- %4.2f   z = %+7.2f",
    printf("found hits %3d, lost hits %3d, chi2 %8.1f ndf %3d outermost pt = %6.1f x = %+7.2f y = %+7.2f z = %+7.2f last-id %12d first-id %12d --- chi2-LxP+FxB = %8.1f",
            t.foundHits(), t.lostHits(), t.chiSquared(), t.foundHits()*2-5, 
            mom.perp(), //mom.perp() * mom.mag() * sqrt(lte.matrix()(0,0)),
            //pos.eta(), float(pos.phi()), 
            pos.x(), /*sqrt(lte.matrix()(3,3)),*/ pos.y(), /*sqrt(lte.matrix()(4,4)),*/
            pos.z(),
            t.lastMeasurement().recHit()->isValid() ? t.lastMeasurement().recHit()->geographicalId()() : 0,
	   t.firstMeasurement().recHit()->isValid() ? t.firstMeasurement().recHit()->geographicalId()() : 0,
	   t.chiSquared() + theLostHitPenalty*t.lostHits() -  theFoundHitBonus*t.foundHits()
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
		    else scores[pair.second.first] += 1; //RA FIXME
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


//to use TempTrajectory methods for the TempCaloTrajectory class
//hugly but seem to work
std::vector<TempTrajectory> HGCTkTrajectoryBuilder::getTempTrajVec(std::vector<TempCaloTrajectory> tmpCT) const{
  std::vector<TempTrajectory> tmpVec;
  for(auto ij : tmpCT) tmpVec.push_back(*(ij.getTempTraj()));
  return tmpVec;
}


std::vector<TempCaloTrajectory> HGCTkTrajectoryBuilder::getTempCaloTrajVec(std::vector<TempCaloTrajectory> tmpCaloTraj, std::vector<TempTrajectory> tmpTraj) const{
  if(tmpCaloTraj.size() != tmpTraj.size()) std::cout << " big size problem!!! HGCTkTrajectoryBuilder::setTempCaloTrajVec " << std::endl;

  std::vector<TempCaloTrajectory> tmpVec;

  for(unsigned int ij=0; ij < tmpTraj.size(); ++ij) 
    tmpVec.push_back(TempCaloTrajectory(std::move(tmpTraj.at(ij)), std::move(tmpCaloTraj.at(ij))) );

  return tmpVec;
}

