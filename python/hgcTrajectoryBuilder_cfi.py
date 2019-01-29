import FWCore.ParameterSet.Config as cms

hgcTrajectoryBuilderPSet = cms.PSet(
    ### RecHit collections
    srcEE = cms.InputTag("HGCalRecHit:HGCEERecHits"),
    srcFH = cms.InputTag("HGCalRecHit:HGCHEFRecHits"),
    #srcBH = cms.InputTag("HGCalRecHit:HGCHEBRecHits"),
    ### Clusters
    srcClusters = cms.InputTag("hgcalLayerClusters"),

    ### Seeding
    # initial uncertainties are rescaled by ( 1 + x * N(lost outer hits) )
    lostHitRescaleErrorFactor = cms.double(0.3), 

    ###### EventSetup services:
    # propagators
    propagator = cms.string("PropagatorWithMaterial"), 
    propagatorOpposite = cms.string("PropagatorWithMaterialOpposite"), 
    #
    # estimator (to select compatible hits)
    estimator = cms.string("hitCollectorHGC"),
    #
    # updator (KF algebra; not really much to configure here)
    updator = cms.string("KFUpdator"),
    #
    # trajectory cleaner for final setup
    trajectoryCleaner = cms.string("TrajectoryCleanerBySharedHits"),
    
    ###### Pattern recognition configuration
    #
    # what to use: 
    #     "singleHit"       = individual calo rechits; 
    #     "hitsAndClusters" = clusters if they are available, hits otherwise
    #     "clusterizing"    = start from calo rechits, but do an on-the-fly clustering the specified radius 
    patternRecoAlgo = cms.string("singleHit"),
    clusterRadius = cms.double(1.2), ## only used for patternRecoAlgo = "clusterizing"; to be tuned
    #
    # number of candidates to propagate each step
    maxCand = cms.uint32(1),
    #
    # during bulding, consider as mutually exclusive two candidates that end on the same hit
    # (meaningful only if maxCand > 1, of course)
    endpointCleaner = cms.bool(False),
    #
    # at the end, pick only the single best trajectory from this track
    # (note: if not, will run the trajectoryCleaner on them)
    fastCleaner = cms.bool(False),
    #
    # parameters to choose the best track: 
    # the figure of merit is  foundHitBouns * hits - lostHitPenalty * (lost hits) - chi2
    foundHitBonus = cms.double(20.0),
    lostHitPenalty = cms.double(5.0),
    #
    # choice of whether to penalize for lost hits in the BH subdetector (True) or not (False)
    lostHitsOnBH = cms.bool(False),
    #
    # allow the track to have at most N layers without hits
    maxStartingEmptyLayers = cms.uint32(2), 
    #
    # start enforcing the maxCand after N layers
    layersBeforeCleaning = cms.uint32(1), 
    #
    # don't consider hits which are this much worse in chi2 with respect to the best hit found
    deltaChiSquareForHits = cms.double(10.),
    #
    # don't consider a lost hit if the chi2 of the best hit is below this threshold
    minChi2ForInvalidHit = cms.double(10.),
    #
    # just use the best hit (whenever the maxCand is enforced; note: use only if maxCand == 1!)
    bestHitOnly = cms.bool(True),

    ##### Trajectory filter configuration
    #
    # determines when to give up on a hit, and whether in the end a candidate
    # should be retained or discarded
    # the most important parameters are marked as # <<== important
    trajectoryFilter = cms.PSet(
        ComponentType = cms.string('CkfBaseTrajectoryFilter'),
        maxConsecLostHits = cms.int32(2), # <<== important, decides when to give up on the track
        maxLostHitsFraction = cms.double(0.3),
        constantValueForLostHitsFractionFilter = cms.double(2.0),
        minimumNumberOfHits = cms.int32(7), # <<== important, decides whether to keep the track or not
        chargeSignificance = cms.double(-1.0),
        maxCCCLostHits = cms.int32(9999),
        maxLostHits = cms.int32(999),
        maxNumberOfHits = cms.int32(100),
        minGoodStripCharge = cms.PSet( refToPSet_ = cms.string('SiStripClusterChargeCutNone') ),
        minHitsMinPt = cms.int32(3),
        minNumberOfHitsForLoopers = cms.int32(13),
        minNumberOfHitsPerLoop = cms.int32(4),
        extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
        minPt = cms.double(0.05),
        nSigmaMinPt = cms.double(5.0),
        pixelSeedExtension = cms.bool(False),
        seedExtension = cms.int32(0),
        seedPairPenalty = cms.int32(0),
        strictSeedExtension = cms.bool(False)
    ),
)

