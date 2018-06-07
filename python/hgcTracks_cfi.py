import FWCore.ParameterSet.Config as cms

from RecoParticleFlow.HGCTracking.hgcTrajectoryBuilder_cfi import *

hgcTracks = cms.EDProducer("HGCTracking",
    hgcTrajectoryBuilderPSet,

    ### Track collection to use for seeding
    srcTk = cms.InputTag("generalTracks"),
    cutTk = cms.string("1.48 < abs(eta) < 3.0 && pt > 0.5 && p > 1 && quality('highPurity')"),
    #srcClusters = cms.InputTag("hgcalLayerClusters"),
    cutSc = cms.string("1.48 < abs(eta) < 3.0"),

    ### Attempt a backwards refit of the HGC hits to get an unbiased track 
    doBackwardsRefit = cms.bool(False),

    ### MC Truth, for debugging
    srcTruth = cms.InputTag("mix","MergedCaloTruth"),

    # verbosity of printout (0 = none)
    debugLevel = cms.untracked.uint32(0)
)

