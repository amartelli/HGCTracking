import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('TRY',eras.Phase2C4)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D28Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('/store/cmst3/user/gpetrucc/l1phase2/081116/DR_PU0/2023D3_guns/ChargedPionGun_PU0_job1.FEVT.root'),
    #fileNames = cms.untracked.vstring('file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/RECO_211Pt10VtxNoSmear.root'),
    fileNames = cms.untracked.vstring('file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/TICLtests/origin/RECO_211Pt10VtxNoSmear.root'),
    #fileNames = cms.untracked.vstring('file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/amartell/HGCAL_Timing_vtxStudies/RECO_211Pt10VtxHLLHC_LowEta/RECO_211Pt10VtxHLLHC_LowEta_10.root'),
    secondaryFileNames = cms.untracked.vstring()
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
# customisation of the process.

process.load('RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff')
process.load('RecoHGCal.HGCTracking.hgcTracks_cff')

process.hgcTracks.debugLevel = 3
process.hgcTracks.patternRecoAlgo = "hitsAndClusters"

process.run = cms.Path(process.hgcalLayerClusters + process.hgcTracks)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("hgctracks_clusters.root"),
    outputCommands = cms.untracked.vstring("drop *", 
        "keep recoTracks_generalTracks_*_*",
        "keep *_HGCalRecHit_*_*", 
        "keep *_particleFlowClusterHGCal_*_*", 
        "keep *_hgcalLayerClusters_*_*",
        "keep *_hgcTracks_*_*",
    ),
)

#add TICL
from RecoHGCal.TICL.ticl_iterations import *
TICL_iterations(process,outputCommands=process.out.outputCommands)
process.myTICL=cms.Path(process.hgcalLayerClusters
                        + process.FilteredLayerClustersMIP
                        + process.TrackstersMIP
                        + process.TrackstersToMultiClusterMIP
                        + process.FilteredLayerClusters
                        + process.Tracksters
                        + process.TrackstersToMultiCluster
                        + process.hgcalMultiClusters)

process.e = cms.EndPath(process.out)
process.schedule=cms.Schedule(process.run,process.myTICL,process.e)


import sys
args = sys.argv[2:] if sys.argv[0] == "cmsRun" else sys.argv[1:]
if args:
    process.source.fileNames = [ args[0] ]
