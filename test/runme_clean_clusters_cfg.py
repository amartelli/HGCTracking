import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('TRY',eras.Phase2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


import FWCore.Utilities.FileUtils as FileUtils
readFiles = cms.untracked.vstring()
#readFiles.extend(FileUtils.loadListFromFile ('/afs/cern.ch/work/a/amartell/HGCAL/clustering/timingStudies/CMSSW_9_4_0_pre1/src/HGCTimingAnalysis/HGCTiming/test/testPU/fileList/PDG_130_pt5_140PU_lowEta.txt') )
readFiles.extend(FileUtils.loadListFromFile ('/afs/cern.ch/work/a/amartell/HGCAL/clustering/timingStudies/CMSSW_9_4_0_pre1/src/HGCTimingAnalysis/HGCTiming/test/fileList_testScope/PDG_130_pt10_0PU_allEta.txt') )
#readFiles.extend(FileUtils.loadListFromFile ('/afs/cern.ch/work/a/amartell/HGCAL/clustering/timingStudies/CMSSW_9_4_0_pre1/src/HGCTimingAnalysis/HGCTiming/test/fileList_testScope/PDG_130_pt5_0PU_allEta.txt') )
#readFiles.extend(FileUtils.loadListFromFile ('/afs/cern.ch/work/a/amartell/HGCAL/clustering/timingStudies/CMSSW_9_4_0_pre1/src/HGCTimingAnalysis/HGCTiming/test/fileList_testScope/PDG_130_pt2_0PU_allEta.txt') )



# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('/store/cmst3/user/gpetrucc/l1phase2/081116/DR_PU0/2023D3_guns/ChargedPionGun_PU0_job1.FEVT.root'),
                            #    fileNames = cms.untracked.vstring('file:/data/g/gpetrucc/ChargedPionGun_PU0_job1.FEVT.root'),

                            fileNames = readFiles,
                            #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/amartell/HGCAL/clustering/CMSSW_10_0_4/src/config/partGun_PDGid211_RECO_multi_pt10.root'),
                            #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/amartell/HGCAL/clustering/CMSSW_10_0_4/src/config/partGun_PDGid211_RECO.root'),
                            #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/amartell/HGCAL/clustering/CMSSW_10_0_4/src/config/partGun_PDGid211_RECO_pt20_PU140.root'),
                            secondaryFileNames = cms.untracked.vstring(),
                            inputCommands=cms.untracked.vstring('keep *',
                                                                'drop EcalEBTriggerPrimitiveDigisSorted_simEcalEBTriggerPrimitiveDigis_*_HLT',
                                                                'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
                                                                'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
                                                                'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
                                                                'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
                                                                'drop l1tEMTFTrack2016s_simEmtfDigis__HLT'
                                                                ), 
                            #eventsToProcess = cms.untracked.VEventRange("1:1:257321", "1:1:339076", "1:1:336444")
                            )

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
# customisation of the process.

process.load('RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cfi')
process.load('RecoParticleFlow.HGCTracking.hgcTracks_cff')
process.load("RecoParticleFlow.HGCTracking.hgcClean2D_cfi")
process.load("RecoParticleFlow.HGCTracking.hgcTrajectoryBuilder_cfi")

process.hgcTracks.debugLevel = 3
process.hgcTracks.patternRecoAlgo = "hitsAndClusters"

##to remove realistic simcluster uncomment
#process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterHGC_cfi")
#process.particleFlowClusterHGCal.initialClusteringStep.algoName = cms.string("GenericSimClusterMapper")


#process.run = cms.Path(process.hgcalLayerClusters + process.hgcTracks)
process.run = cms.Path(process.hgcalLayerClusters + process.hgcClean2D + process.hgcTracks)


process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("hgctracks_clusters_130_pt10_highStat.root"),
                               
#    outputCommands = cms.untracked.vstring("drop *", 
#        "keep recoTracks_generalTracks_*_*",
#        "keep *_HGCalRecHit_*_*", 
#        "keep *_particleFlowClusterHGCal_*_*", 
#        "keep *_hgcalLayerClusters_*_*",
#        "keep *_hgcTracks_*_*",
#    ),
    outputCommands = cms.untracked.vstring("keep *") 
)

process.e = cms.EndPath(process.out)

import sys
args = sys.argv[2:] if sys.argv[0] == "cmsRun" else sys.argv[1:]
if args:
    process.source.fileNames = [ args[0] ]
