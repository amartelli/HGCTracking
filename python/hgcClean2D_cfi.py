import FWCore.ParameterSet.Config as cms

hgcClean2D = cms.EDProducer("Cl2DLayersCleaner",
                            #srcClusters = cms.InputTag("hgcalLayerClusters::RECO2"),
                            srcClusters = cms.InputTag("hgcalLayerClusters::RECO"),
                            )

