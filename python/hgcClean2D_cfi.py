import FWCore.ParameterSet.Config as cms

hgcClean2D = cms.EDProducer("Cl2DLayersCleaner",
                            #srcClusters = cms.InputTag("hgcalLayerClusters::RECO2"), ## test reco with kernel for 2DCl and other cuts
                            srcClusters = cms.InputTag("hgcalLayerClusters::RECO"),
                            )

