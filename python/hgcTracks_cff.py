import FWCore.ParameterSet.Config as cms

import TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi 
hitCollectorHGC = TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi.Chi2MeasurementEstimator.clone(
    ComponentName = cms.string('hitCollectorHGC'),
    MaxChi2 = cms.double(40.0), ## was 30 ## TO BE TUNED
    nSigma  = cms.double(3.),    ## was 3  ## TO BE TUNED 
)

from RecoLocalTracker.SiStripClusterizer.SiStripClusterChargeCut_cfi import SiStripClusterChargeCutNone

from RecoHGCal.HGCTracking.hgcTracks_cfi import *
