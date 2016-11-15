import FWCore.ParameterSet.Config as cms

PulseTree = cms.EDAnalyzer('PulseTree',
                           doAverage = cms.untracked.bool(False), # average over different events with same gain in first sample for each crystal
                           splitByLumi = cms.untracked.bool(False), # split the above average by lumisections (lumis must be in order in the file)
                           nPedestalSamples = cms.untracked.uint32(3), # number of first samples used to measure the pedestal
                           minAmplitudeForAverage = cms.untracked.double(-9e9), # minimum amplitude of any sample (minus pedestal) to be included in the average
                           filterBx = cms.untracked.vuint32(), # list of BX to select
                           invertBxSelection = cms.untracked.bool(False), # invert BX selection logic (veto instead of keep those in the list)
                           processEB = cms.untracked.bool(True),
                           EBDigiCollection = cms.InputTag("ecalDigis","ebDigis"),
                           processEE = cms.untracked.bool(True),
                           EEDigiCollection = cms.InputTag("ecalDigis","eeDigis"),
                           saveRecHits = cms.untracked.bool(False), # whether to save or not associated rechit information (requires rechits to have been produced)
                           minRecHitAmplitudeWeights = cms.untracked.double(-9e9) # filter on minimum amplitude of uncalibrated rechit with weights method
                           )
