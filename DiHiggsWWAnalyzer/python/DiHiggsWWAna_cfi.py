import FWCore.ParameterSet.Config as cms

DiHiggsWWAna = cms.EDAnalyzer('DiHiggsWWAnalyzer',
	verbose = cms.untracked.int32(0),
        runMMC = cms.bool(True)
        )
