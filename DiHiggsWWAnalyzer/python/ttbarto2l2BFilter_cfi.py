import FWCore.ParameterSet.Config as cms

# values tuned also according to slide 3 of :
# https://indico.cern.ch/getFile.py/access?contribId=23&sessionId=2&resId=0&materialId=slides&confId=271548

myttbarto2l2BFilter = cms.EDFilter('ttbarto2l2BFilter',
	moduleLabel = cms.untracked.string("prunedGenParticles"),
	MaxPtLepton = cms.untracked.double(1000),
	MinPtLepton = cms.untracked.double(2.0),
	MaxEtaLepton = cms.untracked.double(5),
	MinEtaLepton = cms.untracked.double(-5),
 )

ttbarFilter = cms.Sequence( myttbarto2l2BFilter )
