import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.maxEvents = cms.untracked.PSet(
	    output = cms.untracked.int32(100)
	)
process.source = cms.Source("PoolSource",
	#skipEvents = cms.untracked.vstring('ProductNotFound'),
	#use voms-proxy-init --voms cms to setup grid
	fileNames = cms.untracked.vstring(
	'file:/eos/uscms/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A65DBF1F-3174-E411-91A9-002481E94B26.root'
	 #'root://cmsxrootd.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/000470E0-3B75-E411-8B90-00266CFFA604.root'
		)
	)


#output files

process.ttbarto2l2BFilter = cms.EDFilter('ttbarto2l2BFilter',
	moduleLabel = cms.untracked.string("prunedGenParticles"),
	MaxPtLepton = cms.untracked.double(1000),
	MinPtLepton = cms.untracked.double(2.0),
	MaxEtaLepton = cms.untracked.double(5),
	MinEtaLepton = cms.untracked.double(-5)
 )

#ttbarFilter = cms.Sequence( ttbarto2l2BFilter )


process.p1 = cms.Path(process.ttbarto2l2BFilter)

process.USER = cms.OutputModule("PoolOutputModule",
	    SelectEvents = cms.untracked.PSet(
		        SelectEvents = cms.vstring( 'p1')
			    ),
	        fileName = cms.untracked.string('test_filtering.root')
	)

#process.p2 = cms.Path(process.filter2)
process.outpath = cms.EndPath(process.USER)

#process.schedule  = cms.Schedule(process.p1)














