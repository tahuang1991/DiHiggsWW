import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.maxEvents = cms.untracked.PSet(
	    output = cms.untracked.int32(100)
	)
process.source = cms.Source("PoolSource",
	#skipEvents = cms.untracked.vstring('ProductNotFound'),
	#use voms-proxy-init --voms cms to setup grid
	fileNames = cms.untracked.vstring(
        'file:out_gen.root'
#	'file:/eos/uscms/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A65DBF1F-3174-E411-91A9-002481E94B26.root'
	 #'root://cmsxrootd.fnal.gov//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/000470E0-3B75-E411-8B90-00266CFFA604.root'
		)
	)


#output files

process.Filter = cms.EDFilter('h2to2W2BFilter',
	moduleLabel = cms.untracked.string("genParticles"),
	Wtotau = cms.untracked.bool(True)
 )

#ttbarFilter = cms.Sequence( ttbarto2l2BFilter )


process.p1 = cms.Path(process.Filter)

process.USER = cms.OutputModule("PoolOutputModule",
	    SelectEvents = cms.untracked.PSet(
		        SelectEvents = cms.vstring( 'p1')
			    ),
	        fileName = cms.untracked.string('test_h2filter.root')
	)

#process.p2 = cms.Path(process.filter2)
process.outpath = cms.EndPath(process.USER)

#process.schedule  = cms.Schedule(process.p1)














