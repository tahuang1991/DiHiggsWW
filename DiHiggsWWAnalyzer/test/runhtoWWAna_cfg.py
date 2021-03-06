import FWCore.ParameterSet.Config as cms
import random
import sys
process = cms.Process("htoWWAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )


process.source = cms.Source("PoolSource",
    # 318 548 618 725 843
    skipEvents = cms.untracked.uint32(58452),
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       # 'file:/fdata/hepx/store/user/taohuang/Hhh/HH-bbWW-B3-Gen-1089680.root'
        #'file:/fdata/hepx/store/user/taohuang/Hhh/HH-bbWW-B3_Gen_100k_0215.root'
        'file:/fdata/hepx/store/user/taohuang/Hhh/out_gen_B3_100k_0714.root'
#	'file:/eos/uscms/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A65DBF1F-3174-E411-91A9-002481E94B26.root'
    )
)



#print "To test array jobs, randint ", random.randint(0,10000000)
#refrootfile = os.getenv( "CMSSW_BASE" ) +"/src/DiHiggsWW/DiHiggsWWAnalyzer/plugins/MMCRefPDF.ROOT"
process.htoWWAna = cms.EDAnalyzer('htoWWAnalyzer',
	verbose = cms.untracked.int32(0),
        finalStates = cms.bool(True),
        Wtotau = cms.bool(False),
	jetLabel = cms.string("ak4GenJetsNoNu"),
	metLabel = cms.string("genMetTrue"),	
	bjetsPt = cms.double(30.0),
	bjetsEta = cms.double(2.40),
	jetsPt = cms.double(20.0),
	jetsEta = cms.double(5.0),
	jetsDeltaR = cms.double(2.00),
	jetleptonDeltaR = cms.double(0.3),
	leptonIso = cms.double(1.0),#not use at generator level
	muonPt1 = cms.double(20.0),
	muonPt2 = cms.double(20.0),
	muonsEta = cms.double(2.40),
	metPt = cms.double(20.0)
)
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService", 
      #fileName = cms.string("/eos/uscms/store/user/tahuang/DiHiggs/ttbar_100k_test_0413_newframe_B3.root"),
      fileName = cms.string("/fdata/hepx/store/user/taohuang/Hhh/htoWWAna_100k_test_0415_B3.root"),
      closeFileFast = cms.untracked.bool(True)
      
  )

process.phtoWWAna = cms.Path(process.htoWWAna)
#process.pdump = cms.Path(process.dump)
#process.schedule = cms.Schedule(process.pdump)
process.schedule = cms.Schedule(process.phtoWWAna)
