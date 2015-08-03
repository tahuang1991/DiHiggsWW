import FWCore.ParameterSet.Config as cms
import random
import sys
process = cms.Process("tt2WWbbAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )


process.source = cms.Source("PoolSource",
    # 318 548 618 725 843
    #skipEvents = cms.untracked.uint32(617),
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
	'file:input.root'
        #'file:/fdata/hepx/store/user/taohuang/Hhh/out_gen_B3_100k_0714.root'
#	'file:/eos/uscms/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A65DBF1F-3174-E411-91A9-002481E94B26.root'
    )
)

#inputfiles
from DiHiggsWW.DiHiggsWWAnalyzer.InputFileHelpers import *
#inputdir = ['/eos/uscms/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/']
inputdir = ['/fdata/hepx/store/user/tahuang/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/TTbarto2l2B_13TeV-madgraph-tauola_PU20bx25_PHYS14_25_V1_brazos_25M/f871079033b0f7f7cadefb46854eddad/']

process = useInputDir(process, inputdir, True)


#print "To test array jobs, randint ", random.randint(0,10000000)
#refrootfile = os.getenv( "CMSSW_BASE" ) +"/src/DiHiggsWW/DiHiggsWWAnalyzer/plugins/MMCRefPDF.ROOT"
process.tt2WWbbAna = cms.EDAnalyzer('tt2WWBBAnalyzer',
	verbose = cms.untracked.int32(0),
        finalStates = cms.bool(True),
        Wtotau = cms.bool(False),
	#jetLabel = cms.string("ak4GenJetsNoNu"),
	jetLabel = cms.string("ak4GenJets"),
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
      fileName = cms.string("/fdata/hepx/store/user/taohuang/Hhh/tt2WWbbAna_100k_test_0415_B3.root"),
      closeFileFast = cms.untracked.bool(True)
      
  )

process.ptt2WWbbAna = cms.Path(process.tt2WWbbAna)
#process.pdump = cms.Path(process.dump)
#process.schedule = cms.Schedule(process.pdump)
process.schedule = cms.Schedule(process.ptt2WWbbAna)
