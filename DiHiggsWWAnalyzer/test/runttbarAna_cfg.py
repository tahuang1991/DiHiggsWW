import FWCore.ParameterSet.Config as cms
import random
import sys
import os
process = cms.Process("ttbarAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


process.source = cms.Source("PoolSource",
    # 318 548 618 725 843
    #skipEvents = cms.untracked.uint32(8),
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       # 'file:/fdata/hepx/store/user/taohuang/Hhh/HH-bbWW-B3-Gen-1089680.root'
        'file:out_filter.root'
	#'file:/eos/uscms/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A65DBF1F-3174-E411-91A9-002481E94B26.root'
    )
)

#inputfiles
from DiHiggsWW.DiHiggsWWAnalyzer.InputFileHelpers import *
#inputdir = ['/eos/uscms/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/']
inputdir = ['/fdata/hepx/store/user/tahuang/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/TTbarto2l2B_13TeV-madgraph-tauola_PU20bx25_PHYS14_25_V1_brazos_25M/f871079033b0f7f7cadefb46854eddad/']

process = useInputDir(process, inputdir, True)


try:
	arg1= sys.argv[1]
	print " seed ",arg1
	#random.seed(arg1)
except:
	print "no input seed"

#print "To test array jobs, randint ", random.randint(0,10000000)
refrootfile = os.getenv( "CMSSW_BASE" ) +"/src/DiHiggsWW/DiHiggsWWAnalyzer/plugins/MMCRefPDF.ROOT"
process.ttbarAna = cms.EDAnalyzer('ttbarAnalyzer',
       # SkipEvent = cms.untracked.vstring('ProductNotFound'),
	verbose = cms.untracked.int32(0),
        finalStates = cms.bool(False),
	genparticleLabel = cms.string("genParticles"),
	jetLabel = cms.string("ak4GenJets"),
	metLabel = cms.string("genMetTrue"),	
	jetsPt = cms.double(30.0),
	jetsEta = cms.double(2.50),
	bjetsPt = cms.double(20.0),
	bjetsEta = cms.double(5.0),
	jetsDeltaR = cms.double(2.00),
	muonPt1 = cms.double(10.0),
	muonPt2 = cms.double(10.0),
	muonsEta = cms.double(2.40),
	metPt = cms.double(20.0),
        runMMC = cms.bool(False),
	simulation = cms.bool(True),
        mmcset = cms.PSet(
	iterations = cms.untracked.int32(100000),
	seed = cms.int32(random.randint(0,100000000)),#may be ignored since we use can take ievent alone as seed
        weightfromonshellnupt_func = cms.bool(False),
        weightfromonshellnupt_hist = cms.bool(True),
        weightfromoffshellWmass_hist = cms.bool(True),
        weightfromonoffshellWmass_hist = cms.bool(True),
	useMET = cms.bool(True),
	RefPDFfile = cms.string("%s"%refrootfile)
        )
)
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService", 
      #fileName = cms.string("/eos/uscms/store/user/tahuang/DiHiggs/ttbar_100k_test_0413_newframe_B3.root"),
      fileName = cms.string("/fdata/hepx/store/user/taohuang/ttbarAna/ttbarAna_100_test_0608_B3.root"),
      closeFileFast = cms.untracked.bool(True)
      
  )

process.pttbarAna = cms.Path(process.ttbarAna)
#process.pdump = cms.Path(process.dump)
#process.schedule = cms.Schedule(process.pdump)
process.schedule = cms.Schedule(process.pttbarAna)
