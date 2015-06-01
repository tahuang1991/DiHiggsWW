import FWCore.ParameterSet.Config as cms
import random
import sys
import os

process = cms.Process("DiHiggsWWAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    # 318 548 618 725 843
    #skipEvents = cms.untracked.uint32(617),
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/fdata/hepx/store/user/taohuang/Hhh/HH-bbWW-B3_Gen_100k_0215.root'
        #'file:/eos/uscms/store/user/tahuang/DiHiggs/HH-bbWW-Gen-10k_0408.root'
    )
)


try:
	arg1= sys.argv[1]
	print " seed ",arg1
	#random.seed(arg1)
except:
	print "no input seed"

print "To test array jobs, randint ", random.randint(0,10000000)
refrootfile = os.getenv( "CMSSW_BASE" ) +"/src/DiHiggsWW/DiHiggsWWAnalyzer/plugins/MMCRefPDF.ROOT"
process.DiHiggsWWAna = cms.EDAnalyzer('DiHiggsWWAnalyzer',
	verbose = cms.untracked.int32(0),
        finalStates = cms.bool(False),
        runMMC = cms.bool(True),
        simulation = cms.bool(True),
        mmcset = cms.PSet(
	iterations = cms.untracked.int32(100000),
	seed = cms.int32(random.randint(0,100000000)),#may be ignored since we use can take ievent alone as seed
        weightfromonshellnupt_func = cms.bool(False),
        weightfromonshellnupt_hist = cms.bool(True),
        weightfromoffshellWmass_hist = cms.bool(True),
        weightfromonoffshellWmass_hist = cms.bool(True),
	RefPDFfile = cms.string("%s"%refrootfile)
        )
)
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService", 
      #fileName = cms.string("/eos/uscms/store/user/tahuang/DiHiggs/DiHiggs_10k_test_0413_newframe_B3.root"),
      fileName = cms.string("/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_10k_test_0413_newframe_B3.root"),
      closeFileFast = cms.untracked.bool(True)
      
  )

process.pDiHiggsWWAna = cms.Path(process.DiHiggsWWAna)
process.pdump = cms.Path(process.dump)
#process.schedule = cms.Schedule(process.pDiHiggsWWAna,process.pdump)
process.schedule = cms.Schedule(process.pDiHiggsWWAna)
