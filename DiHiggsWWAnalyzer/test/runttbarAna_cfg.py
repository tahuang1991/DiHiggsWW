import FWCore.ParameterSet.Config as cms
import random
import sys
process = cms.Process("ttbarAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


process.source = cms.Source("PoolSource",
    # 318 548 618 725 843
    #skipEvents = cms.untracked.uint32(617),
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       # 'file:/fdata/hepx/store/user/taohuang/Hhh/HH-bbWW-B3-Gen-1089680.root'
        'file:/fdata/hepx/store/user/taohuang/Hhh/HH-bbWW-B3_Gen_100k_0215.root'
    )
)


try:
	arg1= sys.argv[1]
	print " seed ",arg1
	#random.seed(arg1)
except:
	print "no input seed"

print "To test array jobs, randint ", random.randint(0,10000000)
process.ttbarAna = cms.EDAnalyzer('ttbarAnalyzer',
	verbose = cms.untracked.int32(0),
        finalStates = cms.bool(False),
        runMMC = cms.bool(False),
        mmcset = cms.PSet(
	iterations = cms.untracked.int32(100000),
	seed = cms.int32(random.randint(0,100000000)),#may be ignored since we use can take ievent alone as seed
        weightfromonshellnupt_func = cms.bool(False),
        weightfromonshellnupt_hist = cms.bool(True),
        weightfromoffshellWmass_hist = cms.bool(True),
	RefPDFfile = cms.string("/home/taohuang/work/CMSSW_7_3_1/src/DiHiggsWW/DiHiggsWWAnalyzer/plugins/MMCRefPDF.ROOT")
        )
)
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService", 
      fileName = cms.string("/fdata/hepx/store/user/taohuang/Hhh/ttbar_100k_test_0410_newframe_B3.root"),
      closeFileFast = cms.untracked.bool(True)
      
  )

process.pttbarAna = cms.Path(process.ttbarAna)
process.pdump = cms.Path(process.dump)
#process.schedule = cms.Schedule(process.pttbarAna,process.pdump)
process.schedule = cms.Schedule(process.pttbarAna)
