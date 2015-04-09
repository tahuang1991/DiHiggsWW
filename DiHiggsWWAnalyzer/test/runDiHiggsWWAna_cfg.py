import FWCore.ParameterSet.Config as cms
import random
import sys
process = cms.Process("DiHiggsWWAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )


process.source = cms.Source("PoolSource",
    # 318 548 618 725 843
    skipEvents = cms.untracked.uint32(617),
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
process.DiHiggsWWAna = cms.EDAnalyzer('DiHiggsWWAnalyzer',
	verbose = cms.untracked.int32(-1),
        runMMC = cms.bool(True),
	iterations = cms.untracked.int32(100000),
	seed = cms.int32(random.randint(0,100000000)),#may be ignored since we use can take ievent alone as seed
        finalStates = cms.bool(False),
        weightfromonshellnupt_func = cms.bool(False),
        weightfromonshellnupt_hist = cms.bool(True),
	RefPDFfile = cms.string("/home/taohuang/work/CMSSW_7_3_1/src/DiHiggsWW/DiHiggsWWAnalyzer/plugins/MMCRefPDF.ROOT")
)
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService", 
      fileName = cms.string("/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_1_618_0406_B3.root"),
      closeFileFast = cms.untracked.bool(True)
      
  )

process.pDiHiggsWWAna = cms.Path(process.DiHiggsWWAna)
process.pdump = cms.Path(process.dump)
#process.schedule = cms.Schedule(process.pDiHiggsWWAna,process.pdump)
process.schedule = cms.Schedule(process.pDiHiggsWWAna)
