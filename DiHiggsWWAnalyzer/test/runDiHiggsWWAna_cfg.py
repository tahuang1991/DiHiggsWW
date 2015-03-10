import FWCore.ParameterSet.Config as cms

process = cms.Process("DiHiggsWWAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )


process.source = cms.Source("PoolSource",
    # 318 548 618 725 843
    skipEvents = cms.untracked.uint32(317),
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       # 'file:/fdata/hepx/store/user/taohuang/Hhh/HH-bbWW-B3-Gen-1089680.root'
        'file:/fdata/hepx/store/user/taohuang/Hhh/HH-bbWW-B3_Gen_100k_0215.root'
    )
)

process.DiHiggsWWAna = cms.EDAnalyzer('DiHiggsWWAnalyzer',
	verbose = cms.untracked.int32(-1),
        runMMC = cms.bool(True)
)
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService", 
      fileName = cms.string("/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_1_318_0309.root"),
      closeFileFast = cms.untracked.bool(True)
      
  )

process.pDiHiggsWWAna = cms.Path(process.DiHiggsWWAna)
process.pdump = cms.Path(process.dump)
#process.schedule = cms.Schedule(process.pDiHiggsWWAna,process.pdump)
process.schedule = cms.Schedule(process.pDiHiggsWWAna)
