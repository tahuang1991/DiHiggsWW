import FWCore.ParameterSet.Config as cms

process = cms.Process("DiHiggsWWAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/fdata/hepx/store/user/taohuang/Hhh/HH-bbWW-B6_Gen_all.root'
    )
)

process.DiHiggsWWAna = cms.EDAnalyzer('DiHiggsWWAnalyzer')
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService", 
      fileName = cms.string("DiHiggs_Test.root"),
      closeFileFast = cms.untracked.bool(True)
  )

process.pDiHiggsWWAna = cms.Path(process.DiHiggsWWAna)
process.pdump = cms.Path(process.dump)
process.schedule = cms.Schedule(process.pDiHiggsWWAna,process.pdump)
