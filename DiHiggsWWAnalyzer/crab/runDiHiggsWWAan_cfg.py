import FWCore.ParameterSet.Config as cms

process = cms.Process("DiHiggsWWAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

#input files
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/uscms_data/d3/tahuang/CMSSW_6_2_0_SLHC16/src/GEMCode/SimMuL1/debug/PU0_10k_Pt20_out_L1.root'
    )
)

#output files
process.TFileService = cms.Service("TFileService",
	    fileName = cms.string("out_Ana.root")
)



process.DiHiggsWWAna = cms.EDAnalyzer('DiHiggsWWAnalyzer')
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.pDiHiggsWWAna = cms.Path(process.DiHiggsWWAna)
process.pdump = cms.Path(process.dump)
process.schedule = cms.Schedule(process.pDiHiggsWWAna,process.pdump)
