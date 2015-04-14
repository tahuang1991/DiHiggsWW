#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("GEN")

process.source = cms.Source("MCFileSource",
	#fileNames = cms.untracked.vstring('file:/fdata/hepx/store/user/taohuang/Hhh/HH-bbWW-test.hepmc')
	fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/tahuang/DiHiggs/HH-bbWW-B3_100k.hepmc')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

process.GEN = cms.OutputModule("PoolOutputModule",
	fileName = cms.untracked.string('/eos/uscms/store/user/tahuang/DiHiggs/HH-bbWW-10k_0408.root')
)

process.outpath = cms.EndPath(process.GEN)
