#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("GEN")

process.source = cms.Source("MCFileSource",
	#fileNames = cms.untracked.vstring('file:/fdata/hepx/store/user/taohuang/Hhh/HH-bbWW-test.hepmc')
	fileNames = cms.untracked.vstring('file:/home/taohuang/work/Sample/SMH-BBWW.hepmc')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

process.GEN = cms.OutputModule("PoolOutputModule",
	fileName = cms.untracked.string('/fdata/hepx/store/user/taohuang/Hhh/SMH-BBWW.root')
)

process.outpath = cms.EndPath(process.GEN)
