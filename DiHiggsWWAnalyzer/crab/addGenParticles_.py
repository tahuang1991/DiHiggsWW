#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("GENP")

process.load('Configuration.StandardSequences.Services_cff')

print "run addGenParticles.py  "
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring('file:HH-bbWW-100_0407.root')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

process.GEN = cms.OutputModule("PoolOutputModule",
	fileName = cms.untracked.string('HH-bbWW-Gen-100_0407_testcrab.root')
)

process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')
process.genParticles.src = cms.InputTag("source")
process.genps = cms.Path(process.genParticles)
process.outpath = cms.EndPath(process.GEN)
