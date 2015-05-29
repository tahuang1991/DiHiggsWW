import FWCore.ParameterSet.Config as cms
from DiHiggsWW.DiHiggsWWAnalyzer.fileNamesPU import fileNamesPU

def addPileup(process):
    process.mix.input.fileNames = fileNamesPU
    return process

