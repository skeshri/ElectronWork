import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
# 
# Only some files are listed here, for debugging purposes. The dataset
# name is  
#  /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM
#

       '/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/029EA4C3-A007-E411-BE0B-D4AE529D9537.root',
       '/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/041B1709-A107-E411-BE7B-001E675A659A.root',
       '/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/061362C4-A007-E411-B6CF-001E67586629.root',
       '/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/088E0822-A207-E411-9941-90B11C1453E1.root',
       '/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0A097C5B-9607-E411-A2A0-0026181D291E.root',
       '/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0C4273AF-A107-E411-9289-D4AE52AAF583.root',
       '/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0CD3EB1A-A007-E411-B601-D4AE52AAF583.root',
       '/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/104CA3E5-A007-E411-A18F-001517E74088.root',
       '/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/10E1138B-9607-E411-9EBF-001517E73B9C.root',
       '/store/mc/Spring14miniaod/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/1AF5235F-8F07-E411-9342-90B11C04FE0C.root',
] );


secFiles.extend( [
               ] )
