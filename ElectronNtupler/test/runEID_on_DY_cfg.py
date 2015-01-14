import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    #
    # A few of the files of the DYJetsToLL PHYS14 files
    #                            
    fileNames = cms.untracked.vstring(
       '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root',
       '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/06C61714-7E6C-E411-9205-002590DB92A8.root',
       '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0EAD09A8-7C6C-E411-B903-0025901D493E.root',
    )
)

process.ntupler = cms.EDAnalyzer('ElectronNtupler',
                                 packed = cms.InputTag("packedGenParticles"),
                                 pruned = cms.InputTag("prunedGenParticles"),
                                 pileup = cms.InputTag("addPileupInfo"),
                                 vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 electrons = cms.InputTag("slimmedElectrons"),
                                 rho = cms.InputTag("fixedGridRhoFastjetAll")
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('electron_ntuple.root')
                                   )


process.p = cms.Path(process.ntupler)
