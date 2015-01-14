import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    #
    # This is a partial PHYS14 TTJets... dataset, for testing
    #
    fileNames = cms.untracked.vstring(
       '/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root',
       '/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00D3EAF1-3174-E411-A5B2-0025904B144E.root',
       '/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02EF3EFC-0475-E411-A9DB-002590DB9166.root',
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
