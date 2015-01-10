import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("ElectronWork/ElectronNtupler/TTJets_CSA14_sc2")

process.ntupler = cms.EDAnalyzer('ElectronNtuplerEventStructure',
                                 packed = cms.InputTag("packedGenParticles"),
                                 pruned = cms.InputTag("prunedGenParticles"),
                                 vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 pileup = cms.InputTag("addPileupInfo"),
                                 electrons = cms.InputTag("slimmedElectrons"),
                                 rho = cms.InputTag("fixedGridRhoFastjetAll")
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('output.root')
                                   )


process.p = cms.Path(process.ntupler)
