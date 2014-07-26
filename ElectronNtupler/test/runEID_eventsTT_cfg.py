import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # 
    # A few files of the TTJets dataset are listed here, for testing
    #
    fileNames = cms.untracked.vstring(
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/00080927-67FD-E311-B049-0025901D4854.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/0025EC82-33FD-E311-B455-0025901D4982.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/002F2689-B7FC-E311-94F1-0025901D46E4.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/00421178-2EFD-E311-ABFC-0025901D4898.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/00428AE4-85FC-E311-97F0-0025901D40A6.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/00678215-05FD-E311-BC19-0025901D40B2.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/006A58D6-CAFC-E311-9ACA-00266CF959BC.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/007E939E-71FD-E311-B43F-0025901D4898.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/008EF5AB-EDFC-E311-9E0B-001EC9B48EB8.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/00B4426C-82FD-E311-A64D-001EC9B48F62.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/0203A95D-70FD-E311-AB62-001EC9B08091.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/02367352-89FC-E311-A988-00266CF94E78.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/028E7CDC-ACFC-E311-8B9D-0025901D40B2.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/02C09CE1-D7FC-E311-A3E2-0025901D46E4.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/02F32F18-76FD-E311-A68B-0025901D4982.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/0426010E-02FD-E311-9412-001EC9B48F26.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/0450282D-3CFD-E311-A4A2-0025901D4892.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/046C2787-35FD-E311-B7E6-0025901D48AA.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU_S14_POSTLS170_V6-v1/00000/048248D5-55FD-E31'
    )
)

process.ntupler = cms.EDAnalyzer('ElectronNtuplerEventStructure',
                                 packed = cms.InputTag("packedGenParticles"),
                                 pruned = cms.InputTag("prunedGenParticles"),
                                 vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 electrons = cms.InputTag("slimmedElectrons"),
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('electron_ntuple_events.root')
                                   )


process.p = cms.Path(process.ntupler)
