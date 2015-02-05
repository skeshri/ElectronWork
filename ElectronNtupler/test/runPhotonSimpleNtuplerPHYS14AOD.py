import FWCore.ParameterSet.Config as cms

process = cms.Process("TestPhotons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for PHYS14 scenario PU4bx50 : global tag is ???
#    for PHYS14 scenario PU20bx25: global tag is PHYS14_25_V1
#  as a rule, find the global tag in the DAS under the Configs for given dataset
process.GlobalTag.globaltag = 'PHYS14_25_V1::All'

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(
   #
   # Just a handful of files from the dataset are listed below, for testing
   #
# The right sample:
       '/store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/000B1591-F96E-E411-8885-00266CFFC198.root',
       '/store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0028736A-346C-E411-BBC8-1CC1DE1CE56C.root',
       '/store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/005A3A97-3B6C-E411-9BDA-1CC1DE1CE01A.root',
# Another dataset example
#       '/store/mc/Phys14DR/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00A074A5-BF72-E411-B455-003048F02CBE.root',
#       '/store/mc/Phys14DR/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04396056-6872-E411-A57B-00266CF32BC4.root',
#       '/store/mc/Phys14DR/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/06A4C6EE-2873-E411-AA26-002590AC4C66.root',
#  A local file
#    'file:/afs/cern.ch/user/i/ikrav/workspace/releases-git/CMSSW_7_2_0/src/ElectronWork/ElectronNtupler/test/00A074A5-BF72-E411-B455-003048F02CBE.root'
 )
)

#
# Run some stuff to produce value maps needed for IDs
#

process.photonIDValueMapProducer = cms.EDProducer('PhotonIDValueMapProducer',
                                          ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                                          eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                                          particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
                                          vertices = cms.InputTag("offlinePrimaryVertices"),
                                          pfCandidates = cms.InputTag("particleFlow"),
                                          esReducedRecHitCollection = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES"),
                                          src = cms.InputTag('gedPhotons'),
                                          dataFormat = cms.string('RECO')
)

#
# Configure an example module for user analysis of electrons
#

process.ntupler = cms.EDAnalyzer('PhotonNtuplerAOD',
                                 photons = cms.InputTag("gedPhotons"),
                                 vertices = cms.InputTag("offlinePrimaryVertices"),
                                 rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                 genParticles = cms.InputTag("genParticles"),
                                 full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
                                 phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
                                 phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
                                 phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation")
                                )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('photon_ntuple.root')
                                   )

process.p = cms.Path(process.photonIDValueMapProducer * process.ntupler)
