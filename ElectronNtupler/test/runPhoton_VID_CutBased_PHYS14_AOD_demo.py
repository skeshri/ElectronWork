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
# A test sample:
       '/store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/000B1591-F96E-E411-8885-00266CFFC198.root',
#       '/store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0028736A-346C-E411-BBC8-1CC1DE1CE56C.root',
#       '/store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/005A3A97-3B6C-E411-9BDA-1CC1DE1CE01A.root',
#  A local file
#    'file:/afs/cern.ch/user/i/ikrav/workspace/releases-git/CMSSW_7_2_0/src/ElectronWork/ElectronNtupler/test/00A074A5-BF72-E411-B455-003048F02CBE.root'
 )
)

#
# Set up photon ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer
switchOnVIDPhotonIdProducer(process)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

#
# Configure an example module for user analysis of photons
#

process.ntupler = cms.EDAnalyzer('PhotonNtuplerVIDDemoAOD',
                                 photons = cms.InputTag("gedPhotons"),
                                 vertices = cms.InputTag("offlinePrimaryVertices"),
                                 rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                 genParticles = cms.InputTag("genParticles"),
                                 phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V1-standalone-loose"),
                                 phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V1-standalone-medium"),
                                 phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V1-standalone-tight"),
                                 full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
                                 phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
                                 phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
                                 phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation")
                                )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('photon_ntuple.root')
                                   )
process.p = cms.Path(process.egmPhotonIDSequence * process.ntupler)
