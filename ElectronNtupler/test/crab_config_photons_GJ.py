from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'photon_MC_analysis_gjet_PU20bx25_4'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runPhotonSimpleNtuplerPHYS14AOD.py'

config.section_("Data")
config.Data.inputDataset = '/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM'
# ??? not sure about the line below. Ken says generally this is in the global
# instance. But anyways this crab config worked
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = False
# config.Data.publishDbsUrl = 'phys03'
# config.Data.publishDataName = 'CRAB3_electron_ntuplizer_something'
config.Data.ignoreLocality = True

config.section_("Site")
# Limit to US Tier2 sites. It appears to give more reliable performance.
config.Site.whitelist = ['T2_US_MIT','T2_US_UCSD','T2_US_Florida','T2_US_Wisconsin','T2_US_Caltech','T2_US_Purdue']
config.Site.storageSite = 'T2_US_Nebraska'
