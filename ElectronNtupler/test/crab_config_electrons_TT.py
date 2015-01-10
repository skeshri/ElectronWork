from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'electron_MC_analysis_ttbar_phys14_pu20bx25_1'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runEID_eventsTT_cfg.py'

config.section_("Data")
config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
# ??? not sure about the line below. Ken says generally this is in the global
# instance. But anyways this crab config worked
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = False
# config.Data.publishDbsUrl = 'phys03'
# config.Data.publishDataName = 'CRAB3_electron_ntuplizer_ttbar_sc2_1'
config.Data.ignoreLocality = True

config.section_("Site")
# Limit to US Tier2 sites. It appears to give more reliable performance.
config.Site.whitelist = ['T2_US_MIT','T2_US_UCSD','T2_US_Florida','T2_US_Wisconsin','T2_US_Caltech','T2_US_Purdue']
config.Site.storageSite = 'T2_US_Nebraska'
