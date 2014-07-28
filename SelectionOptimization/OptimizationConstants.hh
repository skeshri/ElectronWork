#ifndef OPTIMIZATIONCONSTANTS_HH
#define OPTIMIZATIONCONSTANTS_HH

#include "TCut.h"

namespace Opt {

  // Global parameters
  const bool useBarrel = true;
  // Events to test and train. 
  // To use ALL available events, split 50/50, set
  // all of the nTrain and nTest to 0
  const int nTrain_Signal     = 10000;
  const int nTrain_Background = 10000;
  const int nTest_Signal      = 10000;
  const int nTest_Background  = 10000;
  // Cut repository directory
  const TString cutRepositoryDir = "cut_repository";
  
  // TMVA options for MethodCuts
  const TString methodCutsBaseOptions = "!H:!V:FitMethod=GA:EffMethod=EffSel";
  
  //
  // Constants related to working points of interest
  // The sequence needs to be ordered from most loose to most tight
  const int nWP = 4;
  const float eff[nWP]       = {0.95      , 0.90       , 0.80      , 0.70     };
  const TString wpNames[nWP] = {"WP_Veto", "WP_Loose", "WP_Medium", "WP_Tight"};
  enum WorkingPointIndex       {WP_VETO  , WP_LOOSE  , WP_MEDIUM  , WP_TIGHT  };
  
  // 
  // Input ntuples. Note that the signal and background
  // trees can come from the same file and/or can have the same name.
  // Signal and background are separated by cuts.
  //
  TFile *fileSignal = 0;
  TFile *fileBackground = 0;
  TString fnameSignal = "../../ntuples/DYJetsToLL_50ns.root";
  TString signalTreeName = "ntupler/ElectronTree";
  TString fnameBackground = "../../ntuples/TTJets_50ns.root";
  TString backgroundTreeName = "ntupler/ElectronTree";
  
  //
  // Various preselection and truth matching cuts
  //
  // ( 0: fake, 1: true prompt, 2: true from tau, 3: true non-prompt)
  TCut trueEleCut = "isTrueElectron == 1";
  TCut fakeEleCut = "isTrueElectron == 0 || isTrueElectron == 3";
  
  // Kinematics
  TCut ptCut = "pt>=20 && pt<50";
  TCut etaCutBarrel = " abs(etaSC) < 1.4442 " ;
  TCut etaCutEndcap = " abs(etaSC) > 1.566 " ;
  
  // Anything else
  TCut otherPreselectionCuts = "passConversionVeto && abs(dz)<1";
}

#endif
