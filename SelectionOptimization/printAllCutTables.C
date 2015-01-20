#include "TString.h"
#include "TFile.h"
#include "VarCut.hh"
#include "Variables.hh"

const int nWP = 4;
const TString workingPointNames[nWP] = {
  "veto",
  "loose",
  "medium",
  "tight"
};

const TString barrelCutDir = 
  "/home/snow/rkamal/PHYS14/optimizeBarrel/SelectionOptimization/cut_repository/";
const TString barrelCutFiles[nWP] = {
  "cuts_barrel_20141205_205800_WP_Veto.root",
  "cuts_barrel_20141205_205800_WP_Loose.root",
  "cuts_barrel_20141205_205800_WP_Medium.root",
  "cuts_barrel_20141205_205800_WP_Tight.root"
};

const TString endcapCutDir = 
  "/home/snow/rkamal/PHYS14/optimizeENDCAP/SelectionOptimization/cut_repository/";
const TString endcapCutFiles[nWP] = {
  "cuts_endcap_20141206_163000_WP_Veto.root",
  "cuts_endcap_20141206_163000_WP_Loose.root",
  "cuts_endcap_20141206_163000_WP_Medium.root",
  "cuts_endcap_20141206_163000_WP_Tight.root",
};

enum WorkingPointNames {
  WP_VETO = 0,
  WP_LOOSE,
  WP_MEDIUM,
  WP_TIGHT
};

void printAllCutTables(){

  TFile *barrelFiles[nWP];
  TFile *endcapFiles[nWP];
  VarCut *barrelCuts[nWP];
  VarCut *endcapCuts[nWP];
  for(int i=0; i<nWP; i++){

    // Load all barrel cuts
    TString barrelFileName = barrelCutDir + barrelCutFiles[i];
    barrelFiles[i] = new TFile(barrelFileName);
    if( !barrelFiles[i] ){
      printf("Could not open the file %s\n", barrelFileName.Data());
      assert(0);
    }
    barrelCuts[i] = (VarCut*) barrelFiles[i]->Get("cuts");
    if( !barrelCuts[i] ){
      printf("Could not find the cuts object in file %s\n", barrelFileName.Data());
      assert(0);
    }

    // Load all endcap cuts
    TString endcapFileName = endcapCutDir + endcapCutFiles[i];
    endcapFiles[i] = new TFile(endcapFileName);
    if( !endcapFiles[i] ){
      printf("Could not open the file %s\n", endcapFileName.Data());
      assert(0);
    }
    endcapCuts[i] = (VarCut*) endcapFiles[i]->Get("cuts");
    if( !endcapCuts[i] ){
      printf("Could not find the cuts object in file %s\n", endcapFileName.Data());
      assert(0);
    }

  }

  // Print barrel table for the twiki
  printf("\nAll cuts for barrel in a twiki table format:\n\n");
  printf("|                                    |   Veto     |    Loose   |   Medium   |  Tight   |\n");
  for(int i=0; i< Vars::nVariables; i++){
    TString variableName     = Vars::variables[i]->name.Data();
    TString variableTmvaName = Vars::variables[i]->nameTmva.Data();
    printf("|  %30s <  |  %f  |  %f  |  %f  | %f |\n", 
	   variableTmvaName.Data(), 
	   barrelCuts[WP_VETO  ]->getCutValue(variableName),
	   barrelCuts[WP_LOOSE ]->getCutValue(variableName),
	   barrelCuts[WP_MEDIUM]->getCutValue(variableName),
	   barrelCuts[WP_TIGHT ]->getCutValue(variableName)
	   );
  }
  printf("\n");
  
  // Print endcap table for the twiki
  printf("\nAll cuts for endcap in a twiki table format:\n\n");
  printf("|                                    |   Veto     |    Loose   |   Medium   |  Tight   |\n");
  for(int i=0; i< Vars::nVariables; i++){
    TString variableName     = Vars::variables[i]->name.Data();
    TString variableTmvaName = Vars::variables[i]->nameTmva.Data();
    printf("|  %30s <  |  %f  |  %f  |  %f  | %f |\n", 
	   variableTmvaName.Data(), 
	   endcapCuts[WP_VETO  ]->getCutValue(variableName),
	   endcapCuts[WP_LOOSE ]->getCutValue(variableName),
	   endcapCuts[WP_MEDIUM]->getCutValue(variableName),
	   endcapCuts[WP_TIGHT ]->getCutValue(variableName)
	   );
  }
  printf("\n");
  
}
