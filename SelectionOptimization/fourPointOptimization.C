#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "OptimizationConstants.hh"
#include "optimize.hh"

void fourPointOptimization(){
  
  // Define source for the initial cut range
  TString startingCutMaxFileName 
    = "cuts_barrel_eff_0999_20140727_165000.root";
  if( !useBarrel )
    startingCutMaxFileName 
      = "cuts_endcap_eff_0999_20140727_165000.root";

  TString namePrefix = "cuts_barrel_";
  if( !useBarrel )
    namePrefix = "cuts_endcap_";
  TString namePass[nWP] = {"pass1_","pass2_","pass3_","pass4_"};
  TString nameTime = "20140727_182500";

  for( int ipass = 0; ipass<nWP; ipass++){

    // This string is the file name that contains the ROOT file
    // with the VarCut object that defines the range of cut variation.
    // Note: for each subsequence pass, the previous working point
    // is used. For the first pass, the 99.9% efficient cut range is used.
    TString cutMaxFileName = startingCutMaxFileName;
    if( ipass > 0 ){
      cutMaxFileName = namePrefix + namePass[ipass-1] + nameTime + TString("_")
	+ wpNames[ipass-1] + TString(".root");
    }
    
    // The string below is used to construct the file names
    // to save the cut objects
    TString cutOutputBase = namePrefix;
    cutOutputBase += namePass[ipass];
    cutOutputBase += nameTime;

    // This string will be used to construct the dir for the output
    // of TMVA: the dir for weights and the filename for diagnostics
    TString trainingDataOutputBase = "training_results_";
    trainingDataOutputBase += namePass[ipass];
    trainingDataOutputBase += nameTime;

    printf("\n-----------------------------------------------------------------\n");
    printf("\n");
    printf("    Run new optimization pass  \n");
    printf("\n");
    printf(" Input file that defines cut optimization limits:");
    printf(" %s\n", cutMaxFileName.Data());
    printf(" Base for the file name of output cuts          :");
    printf(" %s\n", cutOutputBase.Data());
    printf("------------------------------------------------------------------\n\n");

    optimize(cutMaxFileName, cutOutputBase, trainingDataOutputBase);    

  }
 
  // Finally, create the files containing the working point
  // by copying the appropriate pass files.
  // The first working point is output of pass1, the second of pass2, etc.
  // All other working poitns of all passes are ignored
  printf("\n");
  printf("====================================================\n");
  printf("Final definition of working points\n");
  printf("====================================================\n");
  printf("Copy files with working points info into the final locations\n");
  for(int i=0; i<nWP; i++){

    TString wpPassFileName = cutRepositoryDir + TString("/")
      + namePrefix + namePass[i] + nameTime + TString("_")
      + wpNames[i] + TString(".root");

    TString wpFinalFileName =  cutRepositoryDir + TString("/")
      + namePrefix + nameTime + TString("_")
      + wpNames[i] + TString(".root");

    TString copyCommand = TString::Format("cp %s %s", 
					  wpPassFileName.Data(), 
					  wpFinalFileName.Data());

    gSystem->Exec(copyCommand);

    printf("\nFinal definition for working point %s\n", wpNames[i].Data());
    printf(" file name:   %s\n", wpFinalFileName.Data());
    TFile *fwp = new TFile(wpFinalFileName);
    VarCut *thisWP = (VarCut*)fwp->Get("cuts");
    if( thisWP != 0 )
      thisWP->print();
    else
      printf("???? not found????\n");
    fwp->Close();
    
  }
 
}

