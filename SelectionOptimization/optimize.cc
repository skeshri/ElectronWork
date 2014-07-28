#include "TSystem.h"

#include "optimize.hh"

//
// Main method
//
void optimize(TString cutMaxFileName, TString cutsOutFileNameBase,
	      TString trainingDataOutputBase){

  // Input signal tree
  printf("\n Take true electrons from %s   tree %s\n\n", 
	 fnameSignal.Data(), signalTreeName.Data());
  TTree *signalTree = getTreeFromFile( fnameSignal, signalTreeName, 
				       &fileSignal );

  // Input background tree  
  printf("\n Take background electrons from %s   tree %s\n\n", 
	 fnameBackground.Data(), backgroundTreeName.Data());
  TTree *backgroundTree = getTreeFromFile( fnameBackground, 
					   backgroundTreeName,
					   &fileBackground);
  
  // Configure output details
  TString trainingOutputDir = TString("trainingData/")
    + trainingDataOutputBase;
  printf("The directory where the xml results of the training is:\n");
  printf("         %s\n", trainingOutputDir.Data());
  FileStat_t buf;
  if( gSystem->GetPathInfo(trainingOutputDir.Data(), buf) ){
    printf("     this directory does not exist, creating it.\n");
    gSystem->MakeDirectory(trainingOutputDir.Data());
  }
  TMVA::gConfig().GetIONames().fWeightFileDir = trainingOutputDir;
  
  TString outfileName = trainingOutputDir + TString("/")
    + TString("TMVA_") + trainingDataOutputBase
    + TString(".root");
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  printf("The ROOT file with train/test distributions from TMVA:\n");
  printf("         %s\n", outfileName.Data());

  // Create factory
  TString factoryOptions = "!V:!Silent:Color:DrawProgressBar:Transformations=I";
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification",
					      outputFile, factoryOptions);
  configureFactoryVariables(factory);

  // Define weights and add trees to the factory
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;
  factory->AddSignalTree    ( signalTree,     signalWeight     );
  factory->AddBackgroundTree( backgroundTree, backgroundWeight );
    
  // Configure training and test trees  
  TString trainAndTestOptions = getTrainAndTestOptions();

  TCut signalCuts = "";
  TCut backgroundCuts = "";
  configureCuts(signalCuts, backgroundCuts);

  factory->PrepareTrainingAndTestTree( signalCuts, backgroundCuts,
				       trainAndTestOptions );
  
  // Book the Cuts method with the factory
  TString methodName = "Cuts";
  TString methodOptions = getMethodOptions(cutMaxFileName);
  factory->BookMethod( TMVA::Types::kCuts, methodName,methodOptions);
  
  // Do the work: optimization, testing, and evaluation
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  
  // Save working points into files.
  writeWorkingPoints(factory, cutsOutFileNameBase);

  // Clean up

  outputFile->Close();
  // When the factory is deleted, there appears to be a problem.
  // Things might crash or root does not exit. Commented it out until
  // better understanding.
  // delete factory;

  if( fileSignal != 0 ){
    fileSignal->Close();
  }
  if( fileBackground != 0 ){
    fileBackground->Close();
  }

  return;
}

// Get a given tree from a given file name.
// Note: the **file is the way to return a pointer to a file
// back into the calling method.
TTree *getTreeFromFile(TString fname, TString tname, TFile **file){

  *file = new TFile( fname );
  TTree *tree     = (TTree*) (*file)->Get(tname);
  
  return tree;
}

void configureFactoryVariables(TMVA::Factory *factory){

  // Variables for cut optimization
  printf("Configure factory variables for optimization\n");
  for(int i=0; i<Vars::nVariables; i++){
    TString varName = Vars::variables[i]->nameTmva;
    char varType = Vars::variables[i]->type;
    printf("    add variable %s of the type %c\n", varName.Data(), varType);
    factory->AddVariable( varName, varType );
  }
  
  // Spectator variables
  printf("Configure factory spectator variables\n");
  for(int i=0; i<Vars::nSpectatorVariables; i++){
    TString varName = Vars::spectatorVariables[i]->nameTmva;
    char varType = Vars::spectatorVariables[i]->type;
    printf("    add spectator variable %s of the type %c\n", varName.Data(), varType);
    factory->AddSpectator( varName, varType );
  }
  
}

void configureCuts(TCut &signalCuts, TCut &backgroundCuts){

  // Define all cuts 
 
  TCut etaCut = "";
  if( useBarrel ){
    printf("\n\nTraining for BARREL electrons\n\n");
    etaCut = etaCutBarrel;
  }else{
    printf("\n\nTraining for ENDCAP electrons\n\n");
    etaCut = etaCutEndcap;
  }
  TCut kinematicCuts = ptCut && etaCut;

  TCut preselectionCuts = kinematicCuts && otherPreselectionCuts;
  
  signalCuts = preselectionCuts && trueEleCut;
  backgroundCuts = preselectionCuts && fakeEleCut;  

}

TString getTrainAndTestOptions(){

  TString options = "SplitMode=Random:!V";
  options += ":nTrain_Signal=";
  options += nTrain_Signal;
  options += ":nTrain_Background=";
  options += nTrain_Background;
  options += ":nTest_Signal=";
  options += nTest_Signal;
  options += ":nTest_Background=";
  options += nTest_Background;
 
  printf("INFO: training and test options: %s\n", options.Data());
  return options;
}

TString getMethodOptions(TString cutMaxFileName){

  TString methodOptions = methodCutsBaseOptions;

  // Next, put together cut-specific options
  TString cutsFileName = cutRepositoryDir;
  cutsFileName += "/";
  cutsFileName += cutMaxFileName;

  TFile *cutsFile = new TFile(cutsFileName);
  if( !cutsFile )
    assert(0);
  VarCut *cutMax = (VarCut*)cutsFile->Get("cuts");
  if( !cutMax )
    assert(0);

  // As all cuts are upper cuts, we set the lower cut to -inf
  // Note: we do not have any negative vars, the vars that can be negative
  // are symmetric and enter as abs(XXX).
  for(int i=0; i<Vars::nVariables; i++){
    methodOptions += TString::Format(":VarProp[%d]=FMin",i);
  }
  // Add all cut ranges:
  for(int i=0; i<Vars::nVariables; i++){
    methodOptions += TString::Format(":CutRangeMax[%d]=%.6f", 
				     i, 
				     cutMax->getCutValue(Vars::variables[i]->name));
  }
  
  printf("\nMethod configuration: method options are\n");
  printf("%s\n", methodOptions.Data());
  printf("\n");

  return methodOptions;
}

void writeWorkingPoints(const TMVA::Factory *factory, TString cutsOutFileNameBase){

  TString cutsFileName = cutRepositoryDir;
  cutsFileName += "/";
  cutsFileName += cutsOutFileNameBase;

  // Loop over four working points
  printf("The working points being saved:\n");
  for(int iwp=0; iwp<nWP; iwp++){
    TString cutsFileNameWP = cutsFileName;
    cutsFileNameWP += "_";
    cutsFileNameWP += wpNames[iwp];
    cutsFileNameWP += ".root";
    TFile *cutsFile = new TFile(cutsFileNameWP, "recreate");
    if( !cutsFile )
      assert(0);
    VarCut *cutMax = new VarCut();
    
    const TMVA::MethodCuts *method = dynamic_cast<TMVA::MethodCuts*> (factory->GetMethod("Cuts"));
    if( method == 0 )
      assert(0);

    std::vector <double> cutLo;
    std::vector <double> cutHi;
    method->GetCuts(eff[iwp], cutLo, cutHi);
    // NOTE: this relies on filling the factory with AddVarilables
    // in exactly the same order (using the same loop) 
    // Start with a sanity check:
    if( Vars::nVariables != cutHi.size() )
      assert(0);
    // Now, fill the cut values into the storable object.
    for(uint ivar=0; ivar<cutHi.size(); ivar++){
      cutMax->setCutValue(Vars::variables[ivar]->name, 
			  cutHi.at(ivar));
    }
    printf("   working point %s\n", wpNames[iwp].Data());
    cutMax->print();
    cutMax->Write("cuts");
    cutsFile->Close();
  }

}

