#include "TLatex.h"
#include "TTree.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"

#include "OptimizationConstants.hh"
#include "VarCut.hh"

// Draw barrel or endcap
const bool drawBarrel = true;

// Files with signal and background trees (ideally the ntuples
// that were used for TMVA optimization
const TString fnameSignal = "../../ntuples/DYJetsToLL_50ns.root";
const TString signalTreeName = "ntupler/ElectronTree";
const TString fnameBackground = "../../ntuples/TTJets_50ns.root";
const TString backgroundTreeName = "ntupler/ElectronTree";

// Name TMVA output file that contains the pre-computed ROC, etc
const TString tmvaFileNameBarrel = "trainingData/training_results_barrel_pass1_20140728_194000/TMVA_training_results_barrel_pass1_20140728_194000.root";
const TString tmvaFileNameEndcap = "trainingData/training_results_endcap_pass1_20140728_194000/TMVA_training_results_endcap_pass1_20140728_194000.root";

//
// Names of ROOT files that contain working points
//
const int nWorkingPointSets = 3;
// Set 1
const int markerColorSet1 = kRed;
const int markerStyleSet1 = 20;
const TString legendSet1 = "TMVA optimized";
const int nWP = 4;
const TString cutFileNamesBarrelSet1[nWP] = { 
  "cut_repository/cuts_barrel_20140728_194000_WP_Veto.root",
  "cut_repository/cuts_barrel_20140728_194000_WP_Loose.root",
  "cut_repository/cuts_barrel_20140728_194000_WP_Medium.root",
  "cut_repository/cuts_barrel_20140728_194000_WP_Tight.root"
};
const TString cutFileNamesEndcapSet1[nWP] = {
  "cut_repository/cuts_endcap_20140728_194000_WP_Veto.root",
  "cut_repository/cuts_endcap_20140728_194000_WP_Loose.root",
  "cut_repository/cuts_endcap_20140728_194000_WP_Medium.root",
  "cut_repository/cuts_endcap_20140728_194000_WP_Tight.root"
};

// Set 2
const int markerColorSet2 = kRed;
const int markerStyleSet2 = 24;
const TString legendSet2 = "Preliminary safe";
const TString cutFileNamesBarrelSet2[nWP] = { 
  "cut_repository/cuts_barrel_preliminary_gzdp_50ns_WP_Veto.root",
  "cut_repository/cuts_barrel_preliminary_gzdp_50ns_WP_Loose.root",
  "cut_repository/cuts_barrel_preliminary_gzdp_50ns_WP_Medium.root",
  "cut_repository/cuts_barrel_preliminary_gzdp_50ns_WP_Tight.root"
};
const TString cutFileNamesEndcapSet2[nWP] = {
  "cut_repository/cuts_endcap_preliminary_gzdp_50ns_WP_Veto.root",
  "cut_repository/cuts_endcap_preliminary_gzdp_50ns_WP_Loose.root",
  "cut_repository/cuts_endcap_preliminary_gzdp_50ns_WP_Medium.root",
  "cut_repository/cuts_endcap_preliminary_gzdp_50ns_WP_Tight.root"
};

// Set 3
const int markerColorSet3 = kBlue;
const int markerStyleSet3 = 20;
const TString legendSet3 = "EGM 2012";
const TString cutFileNamesBarrelSet3[nWP] = { 
  "cut_repository/cuts_barrel_EGM2012_WP_Veto.root",
  "cut_repository/cuts_barrel_EGM2012_WP_Loose.root",
  "cut_repository/cuts_barrel_EGM2012_WP_Medium.root",
  "cut_repository/cuts_barrel_EGM2012_WP_Tight.root"
};
const TString cutFileNamesEndcapSet3[nWP] = {
  "cut_repository/cuts_endcap_EGM2012_WP_Veto.root",
  "cut_repository/cuts_endcap_EGM2012_WP_Loose.root",
  "cut_repository/cuts_endcap_EGM2012_WP_Medium.root",
  "cut_repository/cuts_endcap_EGM2012_WP_Tight.root"
};

// Forward declarations
TTree *getTreeFromFile(TString fname, TString tname);
void overlayWorkingPoints(TCanvas *c1, 
			  TTree *signalTree, TTree *backgroundTree, 
			  const TString *cutFileNames,
			  int markerColor, int markerStyle, 
			  TLegend *leg, const TString legendText);
void   findEfficiencies(TTree *signalTree, TTree *backgroundTree,
			float &effSignal, float &effBackground, 
			VarCut *cutObject);
//
// Main function
//
void drawROCandWP(){

  TCanvas *c1 = new TCanvas("c1","",10,10,600,600);
  c1->cd();

  TString tmvaFileName = tmvaFileNameBarrel;
  const TString *cutFileNamesSet1 = cutFileNamesBarrelSet1;
  const TString *cutFileNamesSet2 = cutFileNamesBarrelSet2;
  const TString *cutFileNamesSet3 = cutFileNamesBarrelSet3;
  if( !drawBarrel ){
    tmvaFileName = tmvaFileNameEndcap;
    cutFileNamesSet1 = cutFileNamesEndcapSet1;
    cutFileNamesSet2 = cutFileNamesEndcapSet2;
    cutFileNamesSet3 = cutFileNamesEndcapSet3;
  }

  // Open the file with TMVA output
  TFile *tmvaFile = new TFile(tmvaFileName);
  if( !tmvaFile)
    assert(0);

  // 
  //  Draw the ROC curve
  //
  TH1F *hROC = (TH1F*)tmvaFile->Get("Method_Cuts/Cuts/MVA_Cuts_rejBvsS");
  if( !hROC )
    assert(0);

  // Set histogram attributes and draw the ROC curve
  hROC->SetDirectory(0);
  hROC->SetStats(0);
  hROC->SetLineWidth(2);
  hROC->SetTitle("");
  hROC->GetXaxis()->SetTitle("signal efficiency");
  hROC->GetYaxis()->SetTitle("background rejection");
  hROC->GetYaxis()->SetTitleOffset(1.4);
  if( drawBarrel ){
    hROC->GetXaxis()->SetRangeUser(0.6, 1.0);
    hROC->GetYaxis()->SetRangeUser(0.95, 1.0);
  }else{
    hROC->GetXaxis()->SetRangeUser(0.6, 1.0);
    hROC->GetYaxis()->SetRangeUser(0.8, 1.0);
  }
  hROC->Draw("L");

  TString commentText = "barrel electrons";
  if( !drawBarrel )
    commentText = "endcap electrons";
  TLatex *comment = new TLatex(0.2, 0.2, commentText);
  comment->SetNDC(kTRUE);
  comment->Draw();

  c1->Update();

  // 
  // Overlay the cuts
  //

  // First find the TestTree for measuring efficiency and rejection
  printf("\n Take true electrons from %s   tree %s\n", 
	 fnameSignal.Data(), signalTreeName.Data());
  TTree *signalTree = getTreeFromFile( fnameSignal, signalTreeName);
  // Input background tree  
  printf("\n Take background electrons from %s   tree %s\n", 
	 fnameBackground.Data(), backgroundTreeName.Data());
  TTree *backgroundTree = getTreeFromFile( fnameBackground, backgroundTreeName);
  
  // Next, draw all working point sets
  if( nWorkingPointSets >=1 ){

    TLegend *leg = new TLegend(0.15, 0.45, 0.5, 0.7);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    overlayWorkingPoints(c1, signalTree, backgroundTree, cutFileNamesSet1, 
			 markerColorSet1, markerStyleSet1, leg, legendSet1);
    if( nWorkingPointSets >=2 )
      overlayWorkingPoints(c1, signalTree, backgroundTree, cutFileNamesSet2, 
			   markerColorSet2, markerStyleSet2, leg, legendSet2);
    if( nWorkingPointSets >= 3 )
      overlayWorkingPoints(c1, signalTree, backgroundTree, cutFileNamesSet3, 
			   markerColorSet3, markerStyleSet3, leg, legendSet3);
    leg->Draw("same");

  }


  // Save the figure into a file
  TString outname = "figures/plot_ROCandWP_barrel.png";
  if( !drawBarrel )
    outname = "figures/plot_ROCandWP_endcap.png";
  c1->Print(outname);

  return;
};

// Get a given tree from a given file name.
TTree *getTreeFromFile(TString fname, TString tname){

  TFile *file = new TFile( fname );
  TTree *tree     = (TTree*) file->Get(tname);
  
  return tree;
}

// Draw on a given canvas the full set of working points
void   overlayWorkingPoints(TCanvas *c1, 
			    TTree *signalTree, TTree *backgroundTree, 
			    const TString *cutFileNames,
			    int markerColor, int markerStyle, 
			    TLegend *leg, const TString legendText){


  // Now loop over working points
  for(int iwp = 0; iwp<nWP; iwp++){
    
    // Load the working point from a ROOT file
    TFile *cutFile = new TFile(cutFileNames[iwp]);
    if( !cutFile )
      assert(0);
    VarCut *cutObject = (VarCut*)cutFile->Get("cuts");
    if( !cutObject )
      assert(0);
    
    // Compute the efficiencies
    float effSignal, effBackground;
    findEfficiencies(signalTree, backgroundTree, effSignal, effBackground,
		     cutObject);
    printf("Computed eff for cut from %s, effS= %.3f effB= %.3f\n",
	   cutFileNames[iwp].Data(), effSignal, effBackground);
    
    // Make a marker and draw it.
    TMarker *marker = new TMarker(effSignal, 1.0-effBackground, 20);
    marker->SetMarkerSize(2);
    marker->SetMarkerColor(markerColor);
    marker->SetMarkerStyle(markerStyle);
    marker->Draw("same");

    // Add marker to the legend only once. Do not draw the legend here,
    // it is drawn in the main function later
    if( iwp == 0 ){
      if( !leg )
	assert(0);
      leg->AddEntry(marker, legendText, "p");
    }

    c1->Update();
    
    cutFile->Close();
  }

  
}

// Compute signal and background efficiencies for given cuts
void findEfficiencies(TTree *signalTree, TTree *backgroundTree,
		      float &effSignal, float &effBackground, VarCut *cutObject){

  TCut etaCut = "";
  if( drawBarrel ){
    etaCut = Opt::etaCutBarrel;
  }else{
    etaCut = Opt::etaCutEndcap;
  }
  TCut kinematicCuts = Opt::ptCut && etaCut;

  TCut preselectionCuts = kinematicCuts && Opt::otherPreselectionCuts;
  
  TCut signalCuts = preselectionCuts && Opt::trueEleCut;
  TCut backgroundCuts = preselectionCuts && Opt::fakeEleCut;  
 
  TCut selectionCuts = *(cutObject->getCut());

  // printf("\nSelecton cuts:\n");
  // selectionCuts.Print();
  // printf("\nSignal cuts:\n");
  // signalCuts.Print();
  // printf("\nBackground cuts:\n");
  // backgroundCuts.Print();

  effSignal = (1.0*signalTree->GetEntries(selectionCuts && signalCuts) )
    / signalTree->GetEntries(signalCuts);

  effBackground = (1.0*backgroundTree->GetEntries(selectionCuts 
						  && backgroundCuts) )
    / backgroundTree->GetEntries(backgroundCuts);
  
  return;
}
