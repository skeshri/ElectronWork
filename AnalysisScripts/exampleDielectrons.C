#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TTree.h"
#include "TList.h"
#include "TString.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TLegend.h"

#include <vector>

// For this exercise, we use two MC samples: the signal sample
// and a background-reach sample. Both have the same ntuple structure.
//
// Signal sample: DYToLL
const TString fileNameSignal = "/home/hep/ikrav/releases-git/CMSSW_7_0_6_patch3/src/ElectronWork/ElectronNtupler/test/output.root";
// Background sample: TTbar
const TString fileNameBackground = "/home/hep/ikrav/releases-git/CMSSW_7_0_6_patch3/src/ElectronWork/ElectronNtupler/test/output_TT.root";

const TString treeName = "ntupler/ElectronTree";

const bool verbose = false;
const bool smallEventCount = false;

const float boundaryBarrelEndcap = 1.479;

enum TruthMatching {ELECTRON_FAKE, ELECTRON_TRUE_PROMPT, 
		    ELECTRON_TRUE_FROM_TAU, ELECTRON_TRUE_SECONDARY}; // secondary does not include from tau

// Selection cuts
// Kinematics
const float ptCut = 20; 
// ID/ISO/etc selection
const float cutIsoBarrel = 0.40;
const float cutSeeBarrel = 0.12;

const float cutIsoEndcap = 0.40;
const float cutSeeEndcap = 0.35;

// Other constants
const double electronMass = 0.000511;

// Forward declarations
void setHistogramAttributes(TH1F *hist, TString XTitle);

//
// Main program
//

void exampleDielectrons(){

  // This statement below should not be needed, but in one particular node I had to
  // add it, somehow the vector header was not loaded automatically there.
  gROOT->ProcessLine("#include <vector>"); 

  // Book histograms
  TH1F *hDielectronMass = new TH1F("hDielectronMass","",60, 60, 120);
  TH1F *hDielectronMassCleaned = new TH1F("hDielectronMassCleaned","",
					  60, 60, 120);

  //
  // Open a file and find the tree with electron data
  //
  TFile *fileSignal     = new TFile(fileNameSignal);
  TFile *fileBackground = new TFile(fileNameBackground);
  if( !fileSignal || !fileBackground ){
    printf("Failed to open one of the input files, check\n   %s\n   %s", 
	   fileNameSignal.Data(), fileNameBackground.Data());
    assert(0);
  }
  TTree *treeSignal     = (TTree*)fileSignal->Get(treeName);
  TTree *treeBackground = (TTree*)fileBackground->Get(treeName);
  if( !treeSignal || !treeBackground ){
    printf("Failed to find the tree %s\n", treeName.Data() );
    assert(0);
  }
  // Finall, merge the signal and background trees into a single sample
  // We open a new file: since we create a new tree below,
  // we avoid making the new tree memory-resident by opening this file here.
  TFile *fileOut = new TFile("tmp.root","recreate");
  fileOut->cd();
  TList *list = new TList;
  list->Add(treeSignal);
  list->Add(treeBackground);
  TTree *tree = TTree::MergeTrees(list);
  tree->SetName("mixDYandTT");

  // 
  // Configure reading the information of interest from the TTree,
  // in this example we are working with events containing vectors.
  // A good ROOT reference: 
  //      http://root.cern.ch/root/html/tutorials/tree/hvector.C.html
  //

  // 
  // Set up the branches of interest
  //
  // Declare variables
  //
  // Event-level variables:
  int nEle; // the number of reconstructed electrons in the event
  // Per-eletron variables
  // Kinematics
  std::vector <float> *elePt = 0;         // electron PT
  std::vector <float> *eleEtaSC = 0;      // supercluser eta
  std::vector <float> *elePhiSC = 0;      // supercluser phi
  // ID related variables
  std::vector <float> *eleFull5x5SigmaIEtaIEta = 0;  // sigma_ieta_iet  
  std::vector <float> *eleRelIsoWithDBeta = 0;    // pre-computed (ch+max(0.0,nh+pho-0.5*chPU)/pt from above
  
  // Declare branches
  TBranch *b_nEle = 0;
  TBranch *b_elePt = 0;
  TBranch *b_eleEtaSC = 0;
  TBranch *b_elePhiSC = 0;
  TBranch *b_eleFull5x5SigmaIEtaIEta = 0;
  TBranch *b_eleRelIsoWithDBeta = 0;

  // Connect variables and branches to the tree with the data
  tree->SetBranchAddress("nEle", &nEle, &b_nEle);
  tree->SetBranchAddress("pt", &elePt, &b_elePt);
  tree->SetBranchAddress("etaSC", &eleEtaSC, &b_eleEtaSC);
  tree->SetBranchAddress("phiSC", &elePhiSC, &b_elePhiSC);
  tree->SetBranchAddress("full5x5_sigmaIetaIeta", &eleFull5x5SigmaIEtaIEta, &b_eleFull5x5SigmaIEtaIEta);
  tree->SetBranchAddress("relIsoWithDBeta", &eleRelIsoWithDBeta, &b_eleRelIsoWithDBeta);

  // 
  // Loop over events
  //
  UInt_t maxEvents = tree->GetEntries();
  if( smallEventCount )
    maxEvents = 1000;
  if(verbose)
    printf("Start loop over events, total events = %lld\n", 
	   tree->GetEntries() );
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0)
      printf(".");
    Long64_t tentry = tree->LoadTree(ievent);
    
    // Load the value of the number of the electrons in the event    
    b_nEle->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of electrons %u\n", ievent, nEle);
    
    // To construct dielectrons, we need at least two electrons
    // so skip the event if there are not enough electrons
    if( nEle < 2 ) 
      continue;

    // Get data for all electrons in this event, only vars of interest
    b_elePt->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_elePhiSC->GetEntry(tentry);
    b_eleFull5x5SigmaIEtaIEta->GetEntry(tentry);
    b_eleRelIsoWithDBeta->GetEntry(tentry);

    // Nested loops over the electrons
    for(int iele1 = 0; iele1 < nEle; iele1++){

      // Check kinematics:
      if( !(elePt->at(iele1) > ptCut) )
	continue;

      // Cuts to suppress fakes and non-prompt electrons:
      bool passSelection1 = false;
      bool isBarrel1 = fabs( eleEtaSC->at(iele1) ) < boundaryBarrelEndcap ? true : false;
      float cutSee, cutIso;
      cutIso = cutIsoBarrel; 
      cutSee = cutSeeBarrel; 
      if( !isBarrel1 ){
	cutIso = cutIsoEndcap;
	cutSee = cutSeeEndcap;
      }
      if( eleRelIsoWithDBeta->at(iele1) < cutIso 
	  && eleFull5x5SigmaIEtaIEta->at(iele1) < cutSee )
	passSelection1 = true;

      // Loop over the second electron candidates
      for(int iele2 = iele1 + 1; iele2 < nEle; iele2++){

	// Check kinematics:
	if( !(elePt->at(iele2) > ptCut) )
	  continue;
	
	// Cuts to suppress fakes and non-prompt electrons:
	bool passSelection2 = false;
	bool isBarrel2 = fabs( eleEtaSC->at(iele2) ) < boundaryBarrelEndcap ? true : false;
	cutIso = cutIsoBarrel; 
	cutSee = cutSeeBarrel; 
	if( !isBarrel2 ){
	  cutIso = cutIsoEndcap;
	  cutSee = cutSeeEndcap;
	}
	if( eleRelIsoWithDBeta->at(iele2) < cutIso 
	    && eleFull5x5SigmaIEtaIEta->at(iele2) < cutSee )
	  passSelection2 = true;
	
	// At this point we have a dielectron, compute its invariant 
	// mass and save into histograms
	TLorentzVector e1, e2;
	e1.SetPtEtaPhiM( elePt->at(iele1), eleEtaSC->at(iele1),
			 elePhiSC->at(iele1), electronMass);
	e2.SetPtEtaPhiM( elePt->at(iele2), eleEtaSC->at(iele2),
			 elePhiSC->at(iele2), electronMass);
	float invMass = (e1+e2).M();
	
	hDielectronMass->Fill(invMass);
	if(passSelection1 && passSelection2)
	  hDielectronMassCleaned->Fill(invMass);

      } // end loop over the second electron
    } // end loop over the first electron

  } // end loop over events
  printf("\n");

  // Final clean-up, as recommended by ROOT tutorials
  tree->ResetBranchAddresses();

  fileSignal->Close();
  fileBackground->Close();
  fileOut->Close();
  delete fileOut;
  delete fileSignal;
  delete fileBackground;

  //
  // Draw the invariant mass 
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  c1->cd();
  
  setHistogramAttributes(hDielectronMass, "M_{ee} [GeV]");
  setHistogramAttributes(hDielectronMassCleaned, "M_{ee} [GeV]");

  hDielectronMass->SetLineStyle(kDashed);
  hDielectronMass->Draw("hist");
  
  hDielectronMassCleaned->Draw("same");

  TLegend *leg = new TLegend(0.12, 0.6, 0.55, 0.9);
  leg->AddEntry(hDielectronMass, "all pairs", "l");
  leg->AddEntry(hDielectronMassCleaned, "after selection", "l");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("same");

}

// Set a few histogram attributes
void setHistogramAttributes(TH1F *hist, TString XTitle){

  hist->SetLineWidth(2);
  hist->GetXaxis()->SetTitleOffset(0.8);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitle(XTitle);

  hist->SetStats(0);

  return;
}
