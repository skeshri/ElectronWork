#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TTree.h"
#include "TString.h"
#include "TLatex.h"

#include <vector>

const TString fileName = "/home/hep/ikrav/work/ntuples/CSA14/DYJetsToLL_PU20bx25_event_structure_v3.root";
const TString treeName = "ntupler/ElectronTree";

const bool verbose = false;
const bool smallEventCount = false;

const float boundaryBarrelEndcap = 1.479;

enum TruthMatching {ELECTRON_FAKE, ELECTRON_TRUE_PROMPT, 
		    ELECTRON_TRUE_FROM_TAU, ELECTRON_TRUE_SECONDARY}; // secondary does not include from tau

// Forward declarations
void setHistogramAttributes(TH1F *hist, TString XTitle);

//
// Main program
//
void exampleElectronQuantities(){

  // This statement below should not be needed, but in one particular node I had to
  // add it, somehow the vector header was not loaded automatically there.
  gROOT->ProcessLine("#include <vector>"); 

  // Book histograms for interesting variables
  TH1F *hPt = new TH1F("hPt","",100, 0, 200);
  TH1F *hEta = new TH1F("hEta","",100,-3, 3);
  TH1F *h5x5See_barrel = new TH1F("h5x5See_barrel", "", 100, 0, 0.05); // Full 5x5 sigma_ieta_ieta
  TH1F *h5x5See_endcap = new TH1F("h5x5See_endcap", "", 100, 0, 0.05); // Full 5x5 sigma_ieta_ieta

  //
  // Open a file and find the tree with electron data
  //
  TFile *file = new TFile(fileName);
  if( !file ){
    printf("Failed to open file %s\n", fileName.Data() );
    assert(0);
  }
  TTree *tree = (TTree*)file->Get(treeName);
  if( !tree ){
    printf("Failed to find the tree %s\n", treeName.Data() );
    assert(0);
  }

  // 
  // Configure reading the information of interest from the TTree,
  // in this example we are working with events containing vectors.
  // A good ROOT reference: 
  //      http://root.cern.ch/root/html/tutorials/tree/hvector.C.html
  //

  // 
  // Set up the branches of interest
  //
  // Declare variables (uncomment as needed! )
  //
  // Event-level variables:
  // int nPV;                    // number of reconstructed PVs
  int nEle; // the number of reconstructed electrons in the event
  // Per-eletron variables
  // Kinematics
  std::vector <float> *elePt = 0;         // electron PT
  std::vector <float> *eleEtaSC = 0;      // supercluser eta
  // Matching track-supercluster
  // std::vector <float> *eleDEtaIn = 0;  // deltaEtaIn
  // std::vector <float> *eleDPhiIn = 0;  // deltaPhiIn
  // Misc ID variables
  // std::vector <float> *eleHoverE = 0;  // H/E  
  std::vector <float> *eleFull5x5SigmaIEtaIEta = 0;  // sigma_ieta_ieta from full 5x5 calculation
  // std::vector <float> *eleOOEMOOP = 0; // |1/E - 1/p|
  // Impact parameters
  // std::vector <float> *eleD0 = 0;      // r-phi plane impact parameter
  // std::vector <float> *eleDZ = 0;      // r-z plane impact parameter
  // PF isolation variables, in the cone dR<0.3, electron footprint removed
  // std::vector <float> *eleIsoChargedHadrons = 0;  // total charged hadron energy in the cone
  // std::vector <float> *eleIsoNeutralHadrons = 0;  // total neutral hadron energy in the cone
  // std::vector <float> *eleIsoPhotons = 0;         // total photon energy in the cone
  // std::vector <float> *eleIsoChargedFromPU = 0;   // total charged not from the main PV
  // std::vector <float> *eleRelIsoWithDBeta = 0;    // pre-computed (ch+max(0.0,nh+pho-0.5*chPU)/pt from above
  // Conversion rejection
  // std::vector <float> *eleExpectedMissingInnerHits = 0;
  // std::vector <float> *elePassConversionVeto = 0;
  // Matching to gen truth
  std::vector <int> *eleIsTrueElectron = 0;    // 0-unmatched, 1-true prompt, 2-true from tau, 3-true non-prompt
  
  // Declare branches
  TBranch *b_nEle = 0;
  TBranch *b_elePt = 0;
  TBranch *b_eleEtaSC = 0;
  TBranch *b_eleFull5x5SigmaIEtaIEta = 0;
  TBranch *b_eleIsTrueElectron = 0;
  //
  // .... add here more branch declarations that mirror the variables above.
  //

  // Connect variables and branches to the tree with the data
  tree->SetBranchAddress("nEle", &nEle, &b_nEle);
  tree->SetBranchAddress("pt", &elePt, &b_elePt);
  tree->SetBranchAddress("etaSC", &eleEtaSC, &b_eleEtaSC);
  tree->SetBranchAddress("full5x5_sigmaIetaIeta", &eleFull5x5SigmaIEtaIEta, &b_eleFull5x5SigmaIEtaIEta);
  tree->SetBranchAddress("isTrueElectron", &eleIsTrueElectron, &b_eleIsTrueElectron);
  //
  // ... add here setting more branch addresses if needed
  //

  // 
  // Loop over events
  //
  UInt_t maxEvents = tree->GetEntries();
  if( smallEventCount )
    maxEvents = 1000;
  if(verbose)
    printf("Start loop over events, total events = %lld\n", tree->GetEntries() );
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    Long64_t tentry = tree->LoadTree(ievent);

    // Load the value of the number of the electrons in the event    
    b_nEle->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of electrons %u\n", ievent, nEle);

    // If there are no electrons, skip to the next event
    if( nEle == 0 )
      continue;
    
    // Get data for all electrons in this event, only vars of interest
    b_elePt->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_eleFull5x5SigmaIEtaIEta->GetEntry(tentry);
    b_eleIsTrueElectron->GetEntry(tentry);
    //
    // ...add more get entries for other branches as needed ...
    //

    //
    // Loop over electrons
    //
    for(int iele = 0; iele < nEle; iele++){
      
      // In this example, work only with reconstructed electrons matched
      // to true prompt electrons in MC
      if( eleIsTrueElectron->at(iele) != ELECTRON_TRUE_PROMPT )
	continue;
      
      // Fill histograms
      hPt->Fill( elePt->at(iele) );
      hEta->Fill( eleEtaSC->at(iele) );

      bool isBarrel = fabs( eleEtaSC->at(iele) ) < boundaryBarrelEndcap ? true : false; 
      if( isBarrel ) {
      	h5x5See_barrel->Fill( eleFull5x5SigmaIEtaIEta->at(iele) );
      }else{
      	h5x5See_endcap->Fill( eleFull5x5SigmaIEtaIEta->at(iele) );
      }

    } // end loop over electrons
  
  } // end loop over events
					       
  // Final clean-up, as recommended by ROOT tutorials
  tree->ResetBranchAddresses();

  file->Close();
  delete file;
						
  //
  // Draw the histograms
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  c1->Divide(2,2);
  
  c1->cd(1);
  setHistogramAttributes(hPt, "p_{T} [GeV]");
  hPt->Draw("hist");

  c1->cd(2);
  setHistogramAttributes(hEta, "#eta");
  hEta->Draw("hist");

  c1->cd(3);
  setHistogramAttributes(h5x5See_barrel, "full 5x5 #sigma_{i#etai#eta}");
  h5x5See_barrel->Draw("hist");
  TLatex *comment3 = new TLatex(0.6, 0.6, "BARREL");
  comment3->SetNDC(kTRUE);
  comment3->Draw("same");

  c1->cd(4);
  setHistogramAttributes(h5x5See_endcap, "full 5x5 #sigma_{i#etai#eta}");
  h5x5See_endcap->Draw("hist");
  TLatex *comment4 = new TLatex(0.6, 0.6, "ENDCAP");
  comment4->SetNDC(kTRUE);
  comment4->Draw("same");


}

// Set a few histogram attributes
void setHistogramAttributes(TH1F *hist, TString XTitle){

  hist->SetLineWidth(2);
  hist->GetXaxis()->SetTitleOffset(0.8);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitle(XTitle);

  // hist->GetXaxis()->SetLabelSize(0.05);
  // hist->GetYaxis()->SetLabelSize(0.05);

  hist->SetStats(0);

  return;
}

