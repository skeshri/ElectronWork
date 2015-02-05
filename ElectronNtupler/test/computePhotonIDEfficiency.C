#include "TH2D.h"
#include "TCut.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TTree.h"
#include "TString.h"

#include <vector>

const int nWP = 3;
enum WpType { WP_LOOSE = 0,
	      WP_MEDIUM, 
	      WP_TIGHT};
const TString wpName[nWP] = 
  {"Loose", "Medium", "Tight"};

// Use here one of the WpType values
const WpType wp = WP_LOOSE;

const TString treename = "ntupler/PhotonTree";
//const TString fname1 = "/afs/cern.ch/user/i/ikrav/workspace/ntuples/GJet_Pt40_PU20bx25_photons_event_structure.root";
const TString fname1 = "photon_ntuple.root";

bool verbose = false;
bool smallEventCount = false;

const int nEtaBins = 2;



bool passWorkingPoint(WpType wp, bool isBarrel, float pt,
		      float hOverE, float full5x5_sigmaIetaIeta, 
		      float chIso, float nhIso, float phIso);

void computePhotonIDEfficiency()
{

  // This statement below should not be needed, but in one particular node I had to
  // add it, somehow the vector header was not loaded automatically there.
  gROOT->ProcessLine("#include <vector>"); 

  //
  // Find the tree
  //
  TFile *file1 = new TFile(fname1);
  if( !file1 )
    assert(0);
  TTree *tree = (TTree*)file1->Get(treename);
  if( !tree )
    assert(0);

  // Weight histograms
  TH2D *hSignal = new TH2D("hSignal","",200,-5,5,185,15,200);
  TH2D *hBackground = new TH2D("hBackground","",200,-5,5,185,15,200);
  TCut miscCut = "hasPixelSeed==0";
  TCut sigCut  = "isTrue==1";
  TCut bgCut  = "isTrue!=1";
  if( smallEventCount ){
    tree->Draw("pt:eta>>hSignal",miscCut && sigCut, "colz", 100000);
    tree->Draw("pt:eta>>hBackground",miscCut && bgCut, "colz", 100000);
  }else{
    tree->Draw("pt:eta>>hSignal",miscCut && sigCut, "colz");
    tree->Draw("pt:eta>>hBackground",miscCut && bgCut, "colz");
  }

  TH1D *ptCheckHist = new TH1D("ptCheckHist","",185, 15, 200);
  TH1D *etaCheckHist = new TH1D("etaCheckHist","",200, -5, 5);

  // Event-level variables:
  int nPho;
  // Per-photon variables
  // Kinematics
  std::vector <float> *pt = 0;     
  std::vector <float> *eta = 0;    
  // ID related variables
  std::vector <float> *full5x5_sigmaIetaIeta = 0; 
  std::vector <float> *hOverE = 0;    
  std::vector <float> *isoPhotonsWithEA = 0;    
  std::vector <float> *isoChargedHadronsWithEA = 0;    
  std::vector <float> *isoNeutralHadronsWithEA = 0;    
  std::vector <int> *isTrue = 0;    
  std::vector <int> *hasPixelSeed = 0;    
  
  // Declare branches
  TBranch *b_nPho = 0;
  TBranch *b_pt = 0;
  TBranch *b_eta = 0;
  TBranch *b_full5x5_sigmaIetaIeta = 0;
  TBranch *b_hOverE = 0;
  TBranch *b_isoPhotonsWithEA = 0;
  TBranch *b_isoChargedHadronsWithEA = 0;
  TBranch *b_isoNeutralHadronsWithEA = 0;
  TBranch *b_isTrue = 0;
  TBranch *b_hasPixelSeed = 0;

  // Connect variables and branches to the tree with the data
  tree->SetBranchAddress("nPho", &nPho, &b_nPho);
  tree->SetBranchAddress("pt", &pt, &b_pt);
  tree->SetBranchAddress("eta", &eta, &b_eta);
  tree->SetBranchAddress("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta, &b_full5x5_sigmaIetaIeta);
  tree->SetBranchAddress("hOverE", &hOverE, &b_hOverE);
  tree->SetBranchAddress("isoPhotonsWithEA", &isoPhotonsWithEA, &b_isoPhotonsWithEA);
  tree->SetBranchAddress("isoChargedHadronsWithEA", &isoChargedHadronsWithEA, &b_isoChargedHadronsWithEA);
  tree->SetBranchAddress("isoNeutralHadronsWithEA", &isoNeutralHadronsWithEA, &b_isoNeutralHadronsWithEA);
  tree->SetBranchAddress("isTrue", &isTrue, &b_isTrue);
  tree->SetBranchAddress("hasPixelSeed", &hasPixelSeed, &b_hasPixelSeed);

  //
  double sumSignalDenomEB = 0;
  double sumSignalNumEB   = 0;
  double sumSignalDenomEE = 0;
  double sumSignalNumEE   = 0;
  double sumBackDenomEB = 0;
  double sumBackNumEB   = 0;
  double sumBackDenomEE = 0;
  double sumBackNumEE   = 0;

  double sumSignalDenomEBErr2 = 0;
  double sumSignalNumEBErr2   = 0;
  double sumSignalDenomEEErr2 = 0;
  double sumSignalNumEEErr2   = 0;
  double sumBackDenomEBErr2 = 0;
  double sumBackNumEBErr2   = 0;
  double sumBackDenomEEErr2 = 0;
  double sumBackNumEEErr2   = 0;

  // 
  // Loop over events
  //
  UInt_t maxEvents = tree->GetEntries();
  if( smallEventCount )
    maxEvents = 10000;
  if(verbose)
    printf("Start loop over events, total events = %lld\n", 
           tree->GetEntries() );
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0){
      printf("."); fflush(stdout);
    }
    Long64_t tentry = tree->LoadTree(ievent);

    // Load the value of the number of the photons in the event    
    b_nPho->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of electrons %u\n", ievent, nPho);

    // Get data for all photons in this event, only vars of interest
    b_pt->GetEntry(tentry);
    b_eta->GetEntry(tentry);
    b_full5x5_sigmaIetaIeta->GetEntry(tentry);
    b_hOverE->GetEntry(tentry);
    b_isoPhotonsWithEA->GetEntry(tentry);
    b_isoChargedHadronsWithEA->GetEntry(tentry);
    b_isoNeutralHadronsWithEA->GetEntry(tentry);
    b_isTrue->GetEntry(tentry);
    b_hasPixelSeed->GetEntry(tentry);

    // Loop over photons
    for(int ipho = 0; ipho < nPho; ipho++){
      
      // Preselection
      if( !(pt->at(ipho) > 30 && pt->at(ipho) < 200 ) ) continue;
      if( !( hasPixelSeed->at(ipho) == 0 ) ) continue;

      bool isBarrel = (fabs(eta->at(ipho)) < 1.479);
      bool pass = passWorkingPoint( wp, isBarrel, pt->at(ipho),
				    hOverE->at(ipho),
				    full5x5_sigmaIetaIeta->at(ipho),
				    isoChargedHadronsWithEA->at(ipho),
				    isoNeutralHadronsWithEA->at(ipho),
				    isoPhotonsWithEA->at(ipho) );
      
      // We reweight signal and background (separately) to have 
      // a flat pt and eta distribution. This step is a matter of definition.
      double binContent = 0;
      TH2D *hEvents = hSignal;
      if( isTrue->at(ipho) != 1 )
	hEvents = hBackground;
      binContent = hEvents->GetBinContent
	( hEvents->FindBin( eta->at(ipho), pt->at(ipho) ) );
      double weight = 1;
      if( binContent == 0 ){
	printf("Error! Zero! pt=%f, eta=%f\n", pt->at(ipho), eta->at(ipho));
      }else{
	weight = 1./binContent;
      }
      
      if( isTrue->at(ipho) == 1 ) {
	ptCheckHist->Fill(pt->at(ipho), weight);
	etaCheckHist->Fill(eta->at(ipho), weight);
	  	
	if( isBarrel ) {
	  sumSignalDenomEB += weight;
	  sumSignalDenomEBErr2 += weight*weight;
	  if( pass ) {
	    sumSignalNumEB += weight;
	    sumSignalNumEBErr2 += weight*weight;
	  }
	} else {
	  sumSignalDenomEE += weight;
	  sumSignalDenomEEErr2 += weight*weight;	  
	  if( pass ) {
	    sumSignalNumEE += weight;
	    sumSignalNumEEErr2 += weight*weight;
	  }

	}// end barrel / endcap
      } // end if signal

      if( isTrue->at(ipho) !=1 ) {
	
	if( isBarrel ) {
	  
	  sumBackDenomEB += weight;
	  sumBackDenomEBErr2 += weight*weight;
	  if( pass ) {
	    sumBackNumEB += weight;
	    sumBackNumEBErr2 += weight*weight;
	  }
	} else {
	  sumBackDenomEE += weight;
	  sumBackDenomEEErr2 += weight*weight;	  
	  if( pass ) {
	    sumBackNumEE += weight;
	    sumBackNumEEErr2 += weight*weight;
	  }

	}// end barrel / endcap
      } // end if background
	
    } // end loop over photons
      
      
  }// end loop over events  

  printf("barrel signal a=%f   b=%f\n", sumSignalNumEB, sumSignalDenomEB);
  printf("barrel back a=%f    b=%f\n", sumBackNumEB, sumBackDenomEB);

  printf("\nEfficiencies for the working point %s\n", wpName[wp].Data());

  // Compute signal efficiencies
  double effSignalEB = sumSignalNumEB / sumSignalDenomEB;
  double effSignalEBErr = sqrt( sumSignalDenomEBErr2 
				* effSignalEB*(1-effSignalEB)
				/(sumSignalDenomEB*sumSignalDenomEB) );
  printf("Signal barrel efficiency: %5.1f +- %5.1f %%\n", 
	 effSignalEB*100, effSignalEBErr*100 );

  double effSignalEE = sumSignalNumEE / sumSignalDenomEE;
  double effSignalEEErr = sqrt( sumSignalDenomEEErr2 
				* effSignalEE*(1-effSignalEE)
				/(sumSignalDenomEE*sumSignalDenomEE) );
  printf("Signal endcap efficiency: %5.1f +- %5.1f %%\n", 
	 effSignalEE*100, effSignalEEErr*100 );

  // Compute background efficiencies
  double effBackEB = sumBackNumEB / sumBackDenomEB;
  double effBackEBErr = sqrt( sumBackDenomEBErr2 
				* effBackEB*(1-effBackEB)
				/(sumBackDenomEB*sumBackDenomEB) );
  printf("Background barrel efficiency: %5.1f +- %5.1f %%\n", 
	 effBackEB*100, effBackEBErr*100 );

  double effBackEE = sumBackNumEE / sumBackDenomEE;
  double effBackEEErr = sqrt( sumBackDenomEEErr2 
				* effBackEE*(1-effBackEE)
				/(sumBackDenomEE*sumBackDenomEE) );
  printf("Background endcap efficiency: %5.1f +- %5.1f %%\n", 
	 effBackEE*100, effBackEEErr*100 );
  
  TCanvas *c2 = new TCanvas("c2","c2",10,10,800, 600);
  c2->Divide(2,1);
  c2->cd(1);
  ptCheckHist->Draw();
  c2->cd(2);
  etaCheckHist->Draw();
  
}

const float hOverECut[2][nWP] = 
  { { 0.032, 0.020, 0.012 },
    { 0.023, 0.011, 0.011 } };

const float sieieCut[2][nWP] = 
  { {0.0100, 0.0099, 0.0098},
    {0.0270, 0.0269, 0.0264} };

const float chIsoCut[2][nWP] = 
  { {2.94, 2.62, 1.91},
    {3.07, 1.40, 1.26} };

const float nhIso_A[2][nWP] = 
  { {3.16, 2.69, 2.55},
    {17.16, 4.92, 2.71} };

const float nhIso_B[2][nWP] = 
  { {0.0023, 0.0023, 0.0023},
    {0.0116, 0.0116, 0.0116} };

const float phIso_A[2][nWP] = 
  { {4.43, 1.35, 1.29},
    {2.11, 2.11, 1.91} };

const float phIso_B[2][nWP] = 
  { {0.0004, 0.0004, 0.0004},
    {0.0037, 0.0037, 0.0037} };

bool passWorkingPoint(WpType iwp, bool isBarrel, float pt,
		      float hOverE, float full5x5_sigmaIetaIeta, 
		      float chIso, float nhIso, float phIso){

  int ieta = 0;
  if( !isBarrel )
    ieta = 1;

  bool result = 1
    && hOverE < hOverECut[ieta][iwp]
    && full5x5_sigmaIetaIeta < sieieCut[ieta][iwp]
    && chIso < chIsoCut[ieta][iwp]
    && nhIso < nhIso_A[ieta][iwp] + pt * nhIso_B[ieta][iwp]
    && phIso < phIso_A[ieta][iwp] + pt * phIso_B[ieta][iwp] ;
    ;
  
  return result;
}


