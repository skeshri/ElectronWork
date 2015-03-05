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

const TString treename = "ntupler/PhotonTree";
//const TString fname1 = "/afs/cern.ch/user/i/ikrav/workspace/ntuples/PHYS14/GJet_Pt40_PU20bx25_photons_event_structure_wMVA.root";
const TString fname1 = "photon_ntuple_mva.root";

bool verbose = false;
bool smallEventCount = false;

const int nPtBins  = 2;
const int nEtaBins = 2;

const float mvaCutValue[nEtaBins][nPtBins] = 
  { {0.882, 0.923}, {0.872, 0.909} };

TString binLabel[nEtaBins][nPtBins];

void computePhotonIDEfficiencyMVA(int etaBin, int ptBin)
{

  binLabel[0][0] = TString("pt40to60 barrel");
  binLabel[0][1] = TString("pt60to400 barrel");
  binLabel[1][0] = TString("pt40to60 endcap");
  binLabel[1][1] = TString("pt60to400 endcap");

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
  TH2D *hSignal = new TH2D("hSignal","",100,-5,5,80,0,400);
  TH2D *hBackground = new TH2D("hBackground","",100,-5,5,80,0,400);
  TCut miscCut = "passElectronVeto==1";
  TCut sigCut  = "isTrueAlternative==1";
  TCut bgCut  = "isTrueAlternative!=1";
  if( smallEventCount ){
    tree->Draw("pt:eta>>hSignal",miscCut && sigCut, "colz", 100000);
    tree->Draw("pt:eta>>hBackground",miscCut && bgCut, "colz", 100000);
  }else{
    tree->Draw("pt:eta>>hSignal",miscCut && sigCut, "colz");
    tree->Draw("pt:eta>>hBackground",miscCut && bgCut, "colz");
  }
  if(verbose)
   printf("Finished with weight histograms\n");

  TH1D *ptCheckHistSig = new TH1D("ptCheckHistSig","",80, 0, 400);
  TH1D *etaCheckHistSig = new TH1D("etaCheckHistSig","",100, -5, 5);
  TH1D *ptCheckHistBg = new TH1D("ptCheckHistBg","",80, 0, 400);
  TH1D *etaCheckHistBg = new TH1D("etaCheckHistBg","",100, -5, 5);

  TH1D *weightHist = (TH1D*)hBackground->Clone("weightHist");
  weightHist->Divide(hSignal);

  // Event-level variables:
  int nPho;
  // Per-photon variables
  // Kinematics
  std::vector <float> *pt = 0;     
  std::vector <float> *eta = 0;    
  // ID related variables
  std::vector <int> *isTrueAlternative = 0;    
  std::vector <int> *passElectronVeto = 0;    
  std::vector <float> *mvaValue = 0;    
  
  // Declare branches
  TBranch *b_nPho = 0;
  TBranch *b_pt = 0;
  TBranch *b_eta = 0;
  TBranch *b_isTrueAlternative = 0;
  TBranch *b_passElectronVeto = 0;
  TBranch *b_mvaValue = 0;

  // Connect variables and branches to the tree with the data
  tree->SetBranchAddress("nPho", &nPho, &b_nPho);
  tree->SetBranchAddress("pt", &pt, &b_pt);
  tree->SetBranchAddress("eta", &eta, &b_eta);
  tree->SetBranchAddress("isTrueAlternative", &isTrueAlternative, &b_isTrueAlternative);
  tree->SetBranchAddress("passElectronVeto", &passElectronVeto, &b_passElectronVeto);
  tree->SetBranchAddress("mvaValue", &mvaValue, &b_mvaValue);

  //
  double sumSignalDenom = 0;
  double sumSignalNum   = 0;
  double sumBackDenom = 0;
  double sumBackNum   = 0;

  double sumSignalDenomErr2 = 0;
  double sumSignalNumErr2   = 0;
  double sumBackDenomErr2 = 0;
  double sumBackNumErr2   = 0;

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
      printf("Event %d, number of photons %u\n", ievent, nPho);

    // Get data for all photons in this event, only vars of interest
    b_pt->GetEntry(tentry);
    b_eta->GetEntry(tentry);
    b_isTrueAlternative->GetEntry(tentry);
    b_passElectronVeto->GetEntry(tentry);
    b_mvaValue->GetEntry(tentry);

    // Loop over photons
    for(int ipho = 0; ipho < nPho; ipho++){

      // Preselection
      if( !( pt->at(ipho) > 40 && pt->at(ipho) < 400 ) ) continue;
      if( !( passElectronVeto->at(ipho) == 1 ) ) continue;
      if( fabs(eta->at(ipho))>1.4442 && fabs(eta->at(ipho))<1.566 ) continue;
      if( fabs(eta->at(ipho))>2.5 ) continue;

      // Keep only electrons of the desired pt/eta range
      bool isBarrel = ( fabs(eta->at(ipho)) < 1.479 );
      bool is40to60 = ( pt->at(ipho) >40 && pt->at(ipho) < 60 );
      
      int thisPtBin = 0;
      if( !is40to60 )
	thisPtBin = 1;

      int thisEtaBin = 0;
      if( !isBarrel )
	thisEtaBin = 1;

      if( ! ( ptBin == thisPtBin && etaBin == thisEtaBin) )
	continue;

      // We reweight signal and background (separately) to have 
      // a flat pt and eta distribution. This step is a matter of definition.
      // double binContent = 0;
      // TH2D *hEvents = hSignal;
      // if( isTrue->at(ipho) != 1 )
      // 	hEvents = hBackground;
      // binContent = hEvents->GetBinContent
      // 	( hEvents->FindBin( eta->at(ipho), pt->at(ipho) ) );
      // double weight = 1;
      // if( binContent == 0 ){
      // 	printf("Error! Zero! pt=%f, eta=%f\n", pt->at(ipho), eta->at(ipho));
      // }else{
      // 	weight = 1./binContent;
      // }

      // Reweight signal to fakes
      double weight = 1;
      if( isTrueAlternative->at(ipho) == 1 )
      	weight = weightHist->GetBinContent
      	  ( weightHist->FindBin( eta->at(ipho), pt->at(ipho) ) );

      // Check if mva cut is passed
      bool pass = ( mvaValue->at(ipho) > mvaCutValue[etaBin][ptBin] );
      
      if( isTrueAlternative->at(ipho) == 1 ) {
	ptCheckHistSig->Fill(pt->at(ipho), weight);
	etaCheckHistSig->Fill(eta->at(ipho), weight);
	  	
	sumSignalDenom += weight;
	sumSignalDenomErr2 += weight*weight;
	if( pass ) {
	  sumSignalNum += weight;
	  sumSignalNumErr2 += weight*weight;
	}
      } // end if signal

      if( isTrueAlternative->at(ipho) !=1 ) {
	
	sumBackDenom += weight;
	sumBackDenomErr2 += weight*weight;

	ptCheckHistBg->Fill(pt->at(ipho), weight);
	etaCheckHistBg->Fill(eta->at(ipho), weight);
	if( pass ) {
	  sumBackNum += weight;
	  sumBackNumErr2 += weight*weight;
	}
      } // end if background
	
    } // end loop over photons
      
      
  }// end loop over events  
  printf("\n");

  if(verbose){
    printf("signal a=%f   b=%f\n", sumSignalNum, sumSignalDenom);
    printf("backgr a=%f    b=%f\n", sumBackNum, sumBackDenom);
  }

  printf("\nEfficiencies for the bin %s  and mva cut is %f\n", 
	 (binLabel[etaBin][ptBin]).Data(),
	  mvaCutValue[etaBin][ptBin] );

  // Compute signal efficiencies
  double effSignal = sumSignalNum / sumSignalDenom;
  double effSignalErr = sqrt( sumSignalDenomErr2 
			      * effSignal*(1-effSignal)
			      /(sumSignalDenom*sumSignalDenom) );
  printf("Signal efficiency: %5.1f +- %5.1f %%\n", 
	 effSignal*100, effSignalErr*100 );


  // Compute background efficiencies
  double effBack = sumBackNum / sumBackDenom;
  double effBackErr = sqrt( sumBackDenomErr2 
				* effBack*(1-effBack)
				/(sumBackDenom*sumBackDenom) );
  printf("Background efficiency: %5.1f +- %5.1f %%\n", 
	 effBack*100, effBackErr*100 );

  
  TCanvas *c2 = new TCanvas("c2","c2",10,10,800, 600);
  c2->Divide(2,1);
  c2->cd(1);
  ptCheckHistSig->Draw("pe");
  ptCheckHistBg->Draw("hist,same");
  c2->cd(2);
  etaCheckHistSig->Draw("pe");
  etaCheckHistBg->Draw("hist,same");
  
}

