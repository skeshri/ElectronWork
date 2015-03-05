#include "TVector3.h"
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

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "TMath.h"

#include <vector>

enum Mode {EB_pt40to60 = 0,
	   EB_pt60plus,
	   EE_pt40to60,
	   EE_pt60plus};

const TString modeLabel[4] = {
  "Barrel pt 40-60",
  "Barrel pt 60+",
  "Endcap pt 40-60",
  "Endcap pt 60+"
};

const float mvaCutValue[4] = {
  0.882,
  0.923,
  0.872,
  0.909
};

const float ptMinCut[4] = {40, 60, 40, 60};
const float ptMaxCut[4] = {60, 400, 60, 400};
const float etaMinCut[4] = {-0.1, -0.1, 1.566, 1.566};
const float etaMaxCut[4] = {1.4442, 1.4442, 2.5, 2.5};

const TString treename = "ggNtuplizer/EventTree";
// const TString fname1 = "/afs/cern.ch/user/i/ikrav/workspace/ntuples/PHYS14/ggNtuple_GJets_HT-100to200_singleFile.root";
const TString fname1 = "/afs/cern.ch/work/r/rslu/public/job_phys14_gjet_pt40_20bx25.root";

bool smallEventCount = false;

bool isMatchedRS(float pEt, float pEta, float pPhi,
		 std::vector<int> *mcPID,
		 std::vector<int> *mcMomPID,
		 std::vector<float> *mcPt,
		 std::vector<float> *mcEta,
		 std::vector<float> *mcPhi);

Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) ;

Double_t deltaPhi(Double_t phi1, Double_t phi2) ;

void ggComputePhotonIDEfficiencyMVA(Mode mode){

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

  // Event-level variables:
  int nPho;
  float rho;
  // Per-photon variables
  // Kinematics
  std::vector <float> *phoEt = 0;     
  std::vector <float> *phoSCEta = 0;    
  std::vector <float> *phoPhi = 0;    
  // ID related variables
  std::vector <float> *phoSigmaIEtaIEta_2012 = 0; 
  std::vector <float> *phoHoverE = 0;    
  std::vector <float> *phoPFPhoIso = 0;    
  std::vector <float> *phoPFChIso = 0;    
  std::vector <float> *phoPFNeuIso = 0;    
  std::vector <int> *phohasPixelSeed = 0;    
  std::vector <int> *phoEleVeto = 0;    
  // Other vars
  std::vector <float>* phoR9  = 0;
  std::vector <float>* phoSigmaIEtaIPhi_2012 = 0;
  std::vector <float>* phoE1x3_2012  = 0;
  std::vector <float>* phoE2x2_2012  = 0;
  std::vector <float>* phoE2x5Max_2012  = 0;
  std::vector <float>* phoE5x5_2012  = 0;
  std::vector <float>* phoEtaWidth  = 0;
  std::vector <float>* phoPhiWidth  = 0;
  std::vector <float>* phoRawEnergy  = 0;
  std::vector <float>* phoWorstPFChIso  = 0;
  // For EE
  std::vector<Float_t> *esEffSigmaRR = 0;
  std::vector<Float_t> *esEnergy = 0;
  
  // MC variables  
  std::vector <int> *mcPID = 0;     
  std::vector <int> *mcMomPID = 0;     
  std::vector <int> *mcStatus = 0;     
  std::vector <float> *mcPt  = 0;     
  std::vector <float> *mcEta = 0;     
  std::vector <float> *mcPhi = 0;     

  // Declare branches
  TBranch *b_nPho = 0;
  TBranch *b_rho = 0;
  TBranch *b_phoEt = 0;
  TBranch *b_phoSCEta = 0;
  TBranch *b_phoPhi = 0;
  TBranch *b_phoSigmaIEtaIEta_2012 = 0;
  TBranch *b_phoHoverE = 0;
  TBranch *b_phoPFPhoIso = 0;
  TBranch *b_phoPFChIso = 0;
  TBranch *b_phoPFNeuIso = 0;
  TBranch *b_phohasPixelSeed = 0;
  TBranch *b_phoEleVeto = 0;
  //
  TBranch *b_phoR9  = 0;
  TBranch *b_phoSigmaIEtaIPhi_2012 = 0;
  TBranch *b_phoE1x3_2012  = 0;
  TBranch *b_phoE2x2_2012  = 0;
  TBranch *b_phoE2x5Max_2012  = 0;
  TBranch *b_phoE5x5_2012  = 0;
  TBranch *b_phoEtaWidth  = 0;
  TBranch *b_phoPhiWidth  = 0;
  TBranch *b_phoRawEnergy  = 0;
  TBranch *b_phoWorstPFChIso  = 0;
  // for EE
  TBranch *b_esEffSigmaRR = 0;
  TBranch *b_esEnergy = 0;
  //
  TBranch *b_mcPID = 0;
  TBranch *b_mcMomPID = 0;
  TBranch *b_mcStatus = 0;
  TBranch *b_mcPt = 0;
  TBranch *b_mcEta = 0;
  TBranch *b_mcPhi = 0;
  
  // Connect variables and branches to the tree with the data
  tree->SetBranchAddress("nPho", &nPho, &b_nPho);
  tree->SetBranchAddress("rho", &rho, &b_rho);
  tree->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
  tree->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
  tree->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
  tree->SetBranchAddress("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta_2012, &b_phoSigmaIEtaIEta_2012);
  tree->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
  tree->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
  tree->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
  tree->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
  tree->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
  tree->SetBranchAddress("phoEleVeto"     , &phoEleVeto, &b_phoEleVeto);
  //
  tree->SetBranchAddress("phoR9"          ,&phoR9   	           ,&b_phoR9  );
  tree->SetBranchAddress("phoSigmaIEtaIPhi_2012",&phoSigmaIEtaIPhi_2012      ,&b_phoSigmaIEtaIPhi_2012 );
  tree->SetBranchAddress("phoE1x3_2012"   ,&phoE1x3_2012   	  ,&b_phoE1x3_2012  );
  tree->SetBranchAddress("phoE2x2_2012"   ,&phoE2x2_2012   	  ,&b_phoE2x2_2012  );
  tree->SetBranchAddress("phoE2x5Max_2012",&phoE2x5Max_2012   	  ,&b_phoE2x5Max_2012  );
  tree->SetBranchAddress("phoE5x5_2012"   ,&phoE5x5_2012    	  ,&b_phoE5x5_2012  );
  tree->SetBranchAddress("phoSCEtaWidth"  ,&phoEtaWidth   	  ,&b_phoEtaWidth  );
  tree->SetBranchAddress("phoSCPhiWidth"  ,&phoPhiWidth   	  ,&b_phoPhiWidth  );
  tree->SetBranchAddress("phoSCRawE"      ,&phoRawEnergy   	  ,&b_phoRawEnergy  );
  tree->SetBranchAddress("phoPFChWorstIso",&phoWorstPFChIso       ,&b_phoWorstPFChIso  );
  // for ee
  tree->SetBranchAddress("phoESEffSigmaRR" ,&esEffSigmaRR  ,&b_esEffSigmaRR  );
  tree->SetBranchAddress("phoESEn"         ,&esEnergy      ,&b_esEnergy  );
  //
  tree->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
  tree->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
  tree->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
  tree->SetBranchAddress("mcPt", &mcPt , &b_mcPt);
  tree->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
  tree->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);

  //
  // Set up MVA 
  //
  // Variables that will be containers on which TMVA Reader works
  // The variables
  float varPhi_;
  float varR9_; 
  float varSieie_;
  float varSieip_; 
  float varE1x3overE5x5_; 
  float varE2x2overE5x5_; 
  float varE2x5overE5x5_; 
  float varSCEta_; 
  float varRawE_; 
  float varSCEtaWidth_; 
  float varSCPhiWidth_; 
  float varESEnOverRawE_; // only EE
  float varESEffSigmaRR_; // only EE
  float varRho_;
  float varPhoIsoRaw_;
  float varChIsoRaw_; 
  float varWorstChRaw_;
  // The spectators
  float varPt_; 
  float varEta_;

  TMVA::Reader *tmvaReader  = new TMVA::Reader( "!Color:!Silent:Error" );  
  tmvaReader ->SetVerbose(kFALSE);
  // Add all the vars, we take the string with variable name from the weights file (the Expression field)
  tmvaReader ->AddVariable("recoPhi"   , &varPhi_);
  tmvaReader ->AddVariable("r9"        , &varR9_);
  tmvaReader ->AddVariable("sieie_2012", &varSieie_);
  tmvaReader ->AddVariable("sieip_2012", &varSieip_);
  tmvaReader ->AddVariable("e1x3_2012/e5x5_2012"        , &varE1x3overE5x5_);
  tmvaReader ->AddVariable("e2x2_2012/e5x5_2012"        , &varE2x2overE5x5_);
  tmvaReader ->AddVariable("e2x5_2012/e5x5_2012"        , &varE2x5overE5x5_);
  tmvaReader ->AddVariable("recoSCEta" , &varSCEta_);
  tmvaReader ->AddVariable("rawE"      , &varRawE_);
  tmvaReader ->AddVariable("scEtaWidth", &varSCEtaWidth_);
  tmvaReader ->AddVariable("scPhiWidth", &varSCPhiWidth_);
  if( mode == EE_pt40to60 || mode == EE_pt60plus ){
    tmvaReader->AddVariable("esEn/rawE" , &varESEnOverRawE_);
    tmvaReader->AddVariable("esRR"      , &varESEffSigmaRR_);
  }
  tmvaReader ->AddVariable("rho"       , &varRho_);
  tmvaReader ->AddVariable("phoIsoRaw" , &varPhoIsoRaw_);
  tmvaReader ->AddVariable("chIsoRaw"  , &varChIsoRaw_);
  tmvaReader ->AddVariable("chWorstRaw", &varWorstChRaw_);
  // Add spectators
  tmvaReader ->AddSpectator("recoPt" , &varPt_);
  tmvaReader ->AddSpectator("recoEta", &varEta_);
  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
  //
  TString localFileName = "EgammaAnalysis/PhotonTools/data/PHYS14/photon_general_MVA_phys14_pu20bx25_EB_V1.weights.xml";
  if( mode == EE_pt40to60 || mode == EE_pt60plus ){
    localFileName = "EgammaAnalysis/PhotonTools/data/PHYS14/photon_general_MVA_phys14_pu20bx25_EE_V1.weights.xml";
  }  
  TString weightsFileName = TString(cmssw_base_src) + localFileName;
  TString methodName = "BDTG photons";
  tmvaReader->BookMVA(methodName, weightsFileName);

  // Book some useful histograms
  TH1D *hSignalMVA = new TH1D("hSignalMVA","",260,-1.3, 1.3); 
  TH1D *hBackgroundMVA = new TH1D("hBackgroundMVA","",260,-1.3, 1.3); 

  // Prepare for efficiency caclulation:
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
    maxEvents = 200;

  //
  // First loop over events, fill weight histograms
  //
  // Weight histograms
  TH2D *hSignal = new TH2D("hSignal","",100,-5,5,80,0,400);
  TH2D *hBackground = new TH2D("hBackground","",100,-5,5,80,0,400);
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){
    
    if( ievent%100000 == 0){
      printf("."); fflush(stdout);
    }
    Long64_t tentry = tree->LoadTree(ievent);
    // Load the value of the number of the photons in the event    
    b_nPho->GetEntry(tentry);

    b_phoEt->GetEntry(tentry);
    b_phoSCEta->GetEntry(tentry);
    b_phoPhi->GetEntry(tentry);
    b_phoEleVeto->GetEntry(tentry);

    b_mcPID->GetEntry(tentry);
    b_mcMomPID->GetEntry(tentry);
    b_mcStatus->GetEntry(tentry);
    b_mcPt->GetEntry(tentry);
    b_mcEta->GetEntry(tentry);
    b_mcPhi->GetEntry(tentry);

    // Loop over photons
    for(int ipho = 0; ipho < nPho; ipho++){
      
      // Preselection, keep only true barrel photons, with electron veto applied
      // if( !( phoEt->at(ipho) > ptMinCut[mode] && phoEt->at(ipho) < ptMaxCut[mode] ) ) continue;
      // if( !( fabs( phoSCEta->at(ipho) ) > etaMinCut[mode] &&  fabs( phoSCEta->at(ipho) ) < etaMaxCut[mode] ) ) continue;
      if( !( phoEleVeto->at(ipho) == 1 ) ) continue;

      // Match to MC truth
      bool isTrue = isMatchedRS( phoEt->at(ipho), phoSCEta->at(ipho), phoPhi->at(ipho),
				 mcPID, mcMomPID,
				 mcPt, mcEta, mcPhi);
      if( isTrue )
	hSignal->Fill(phoSCEta->at(ipho), phoEt->at(ipho));
      else
	hBackground->Fill(phoSCEta->at(ipho), phoEt->at(ipho));

    } // end loop over photons

  } // end first loop over events
  TH1D *weightHist = (TH1D*)hBackground->Clone("weightHist");
  weightHist->Divide(hSignal);
  
  //
  // Second loop over events, compute efficiencies
  //
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){
    
    if( ievent%100000 == 0){
      printf("."); fflush(stdout);
    }
    Long64_t tentry = tree->LoadTree(ievent);
    // Load the value of the number of the photons in the event    
    b_nPho->GetEntry(tentry);
    b_rho->GetEntry(tentry);

    // Get data for all photons in this event, only vars of interest
    b_phoEt->GetEntry(tentry);
    b_phoSCEta->GetEntry(tentry);
    b_phoPhi->GetEntry(tentry);
    b_phoSigmaIEtaIEta_2012->GetEntry(tentry);
    b_phoHoverE->GetEntry(tentry);
    b_phoPFPhoIso->GetEntry(tentry);
    b_phoPFChIso->GetEntry(tentry);
    b_phoPFNeuIso->GetEntry(tentry);
    b_phohasPixelSeed->GetEntry(tentry);
    b_phoEleVeto->GetEntry(tentry);

    b_phoR9 ->GetEntry(tentry);
    b_phoSigmaIEtaIPhi_2012->GetEntry(tentry);
    b_phoE1x3_2012 ->GetEntry(tentry);
    b_phoE2x2_2012 ->GetEntry(tentry);
    b_phoE2x5Max_2012 ->GetEntry(tentry);
    b_phoE5x5_2012 ->GetEntry(tentry);
    b_phoEtaWidth ->GetEntry(tentry);
    b_phoPhiWidth ->GetEntry(tentry);
    b_phoRawEnergy ->GetEntry(tentry);
    b_phoWorstPFChIso ->GetEntry(tentry);
    // For EE only:
    b_esEnergy     ->GetEntry(tentry);
    b_esEffSigmaRR ->GetEntry(tentry);

    b_mcPID->GetEntry(tentry);
    b_mcMomPID->GetEntry(tentry);
    b_mcStatus->GetEntry(tentry);
    b_mcPt->GetEntry(tentry);
    b_mcEta->GetEntry(tentry);
    b_mcPhi->GetEntry(tentry);

    // Loop over photons
    for(int ipho = 0; ipho < nPho; ipho++){
      
      // Preselection, keep only true barrel photons, with electron veto applied
      if( !( phoEt->at(ipho) > ptMinCut[mode] && phoEt->at(ipho) < ptMaxCut[mode] ) ) continue;
      if( !( fabs( phoSCEta->at(ipho) ) > etaMinCut[mode] &&  fabs( phoSCEta->at(ipho) ) < etaMaxCut[mode] ) ) continue;
      if( !( phoEleVeto->at(ipho) == 1 ) ) continue;

      // Match to MC truth
      bool isTrue = isMatchedRS( phoEt->at(ipho), phoSCEta->at(ipho), phoPhi->at(ipho),
				 mcPID, mcMomPID,
				 mcPt, mcEta, mcPhi);

      // if( !isTrue ) continue;

      varPhi_          = phoPhi->at(ipho);
      varR9_           = phoR9->at(ipho); 
      varSieie_        = phoSigmaIEtaIEta_2012->at(ipho);
      varSieip_        = phoSigmaIEtaIPhi_2012->at(ipho); 
      varE1x3overE5x5_ = phoE1x3_2012->at(ipho) / phoE5x5_2012->at(ipho); 
      varE2x2overE5x5_ = phoE2x2_2012->at(ipho) / phoE5x5_2012->at(ipho); 
      varE2x5overE5x5_ = phoE2x5Max_2012->at(ipho) / phoE5x5_2012->at(ipho); 
      varSCEta_        = phoSCEta->at(ipho); 
      varRawE_         = phoRawEnergy->at(ipho); 
      varSCEtaWidth_   = phoEtaWidth->at(ipho); 
      varSCPhiWidth_   = phoPhiWidth->at(ipho); 
      varESEnOverRawE_ = esEnergy->at(ipho) / phoRawEnergy->at(ipho);
      varESEffSigmaRR_ = esEffSigmaRR->at(ipho);
      varRho_          = rho;
      varPhoIsoRaw_    = phoPFPhoIso->at(ipho);
      varChIsoRaw_     = phoPFChIso->at(ipho); 
      varWorstChRaw_   = phoWorstPFChIso->at(ipho);
      // The spectators
      varPt_           = phoEt->at(ipho); 
      varEta_          = phoSCEta->at(ipho);

      // printf("--->\n");
      // printf("phi=%f  R9=%f  sieie=%f  sieip=%f  e1x3=%f    e2x2=%f   e2x5=%f    e5x5=%f\n",
      // 	     );   
      // printf("etaSC=%f   rawE=%f   esEnergy=%f   esEffSigmaRR=%f\n");

      float mvaValue =  tmvaReader->EvaluateMVA(methodName) ;
      bool pass = (mvaValue > mvaCutValue[mode]);

      if ( isTrue) 
	hSignalMVA->Fill(mvaValue);
      else
	hBackgroundMVA->Fill(mvaValue);

      double weight = 1;
      if( isTrue )
        weight = weightHist->GetBinContent
          ( weightHist->FindBin( phoSCEta->at(ipho), phoEt->at(ipho) ) );

      // printf("pt=%f  eta=%f  mva=%f    isTrue=%d\n", phoEt->at(ipho), phoSCEta->at(ipho), mvaValue, isTrue);

      if( isTrue ) {
                
        sumSignalDenom += weight;
        sumSignalDenomErr2 += weight*weight;
        if( pass ) {
          sumSignalNum += weight;
          sumSignalNumErr2 += weight*weight;
        }
      } // end if signal
      else {
        sumBackDenom += weight;
        sumBackDenomErr2 += weight*weight;

        if( pass ) {
          sumBackNum += weight;
          sumBackNumErr2 += weight*weight;
        }
      } // end if background


      //printf("Found good photon\n");

    }// end loop over photons
  } // end loop over events

  // 
  // Finalize efficiencies and print them
  //
  printf("signal a=%f   b=%f\n", sumSignalNum, sumSignalDenom);
  printf("backgr a=%f    b=%f\n", sumBackNum, sumBackDenom);

  printf("\nEfficiencies for the %s  and mva cut is %f\n", 
         modeLabel[mode].Data(),
	 mvaCutValue[mode] );

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


  //
  // Draw MVA values
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  c1->cd();

  hSignalMVA->SetLineColor(kBlue);
  hSignalMVA->SetLineWidth(2);
  hBackgroundMVA->SetLineColor(kRed);
  hBackgroundMVA->SetLineWidth(2);

  hBackgroundMVA->Scale(hSignalMVA->GetSumOfWeights() / hBackgroundMVA->GetSumOfWeights());

  hSignalMVA->Draw("hist");
  hBackgroundMVA->Draw("hist,same");

}

//
// MC truth matching
//


bool isMatchedRS(float pEt, float pEta, float pPhi,
               std::vector<int> *mcPID,
               std::vector<int> *mcMomPID,
               std::vector<float> *mcPt,
               std::vector<float> *mcEta,
               std::vector<float> *mcPhi){

  bool isMatched    = false;
  int nMC = mcPID->size();
  // printf("pho Et %.2f, eta %.2f, phi %.2f  \n", pEt, pEta, pPhi);
  for (Int_t k=0; k<nMC; ++k) {     
    // printf("   mc  pt %.2f, eta %.2f, phi %.2f  pid= %d   mpid= %d\n", 
    // 	   (*mcPt)[k], (*mcEta)[k], (*mcPhi)[k], (*mcPID)[k], (*mcMomPID)[k]);
    if ((*mcPID)[k] == 22 && TMath::Abs( (*mcMomPID)[k] ) <= 22 ) {
      // printf("      matched id\n");
      float dr = deltaR(pEta, pPhi, (*mcEta)[k], (*mcPhi)[k]);
      float dpt = fabs((pEt - (*mcPt)[k])/(*mcPt)[k]);
      if (dr < 0.2 && dpt < 0.2){
	isMatched = true;
	// printf("         matched kine\n");
      }
    }
  }
  return isMatched;
 }


Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);
  return sqrt(dEta*dEta+dPhi*dPhi);
}
Double_t deltaPhi(Double_t phi1, Double_t phi2) {
  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
  return dPhi;
}

