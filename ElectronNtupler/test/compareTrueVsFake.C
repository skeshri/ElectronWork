#include "TCut.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TTree.h"

const TString treename = "ntupler/PhotonTree";
// const TString fname1 = "photon_ntuple_wfoot_mini.root";
// const TString fname1 = "photon_ntuple_nofoot_keymatching.root";
// const TString fname1 = "photon_ntuple_wfoot_keymatching_v2.root";
// const TString fname2 = "photon_ntuple_wfoot_keymatching_v3.root";

// const TString fname1 = "photon_ntuple_nofoot_mini_v2.root";
// const TString fname2 = "photon_ntuple_wfoot_mini_v2.root";

// const TString fname1 = "photon_ntuple_nofoot_v4.root";
// const TString fname1 = "photon_ntuple_v5.root";
// const TString fname1 = "photon_ntuple.root";
// const TString fname2 = "photon_ntuple_wfoot_v4.root";

//const TString fname1 = "/afs/cern.ch/user/i/ikrav/workspace/ntuples/GJet_Pt40_PU20bx25_photons_event_structure.root";

// const TString fname1 = "photon_ntuple_v7.root";
// const TString label  = "w_footprint_removal";

const TString fname1 = "photon_ntuple_nofoot_v7.root";
const TString label  = "NO_footprint_removal";

const Bool_t isLog = kTRUE;
const Bool_t toNormalize = kTRUE;

void compareTrueVsFake(TString varname, float xlow, float xhigh, 
		       TString extraCuts1, TString extraCuts2, bool isBarrel)
{

  //
  // Find the two trees
  //
  TFile *file1 = new TFile(fname1);
  if( !file1 )
    assert(0);
  TTree *tree1 = (TTree*)file1->Get(treename);
  if( !tree1 )
    assert(0);

  // 
  // Make histograms
  //
  
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  c1->cd();
  c1->SetLogy(isLog);

  TH1F *hist1 = new TH1F("hist1","", 100, xlow, xhigh);
  TH1F *hist2 = new TH1F("hist2","", 100, xlow, xhigh);

  TString drawCommand1 = TString::Format("%s>>hist1",varname.Data());
  TString drawCommand2 = TString::Format("%s>>hist2",varname.Data());

  TString etaCutsString = " abs(eta)<1.479 ";
  TString etaLabel = "barrel";
  if( !isBarrel){
    etaCutsString = " abs(eta)>1.479 ";
    etaLabel = "endcap";
  }
  TCut etaCut = etaCutsString.Data();

  tree1->Draw(drawCommand1,extraCuts1 && etaCut);
  tree1->Draw(drawCommand2,extraCuts2 && etaCut);

  //
  // Draw everything
  //

  hist1->SetLineWidth(2);
  hist1->SetLineColor(kRed);
  hist1->SetMarkerStyle(20);
  hist1->SetMarkerSize(1);
  hist1->GetXaxis()->SetTitle(varname);

  hist2->SetLineWidth(2);
  hist2->SetLineColor(kBlue);

  hist1->SetStats(0);
  hist2->SetStats(0);

  if( toNormalize )
    hist2->Scale(hist1->GetSumOfWeights()/hist2->GetSumOfWeights());

  hist1->Draw("hist");
  hist2->Draw("hist,same");

  TLegend *leg = new TLegend(0.3,0.75,0.91,0.88);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hist1, extraCuts1, "p");
  leg->AddEntry(hist2, extraCuts2, "l");
  leg->Draw("same");

  TLatex *comment0 = new TLatex(0.6, 0.7, etaLabel.Data());
  comment0->SetNDC(kTRUE);
  comment0->Draw();

  // TString label = TString::Format("source: %s\n", fname1.Data());
  TLatex *comment = new TLatex(0.1, 0.95, label.Data());
  comment->SetNDC(kTRUE);
  comment->Draw();

  TString varname_adjusted = varname;
  varname_adjusted.ReplaceAll("/","_relativeTo_");
  TString outFileName = TString::Format("figures/pho_%s_%s_%s_true_vs_fake.pdf",
					varname_adjusted.Data(), label.Data(), etaLabel.Data()); 

  c1->Print(outFileName);
}

