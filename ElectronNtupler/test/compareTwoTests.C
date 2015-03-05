#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TTree.h"

const TString treename = "ntupler/PhotonTree";

const TString fname1 = "photon_ntuple_mva.root";
const TString fname2 = "photon_ntuple_mva_mini.root";

const Bool_t isLog = kTRUE;

void compareTwoTests(TString varname, float xlow, float xhigh, TString extraCuts)
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

  TFile *file2 = new TFile(fname2);
  if( !file2 )
    assert(0);
  TTree *tree2 = (TTree*)file2->Get(treename);
  if( !tree2 )
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

  tree1->Draw(drawCommand1,extraCuts);
  tree2->Draw(drawCommand2,extraCuts);

  //
  // Draw everything
  //

  hist1->SetLineWidth(2);
  hist1->SetMarkerStyle(20);
  hist1->SetMarkerSize(1);
  hist1->GetXaxis()->SetTitle(varname);

  hist2->SetLineWidth(2);
  hist2->SetLineColor(kBlue);

  hist1->SetStats(0);
  hist2->SetStats(0);

  hist1->Draw("PE");
  hist2->Draw("hist,same");

  TLegend *leg = new TLegend(0.3,0.8,0.91,0.92);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hist1, fname1, "p");
  leg->AddEntry(hist2, fname2, "l");
  leg->Draw("same");
}

