#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include <iostream>
#include "TH2F.h"
#include "TH1F.h"
#include "TString.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TRint.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TCut.h"
#include "THStack.h"

void drawOverlay(TString fileName, TString mvatype, int njet, int mH, TString plotDir, TString flavor);
 
//###################
//# main function
//###################
void mvaoverlay(TString inputFDir = "../cards/LP2011_ME/", int mH = 115, int njet = 0,  TString flavor = "sf", TString plotDir = "plots/", TString mvatype="ME") {
  
  TString fileName =  inputFDir+Form("/%i/hww%s_%ij.input_8TeV.root", mH, flavor.Data(), njet);
  if ( mvatype == "ME") 
    fileName =  inputFDir+Form("/%i/hww%s_%ij_me.input_8TeV.root", mH, flavor.Data(), njet);

  drawOverlay(fileName, mvatype, njet, mH, plotDir, flavor);
}  

void drawOverlay(TString fileName, TString mvatype, int njet, int mH, TString plotDir, TString flavor) {
  // make it look nice
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  gROOT->ForceStyle();
  TFile *mvaFile = TFile::Open(fileName, "READ");
  if ( mvaFile == 0x0 ) {
    std::cout << "Error, " << fileName << " does not exist, exitting...\n";
    return;
  }
  std::cout << "Opening..." << fileName << "\n";
  gROOT->cd();
  
  TH1F *hist_data = (TH1F*) mvaFile->Get("histo_Data");
  hist_data->SetMarkerStyle(20);
  
  // get the signal shape by adding 4 sub-channels
  TH1F *hist_ggH = (TH1F*) mvaFile->Get("histo_ggH");
  TH1F *hist_qqH = (TH1F*) mvaFile->Get("histo_qqH");
  TH1F *hist_WH = (TH1F*) mvaFile->Get("histo_WH");
  TH1F *hist_ZH = (TH1F*) mvaFile->Get("histo_ZH");
  TH1F *hist_sig = (TH1F*) mvaFile->Get("histo_ggH")->Clone("histo_Higgs");
  hist_sig->Add(hist_sig, hist_qqH);
  hist_sig->Add(hist_sig, hist_WH);
  hist_sig->Add(hist_sig, hist_ZH);
  hist_sig->SetLineColor(kRed+1);
  hist_sig->SetMarkerColor(kRed+1);
  hist_sig->SetLineWidth(3);

  // get the backgrounds
  std::vector<TH1F*> hist_bkg;
  THStack *stack_bkg = new THStack(); 

  // define the stack and overlay legends..
  TLegend *stacklg = new TLegend(0.0, 0.4, 1.0, 0.95);
  stacklg->SetBorderSize(0);
  stacklg->SetFillStyle(0);
  stacklg->SetShadowColor(0);
  stacklg->SetTextSize(0.16);

  TLegend *overlaylg = new TLegend(0.0, 0.4, 1.0, 0.95);
  overlaylg->SetBorderSize(0);
  overlaylg->SetFillStyle(0);
  overlaylg->SetShadowColor(0);
  overlaylg->SetTextSize(0.16);
  
  // qqWW
  TH1F *hist_qqWW = (TH1F*) mvaFile->Get("histo_qqWW");
  hist_qqWW->SetLineColor(kAzure-9);
  hist_qqWW->SetMarkerColor(kAzure-9);
  hist_qqWW->SetFillColor(kAzure-9);
  hist_bkg.push_back(hist_qqWW);
  stack_bkg->Add(hist_qqWW);
  stacklg->AddEntry(hist_qqWW, "qqWW", "f");
  overlaylg->AddEntry(hist_qqWW, "qqWW", "l");

  // ggWW
  TH1F *hist_ggWW = (TH1F*) mvaFile->Get("histo_ggWW");
  hist_ggWW->SetLineColor(kAzure-9);
  hist_ggWW->SetMarkerColor(kAzure-9);
  hist_ggWW->SetFillColor(kAzure-9);
  hist_bkg.push_back(hist_ggWW);
  stack_bkg->Add(hist_ggWW);
  stacklg->AddEntry(hist_ggWW, "ggWW", "f");
  overlaylg->AddEntry(hist_ggWW, "ggWW", "l");
  
  // Wjets
  TH1F *hist_Wjets = (TH1F*) mvaFile->Get("histo_Wjets");
  hist_Wjets->SetLineColor(kGray+1);
  hist_Wjets->SetMarkerColor(kGray+1);
  hist_Wjets->SetFillColor(kGray+1);
  hist_bkg.push_back(hist_Wjets);
  stack_bkg->Add(hist_Wjets);
  stacklg->AddEntry(hist_Wjets, "Wjets", "f");
  overlaylg->AddEntry(hist_Wjets, "Wjets", "l");

  // Top
  TH1F *hist_Top = (TH1F*) mvaFile->Get("histo_Top");
  hist_Top->SetLineColor(kYellow);
  hist_Top->SetMarkerColor(kYellow);
  hist_Top->SetFillColor(kYellow);
  hist_bkg.push_back(hist_Top);
  stack_bkg->Add(hist_Top);
  stacklg->AddEntry(hist_Top, "Top", "f");
  overlaylg->AddEntry(hist_Top, "Top", "l");

  // VV
  TH1F *hist_VV = (TH1F*) mvaFile->Get("histo_VV");
  hist_VV->SetLineColor(kAzure-2);
  hist_VV->SetMarkerColor(kAzure-2);
  hist_VV->SetFillColor(kAzure-2);
  hist_bkg.push_back(hist_VV);
  stack_bkg->Add(hist_VV);
  stacklg->AddEntry(hist_VV, "VV", "f");
  overlaylg->AddEntry(hist_VV, "VV", "l");

  // Zjets
  TH1F *hist_Zjets = (TH1F*) mvaFile->Get("histo_Zjets");
  hist_Zjets->SetLineColor(kGreen+2);
  hist_Zjets->SetMarkerColor(kGreen+2);
  hist_Zjets->SetFillColor(kGreen+2);
  hist_bkg.push_back(hist_Zjets);
  stack_bkg->Add(hist_Zjets);
  stacklg->AddEntry(hist_Zjets, "Zjets", "f");
  overlaylg->AddEntry(hist_Zjets, "Zjets", "l");

  // Wgamma
  TH1F *hist_Wgamma = (TH1F*) mvaFile->Get("histo_Wgamma");
  hist_Wgamma->SetLineColor(kYellow+2);
  hist_Wgamma->SetMarkerColor(kYellow+2);
  hist_Wgamma->SetFillColor(kYellow+2);
  hist_bkg.push_back(hist_Wgamma);
  stack_bkg->Add(hist_Wgamma);
  stacklg->AddEntry(hist_Wgamma, "Wgamma", "f");
  overlaylg->AddEntry(hist_Wgamma, "Wgamma", "l");
  
  // Ztt
  TH1F *hist_Ztt = (TH1F*) mvaFile->Get("histo_Ztt");
  hist_Ztt->SetLineColor(kGreen+2);
  hist_Ztt->SetMarkerColor(kGreen+2);
  hist_Ztt->SetFillColor(kGreen+2);
  hist_bkg.push_back(hist_Ztt);
  stack_bkg->Add(hist_Ztt);
  stacklg->AddEntry(hist_Ztt, "Ztt", "f");
  overlaylg->AddEntry(hist_Ztt, "Ztt", "l");


  stacklg->AddEntry(hist_sig, Form("HWW%i", mH), "l");
  overlaylg->AddEntry(hist_sig, Form("HWW%i", mH), "l");

  stacklg->AddEntry(hist_data, "Data", "lp");
  //overlaylg->AddEntry(hist_data, "Data", "lp");
  
  
  float yMax = hist_sig->GetMaximum();
  for (unsigned int i = 0; i < hist_bkg.size(); i++) {
    hist_bkg[i]->SetLineWidth(3);
    hist_bkg[i]->SetFillStyle(1001);
    hist_bkg[i]->Draw("SAMEHIST");
    yMax = yMax > hist_bkg[i]->GetMaximum() ? yMax : hist_bkg[i]->GetMaximum();
  }
  hist_sig->SetMaximum(yMax  + 1.0 *sqrt(yMax));
  yMax = yMax > hist_data->GetMaximum()  ? yMax : hist_data->GetMaximum();
  stack_bkg->SetMaximum(yMax  + 1.5 *sqrt(yMax));
  hist_sig->SetMinimum(0.1);
  stack_bkg->SetMinimum(0.1);


  // == Draw histograms
  // draw stacked plots
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
  c1->cd();
  
  TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.0, 0.82, 1.0);
  pad1->SetBottomMargin(0.13);
  pad1->SetRightMargin(0.07);
  pad1->Draw();
  c1->cd();
  
  TPad *pad2 = new TPad("p_leg", "p_leg", 0.82, 0.0, 1.0, 1.0);
  pad2->SetTopMargin(0.01);
  pad2->SetRightMargin(0.01);
  pad2->SetBottomMargin(0.13);
  pad2->Draw();
  
  pad1->cd();
  pad1->SetLogy(0);
  stack_bkg->Draw("HIST");
  stack_bkg->GetXaxis()->SetTitle(mvatype);
  stack_bkg->GetYaxis()->SetTitle(Form("Number of Events / %.2f", stack_bkg->GetXaxis()->GetBinWidth(1)));
  hist_sig->Draw("SAMEHIST");
  hist_data->Draw("SAMEE1");
  pad2->cd();
  stacklg->Draw();
  
  // draw stacked plots - linear
  c1->SaveAs( plotDir + Form("/%s_mH%i_%ij_%s_stack_lin.png",mvatype.Data(), mH, njet, flavor.Data()));
  c1->SaveAs( plotDir + Form("/%s_mH%i_%ij_%s_stack_lin.eps",mvatype.Data(), mH, njet, flavor.Data()));

  // draw stacked plots - log
  pad1->SetLogy();
  c1->SaveAs( plotDir + Form("/%s_mH%i_%ij_%s_stack_log.png",mvatype.Data(), mH, njet, flavor.Data()));
  c1->SaveAs( plotDir + Form("/%s_mH%i_%ij_%s_stack_log.eps",mvatype.Data(), mH, njet, flavor.Data()));

  // === draw overlay plots
  c1->cd();
  pad1->Clear();
  pad1->cd();
  pad1->SetLogy(0);
  hist_sig->Draw("HIST"); 
  for ( int s = 0; s < hist_bkg.size(); s++) {
    hist_bkg[s]->SetFillColor(0);
    hist_bkg[s]->SetLineWidth(2);
    hist_bkg[s]->Draw("SAMEHIST");
  }
  pad2->Clear();
  pad2->cd();
  overlaylg->Draw();

  // draw overlay plots - linear

  c1->SaveAs( plotDir + Form("/%s_mH%i_%ij_%s_overlay_lin.png",mvatype.Data(), mH, njet, flavor.Data()));
  c1->SaveAs( plotDir + Form("/%s_mH%i_%ij_%s_overlay_lin.eps",mvatype.Data(), mH, njet, flavor.Data()));
  
  // draw overlay plots - log
  pad1->SetLogy(1);
  c1->SaveAs( plotDir + Form("/%s_mH%i_%ij_%s_overlay_log.png",mvatype.Data(), mH, njet, flavor.Data()));
  c1->SaveAs( plotDir + Form("/%s_mH%i_%ij_%s_overlay_log.eps",mvatype.Data(), mH, njet, flavor.Data()));
  
  delete pad1;
  delete pad2;
  // end of drawing histograms
  
  
  // tidy up
  delete hist_sig; 
  delete hist_data;
  for (unsigned int s = 0; s < hist_bkg.size(); s++) delete hist_bkg[s];
  delete stack_bkg;
  delete stacklg;
  delete overlaylg;
  delete c1;
  
  // tidy up
  mvaFile->Close();
}

