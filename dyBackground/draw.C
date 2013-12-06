#include "TH1F.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TArrow.h"
#include "TStyle.h"
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include "TH2F.h"
#include "TF1.h"


void setStyle(TH1F *& hist, int rebin, bool scale, Color_t color, Style_t marker, TString x_title, TString y_title);

void draw(int mH = 0, int njet=0, float lumi = 1.092, TString dysample="dy" ) {

  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  gStyle->SetOptFit(1);

  TFile *file = TFile::Open(Form("Routin_mH%i_%.0fpb.root", mH, lumi ), "READ");
  assert(file);
  gROOT->cd();

  TString x_title = "min(proj_{tk-MET}, proj_{pfMET}) (GeV)";
  if ( mH > 0 && mH <= 140 ) x_title = "DY MVA output";
  if ( njet == 2 ) 
    x_title = "pMet (GeV)";
  
  TH1F *Ree_vs_met_mc = (TH1F*) file->Get(Form("Ree_vs_met_mc_%iJet", njet));
  TH1F *Rmm_vs_met_mc = (TH1F*) file->Get(Form("Rmm_vs_met_mc_%iJet", njet));
  TH1F *R_vs_met_mc = (TH1F*) file->Get(Form("R_vs_met_mc_%iJet", njet));
  
  
  setStyle(Ree_vs_met_mc, 1, false, kRed, 21, x_title, "R_{out/in}");
  setStyle(Rmm_vs_met_mc, 1, false, kRed, 21, x_title, "R_{out/in}");
  setStyle(R_vs_met_mc, 1, false, kRed, 24, x_title, "R_{out/in}");

  Ree_vs_met_mc->SetMarkerSize(1.5);
  Rmm_vs_met_mc->SetMarkerSize(1.5);
  R_vs_met_mc->SetMarkerSize(1.5);
  
  TH1F *Ree_vs_met_data = (TH1F*) file->Get(Form("Ree_vs_met_data_%iJet", njet));
  TH1F *Rmm_vs_met_data = (TH1F*) file->Get(Form("Rmm_vs_met_data_%iJet", njet));
  TH1F *R_vs_met_data = (TH1F*) file->Get(Form("R_vs_met_data_%iJet", njet));
  
  Ree_vs_met_data->SetBinContent(Ree_vs_met_data->GetNbinsX(), 0);
  Ree_vs_met_data->SetBinError(Ree_vs_met_data->GetNbinsX(), 0);

  Rmm_vs_met_data->SetBinContent(Ree_vs_met_data->GetNbinsX(), 0);
  Rmm_vs_met_data->SetBinError(Ree_vs_met_data->GetNbinsX(), 0);

  R_vs_met_data->SetBinContent(Ree_vs_met_data->GetNbinsX(), 0);
  R_vs_met_data->SetBinError(Ree_vs_met_data->GetNbinsX(), 0);

  setStyle(Ree_vs_met_data, 1, false, kBlack, 20, x_title, "R_{out/in}");
  setStyle(Rmm_vs_met_data, 1, false, kBlack, 20, x_title, "R_{out/in}");
  setStyle(R_vs_met_data, 1, false, kBlack, 20, x_title, "R_{out/in}");

  Ree_vs_met_data->SetMarkerSize(1.5);
  Rmm_vs_met_data->SetMarkerSize(1.5);
  R_vs_met_data->SetMarkerSize(1.5);

  float yMax = 1.0;
  if (njet == 1) {
    if (mH <= 400) yMax = 1.0;
    else yMax = 5.0;
  }
  
  if (njet == 0) {
    switch (mH) {
    case 115:
    case 120:
    case 125:
    case 130:
    case 135:
    case 140:
    case 145:
      yMax = 3.0;
      break;
    case 150:
    case 160:
    case 170:
    case 180:
      yMax = 5.0;
      break;
    case 190:
    case 200:
    case 250:
    case 300:
      yMax = 1.0;
      break;
    case 350:
    case 400:
    case 450:
    case 500:
    case 550:
    case 600:
      yMax = 5.0;
      break;
    default:
      yMax = 1.0;
      break;
    }
  }
  
  
  Ree_vs_met_mc->GetYaxis()->SetRangeUser(0,yMax);
  Rmm_vs_met_mc->GetYaxis()->SetRangeUser(0,yMax);
  R_vs_met_mc->GetYaxis()->SetRangeUser(0,yMax);

  Rmm_vs_met_data->GetYaxis()->SetRangeUser(0,yMax);
  Ree_vs_met_data->GetYaxis()->SetRangeUser(0,yMax);
  R_vs_met_data->GetYaxis()->SetRangeUser(0,yMax);


  TCanvas *c1 = new TCanvas();
  c1->cd();  
  
  TLegend *leg;
  if (mH == 0)
    leg = new TLegend(0.2, 0.75, 0.6, 0.90, Form("R_{out/in} with WW Selection",mH), "brNDC");
  else
    leg = new TLegend(0.2, 0.75, 0.6, 0.90, Form("R_{out/in} with mH(%i) Selection",mH), "brNDC");
  
  leg->SetFillColor(0);
  leg->SetTextSize(0.035);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetShadowColor(0);

  // do basic fit to a straight line
  //R_vs_met_data->Fit("pol0");
  //R_vs_met_data->SetStats(1);
  R_vs_met_data->Draw("HISTE1FUNC");
  R_vs_met_mc->Draw("SAMEHISTE1");
  
  TLine *metsig;
  
  if ( mH == 0 || mH > 140 ) 
    metsig = new TLine(45, 0., 45, yMax * 0.7 );
  else {
    if ( njet == 0 ) 
      metsig = new TLine(0.6, 0., 0.6, yMax * 0.7 );
    if ( njet == 1 ) 
      metsig = new TLine(0.3, 0., 0.3, yMax * 0.7 );
  }
    
  metsig->SetLineStyle(kDashed);
  metsig->SetLineWidth(2);
  metsig->Draw("same");

  leg->AddEntry(R_vs_met_data, "OF-Subt in data", "lp");
  if (dysample == "dysoup") 
    leg->AddEntry(R_vs_met_mc, "Z/ZZ/WZ MC Prediction", "lp");
  if (dysample == "dy")
    leg->AddEntry(R_vs_met_mc, "Drell-Yan MC Prediction", "lp");
  if (dysample == "vz")
    leg->AddEntry(R_vs_met_mc, "ZZ/WZ MC Prediction", "lp");

  leg->Draw();
  c1->Print(Form("Routin_%iJet_mH%i_%.0fpb_%s.eps", njet, mH,lumi,dysample.Data()));
  c1->Print(Form("Routin_%iJet_mH%i_%.0fpb_%s.png", njet, mH,lumi,dysample.Data()));
  
  // draw ee
  c1->Clear();
  Ree_vs_met_data->Draw("HISTE1FUNC");
  Ree_vs_met_mc->Draw("SAMEHISTE1");
  leg->Draw();
  c1->Print(Form("Routin_ee_%iJet_mH%i_%.0fpb_%s.eps", njet, mH,lumi,dysample.Data()));
  c1->Print(Form("Routin_ee_%iJet_mH%i_%.0fpb_%s.png", njet, mH,lumi,dysample.Data()));
  
  // draw mm
  c1->Clear();
  Rmm_vs_met_data->Draw("HISTE1FUNC");
  Rmm_vs_met_mc->Draw("SAMEHISTE1");
  leg->Draw();
  c1->Print(Form("Routin_mm_%iJet_mH%i_%.0fpb_%s.eps", njet, mH,lumi,dysample.Data()));
  c1->Print(Form("Routin_mm_%iJet_mH%i_%.0fpb_%s.png", njet, mH,lumi,dysample.Data()));


  delete c1;
  delete Ree_vs_met_data;
  delete Rmm_vs_met_data;
  delete R_vs_met_data;
  delete Ree_vs_met_mc;
  delete Rmm_vs_met_mc;
  delete R_vs_met_mc;
  delete metsig;
  
  file->Close();
  
}

void setStyle(TH1F *& hist, int rebin, bool scale, Color_t color, Style_t marker, TString x_title, TString y_title)
{
  hist->Rebin(rebin);
  if(scale)
    hist->Scale(1.0/rebin);
  hist->SetLineColor(color);
  hist->SetLineWidth(3);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(1.0);

  hist->GetXaxis()->SetTitle(x_title);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.04);

  hist->GetYaxis()->SetTitle(y_title);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.04);

  hist->SetTitleOffset(1.5,"Y");
  hist->SetTitleOffset(1.1,"X");
  hist->SetTitleFont(42);

  hist->SetStats(0);
  hist->SetFillColor(0);
  //hist->SetBinContent(1, hist->GetBinContent(0)+hist->GetBinContent(1));
  // hist->GetXaxis()->SetNdivisions(505);
  // hist->GetYaxis()->SetNdivisions(505);
  
}
