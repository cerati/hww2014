
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include <iostream>
#include "TH2F.h"
#include "TH1F.h"
#include "TString.h"
#include "TRint.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"


void drawsingle2d(TString fileName, TString histName, TString appendix); 
void drawsingle1d(TString fileName, bool normalize);
TH2F* Roll1DTo2D(TH1* h1, const char* hname);
bool signalregion = true;
  int njet = 1; 

void drawtempfromcard() {

  TString fileName = Form("../cards/hwwjcp_19p5fb_hypsep/125/hwwof_%ij.input_8TeV.root", njet);
  drawsingle2d(fileName, "histo_ggH", "hww");
  drawsingle2d(fileName, "histo_qqWW", "qqWW");
  drawsingle2d(fileName, "histo_Top", "Top");
  drawsingle2d(fileName, "histo_WjetsE", "WjetsE");
  drawsingle2d(fileName, "histo_WjetsM", "WjetsM");
  drawsingle2d(fileName, "histo_Wgamma", "Wgamma");
  drawsingle2d(fileName, "histo_Wg3l", "Wgstar");
  // drawsingle2d(fileName, "histo_VV", "VV");
  
  TString fileName = Form("../cards/hwwjcp_19p5fb_hypsep/125/xwwof_%ij.input_8TeV.root", njet);
  drawsingle2d(fileName, "histo_ggH", "xww");
  
}

void drawsingle2d(TString fileName, TString histName, TString appendix) {

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleXOffset(1.1);                                                                                   
  TGaxis::SetMaxDigits(3);
  gStyle->SetNdivisions(505, "XYZ");                                               
  gROOT->ForceStyle();



  TFile *f = TFile::Open(fileName, "READ");
  gROOT->cd();
  TH1F *h1 = (TH1F*)f->Get(histName);
  TH2F *h2 =  Roll1DTo2D(h1,""); 

  if ( signalregion ) {
    h2->GetXaxis()->SetRangeUser(60, 119);     
    h2->GetYaxis()->SetRangeUser(12, 70); 
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  if ( ! signalregion ) 
    c1->SetWindowSize(1000, 600);

  h2->SetXTitle("m_{T} [GeV]");
  h2->SetYTitle("m_{ll} [GeV]");
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();
  h2->Draw("colz");
 
  if ( signalregion ) {
    c1->SaveAs(Form("tempplots/mtvsmll_%s_%ij.eps", appendix.Data(), njet));
    c1->SaveAs(Form("tempplots/mtvsmll_%s_%ij.png", appendix.Data(), njet));
    c1->SaveAs(Form("tempplots/mtvsmll_%s_%ij.pdf", appendix.Data(), njet));
  } else {
    c1->SaveAs(Form("tempplots/mtvsmll_%s_%ij_fullrange.eps", appendix.Data(), njet));
    c1->SaveAs(Form("tempplots/mtvsmll_%s_%ij_fullrange.png", appendix.Data(), njet));
    c1->SaveAs(Form("tempplots/mtvsmll_%s_%ij_fullrange.pdf", appendix.Data(), njet));
  }
 
  delete c1;
  f->Close();

  
}


void drawsingle1d(TString fileName, bool normalize) {

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleXOffset(1.1);                                                                                   
  TGaxis::SetMaxDigits(3);
  gStyle->SetNdivisions(505, "XYZ");                                               
  gROOT->ForceStyle();



  TFile *f = TFile::Open(fileName, "READ");
  gROOT->cd();

  TH1F *h1_ggH = (TH1F*)f->Get("histo_ggH");
  h1_ggH->SetLineColor(kRed);
  h1_ggH->SetLineWidth(3);
  
  TH1F *h1_qqWW = (TH1F*)f->Get("histo_qqWW");
  h1_qqWW->SetLineColor(kBlack);
  h1_qqWW->SetLineWidth(3);

  TH1F *h1_Wjets = (TH1F*)f->Get("histo_Wjets");
  h1_Wjets->SetLineColor(kBlue);
  h1_Wjets->SetLineWidth(3);
  
  
  if ( normalize ) {
    h1_ggH->Scale(1./h1_ggH->Integral(0,1000));
    h1_qqWW->Scale(1./h1_qqWW->Integral(0,1000));
    h1_Wjets->Scale(1./h1_Wjets->Integral(0,1000));
  }

  double yMax = h1_ggH->GetMaximum();
  yMax = yMax > h1_qqWW->GetMaximum() ? yMax : h1_qqWW->GetMaximum();
  yMax = yMax > h1_Wjets->GetMaximum() ? yMax : h1_Wjets->GetMaximum();


  TLegend *leg = new TLegend(0.55, 0.7, 0.82, 0.9);
  leg->SetFillColor(0);                                                                                                              
  leg->AddEntry(h1_ggH, "HWW(125)", "l");
  leg->AddEntry(h1_qqWW, "qq#rightarrow WW", "l");
  leg->AddEntry(h1_Wjets, "Wjets", "l");

  TCanvas *c1 = new TCanvas();
  gPad->SetLogy(1);
  h1_ggH->SetXTitle("BDT output");
  h1_ggH->SetYTitle(Form("Number of Events Per 0.1"));
  h1_ggH->GetXaxis()->CenterTitle();
  h1_ggH->SetMaximum(yMax * 1.1);
  h1_ggH->Draw("hist");
  h1_qqWW->Draw("samehist");
  h1_Wjets->Draw("samehist");
  leg->Draw("same");
   
  c1->SaveAs(Form("tempplots/bdt_mH125.eps"));
  c1->SaveAs(Form("tempplots/bdt_mH125.png"));
  
  delete c1;
  f->Close();

  
}


// rolling 
TH2F* Roll1DTo2D(TH1* h1, const char* hname) {

	unsigned int nbins  = h1->GetXaxis()->GetNbins();
	float mtbins[15]	=	{60,70,80,90,100,110,120,140,160,180,200,220,240,260,280};
	float mllbins[10]	=	{12,30,45,60,75,100,125,150,175,200}; 
	TH2F* h2 = new TH2F(hname, hname, 14, mtbins, 9, mllbins);
	
	for(unsigned int ibin=1; ibin <= nbins; ++ibin) {
		int ibinX = (int) (ibin-1)/9+1;
		int ibinY = (int) (ibin-1)%9+1;
		h2->SetBinContent(  ibinX, ibinY, h1->GetBinContent(ibin) );
		h2->SetBinError(    ibinX, ibinY, h1->GetBinError(ibin) );
	}
	h2->SetStats(0);
	return h2;
}
