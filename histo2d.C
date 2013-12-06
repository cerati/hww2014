#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include <iostream>
#include "TH2F.h"
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
#include "TPaletteAxis.h"
#include "TStyle.h"
#include "THStack.h"

//#define ERRORMIN 0.001
#define ERRORMIN 0.0

void draw(TString fileName, TString mvatype, int njet, int mH, TString plotDir, TString flavor);

// make uncertainty histograms
TH2F* filluncer(TH2F* h2, TH2F* h2_total) {
	
	TH2F* output = (TH2F*) h2->Clone();
	
	for (unsigned int x = 1; x <= h2->GetXaxis()->GetNbins(); ++x) {
		for (unsigned int y = 1; y <= h2->GetYaxis()->GetNbins(); ++y) {
			output->SetBinContent(x, y, h2->GetBinError(x, y)/h2_total->GetBinContent(x, y));
		}
	}
	return output;
}

float getHistMinumum(TH2F* h2) {
	
	float min = 999999.;

	for (unsigned int x = 1; x <= h2->GetXaxis()->GetNbins(); ++x) {
		for (unsigned int y = 1; y <= h2->GetYaxis()->GetNbins(); ++y) {
			if( h2->GetBinContent(x,y) > 0 && h2->GetBinContent(x,y) < min)	 min = h2->GetBinContent(x,y);
		}
	}
	
	return min;
}

//###################
//# main function
//###################
void histo2d(TString inputFDir = "../cards/2D", int mH = 125, int njet = 0,  TString flavor = "sf", TString plotDir = "plots/2D", TString mvatype="2D") {
 
  TString fileName = "";
  fileName =  inputFDir+Form("/%i/hww%s_%ij_2d.input_8TeV.root", mH, flavor.Data(), njet);
  //fileName =  inputFDir+Form("/%i/hww%s_%ij_2d.input_8TeV_7by10bin.root", mH, flavor.Data(), njet);

  draw(fileName, mvatype, njet, mH, plotDir, flavor);
}  

void draw(TString fileName, TString mvatype, int njet, int mH, TString plotDir, TString flavor) {
  // make it look nice
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  gROOT->ForceStyle(); 
  gStyle->SetPaintTextFormat("6.2f");
  gStyle->SetMarkerSize(1.5);
  TFile *File = TFile::Open(fileName, "READ");
  if ( File == 0x0 ) {
    std::cout << "Error, " << fileName << " does not exist, exitting...\n";
    return;
  }
  std::cout << "Opening..." << fileName << "\n";
  gROOT->cd();
 

  TString xvar = "M_{T}";
  TString yvar = "M_{ll}";


  // ------------------------ 
  // 	Total background MC 
  // ------------------------ 

  // First get a histogram with all backgrounds added   
  TH2F *hist_bkg_qqWW 	= (TH2F*) File->Get("histo2_qqWW");
  TH2F *hist_bkg_ggWW 	= (TH2F*) File->Get("histo2_ggWW");
  TH2F *hist_bkg_Wjets 	= (TH2F*) File->Get("histo2_Wjets");
  TH2F *hist_bkg_Top 	= (TH2F*) File->Get("histo2_Top");
  TH2F *hist_bkg_VV 	= (TH2F*) File->Get("histo2_VV");
  TH2F *hist_bkg_Zjets 	= (TH2F*) File->Get("histo2_Zjets");
  TH2F *hist_bkg_Wgamma = (TH2F*) File->Get("histo2_Wgamma");
  TH2F *hist_bkg_Ztt 	= (TH2F*) File->Get("histo2_Ztt");
  
  TH2F *hist_bkg	= (TH2F*) hist_bkg_qqWW->Clone("hist_bkg");
  hist_bkg->Add(hist_bkg_ggWW);
  hist_bkg->Add(hist_bkg_Wjets);
  hist_bkg->Add(hist_bkg_Top);
  hist_bkg->Add(hist_bkg_VV);
  hist_bkg->Add(hist_bkg_Zjets);
  hist_bkg->Add(hist_bkg_Wgamma);
  hist_bkg->Add(hist_bkg_Ztt);

  float histmax 	= hist_bkg->GetMaximum();
  
   
  float histmin 	= 999.;
  if( getHistMinumum(hist_bkg_qqWW)<histmin ) 	histmin = getHistMinumum(hist_bkg_qqWW);
  if( getHistMinumum(hist_bkg_ggWW)<histmin ) 	histmin = getHistMinumum(hist_bkg_ggWW);
  if( getHistMinumum(hist_bkg_Wjets)<histmin ) 	histmin = getHistMinumum(hist_bkg_Wjets);
  if( getHistMinumum(hist_bkg_Top)<histmin ) 	histmin = getHistMinumum(hist_bkg_Top);
  if( getHistMinumum(hist_bkg_VV)<histmin ) 	histmin = getHistMinumum(hist_bkg_VV);
  if( getHistMinumum(hist_bkg_Zjets)<histmin ) 	histmin = getHistMinumum(hist_bkg_Zjets);
  if( getHistMinumum(hist_bkg_Wgamma)<histmin ) 	histmin = getHistMinumum(hist_bkg_Wgamma);
  if( getHistMinumum(hist_bkg_Ztt)<histmin ) 	histmin = getHistMinumum(hist_bkg_Ztt);

  histmin = 0.02;

  // ------------------------ 
  // 	DATA 
  // ------------------------ 

  TCanvas *c_data = new TCanvas("c_data", "c_data", 600, 500);
  c_data->cd(1);
  //c_data->SetLogz(1);
  c_data->SetRightMargin(0.15);
  TH2F *hist_data = (TH2F*) File->Get("histo2_Data");
  hist_data->SetTitle("DATA");
  hist_data->SetXTitle(xvar);
  hist_data->SetYTitle(yvar);
  hist_data->SetMaximum(histmax);
  hist_data->SetMinimum(histmin*0.5); 
  hist_data->Draw("colz");
  c_data->SaveAs( plotDir + Form("/data_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_data->SaveAs( plotDir + Form("/data_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));
  // uncertainty histogram  
  TH2F* hist_data_error = filluncer(hist_data, hist_data);
  hist_data_error->Draw("colz");
  hist_data_error->SetTitle("SIGNAL");
  hist_data_error->SetXTitle(xvar);
  hist_data_error->SetYTitle(yvar);
  c_data->SaveAs( plotDir + Form("/dataerr_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_data->SaveAs( plotDir + Form("/dataerr_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));


  // ------------------------ 
  // 	MC 
  // ------------------------ 

  // 
  // get the signal shape by adding 4 sub-channels
  // 
  TCanvas *c_sig = new TCanvas("c_sig", "c_sig", 600, 500);
  c_sig->cd(1);
  //c_sig->SetLogz(1);
  c_sig->SetRightMargin(0.15);
  TH2F *hist_ggH = (TH2F*) File->Get("histo2_ggH");
  TH2F *hist_qqH = (TH2F*) File->Get("histo2_qqH");
  TH2F *hist_WH = (TH2F*) File->Get("histo2_WH");
  TH2F *hist_ZH = (TH2F*) File->Get("histo2_ZH");
  TH2F *hist_sig = (TH2F*) File->Get("histo2_ggH")->Clone("histo2_Higgs");
  hist_sig->Add(hist_sig, hist_qqH);
  hist_sig->Add(hist_sig, hist_WH);
  hist_sig->Add(hist_sig, hist_ZH);
  hist_sig->SetTitle("SIGNAL");
  hist_sig->SetXTitle(xvar);
  hist_sig->SetYTitle(yvar);
//  hist_sig->SetMinimum(histmin*0.5); 
//  hist_sig->SetMaximum(histmax);
  hist_sig->Draw("colz");
  c_sig->SaveAs( plotDir + Form("/sig_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_sig->SaveAs( plotDir + Form("/sig_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));
  // uncertainty histogram  
  TH2F* hist_sig_error = filluncer(hist_sig, hist_bkg);
  hist_sig_error->Draw("colz");
  hist_sig_error->SetTitle("SIGNAL");
  hist_sig_error->SetXTitle(xvar);
  hist_sig_error->SetYTitle(yvar);
  hist_sig_error->SetMinimum(ERRORMIN);	
  //hist_sig_error->SetMinimum(getHistMinumum(hist_sig_error));	
  hist_sig_error->SetMaximum(1);
  c_sig->SaveAs( plotDir + Form("/sigerr_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_sig->SaveAs( plotDir + Form("/sigerr_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));

  // S / B plots
  TCanvas *c_SoverB = new TCanvas("c_SoverB", "c_SoverB", 600, 500);
  c_SoverB->cd(1);
  //c_SoverB->SetLogz(1);
  c_SoverB->SetRightMargin(0.15);
  TH2F *hist_SoverB = (TH2F*)hist_sig->Clone("histo2_SoverB");
  hist_SoverB->Divide(hist_bkg);
  hist_SoverB->SetTitle("SoverB");
  hist_SoverB->SetXTitle(xvar);
  hist_SoverB->SetYTitle(yvar);
//  hist_SoverB->SetMinimum(0.0001); 
//  hist_SoverB->SetMaximum(0.35);
  hist_SoverB->Draw("colz");
  c_SoverB->SaveAs( plotDir + Form("/SoverB_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_SoverB->SaveAs( plotDir + Form("/SoverB_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));



  // qqWW
  TCanvas *c_qqWW = new TCanvas("c_qqWW", "c_qqWW", 600, 500);
  c_qqWW->cd(1);
  //c_qqWW->SetLogz(1); 
  c_qqWW->SetRightMargin(0.15);
  TH2F *hist_qqWW = (TH2F*) File->Get("histo2_qqWW");
  hist_qqWW->SetTitle("QQWW");
  hist_qqWW->SetXTitle(xvar);
  hist_qqWW->SetYTitle(yvar);
//  hist_qqWW->SetMinimum(histmin*0.5); 
//  hist_qqWW->SetMaximum(histmax); 
  hist_qqWW->Draw("colz");
  c_qqWW->SaveAs( plotDir + Form("/qqWW_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_qqWW->SaveAs( plotDir + Form("/qqWW_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));
  // uncertainty histogram  
  TH2F* hist_qqWW_error = filluncer(hist_qqWW, hist_bkg);
  hist_qqWW_error->Draw("colz");
  hist_qqWW_error->SetTitle("SIGNAL");
  hist_qqWW_error->SetXTitle(xvar);
  hist_qqWW_error->SetYTitle(yvar);
  hist_qqWW_error->SetMinimum(ERRORMIN);	
  //hist_qqWW_error->SetMinimum(getHistMinumum(hist_qqWW_error));	
  hist_qqWW_error->SetMaximum(1);
  c_qqWW->SaveAs( plotDir + Form("/qqWWerr_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_qqWW->SaveAs( plotDir + Form("/qqWWerr_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));


  // qqWW : mc@nlo
  TCanvas *c_qqWWmcnlo = new TCanvas("c_qqWWmcnlo", "c_qqWWmcnlo", 600, 500);
  c_qqWWmcnlo->cd(1);
  //c_qqWWmcnlo->SetLogz(1); 
  c_qqWWmcnlo->SetRightMargin(0.15);
  TH2F *hist_qqWWmcnlo = (TH2F*) File->Get("histo2_qqWW");
  hist_qqWWmcnlo->SetTitle("QQWW");
  hist_qqWWmcnlo->SetXTitle(xvar);
  hist_qqWWmcnlo->SetYTitle(yvar);
//  hist_qqWW->SetMinimum(histmin*0.5); 
//  hist_qqWW->SetMaximum(histmax); 
  hist_qqWW->Draw("colz");
  c_qqWW->SaveAs( plotDir + Form("/qqWW_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_qqWW->SaveAs( plotDir + Form("/qqWW_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));



  // ggWW
  TCanvas *c_ggWW = new TCanvas("c_ggWW", "c_ggWW", 600, 500);
  c_ggWW->cd(1);
  //c_ggWW->SetLogz(1);
  c_ggWW->SetRightMargin(0.15);
  TH2F *hist_ggWW = (TH2F*) File->Get("histo2_ggWW");
  hist_ggWW->SetTitle("GGWW");
  hist_ggWW->SetXTitle(xvar);
  hist_ggWW->SetYTitle(yvar);
//  hist_ggWW->SetMinimum(histmin*0.5); 
//  hist_ggWW->SetMaximum(histmax);
  hist_ggWW->Draw("colz");
  c_ggWW->SaveAs( plotDir + Form("/ggWW_%s_mH%i_%ij_%s.eps",mvatype.Data(), mH, njet, flavor.Data()));
  c_ggWW->SaveAs( plotDir + Form("/ggWW_%s_mH%i_%ij_%s.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  // uncertainty histogram  
  TH2F* hist_ggWW_error = filluncer(hist_ggWW, hist_bkg);
  hist_ggWW_error->Draw("colz");
  hist_ggWW_error->SetTitle("ggWWNAL");
  hist_ggWW_error->SetXTitle(xvar);
  hist_ggWW_error->SetYTitle(yvar);
  hist_ggWW_error->SetMinimum(ERRORMIN);	
  hist_ggWW_error->SetMaximum(1);
  c_ggWW->SaveAs( plotDir + Form("/ggWWerr_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_ggWW->SaveAs( plotDir + Form("/ggWWerr_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));

  
  // Wjets
  TCanvas *c_Wjets = new TCanvas("c_Wjets", "c_Wjets", 600, 500);
  c_Wjets->cd(1);
  //c_Wjets->SetLogz(1);
  c_Wjets->SetRightMargin(0.15);
  TH2F *hist_Wjets = (TH2F*) File->Get("histo2_Wjets");
  hist_Wjets->SetTitle("Wjets");
  hist_Wjets->SetXTitle(xvar);
  hist_Wjets->SetYTitle(yvar);
//  hist_Wjets->SetMinimum(histmin*0.5); 
//  hist_Wjets->SetMaximum(histmax);
  hist_Wjets->Draw("colz");
  c_Wjets->SaveAs( plotDir + Form("/Wjets_%s_mH%i_%ij_%s.eps",mvatype.Data(), mH, njet, flavor.Data()));
  c_Wjets->SaveAs( plotDir + Form("/Wjets_%s_mH%i_%ij_%s.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  // uncertainty histogram  
  TH2F* hist_Wjets_error = filluncer(hist_Wjets, hist_bkg);
  hist_Wjets_error->Draw("colz");
  hist_Wjets_error->SetTitle("WjetsNAL");
  hist_Wjets_error->SetXTitle(xvar);
  hist_Wjets_error->SetYTitle(yvar);
  hist_Wjets_error->SetMinimum(ERRORMIN);	
  hist_Wjets_error->SetMaximum(1);
  c_Wjets->SaveAs( plotDir + Form("/Wjetserr_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_Wjets->SaveAs( plotDir + Form("/Wjetserr_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));


  // Top
  TCanvas *c_Top = new TCanvas("c_Top", "c_Top", 600, 500);
  c_Top->cd(1);
  //c_Top->SetLogz(1);
  c_Top->SetRightMargin(0.15);
  TH2F *hist_Top = (TH2F*) File->Get("histo2_Top");
  hist_Top->SetTitle("TOP");
  hist_Top->SetXTitle(xvar);
  hist_Top->SetYTitle(yvar);
//  hist_Top->SetMinimum(histmin*0.5); 
//  hist_Top->SetMaximum(histmax);
  hist_Top->Draw("colz");
  c_Top->SaveAs( plotDir + Form("/Top_%s_mH%i_%ij_%s.eps",mvatype.Data(), mH, njet, flavor.Data()));
  c_Top->SaveAs( plotDir + Form("/Top_%s_mH%i_%ij_%s.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  // uncertainty histogram  
  TH2F* hist_Top_error = filluncer(hist_Top, hist_bkg);
  hist_Top_error->Draw("colz");
  hist_Top_error->SetTitle("TopNAL");
  hist_Top_error->SetXTitle(xvar);
  hist_Top_error->SetYTitle(yvar);
  hist_Top_error->SetMinimum(ERRORMIN);	
  hist_Top_error->SetMaximum(1);
  c_Top->SaveAs( plotDir + Form("/Toperr_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_Top->SaveAs( plotDir + Form("/Toperr_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));


  // VV
  TCanvas *c_VV = new TCanvas("c_VV", "c_VV", 600, 500);
  c_VV->cd(1);
  //c_VV->SetLogz(1);
  c_VV->SetRightMargin(0.15);
  TH2F *hist_VV = (TH2F*) File->Get("histo2_VV");
  hist_VV->SetTitle("VV");
  hist_VV->SetXTitle(xvar);
  hist_VV->SetYTitle(yvar);
//  hist_VV->SetMinimum(histmin*0.5); 
//  hist_VV->SetMaximum(histmax);
  hist_VV->Draw("colz");
  c_VV->SaveAs( plotDir + Form("/VV_%s_mH%i_%ij_%s.eps",mvatype.Data(), mH, njet, flavor.Data()));
  c_VV->SaveAs( plotDir + Form("/VV_%s_mH%i_%ij_%s.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  // uncertainty histogram  
  TH2F* hist_VV_error = filluncer(hist_VV, hist_bkg);
  hist_VV_error->Draw("colz");
  hist_VV_error->SetTitle("VVNAL");
  hist_VV_error->SetXTitle(xvar);
  hist_VV_error->SetYTitle(yvar);
  hist_VV_error->SetMinimum(ERRORMIN);	
  hist_VV_error->SetMaximum(1);
  c_VV->SaveAs( plotDir + Form("/VVerr_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_VV->SaveAs( plotDir + Form("/VVerr_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));


  // Zjets
  TCanvas *c_Zjets = new TCanvas("c_Zjets", "c_Zjets", 600, 500);
  c_Zjets->cd(1);
  //c_Zjets->SetLogz(1);
  c_Zjets->SetRightMargin(0.15);
  TH2F *hist_Zjets = (TH2F*) File->Get("histo2_Zjets");
  hist_Zjets->SetTitle("ZJETS");
  hist_Zjets->SetXTitle(xvar);
  hist_Zjets->SetYTitle(yvar);
//  hist_Zjets->SetMinimum(histmin*0.5); 
//  hist_Zjets->SetMaximum(histmax);
  hist_Zjets->Draw("colz");
  c_Zjets->SaveAs( plotDir + Form("/Zjets_%s_mH%i_%ij_%s.eps",mvatype.Data(), mH, njet, flavor.Data()));
  c_Zjets->SaveAs( plotDir + Form("/Zjets_%s_mH%i_%ij_%s.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  // uncertainty histogram  
  TH2F* hist_Zjets_error = filluncer(hist_Zjets, hist_bkg);
  hist_Zjets_error->Draw("colz");
  hist_Zjets_error->SetTitle("ZjetsNAL");
  hist_Zjets_error->SetXTitle(xvar);
  hist_Zjets_error->SetYTitle(yvar);
  hist_Zjets_error->SetMinimum(ERRORMIN);	
  hist_Zjets_error->SetMaximum(1);
  c_Zjets->SaveAs( plotDir + Form("/Zjetserr_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_Zjets->SaveAs( plotDir + Form("/Zjetserr_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));


  // Wgamma
  TCanvas *c_Wgamma = new TCanvas("c_Wgamma", "c_Wgamma", 600, 500);
  c_Wgamma->cd(1);
  //c_Wgamma->SetLogz(1);
  c_Wgamma->SetRightMargin(0.15);
  TH2F *hist_Wgamma = (TH2F*) File->Get("histo2_Wgamma");
  hist_Wgamma->SetTitle("WGAMMA");
  hist_Wgamma->SetXTitle(xvar);
  hist_Wgamma->SetYTitle(yvar);
//  hist_Wgamma->SetMinimum(histmin*0.5); 
//  hist_Wgamma->SetMaximum(histmax);
  hist_Wgamma->Draw("colz");
  c_Wgamma->SaveAs( plotDir + Form("/Wgamma_%s_mH%i_%ij_%s.eps",mvatype.Data(), mH, njet, flavor.Data()));
  c_Wgamma->SaveAs( plotDir + Form("/Wgamma_%s_mH%i_%ij_%s.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  // uncertainty histogram  
  TH2F* hist_Wgamma_error = filluncer(hist_Wgamma, hist_bkg);
  hist_Wgamma_error->Draw("colz");
  hist_Wgamma_error->SetTitle("WgammaNAL");
  hist_Wgamma_error->SetXTitle(xvar);
  hist_Wgamma_error->SetYTitle(yvar);
  hist_Wgamma_error->SetMinimum(ERRORMIN);	
  hist_Wgamma_error->SetMaximum(1);
  c_Wgamma->SaveAs( plotDir + Form("/Wgammaerr_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_Wgamma->SaveAs( plotDir + Form("/Wgammaerr_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));

  
  // Ztt
  TCanvas *c_Ztt = new TCanvas("c_Ztt", "c_Ztt", 600, 500);
  c_Ztt->cd(1);
  //c_Ztt->SetLogz(1);
  c_Ztt->SetRightMargin(0.15);
  TH2F *hist_Ztt = (TH2F*) File->Get("histo2_Ztt");
  hist_Ztt->SetTitle("Ztt");
  hist_Ztt->SetXTitle(xvar);
  hist_Ztt->SetYTitle(yvar);
//  hist_Ztt->SetMinimum(histmin*0.5); 
//  hist_Ztt->SetMaximum(histmax);
  hist_Ztt->Draw("colz");
  c_Ztt->SaveAs( plotDir + Form("/Ztt_%s_mH%i_%ij_%s.eps",mvatype.Data(), mH, njet, flavor.Data()));
  c_Ztt->SaveAs( plotDir + Form("/Ztt_%s_mH%i_%ij_%s.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  // uncertainty histogram  
  TH2F* hist_Ztt_error = filluncer(hist_Ztt, hist_bkg);
  hist_Ztt_error->Draw("colz");
  hist_Ztt_error->SetTitle("ZttNAL");
  hist_Ztt_error->SetXTitle(xvar);
  hist_Ztt_error->SetYTitle(yvar);
  hist_Ztt_error->SetMinimum(ERRORMIN);	
  hist_Ztt_error->SetMaximum(1);
  c_Ztt->SaveAs( plotDir + Form("/Ztterr_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_Ztt->SaveAs( plotDir + Form("/Ztterr_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));


  // All bkg 
  TCanvas *c_bkg = new TCanvas("c_bkg", "c_bkg", 600, 500);
  c_bkg->cd(1);
  //c_bkg->SetLogz(1);
  c_bkg->SetRightMargin(0.15);
  hist_bkg->SetTitle("bkg");
  hist_bkg->SetXTitle(xvar);
  hist_bkg->SetYTitle(yvar);
//  hist_bkg->SetMinimum(histmin*0.5); 
//  hist_bkg->SetMaximum(histmax);
  hist_bkg->Draw("colz");
  c_bkg->SaveAs( plotDir + Form("/bkg_%s_mH%i_%ij_%s.eps",mvatype.Data(), mH, njet, flavor.Data()));
  c_bkg->SaveAs( plotDir + Form("/bkg_%s_mH%i_%ij_%s.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  // uncertainty histogram  
  TH2F* hist_bkg_error = filluncer(hist_bkg, hist_bkg);
  hist_bkg_error->Draw("colz");
  hist_bkg_error->SetTitle("bkgNAL");
  hist_bkg_error->SetXTitle(xvar);
  hist_bkg_error->SetYTitle(yvar);
  hist_bkg_error->SetMinimum(ERRORMIN);	
  hist_bkg_error->SetMaximum(1);
  c_bkg->SaveAs( plotDir + Form("/bkgerr_%s_mH%i_%ij_%s.eps", mvatype.Data(), mH, njet, flavor.Data()));
  c_bkg->SaveAs( plotDir + Form("/bkgerr_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));

  // ------------------------
  //  Draw 1D projection 
  // ------------------------
  
  // mT
  std::vector<TH1F*> hist_bkg_ProjX;
  THStack *stack_bkg_ProjX = new THStack(); 
  
  // define the stack and overlay legends..
  //TLegend *stacklg = new TLegend(0.0, 0.4, 1.0, 0.95);
  TLegend *stacklg = new TLegend(0.0, 0.1, 1.0, 0.95);
  stacklg->SetBorderSize(0);
  stacklg->SetFillStyle(0);
  stacklg->SetShadowColor(0);
  stacklg->SetTextSize(0.16);

  //TLegend *overlaylg = new TLegend(0.0, 0.4, 1.0, 0.95);
  TLegend *overlaylg = new TLegend(0.0, 0.1, 1.0, 0.95);
  overlaylg->SetBorderSize(0);
  overlaylg->SetFillStyle(0);
  overlaylg->SetShadowColor(0);
  overlaylg->SetTextSize(0.16);

  TCanvas *c_mt = new TCanvas("c_mt", "c_mt", 600, 500);
  c_mt->cd(1);
  c_mt->SaveAs( plotDir + Form("/mt_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));

  // signal 
  TH1F* hist_sig_ProjX = (TH1F*)hist_sig->ProjectionX(); 
  hist_sig_ProjX->SetLineColor(kRed+1);
  hist_sig_ProjX->SetMarkerColor(kRed+1);
  hist_sig_ProjX->SetLineWidth(3);

  // qqWW
  TH1F* hist_qqWW_ProjX = (TH1F*)hist_qqWW->ProjectionX(); 
  hist_qqWW_ProjX->SetLineColor(kAzure-9);
  hist_qqWW_ProjX->SetMarkerColor(kAzure-9);
  hist_qqWW_ProjX->SetFillColor(kAzure-9);
  hist_bkg_ProjX.push_back(hist_qqWW_ProjX);
  stack_bkg_ProjX->Add(hist_qqWW_ProjX);
  stacklg->AddEntry(hist_qqWW_ProjX, "qqWW", "f");
  overlaylg->AddEntry(hist_qqWW_ProjX, "qqWW", "l");

  // ggWW
  TH1F* hist_ggWW_ProjX = (TH1F*)hist_ggWW->ProjectionX(); 
  hist_ggWW_ProjX->SetLineColor(kAzure-9);
  hist_ggWW_ProjX->SetMarkerColor(kAzure-9);
  hist_ggWW_ProjX->SetFillColor(kAzure-9);
  hist_bkg_ProjX.push_back(hist_ggWW_ProjX);
  stack_bkg_ProjX->Add(hist_ggWW_ProjX);
  stacklg->AddEntry(hist_ggWW_ProjX, "ggWW", "f");
  overlaylg->AddEntry(hist_ggWW_ProjX, "ggWW", "l");

  // Wjets
  TH1F* hist_Wjets_ProjX = (TH1F*)hist_Wjets->ProjectionX(); 
  hist_Wjets_ProjX->SetLineColor(kGray+1);
  hist_Wjets_ProjX->SetMarkerColor(kGray+1);
  hist_Wjets_ProjX->SetFillColor(kGray+1);
  hist_bkg_ProjX.push_back(hist_Wjets_ProjX);
  stack_bkg_ProjX->Add(hist_Wjets_ProjX);
  stacklg->AddEntry(hist_Wjets_ProjX, "Wjets", "f");
  overlaylg->AddEntry(hist_Wjets_ProjX, "Wjets", "l");

  // Top
  TH1F* hist_Top_ProjX = (TH1F*)hist_Top->ProjectionX(); 
  hist_Top_ProjX->SetLineColor(kYellow);
  hist_Top_ProjX->SetMarkerColor(kYellow);
  hist_Top_ProjX->SetFillColor(kYellow);
  hist_bkg_ProjX.push_back(hist_Top_ProjX);
  stack_bkg_ProjX->Add(hist_Top_ProjX);
  stacklg->AddEntry(hist_Top_ProjX, "Top", "f");
  overlaylg->AddEntry(hist_Top_ProjX, "Top", "l");

  // VV
  TH1F* hist_VV_ProjX = (TH1F*)hist_VV->ProjectionX(); 
  hist_VV_ProjX->SetLineColor(kAzure-2);
  hist_VV_ProjX->SetMarkerColor(kAzure-2);
  hist_VV_ProjX->SetFillColor(kAzure-2);
  hist_bkg_ProjX.push_back(hist_VV_ProjX);
  stack_bkg_ProjX->Add(hist_VV_ProjX);
  stacklg->AddEntry(hist_VV_ProjX, "VV", "f");
  overlaylg->AddEntry(hist_VV_ProjX, "VV", "l");

  // Zjets
  TH1F* hist_Zjets_ProjX = (TH1F*)hist_Zjets->ProjectionX(); 
  hist_Zjets_ProjX->SetLineColor(kGreen+2);
  hist_Zjets_ProjX->SetMarkerColor(kGreen+2);
  hist_Zjets_ProjX->SetFillColor(kGreen+2);
  hist_bkg_ProjX.push_back(hist_Zjets_ProjX);
  stack_bkg_ProjX->Add(hist_Zjets_ProjX);
  stacklg->AddEntry(hist_Zjets_ProjX, "Zjets", "f");
  overlaylg->AddEntry(hist_Zjets_ProjX, "Zjets", "l");

  // Wgamma
  TH1F* hist_Wgamma_ProjX = (TH1F*)hist_Wgamma->ProjectionX(); 
  hist_Wgamma_ProjX->SetLineColor(kYellow+2);
  hist_Wgamma_ProjX->SetMarkerColor(kYellow+2);
  hist_Wgamma_ProjX->SetFillColor(kYellow+2);
  hist_bkg_ProjX.push_back(hist_Wgamma_ProjX);
  stack_bkg_ProjX->Add(hist_Wgamma_ProjX);
  stacklg->AddEntry(hist_Wgamma_ProjX, "Wgamma", "f");
  overlaylg->AddEntry(hist_Wgamma_ProjX, "Wgamma", "l");

  // Ztt
  TH1F* hist_Ztt_ProjX = (TH1F*)hist_Ztt->ProjectionX(); 
  hist_Ztt_ProjX->SetLineColor(kGreen+2);
  hist_Ztt_ProjX->SetMarkerColor(kGreen+2);
  hist_Ztt_ProjX->SetFillColor(kGreen+2);
  hist_bkg_ProjX.push_back(hist_Ztt_ProjX);
  stack_bkg_ProjX->Add(hist_Ztt_ProjX);
  stacklg->AddEntry(hist_Ztt_ProjX, "Ztt", "f");
  overlaylg->AddEntry(hist_Ztt_ProjX, "Ztt", "l");

  stacklg->AddEntry(hist_sig_ProjX, Form("HWW%i", mH), "l");
  overlaylg->AddEntry(hist_sig_ProjX, Form("HWW%i", mH), "l");

  //stacklg->AddEntry(hist_data, "Data", "lp");
  
  TH1F* hist_bkg_ProjX_tmp = (TH1F*)hist_bkg->ProjectionX(); 
  float yMax = hist_bkg_ProjX_tmp->GetMaximum();
/*
  float yMax = hist_sig_ProjX->GetMaximum();
  for (unsigned int i = 0; i < hist_bkg_ProjX.size(); i++) {
	  hist_bkg_ProjX[i]->SetLineWidth(3);
	  hist_bkg_ProjX[i]->SetFillStyle(1001);
	  hist_bkg_ProjX[i]->Draw("SAMEHIST");
	  yMax = yMax > hist_bkg_ProjX[i]->GetMaximum() ? yMax : hist_bkg_ProjX[i]->GetMaximum();
  }
*/
  hist_sig_ProjX->SetMaximum(1.4*yMax);
  stack_bkg_ProjX->SetMaximum(1.4*yMax);
  hist_sig_ProjX->SetMinimum(0.1);
  stack_bkg_ProjX->SetMinimum(0.1);

  // == Draw histograms
  // draw stacked plots
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 300);
  c1->cd();

  TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.0, 0.82, 1.0);
  pad1->SetBottomMargin(0.13);
  pad1->SetRightMargin(0.07);
  pad1->Draw();
  c1->cd();

  TPad *pad2 = new TPad("p_leg", "p_leg", 0.82, 0.0, 1.0, 1.0);
  pad2->SetTopMargin(0.01);
  pad2->SetRightMargin(0.01);
  pad2->SetBottomMargin(0.01);
  pad2->Draw();

  pad1->cd();
  pad1->SetLogy(0);
  stack_bkg_ProjX->Draw("HIST");
  //stack_bkg->GetXaxis()->SetTitle(mvatype);
  stack_bkg_ProjX->GetYaxis()->SetTitle(Form("Number of Events / %.2f", stack_bkg_ProjX->GetXaxis()->GetBinWidth(1)));
  hist_sig_ProjX->Draw("SAMEHIST");
  //hist_data->Draw("SAMEE1");
  pad2->cd();
  stacklg->Draw();

  // draw stacked plots - linear
  c1->SaveAs( plotDir + Form("/mtproj_%s_mH%i_%ij_%s_stack_lin.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  c1->SaveAs( plotDir + Form("/mtproj_%s_mH%i_%ij_%s_stack_lin.eps",mvatype.Data(), mH, njet, flavor.Data()));

  // draw stacked plots - log
  pad1->SetLogy();
  c1->SaveAs( plotDir + Form("/mtproj_%s_mH%i_%ij_%s_stack_log.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  c1->SaveAs( plotDir + Form("/mtproj_%s_mH%i_%ij_%s_stack_log.eps",mvatype.Data(), mH, njet, flavor.Data()));


  // === draw overlay plots
  c1->cd();
  pad1->Clear();
  pad1->cd();
  pad1->SetLogy(0);
  hist_sig_ProjX->Draw("HIST");
  for ( int s = 0; s < hist_bkg_ProjX.size(); s++) {
	  hist_bkg_ProjX[s]->SetFillColor(0);
	  hist_bkg_ProjX[s]->SetLineWidth(2);
	  hist_bkg_ProjX[s]->Draw("SAMEHIST");
  }
  pad2->Clear();
  pad2->cd();
  overlaylg->Draw();

  // draw overlay plots - linear

  c1->SaveAs( plotDir + Form("/mtproj_%s_mH%i_%ij_%s_overlay_lin.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  c1->SaveAs( plotDir + Form("/mtproj_%s_mH%i_%ij_%s_overlay_lin.eps",mvatype.Data(), mH, njet, flavor.Data()));

  // draw overlay plots - log
  pad1->SetLogy(1);
  c1->SaveAs( plotDir + Form("/mtproj_%s_mH%i_%ij_%s_overlay_log.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  c1->SaveAs( plotDir + Form("/mtproj_%s_mH%i_%ij_%s_overlay_log.eps",mvatype.Data(), mH, njet, flavor.Data()));

  delete pad1;
  delete pad2;
  // end of drawing histograms




  // mll 
  std::vector<TH1F*> hist_bkg_ProjY;
  THStack *stack_bkg_ProjY = new THStack(); 
  
  // define the stack and overlay legends..
  //TLegend *stacklg2 = new TLegend(0.0, 0.4, 1.0, 0.95);
  TLegend *stacklg2 = new TLegend(0.0, 0.1, 1.0, 0.95);
  stacklg2->SetBorderSize(0);
  stacklg2->SetFillStyle(0);
  stacklg2->SetShadowColor(0);
  stacklg2->SetTextSize(0.16);

  //TLegend *overlaylg2 = new TLegend(0.0, 0.4, 1.0, 0.95);
  TLegend *overlaylg2 = new TLegend(0.0, 0.1, 1.0, 0.95);
  overlaylg2->SetBorderSize(0);
  overlaylg2->SetFillStyle(0);
  overlaylg2->SetShadowColor(0);
  overlaylg2->SetTextSize(0.16);

  TCanvas *c_mll = new TCanvas("c_mll", "c_mll", 600, 500);
  c_mll->cd(1);
  c_mll->SaveAs( plotDir + Form("/mll_%s_mH%i_%ij_%s.pdf", mvatype.Data(), mH, njet, flavor.Data()));

  // signal 
  TH1F* hist_sig_ProjY = (TH1F*)hist_sig->ProjectionY(); 
  hist_sig_ProjY->SetLineColor(kRed+1);
  hist_sig_ProjY->SetMarkerColor(kRed+1);
  hist_sig_ProjY->SetLineWidth(3);

  // qqWW
  TH1F* hist_qqWW_ProjY = (TH1F*)hist_qqWW->ProjectionY(); 
  hist_qqWW_ProjY->SetLineColor(kAzure-9);
  hist_qqWW_ProjY->SetMarkerColor(kAzure-9);
  hist_qqWW_ProjY->SetFillColor(kAzure-9);
  hist_bkg_ProjY.push_back(hist_qqWW_ProjY);
  stack_bkg_ProjY->Add(hist_qqWW_ProjY);
  stacklg2->AddEntry(hist_qqWW_ProjY, "qqWW", "f");
  overlaylg2->AddEntry(hist_qqWW_ProjY, "qqWW", "l");

  // ggWW
  TH1F* hist_ggWW_ProjY = (TH1F*)hist_ggWW->ProjectionY(); 
  hist_ggWW_ProjY->SetLineColor(kAzure-9);
  hist_ggWW_ProjY->SetMarkerColor(kAzure-9);
  hist_ggWW_ProjY->SetFillColor(kAzure-9);
  hist_bkg_ProjY.push_back(hist_ggWW_ProjY);
  stack_bkg_ProjY->Add(hist_ggWW_ProjY);
  stacklg2->AddEntry(hist_ggWW_ProjY, "ggWW", "f");
  overlaylg2->AddEntry(hist_ggWW_ProjY, "ggWW", "l");

  // Wjets
  TH1F* hist_Wjets_ProjY = (TH1F*)hist_Wjets->ProjectionY(); 
  hist_Wjets_ProjY->SetLineColor(kGray+1);
  hist_Wjets_ProjY->SetMarkerColor(kGray+1);
  hist_Wjets_ProjY->SetFillColor(kGray+1);
  hist_bkg_ProjY.push_back(hist_Wjets_ProjY);
  stack_bkg_ProjY->Add(hist_Wjets_ProjY);
  stacklg2->AddEntry(hist_Wjets_ProjY, "Wjets", "f");
  overlaylg2->AddEntry(hist_Wjets_ProjY, "Wjets", "l");

  // Top
  TH1F* hist_Top_ProjY = (TH1F*)hist_Top->ProjectionY(); 
  hist_Top_ProjY->SetLineColor(kYellow);
  hist_Top_ProjY->SetMarkerColor(kYellow);
  hist_Top_ProjY->SetFillColor(kYellow);
  hist_bkg_ProjY.push_back(hist_Top_ProjY);
  stack_bkg_ProjY->Add(hist_Top_ProjY);
  stacklg2->AddEntry(hist_Top_ProjY, "Top", "f");
  overlaylg2->AddEntry(hist_Top_ProjY, "Top", "l");

  // VV
  TH1F* hist_VV_ProjY = (TH1F*)hist_VV->ProjectionY(); 
  hist_VV_ProjY->SetLineColor(kAzure-2);
  hist_VV_ProjY->SetMarkerColor(kAzure-2);
  hist_VV_ProjY->SetFillColor(kAzure-2);
  hist_bkg_ProjY.push_back(hist_VV_ProjY);
  stack_bkg_ProjY->Add(hist_VV_ProjY);
  stacklg2->AddEntry(hist_VV_ProjY, "VV", "f");
  overlaylg2->AddEntry(hist_VV_ProjY, "VV", "l");

  // Zjets
  TH1F* hist_Zjets_ProjY = (TH1F*)hist_Zjets->ProjectionY(); 
  hist_Zjets_ProjY->SetLineColor(kGreen+2);
  hist_Zjets_ProjY->SetMarkerColor(kGreen+2);
  hist_Zjets_ProjY->SetFillColor(kGreen+2);
  hist_bkg_ProjY.push_back(hist_Zjets_ProjY);
  stack_bkg_ProjY->Add(hist_Zjets_ProjY);
  stacklg2->AddEntry(hist_Zjets_ProjY, "Zjets", "f");
  overlaylg2->AddEntry(hist_Zjets_ProjY, "Zjets", "l");

  // Wgamma
  TH1F* hist_Wgamma_ProjY = (TH1F*)hist_Wgamma->ProjectionY(); 
  hist_Wgamma_ProjY->SetLineColor(kYellow+2);
  hist_Wgamma_ProjY->SetMarkerColor(kYellow+2);
  hist_Wgamma_ProjY->SetFillColor(kYellow+2);
  hist_bkg_ProjY.push_back(hist_Wgamma_ProjY);
  stack_bkg_ProjY->Add(hist_Wgamma_ProjY);
  stacklg2->AddEntry(hist_Wgamma_ProjY, "Wgamma", "f");
  overlaylg2->AddEntry(hist_Wgamma_ProjY, "Wgamma", "l");

  // Ztt
  TH1F* hist_Ztt_ProjY = (TH1F*)hist_Ztt->ProjectionY(); 
  hist_Ztt_ProjY->SetLineColor(kGreen+2);
  hist_Ztt_ProjY->SetMarkerColor(kGreen+2);
  hist_Ztt_ProjY->SetFillColor(kGreen+2);
  hist_bkg_ProjY.push_back(hist_Ztt_ProjY);
  stack_bkg_ProjY->Add(hist_Ztt_ProjY);
  stacklg2->AddEntry(hist_Ztt_ProjY, "Ztt", "f");
  overlaylg2->AddEntry(hist_Ztt_ProjY, "Ztt", "l");

  stacklg2->AddEntry(hist_sig_ProjY, Form("HWW%i", mH), "l");
  overlaylg2->AddEntry(hist_sig_ProjY, Form("HWW%i", mH), "l");

  //stacklg2->AddEntry(hist_data, "Data", "lp");

  TH1F* hist_bkg_ProjY_tmp = (TH1F*)hist_bkg->ProjectionY(); 
  yMax = hist_bkg_ProjY_tmp->GetMaximum(); 
  /*
  for (unsigned int i = 0; i < hist_bkg_ProjY.size(); i++) {
	  hist_bkg_ProjY[i]->SetLineWidth(3);
	  hist_bkg_ProjY[i]->SetFillStyle(1001);
	  hist_bkg_ProjY[i]->Draw("SAMEHIST");
	  yMax = yMax > hist_bkg_ProjY[i]->GetMaximum() ? yMax : hist_bkg_ProjY[i]->GetMaximum();
  }
  */
  hist_sig_ProjY->SetMaximum(1.4*yMax);
  stack_bkg_ProjY->SetMaximum(1.4*yMax);
  hist_sig_ProjY->SetMinimum(0.1);
  stack_bkg_ProjY->SetMinimum(0.1);

  // == Draw histograms
  // draw stacked plots
  TCanvas *c2 = new TCanvas("c2", "c2", 600, 300);
  c2->cd();

  TPad *pad3 = new TPad("p_main", "p_main", 0.0, 0.0, 0.82, 1.0);
  pad3->SetBottomMargin(0.13);
  pad3->SetRightMargin(0.07);
  pad3->Draw();
  c2->cd();

  TPad *pad4 = new TPad("p_leg", "p_leg", 0.82, 0.0, 1.0, 1.0);
  pad4->SetTopMargin(0.01);
  pad4->SetRightMargin(0.01);
  pad4->SetBottomMargin(0.01);
  pad4->Draw();

  pad3->cd();
  pad3->SetLogy(0);
  stack_bkg_ProjY->Draw("HIST");
  //stack_bkg->GetXaxis()->SetTitle(mvatype);
  stack_bkg_ProjY->GetYaxis()->SetTitle(Form("Number of Events / %.2f", stack_bkg_ProjY->GetXaxis()->GetBinWidth(1)));
  hist_sig_ProjY->Draw("SAMEHIST");
  //hist_data->Draw("SAMEE1");
  pad4->cd();
  stacklg2->Draw();

  // draw stacked plots - linear
  c2->SaveAs( plotDir + Form("/mllproj_%s_mH%i_%ij_%s_stack_lin.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  c2->SaveAs( plotDir + Form("/mllproj_%s_mH%i_%ij_%s_stack_lin.eps",mvatype.Data(), mH, njet, flavor.Data()));

  // draw stacked plots - log
  pad3->SetLogy();
  c2->SaveAs( plotDir + Form("/mllproj_%s_mH%i_%ij_%s_stack_log.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  c2->SaveAs( plotDir + Form("/mllproj_%s_mH%i_%ij_%s_stack_log.eps",mvatype.Data(), mH, njet, flavor.Data()));


  // === draw overlay plots
  c2->cd();
  pad3->Clear();
  pad3->cd();
  pad3->SetLogy(0);
  hist_sig_ProjY->Draw("HIST");
  for ( int s = 0; s < hist_bkg_ProjY.size(); s++) {
	  hist_bkg_ProjY[s]->SetFillColor(0);
	  hist_bkg_ProjY[s]->SetLineWidth(2);
	  hist_bkg_ProjY[s]->Draw("SAMEHIST");
  }
  pad4->Clear();
  pad4->cd();
  overlaylg2->Draw();

  // draw overlay plots - linear

  c2->SaveAs( plotDir + Form("/mllproj_%s_mH%i_%ij_%s_overlay_lin.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  c2->SaveAs( plotDir + Form("/mllproj_%s_mH%i_%ij_%s_overlay_lin.eps",mvatype.Data(), mH, njet, flavor.Data()));

  // draw overlay plots - log
  pad3->SetLogy(1);
  c2->SaveAs( plotDir + Form("/mllproj_%s_mH%i_%ij_%s_overlay_log.pdf",mvatype.Data(), mH, njet, flavor.Data()));
  c2->SaveAs( plotDir + Form("/mllproj_%s_mH%i_%ij_%s_overlay_log.eps",mvatype.Data(), mH, njet, flavor.Data()));

  delete pad3;
  delete pad4;
  // end of drawing histograms



  //============================================
  
  
  // tidy up
  delete hist_data;
  delete hist_bkg;
  delete hist_sig; 
  delete hist_ggH; 
  delete hist_qqH; 
  delete hist_WH; 
  delete hist_ZH; 
  delete hist_qqWW; 
  delete hist_ggWW; 
  delete hist_Wjets; 
  delete hist_Top; 
  delete hist_VV; 
  delete hist_Zjets; 
  delete hist_Wgamma; 
  delete hist_Ztt; 
  
  delete hist_data_error;
  delete hist_sig_error; 
  delete hist_qqWW_error; 
  delete hist_ggWW_error; 
  delete hist_Wjets_error; 
  delete hist_Top_error; 
  delete hist_VV_error; 
  delete hist_Zjets_error; 
  delete hist_Wgamma_error; 
  delete hist_Ztt_error; 
  
  delete c_data;
  delete c_sig; 
  delete c_qqWW; 
  delete c_ggWW; 
  delete c_Wjets; 
  delete c_Top; 
  delete c_VV; 
  delete c_Zjets; 
  delete c_Wgamma; 
  delete c_Ztt; 
  
  File->Close();

}

