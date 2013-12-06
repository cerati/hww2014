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

void drawOverlay(TString fileName1, TString fileName2, TString mvatype, int njet, int mH, TString plotDir, TString flavor, TString histName, TString input1Name, TString input2Name);

//###################
//# main function
//###################
void drawhwwcarddiff(TString inputFDir1 = "../cards/ana_v7_Full2011/", TString inputFDir2 = "ana_v7_Full2011_Guil", 
		     int mH = 115, int njet = 0,  TString flavor = "of", TString plotDir = "plots/carddiff/", TString mvatype="BDT", 
		     TString input1Name = "Yanyan", TString input2Name = "Guillelmo") {
  

  TString fileName1 =  inputFDir1+Form("/%i/hww%s_%ij.input_8TeV.root", mH, flavor.Data(), njet);
  if ( mvatype == "ME") 
    fileName1 =  inputFDir1+Form("/%i/hww%s_%ij_me.input_8TeV.root", mH, flavor.Data(), njet);

  TString fileName2 =  inputFDir2+Form("/%i/hww%s_%ij.input_8TeV.root", mH, flavor.Data(), njet);
  if ( mvatype == "ME") 
    fileName2 =  inputFDir2+Form("/%i/hww%s_%ij_me.input_8TeV.root", mH, flavor.Data(), njet);  

  bool blind = false;
  if ( !blind) 
    drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_Data", input1Name, input2Name);
  drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_ggH", input1Name, input2Name);
  drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_qqH", input1Name, input2Name);
  drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_qqWW", input1Name, input2Name);
  drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_ggWW", input1Name, input2Name);
  drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_Top", input1Name, input2Name);
  drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_Wjets", input1Name, input2Name);
  drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_VV", input1Name, input2Name);
  if ( flavor == "sf" ) {
    drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_Zjets", input1Name, input2Name);
    drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, Form("histo_Zjets_CMS_hwwsf_%ij_MVAZBoundingUp", njet), input1Name, input2Name);
  }
  drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_Wgamma", input1Name, input2Name);
  drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_qqWW_CMS_hww_MVAWWBoundingUp", input1Name, input2Name);
  drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_qqWW_CMS_hww_MVAWWBoundingDown", input1Name, input2Name);
  drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_Top_CMS_hww_MVATopBoundingUp", input1Name, input2Name);
  drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_Top_CMS_hww_MVATopBoundingDown", input1Name, input2Name);
  
  // drawOverlay(fileName1, fileName2, mvatype, njet, mH, plotDir, flavor, "histo_Ztt", input1Name, input2Name);

}  

void drawOverlay(TString fileName1, TString fileName2, TString mvatype, int njet, int mH, TString plotDir, TString flavor, TString histName,
		 TString input1Name, TString input2Name) {
  // make it look nice
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  gROOT->ForceStyle();

  TFile *mvaFile1 = TFile::Open(fileName1, "READ");
  if ( mvaFile1 == 0x0 ) {
    std::cout << "Error, " << fileName1 << " does not exist, exitting...\n";
    return;
  }
  std::cout << "Opening..." << fileName1 << "\n";
  gROOT->cd();

  
  TH1F *hist1 = (TH1F*) mvaFile1->Get(histName);
  
 TFile *mvaFile2 = TFile::Open(fileName2, "READ");
  if ( mvaFile2 == 0x0 ) {
    std::cout << "Error, " << fileName2 << " does not exist, exitting...\n";
    return;
  }
  std::cout << "Opening..." << fileName2 << "\n";
  gROOT->cd();

  
  TH1F *hist2 = (TH1F*) mvaFile2->Get(histName);
  
  
  TLegend *leg = new TLegend(0.2, 0.75, 0.6, 0.92, histName);                                                                                  
  leg->SetFillColor(0);                                                                                                              
  leg->SetTextSize(0.04);
  leg->AddEntry(hist1, input1Name, "l");                                                                                            
  leg->AddEntry(hist2, input2Name, "pe");                                                                                                      

  hist1->SetXTitle(Form("%s output", mvatype.Data()) );
  hist2->SetXTitle(Form("%s output", mvatype.Data()) );
  hist1->SetLineColor(kRed);
  hist1->SetMarkerColor(kRed);
  hist1->SetMarkerSize(0);

  hist2->SetLineColor(kBlue);
  hist2->SetMarkerColor(kBlue);
  hist2->SetMarkerSize(1);

  hist1->SetMaximum(hist1->GetMaximum()*1.5);

  // == Draw

  // == Draw histograms
  // draw stacked plots
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
  c1->cd();
  

  hist1->Draw("hist");
  hist2->Draw("samee0");
  leg->Draw("same");

  c1->SaveAs( plotDir + Form("/%s_%s_mH%i_%ij_%s.png",mvatype.Data(), histName.Data(), mH, njet, flavor.Data()));
  c1->SaveAs( plotDir + Form("/%s_%s_mH%i_%ij_%s.eps",mvatype.Data(), histName.Data(), mH, njet, flavor.Data()));

  
  // tidy up
  delete c1;
  mvaFile1->Close();
  mvaFile2->Close();

}

