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
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void drawsingle(TString mvaType, TString cardDir, TString outputDir, ofstream &text,  int mH, int njet, const char*flavor, const char* process, const char* systvar, bool flavordep, bool hwwonly);
void drawshapesystsinglemass(TString mvaType, ofstream & text, int mH, int njet, TString cardDir, TString outputDir);


void drawhwwshapevar(TString mvaType = "2D", int mH = 125, TString cardDir="../cards/2D/", TString outputDir="plots/2D/Unrolled") 
{
  ofstream text;
  text.open(outputDir+Form("/shapevar_mH%i_%s.tex",  mH, mvaType.Data()));
  
  text << "\\documentclass{article} \n";
  text << "\\usepackage{times}\n"; 
  text << "\\usepackage{epsfig}\n"; 
  text << "\\begin{document}\n";
  
  drawshapesystsinglemass(mvaType, text, mH, 0, cardDir, outputDir);	text << "\\newpage\n";
  //drawshapesystsinglemass(mvaType, text, mH, 1, cardDir, outputDir);	text << "\\newpage\n";
  //drawshapesystsinglemass(mvaType, text, mH, 2, cardDir, outputDir);

  text << "\\clearpage\n\n";
  text << "\n\n";

  text << "\\end{document}\n"; 
  
}




void drawshapesystsinglemass(TString mvaType, ofstream & text, int mH, int njet, TString cardDir, TString outputDir)
{

  const char *table_head = "\\begin{figure}[!htb]\n\\begin{center}\n\\begin{tabular}{cc}\n";
  const char *table_end = "\\end{tabular}\n\\end{center}\n";
  const char *table_caption = "The central shape is shown with black lines with error bars. The up and down shape variations are shown with red and blue dashed lines respectively. The ratio between the up/down shapes to the central one is shown as well. }\n\\end{figure}\n";
  
  // write down the plots in latex format

  // ggH variations
  // text << table_head;
  //  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "sf", "ggH", "ggHBounding", false, false);
  // drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "ggH", "ggHBounding", false, false);
  // text << table_end;
  // text << Form("\\caption{%s output shape variations for ggH due to higher order effects in sf (left) and of (right) for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // WW variations
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "sf", "qqWW", "WWBounding", false, true);
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "qqWW", "WWBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for WW due to higher order effects from mc@nlo and madgraph differences in sf (left) and of (right) for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // WW variations
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "sf", "qqWW", "WWNLOBounding", false, true);
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "qqWW", "WWNLOBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for WW due to qcd scale variations from mc@nlo in sf (left) and of (right) for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;
/*
  // Wjets variations
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "sf", "Wjets", "WBounding", false, true);
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "Wjets", "WBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for Wjets due to FR variations in data, in sf (left) and of (right) for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;
*/
 /* 
  // Wjets variations
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "sf", "Wjets", "WMCBounding", false, true);
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "Wjets", "WMCBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for Wjets due to data/MC differences, in sf (left) and of (right) for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;
*/ 
/*
  // Zjets variations
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "sf", "Zjets", "ZBounding", true, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for Zjets due to the difference between the loose met region in data and the signal region in MC in the sf for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;
  
 
  // Top variation
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "sf", "Top", "TopBounding", false, true);
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "Top", "TopBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for Top due to powheg/madgraph differences, in sf (left) and of (right) for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // MET variation for qqWW
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "sf", "qqWW", "METResBounding", false, true);
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "qqWW", "METResBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for qqWW due to MET resolution, in sf (left) and of (right) for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // LEPRES variation for qqWW
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "sf", "qqWW", "LepResBounding", false, true);
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "qqWW", "LepResBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for qqWW due to lepton resolution, in sf (left) and of (right) for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // LepEff variation for qqWW
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "sf", "qqWW", "LepEffBounding", false, true);
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "qqWW", "LepEffBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for qqWW due to lepton efficiency, in sf (left) and of (right) for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // LEPRES variation for qqWW
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "sf", "qqWW", "JESBounding", false, true);
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "qqWW", "JESBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for qqWW due to Jet Energy Scale, in sf (left) and of (right) for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;
*/
}
    

void drawsingle(TString mvaType, TString cardDir, TString outputDir, ofstream & text, int mH, int njet, const char* flavor, const char* process, const char* systvar, bool flavordep, bool hwwonly)
{  
  
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  gROOT->ForceStyle();

  TString runera = "_8TeV";
  
  TString fileName;
  
  if ( mvaType == "BDT") 
    fileName = Form("%s/%i/hww%s_%ij.input%s.root", cardDir.Data(), mH, flavor, njet, runera.Data());
  if ( mvaType == "ME") 
    fileName = Form("%s/%i/hww%s_%ij_me.input%s.root", cardDir.Data(), mH, flavor, njet, runera.Data());
  if ( mvaType == "2D") 
    fileName = Form("%s/%i/hww%s_%ij.input%s.root", cardDir.Data(), mH, flavor, njet, runera.Data());
  
  TFile *file = TFile::Open(fileName, "READ");
  if ( file == 0x0 ) {
    std::cout << "Error " << fileName << " does not exist, exitting..\n";
  }
  
  TString centralHistName = Form("histo_%s", process);
  TH1F *hist_central = (TH1F*) file->Get(centralHistName);
  if ( hist_central == 0x0 || hist_central->Integral() <= 0.0) {
    std::cout << "Error " << centralHistName << " is not in " << fileName << "; exitting \n";
    return;
  }

  TString upHistName; 
  TString downHistName; 
  
  if ( flavordep) {
    if ( ! hwwonly ) {
      upHistName = Form("histo_%s_CMS_hww_MVA%s_%sUp", process, systvar, flavor);
      downHistName = Form("histo_%s_CMS_hww_MVA%s_%sDown", process, systvar, flavor);
    } else {
      upHistName = Form("histo_%s_CMS_hww_MVA%s%sUp", process, systvar, flavor);
      downHistName = Form("histo_%s_CMS_hww_MVA%sDown", process, systvar, flavor);
    }
  }
  else {
    if ( ! hwwonly) {
      upHistName = Form("histo_%s_CMS_hww_MVA%sUp", process, systvar);
      downHistName = Form("histo_%s_CMS_hww_MVA%sDown", process, systvar);
    }
    else {
      upHistName = Form("histo_%s_CMS_hww_MVA%sUp", process, systvar);
      downHistName = Form("histo_%s_CMS_hww_MVA%sDown", process, systvar);
    }
  }
  
  // special cases

  if ( TString(systvar) == "ZBounding") {
    upHistName = Form("histo_%s_CMS_hww%s_%ij_MVA%sUp", process, flavor, njet, systvar);
    downHistName = Form("histo_%s_CMS_hww%s_%ij_MVA%sDown", process, flavor, njet, systvar);

  }
    


  TH1F *hist_up = (TH1F*) file->Get(upHistName);
  if ( hist_up == 0x0 || hist_up->Integral() <= 0.0 ) {
    std::cout << "Error " << upHistName << " is not in " << fileName << "; exitting \n";
    return;
  }

  TH1F *hist_down = (TH1F*) file->Get(downHistName);
  if ( hist_down == 0x0 || hist_down->Integral() <= 0.0 ) {
    std::cout << "Error " << downHistName << " is not in " << fileName << "; exitting \n";
    return;
  }
  
  // draw overlay 
  float yMax = hist_up->GetMaximum() > hist_central->GetMaximum() ? hist_up->GetMaximum() : hist_central->GetMaximum();
  yMax = yMax > hist_down->GetMaximum() ? yMax : hist_down->GetMaximum();

  hist_central->SetMaximum(yMax * 1.2);
  hist_up->SetMaximum(yMax * 1.2);
  hist_down->SetMaximum(yMax * 1.2);
  
  TString xTitle = Form("%s Output", mvaType.Data());

  hist_central->SetXTitle(xTitle);
  hist_up->SetXTitle(xTitle);
  hist_down->SetXTitle(xTitle);

  hist_central->SetLineWidth(3);

  hist_central->SetMinimum(0.);
  hist_up->SetMinimum(0.);
  hist_down->SetMinimum(0.);

  //TString yTitle = Form("Number of Events / %.1f GeV", hist_up->GetBinWidth(1));
  TString yTitle = "";
  hist_central->SetYTitle(yTitle);
  hist_up->SetYTitle(yTitle);
  hist_down->SetYTitle(yTitle);
  
  hist_up->SetLineColor(kRed);
  hist_up->SetMarkerColor(kRed);
  hist_up->SetLineStyle(2);
  hist_up->SetLineWidth(3);
    
  hist_down->SetLineColor(kBlue);
  hist_down->SetMarkerColor(kBlue);
  hist_down->SetLineStyle(2);
  hist_down->SetLineWidth(3);

  // get legend
  TLegend *leg = new TLegend(0.4, 0.75, 0.6, 0.92);
  leg->SetFillColor(0);
  leg->AddEntry(hist_central, "Central");
  leg->AddEntry(hist_up, "Up");
  leg->AddEntry(hist_down, "Down");

  TH1F* hist_ratio_up = (TH1F*) hist_up->Clone("hist_ratio_up");
  hist_ratio_up->Divide(hist_ratio_up, hist_central, 1, 1);
  hist_ratio_up->GetXaxis()->SetTitleSize(0.15);
  hist_ratio_up->GetYaxis()->SetTitleSize(0.12);
  hist_ratio_up->GetYaxis()->SetTitleOffset(0.5);
  hist_ratio_up->GetXaxis()->SetLabelSize(0.1);
  hist_ratio_up->GetYaxis()->SetLabelSize(0.1);
  hist_ratio_up->GetYaxis()->SetRangeUser(0.0, 2.0);
  hist_ratio_up->GetYaxis()->SetNdivisions(504);
  

  TH1F* hist_ratio_down = (TH1F*) hist_down->Clone("hist_ratio_down");
  hist_ratio_down->Divide(hist_ratio_down, hist_central, 1, 1);
  hist_ratio_down->GetXaxis()->SetTitleSize(0.15);
  hist_ratio_down->GetYaxis()->SetTitleSize(0.12);
  hist_ratio_down->GetYaxis()->SetTitleOffset(0.5);
  hist_ratio_down->GetXaxis()->SetLabelSize(0.1);
  hist_ratio_down->GetYaxis()->SetLabelSize(0.1);
  hist_ratio_down->GetYaxis()->SetRangeUser(0.0, 2.0);
  hist_ratio_down->GetYaxis()->SetNdivisions(504);


  TCanvas *c1 = new TCanvas("c1", "c1", 600, 800);
  c1->cd();
  
  //TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.3, 1.0, 1.0);
  TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.5, 1.0, 1.0);
  pad1->SetBottomMargin(0.12);
  pad1->SetRightMargin(0.07);
  pad1->SetLeftMargin(0.18);
  pad1->Draw();
  
  c1->cd();
  //TPad *pad2 = new TPad("p_leg", "p_leg", 0.0, 0.0, 1.0, 0.3);
  TPad *pad2 = new TPad("p_leg", "p_leg", 0.0, 0.0, 1.0, 0.5);
  pad2->SetTopMargin(0.01);
  pad2->SetLeftMargin(0.18);
  pad2->SetRightMargin(0.07);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();

  pad1->cd();
  
  hist_central->Draw("HISTE0");
  hist_up->Draw("SAMEHIST");
  hist_down->Draw("SAMEHIST");
  
  leg->Draw("same");
  
  pad2->cd();  
  float min = hist_central->GetXaxis()->GetXmin();
  float max = hist_central->GetXaxis()->GetXmax();
  TLine *line = new TLine(min, 1, max, 1);
  line->SetLineWidth(2);
  TLine *line_up = new TLine(min, 1.25, max, 1.25);
  line_up->SetLineWidth(2);
  line_up->SetLineStyle(2);
  TLine *line_down = new TLine(min, 0.75, max, 0.75);
  line_down->SetLineWidth(2);
  line_down->SetLineStyle(2);
  
  hist_ratio_up->SetMinimum(0.75); 
  hist_ratio_up->SetMaximum(1.25); 

  hist_ratio_up->Draw("e0");
  hist_ratio_down->Draw("samee0");
  line->Draw("same");
  line_up->Draw("same");
  line_down->Draw("same"); 
  
  c1->cd()->SetLogy(1);

  c1->SaveAs(Form("%s/%s_%s_mT_mH%i_%ij_%s_lin.png", outputDir.Data(), process, systvar, mH, njet, flavor));
  c1->SaveAs(Form("%s/%s_%s_mT_mH%i_%ij_%s_lin.eps", outputDir.Data(), process, systvar, mH, njet, flavor));
  
  // pad1->SetLogy();

  // tidy up
  file->Close();
  delete c1;

  text << Form("\\epsfig{figure=%s/%s_%s_mT_mH%i_%ij_%s_lin.eps, width=2.0in}\n", outputDir.Data(), process, systvar, mH, njet, flavor);

}


