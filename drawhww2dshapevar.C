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
#include "TH2.h"
#include "TF1.h"
#include "TAxis.h"
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void drawsingle(TString mvaType, TString cardDir, TString outputDir, ofstream &text,  int mH, int njet, const char*flavor, const char* process, const char* systvar, bool flavordep, bool hwwonly);
void drawshapesystsinglemass(TString mvaType, ofstream & text, int mH, int njet, TString cardDir, TString outputDir);
bool zoom = true; // global variable for whether to draw it in the signal region

// rolling 
TH2F* Roll1DTo2D(TH1* h1, const char* hname="") {

	unsigned int nbins  = h1->GetXaxis()->GetNbins();
	//TH2F* h2 = new TH2F(hname, hname, 10, 80, 280, 8, 0, 200);     
        float mtbins[15]        =       {60,70,80,90,100,110,120,140,160,180,200,220,240,260,280};
        float mllbins[10]       =       {12,30,45,60,75,100,125,150,175,200}; 
        TH2F* h2 = new TH2F(hname, hname, 14, mtbins, 9, mllbins);

	for(unsigned int ibin=1; ibin <= nbins; ++ibin) {
	  //		int ibinX = (int) (ibin-1)/8+1;
	  //		int ibinY = (int) (ibin-1)%8+1;
	  	int ibinX = (int) (ibin-1)/9+1;
	  	int ibinY = (int) (ibin-1)%9+1;
		h2->SetBinContent(  ibinX, ibinY, h1->GetBinContent(ibin) );
		h2->SetBinError(    ibinX, ibinY, h1->GetBinError(ibin) );
	}
	h2->SetStats(0);
	return h2;
}

void drawhww2dshapevar(TString mvaType = "2D", int mH = 125, TString cardDir="../cards/hwwjcp_19fb/", TString outputDir="gghshapevar/") 
{
  ofstream text;
  text.open(outputDir+Form("/shapevar_mH%i_%s.tex",  mH, mvaType.Data()));
  
  text << "\\documentclass{article} \n";
  text << "\\usepackage{times}\n"; 
  text << "\\usepackage{epsfig}\n"; 
  text << "\\begin{document}\n";
  

  drawshapesystsinglemass(mvaType, text, mH, 0, cardDir, outputDir);
  // drawshapesystsinglemass(mvaType, text, mH, 1, cardDir, outputDir);

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

  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "ggH", "ggHBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for ggH due to higher order effects from jhugen and powheg differences for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;


 // MET variation for ggH
  text << table_head;
   drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "ggH", "METResBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for ggH due to MET resolution, in of for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // LEPRES variation for ggH
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "ggH", "LepResBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for ggH due to lepton resolution, in of for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // LepEff variation for ggH
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "ggH", "LepEffBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for ggH due to lepton efficiency, in of for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // JES variation for ggH
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "ggH", "JESBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for ggH due to Jet Energy Scale, in of for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // WW variations
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "qqWW", "WWBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for WW due to higher order effects from mc@nlo and madgraph differences in of for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // WW variations
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "qqWW", "WWNLOBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for WW due to qcd scale variations from mc@nlo in of for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;


  // Wjets variations
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "Wjets", "WBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for Wjets due to FR variations in data, in of for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // Top variation
  text << table_head;

  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "Top", "TopBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for Top due to powheg/madgraph differences, in of for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // MET variation for qqWW
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "qqWW", "METResBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for qqWW due to MET resolution, in of for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // LEPRES variation for qqWW
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "qqWW", "LepResBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for qqWW due to lepton resolution, in of for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // LepEff variation for qqWW
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "qqWW", "LepEffBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for qqWW due to lepton efficiency, in of for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

  // JES variation for qqWW
  text << table_head;
  drawsingle(mvaType, cardDir, outputDir, text, mH, njet, "of", "qqWW", "JESBounding", false, true);
  text << table_end;
  text << Form("\\caption{%s output shape variations for qqWW due to Jet Energy Scale, of (right) for mH=%i GeV in the %i jet bin.", mvaType.Data(), mH, njet) << table_caption;

}
    

void drawsingle(TString mvaType, TString cardDir, TString outputDir, ofstream & text, int mH, int njet, const char* flavor, const char* process, const char* systvar, bool flavordep, bool hwwonly)
{  
  
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle();");
  gROOT->ForceStyle();
  gStyle->SetPaintTextFormat("4.0f");
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleXOffset(1.1);                                                                                   
  gStyle->SetNdivisions(506, "XYZ");
  gStyle->SetPaintTextFormat(".1f")                                               ;
  gROOT->ForceStyle();
  gStyle->SetOptTitle(1);
  
  TString runera = "_8TeV";
  
  TString fileName;
  
  if ( mvaType == "BDT") 
    fileName = Form("%s/%i/hww%s_%ij.input%s.root", cardDir.Data(), mH, flavor, njet, runera.Data());
  if ( mvaType == "ME") 
    fileName = Form("%s/%i/hww%s_%ij_me.input%s.root", cardDir.Data(), mH, flavor, njet, runera.Data());
  if ( mvaType == "2D") 
    fileName = Form("%s/%i/xww%s_%ij.input%s.root", cardDir.Data(), mH, flavor, njet, runera.Data());
  
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
      downHistName = Form("histo_%s_CMS_hww_MVA%s%sDown", process, systvar, flavor);
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
 
  // rolling central / up / down
  TH2F* h2_central 	= Roll1DTo2D(hist_central); 
  h2_central->SetTitle("Central");
  if ( zoom ) {
    h2_central->GetXaxis()->SetRangeUser(59, 119.);    
    h2_central->GetYaxis()->SetRangeUser(0, 70.);
  }

  TH2F* h2_up 		= Roll1DTo2D(hist_up);
  h2_up->SetTitle("Up");
  TH2F* h2_down 	= Roll1DTo2D(hist_down);
  h2_down->SetTitle("Down");

  // ratio = ( up(down) - central ) / central
  TH2F* h2_up_minus_central 		= (TH2F*) h2_up->Clone(); 
  h2_up_minus_central->SetTitle("(Up - central)/central (%)"); 
  h2_up_minus_central->Add(h2_central, -1); 
  h2_up_minus_central->Divide(h2_central); 
  h2_up_minus_central->Scale(100.); 
  h2_up_minus_central->GetZaxis()->SetRangeUser(-50, 50.);
  if ( zoom ) {
    h2_up_minus_central->GetXaxis()->SetRangeUser(59, 119.);    
    h2_up_minus_central->GetYaxis()->SetRangeUser(0, 70.);
    h2_up_minus_central->GetZaxis()->SetRangeUser(-50, 50.);
    h2_up_minus_central->SetMarkerSize(3);
  }
  
  TH2F* h2_down_minus_central 	= (TH2F*) h2_down->Clone(); 
  h2_down_minus_central->SetTitle("(Down - central)/central (%)"); 
  h2_down_minus_central->Add(h2_central, -1); 
  h2_down_minus_central->Divide(h2_central); 
  h2_down_minus_central->Scale(100.);
  h2_down_minus_central->GetZaxis()->SetRangeUser(-50, 50.);
  if ( zoom ) {
    h2_down_minus_central->GetXaxis()->SetRangeUser(59, 119.);    
    h2_down_minus_central->GetYaxis()->SetRangeUser(0, 70.);
    h2_down_minus_central->GetZaxis()->SetRangeUser(-50, 50.);
    h2_down_minus_central->SetMarkerSize(3);
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 2000, 500);
  c1->Divide(3,1);
  c1->cd(1);
  // h2_central->Draw("colz text"); 
  h2_central->SetXTitle("m_{T} [GeV]");
  h2_central->SetYTitle("m_{ll} [GeV]");
  h2_central->GetXaxis()->CenterTitle();
  h2_central->GetYaxis()->CenterTitle();
  h2_central->Draw("COLZ"); 
  c1->cd(2);
  h2_up_minus_central->SetXTitle("m_{T} [GeV]");
  h2_up_minus_central->SetYTitle("m_{ll} [GeV]");
  h2_up_minus_central->GetXaxis()->CenterTitle();
  h2_up_minus_central->GetYaxis()->CenterTitle();
  if ( zoom ) 
    h2_up_minus_central->Draw("COLZ text"); 
  else
    h2_up_minus_central->Draw("COLZ"); 


  c1->cd(3);
  h2_down_minus_central->SetXTitle("m_{T} [GeV]");
  h2_down_minus_central->SetYTitle("m_{ll} [GeV]");
  h2_down_minus_central->GetXaxis()->CenterTitle();
  h2_down_minus_central->GetYaxis()->CenterTitle();
  if ( zoom )
      h2_down_minus_central->Draw("COLZ text"); 
  else 
    h2_down_minus_central->Draw("COLZ"); 
  /*
  c1->cd(4);
  h2_up->Draw("colz text"); 
  c1->cd(5);
  h2_down->Draw("colz text"); 
  */
  c1->SaveAs(Form("%s/%s_%s_2D_mH%i_%ij_%s.png", outputDir.Data(), process, systvar, mH, njet, flavor));
  c1->SaveAs(Form("%s/%s_%s_2D_mH%i_%ij_%s.eps", outputDir.Data(), process, systvar, mH, njet, flavor));

  // tidy up
  file->Close();
  delete c1;

  text << Form("\\epsfig{figure=%s/%s_%s_mT_mH%i_%ij_%s_lin.eps, width=2.0in}\n", outputDir.Data(), process, systvar, mH, njet, flavor);

}


