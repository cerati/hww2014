#include <cmath>
#include <set>
#include <iostream>
#include <stdio.h>
#include "TH1F.h"
#include "TH2F.h"

//
// This code replaces the histo_Data in the card by the 
// S+B expectation and changes the corresponding card text
// 

double replaceinputroot(TString inputFileName);
void  replacecard(TString cardName, double yield);


void replacecarddata(TString carddir, int mH, int njet) {

  // replace hww cards
  const char *hwwinputFileName = Form("%s/%i/hwwof_%ij.input_8TeV.root", carddir.Data(), mH, njet); 
  const char *hwwcardFileName = Form("%s/%i/hwwof_%ij_shape_8TeV.txt", carddir.Data(), mH, njet); 
  double yield = replaceinputroot(hwwinputFileName);
  replacecard(hwwcardFileName, yield);

  // replace xwwcard
  const char *xwwinputFileName = Form("%s/%i/xwwof_%ij.input_8TeV.root", carddir.Data(), mH, njet); 
  const char *xwwcardFileName = Form("%s/%i/xwwof_%ij_shape_8TeV.txt", carddir.Data(), mH, njet); 
  yield = replaceinputroot(xwwinputFileName);
  replacecard(xwwcardFileName, yield);
}


double replaceinputroot(TString inputFileName) {

  std::cout << "Replacing " << inputFileName << "\n";

  TFile *inputFile = new TFile(inputFileName, "UPDATE");
  gROOT->cd();
  inputFile->cd();
  
  // replace the histo_Data
  TH1F *h1_Data = (TH1F*)inputFile->Get("histo_ggH")->Clone("histo_Data");
  cout << "ggH : " << h1_Data->Integral() <<  endl;   
  
  TH1F *h1_qqH = (TH1F*)inputFile->Get("histo_qqH");
  h1_Data->Add(h1_qqH);

  TH1F *h1_WH = (TH1F*)inputFile->Get("histo_WH");
  h1_Data->Add(h1_WH);

  TH1F *h1_ZH = (TH1F*)inputFile->Get("histo_ZH");
  h1_Data->Add(h1_ZH);
  
  TH1F *h1_qqWW = (TH1F*)inputFile->Get("histo_qqWW");
  h1_Data->Add(h1_qqWW);
  
  TH1F *h1_ggWW = (TH1F*)inputFile->Get("histo_ggWW");
  h1_Data->Add(h1_ggWW);

  TH1F *h1_VV = (TH1F*)inputFile->Get("histo_VV");
  h1_Data->Add(h1_VV);
  
  TH1F *h1_Top = (TH1F*)inputFile->Get("histo_Top");
  h1_Data->Add(h1_Top);
  
  TH1F *h1_Zjets = (TH1F*)inputFile->Get("histo_Zjets");
  h1_Data->Add(h1_Zjets);
  
  TH1F *h1_WjetsE = (TH1F*)inputFile->Get("histo_WjetsE");
  h1_Data->Add(h1_WjetsE);

  TH1F *h1_WjetsM = (TH1F*)inputFile->Get("histo_WjetsM");
  h1_Data->Add(h1_WjetsM);

  TH1F *h1_Wgamma = (TH1F*)inputFile->Get("histo_Wgamma");
  h1_Data->Add(h1_Wgamma);

  TH1F *h1_Wg3l = (TH1F*)inputFile->Get("histo_Wg3l");
  h1_Data->Add(h1_Wg3l);

  TH1F *h1_Ztt = (TH1F*)inputFile->Get("histo_Ztt");
  h1_Data->Add(h1_Ztt);
  

   cout << "qqH : " << h1_qqH->Integral() <<  endl;   
   cout << "WH  : " << h1_WH->Integral() <<  endl;   
   cout << "ZH  : " << h1_ZH->Integral() <<  endl;   
   cout << "qqWW  : " << h1_qqWW->Integral() <<  endl;   
   cout << "ggWW  : " << h1_ggWW->Integral() <<  endl;   
   cout << "VV    : " << h1_VV->Integral() <<  endl;   
   cout << "Top   : " << h1_Top->Integral() <<  endl;   
   cout << "Zjets : " << h1_Zjets->Integral() <<  endl;   
   cout << "WjetsE : " << h1_WjetsE->Integral() <<  endl;   
   cout << "WjetsM : " << h1_WjetsM->Integral() <<  endl;   
   cout << "Wgamma: " << h1_Wgamma->Integral() <<  endl;   
   cout << "Wg3l  : " << h1_Wg3l->Integral() <<  endl;   
   cout << "Ztt   : " << h1_Ztt->Integral() <<  endl;   

  double yield = histo_Data->Integral(0, 1000);
  std::cout << "yield = " << yield << "\n";
  

  h1_Data->SetDirectory(0);
  h1_Data->Write();
  
  inputFile->Close();
  return yield;
}


void replacecard(TString cardFileName, double yield) {

  ifstream inputfile(cardFileName);
  TString tempfileName = cardFileName;
  tempfileName.ReplaceAll(".txt", "_temp.txt");
  ofstream tempfile;
  tempfile.open(tempfileName);
  
  string line;
  if (inputfile.is_open())
    {
      while ( inputfile.good() )
	{
	  getline (inputfile,line);
	  if (TString(line).Contains("Observation", TString::kExact)) 
	    tempfile << Form("Observation %.0f\n", yield);
	  else 
	    tempfile << line << endl;
	}
      inputfile.close();
      tempfile.close();
    } else {
    std::cout << "Unable to open " << cardFileName << "\n";
  }
}
