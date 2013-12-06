//
// This code does the Drell-Yan Estimations
// It produces the DYBkgScaleFactor_8TeV.h
// 

#include "../core/Enums.h"

#include "TROOT.h"

#include <vector>
#include <iostream>

void doDYEstimation(RunEra runEra = RUN2012)
{

    //
    // load the libraries
    //

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gROOT->ProcessLine(".L ../../../../../Smurf/Core/SmurfTree.h+");
    gROOT->ProcessLine(".L ../../../../../Smurf/Core/LeptonScaleLookup.cc+");
    gROOT->ProcessLine(".L ../../../../NtupleMacros/Tools/goodrun.cc+");
    gROOT->ProcessLine(".L libSmurfDYLooper.so");
    
    // 
    // define the  complete mass points considered and dy background file
    // 

    int mHiggs[26] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
    double DYBkgScaleFactorWWPreselection[3];
    double DYBkgScaleFactorWWPreselectionKappa[3]; // this is 1 + relative systematic errors
    double DYBkgScaleFactorHiggsSelection[3][26]; 
    double DYBkgScaleFactorHiggsSelectionKappa[3][26]; // this is 1 + relative systematic errors
    double DYBkgScaleFactorHiggsSelectionMVA[3][26]; 
    double DYBkgScaleFactorHiggsSelectionKappaMVA[3][26]; // this is 1 + relative systematic errors

	// initialize the values 
	for (int njet = 0; njet < 3; njet++) {
		DYBkgScaleFactorWWPreselection[njet] = 1.0;
		DYBkgScaleFactorWWPreselectionKappa[njet] = 1.0;
		for ( int mH = 0; mH < 26; mH++) {
			DYBkgScaleFactorHiggsSelection[njet][mH] = 1.0;
			DYBkgScaleFactorHiggsSelectionKappa[njet][mH] = 1.0;
			DYBkgScaleFactorHiggsSelectionMVA[njet][mH] = 1.0;
			DYBkgScaleFactorHiggsSelectionKappaMVA[njet][mH] = 1.0;
		}
	}
    

    //
    // configuration  
    //
    bool runWWXSec = 0; 
    bool runWWPresel = 1;
    bool run115 = 1;
    bool run120 = 1;
    bool run125 = 1;
    bool run130 = 1; 
    bool run135 = 1; 
    bool run140 = 1; 
    bool run145 = 1; 
    bool run150 = 1; 
    bool run160 = 1; 
    bool run170 = 1; 
    bool run180 = 1; 
    bool run190 = 1; 
    bool run200 = 1; 
    bool run250 = 1; 
    bool run300 = 1;
    bool run350 = 0;
    bool run400 = 0;
    bool run450 = 0;
    bool run500 = 0;
    bool run550 = 0;
    bool run600 = 0;
    
    //
    // run Drell-Yan estimation for masses considered above
    //


    if ( runWWXSec )
      doMassPoint(0, WW_OPT_SMURFXSECSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );

    if ( runWWPresel )
      doMassPoint(0, HWW_OPT_SMURFPRESEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    
    if ( run115) {
      doMassPoint(115.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(115.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run120) {
      doMassPoint(120.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(120.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }

    if ( run125) {
      doMassPoint(125.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(125.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run130) {
      doMassPoint(130.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
      
//      doMassPoint(130.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }

    if ( run135) {
      doMassPoint(135.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(135.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run140) {
      doMassPoint(140.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(140.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run145) {
      doMassPoint(145.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(145.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run150) {
      doMassPoint(150.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(150.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }

    if ( run160) {
      doMassPoint(160.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(160.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run170) {
      doMassPoint(170.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(170.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run180) {
      doMassPoint(180.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(180.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run190) {
      doMassPoint(190.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(190.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run200) {
      doMassPoint(200.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(200.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run250) {
      doMassPoint(250.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
//      doMassPoint(250.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run300) {
      doMassPoint(300.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );

//      doMassPoint(300.0, HWW_OPT_SMURFMVASEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
//		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
//		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }

    if ( run350) {
      doMassPoint(350.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run400) {
      doMassPoint(400.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run450) {
      doMassPoint(450.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run500) {
      doMassPoint(500.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }

    if ( run550) {
      doMassPoint(550.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }
    if ( run600) {
      doMassPoint(600.0, HWW_OPT_SMURFCUTSEL, runEra, mHiggs, DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
		  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
		  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    }


    //
    // write the data/MC scale factors into the code
    //

    const std::string dyestFileName = "DYBkgScaleFactors_8TeV.h";
    FILE *dysftext = fopen(dyestFileName.c_str(), "w"); 

    fputs("static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {\n", dysftext);
    fputs("  Int_t mHiggs[20] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190,200,250,300};\n",dysftext);
    fputs(Form("  Double_t DYBkgScaleFactorWWPreselection[3] = { %.5f, %.5f, %.5f};\n", 
	       DYBkgScaleFactorWWPreselection[0], DYBkgScaleFactorWWPreselection[1], DYBkgScaleFactorWWPreselection[2]), dysftext);
    fputs("  Double_t DYBkgScaleFactorHiggsSelection[3][20] = { \n", dysftext);
    for ( int njet = 0; njet < 2; njet++) {
      fputs("    {", dysftext);
      for (int i = 0; i < 19; i++) 
	fputs(Form("%.5f,",  DYBkgScaleFactorHiggsSelection[njet][i]), dysftext);
      fputs(Form("%.5f},\n",  DYBkgScaleFactorHiggsSelection[njet][19]), dysftext);
    }
    fputs("    {", dysftext);
    for (int i = 0; i < 19; i++) 
      // fputs(Form("%.5f,",  DYBkgScaleFactorWWPreselection[2]), dysftext);
      fputs(Form("%.5f,",  DYBkgScaleFactorHiggsSelection[2][i]), dysftext);
    fputs(Form("%.5f} };\n",  DYBkgScaleFactorHiggsSelection[2][19]), dysftext);
    
    fputs("  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];\n", dysftext);
    fputs("  Int_t massIndex = -1;\n", dysftext);
    fputs("  for (UInt_t m=0; m < 20 ; ++m) {\n", dysftext);
    fputs("    if (mH == mHiggs[m]) massIndex = m;\n", dysftext);
    fputs("  }\n", dysftext);
    fputs("  if (massIndex >= 0) {\n", dysftext);
    fputs("    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];\n", dysftext);
    fputs("  } else {\n", dysftext);
    fputs("    return DYBkgScaleFactorWWPreselection[jetBin];\n", dysftext);
    fputs("  }\n", dysftext);
    fputs("}\n\n", dysftext);



    fputs("Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {\n", dysftext);
    fputs("  Int_t mHiggs[20] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190,200,250,300};\n",dysftext);
    fputs(Form("  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { %.5f, %.5f, %.5f};\n", 
	       DYBkgScaleFactorWWPreselectionKappa[0], DYBkgScaleFactorWWPreselectionKappa[1], DYBkgScaleFactorWWPreselectionKappa[2]), dysftext);    
    fputs("  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][20] = { \n", dysftext);   
    for ( int njet = 0; njet < 2; njet++) {
      fputs("    {", dysftext);
      for (int i = 0; i < 19; i++) 
	fputs(Form("%.5f,",  DYBkgScaleFactorHiggsSelectionKappa[njet][i]), dysftext);
      fputs(Form("%.5f},\n",  DYBkgScaleFactorHiggsSelectionKappa[njet][19]), dysftext);
    }
    fputs("    {", dysftext);
    for (int i = 0; i < 19; i++) 
      fputs(Form("%.5f,",  DYBkgScaleFactorHiggsSelectionKappa[2][i]), dysftext);
    fputs(Form("%.5f} };\n",  DYBkgScaleFactorHiggsSelectionKappa[2][19]), dysftext);
    fputs("  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];\n", dysftext);
    fputs("  Int_t massIndex = -1;\n", dysftext);
    fputs("  for (UInt_t m=0; m < 20 ; ++m) {\n", dysftext);
    fputs("    if (mH == mHiggs[m]) massIndex = m;\n", dysftext);
    fputs("  }\n", dysftext);
    fputs("  if (massIndex >= 0) {\n", dysftext);
    fputs("    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];\n", dysftext);
    fputs("  } else {\n", dysftext);
    fputs("    return DYBkgScaleFactorWWPreselectionKappa[jetBin];\n", dysftext);
    fputs("  }\n", dysftext);
    fputs("}\n", dysftext);

    
    fputs("static Double_t DYBkgScaleFactorBDT(Int_t mH, Int_t jetBin) {\n", dysftext);
    fputs("  Int_t mHiggs[20] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190,200,250,300};\n",dysftext);
    fputs(Form("  Double_t DYBkgScaleFactorWWPreselection[3] = { %.5f, %.5f, %.5f};\n", 
	       DYBkgScaleFactorWWPreselection[0], DYBkgScaleFactorWWPreselection[1], DYBkgScaleFactorWWPreselection[2]), dysftext);
    fputs("  Double_t DYBkgScaleFactorHiggsSelection[3][20] = { \n", dysftext);
    for ( int njet = 0; njet < 2; njet++) {
      fputs("    {", dysftext);
      for (int i = 0; i < 19; i++) 
	fputs(Form("%.5f,",  DYBkgScaleFactorHiggsSelectionMVA[njet][i]), dysftext);
      fputs(Form("%.5f},\n",  DYBkgScaleFactorHiggsSelectionMVA[njet][19]), dysftext);
    }
    fputs("    {", dysftext);
    for (int i = 0; i < 19; i++) 
      fputs(Form("%.5f,",  DYBkgScaleFactorWWPreselection[2]), dysftext);
    fputs(Form("%.5f} };\n",  DYBkgScaleFactorWWPreselection[2]), dysftext);
    
    fputs("  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];\n", dysftext);
    fputs("  Int_t massIndex = -1;\n", dysftext);
    fputs("  for (UInt_t m=0; m < 20 ; ++m) {\n", dysftext);
    fputs("    if (mH == mHiggs[m]) massIndex = m;\n", dysftext);
    fputs("  }\n", dysftext);
    fputs("  if (massIndex >= 0) {\n", dysftext);
    fputs("    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];\n", dysftext);
    fputs("  } else {\n", dysftext);
    fputs("    return DYBkgScaleFactorWWPreselection[jetBin];\n", dysftext);
    fputs("  }\n", dysftext);
    fputs("}\n\n", dysftext);




    fputs("Double_t DYBkgScaleFactorBDTKappa(Int_t mH, Int_t jetBin) {\n", dysftext);
    fputs("  Int_t mHiggs[20] = {115,118,120,122,124,125,126,128,130,135,140,145,150,160,170,180,190,200,250,300};\n",dysftext);
    fputs(Form("  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { %.5f, %.5f, %.5f};\n", 
	       DYBkgScaleFactorWWPreselectionKappa[0], DYBkgScaleFactorWWPreselectionKappa[1], DYBkgScaleFactorWWPreselectionKappa[2]), dysftext);    
    fputs("  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][20] = { \n", dysftext);   
    for ( int njet = 0; njet < 2; njet++) {
      fputs("    {", dysftext);
      for (int i = 0; i < 19; i++) 
	fputs(Form("%.5f,",  DYBkgScaleFactorHiggsSelectionKappaMVA[njet][i]), dysftext);
      fputs(Form("%.5f},\n",  DYBkgScaleFactorHiggsSelectionKappaMVA[njet][19]), dysftext);
    }
    fputs("    {", dysftext);
    for (int i = 0; i < 19; i++) 
      fputs(Form("%.5f,",  DYBkgScaleFactorWWPreselectionKappa[2]), dysftext);
    fputs(Form("%.5f} };\n",  DYBkgScaleFactorWWPreselectionKappa[2]), dysftext);
    fputs("  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];\n", dysftext);
    fputs("  Int_t massIndex = -1;\n", dysftext);
    fputs("  for (UInt_t m=0; m < 20 ; ++m) {\n", dysftext);
    fputs("    if (mH == mHiggs[m]) massIndex = m;\n", dysftext);
    fputs("  }\n", dysftext);
    fputs("  if (massIndex >= 0) {\n", dysftext);
    fputs("    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];\n", dysftext);
    fputs("  } else {\n", dysftext);
    fputs("    return DYBkgScaleFactorWWPreselectionKappa[jetBin];\n", dysftext);
    fputs("  }\n", dysftext);
    fputs("}\n", dysftext);

    fclose(dysftext);    



    

}


void doMassPoint(float analysis, Option option, RunEra runEra, int mHiggs[20],  
		 double DYBkgScaleFactorWWPreselection[3], double DYBkgScaleFactorWWPreselectionKappa[3], 
		 double DYBkgScaleFactorHiggsSelection[3][20], double DYBkgScaleFactorHiggsSelectionKappa[3][20],
		 double DYBkgScaleFactorHiggsSelectionMVA[3][20], double DYBkgScaleFactorHiggsSelectionKappaMVA[3][20])
{

    float lumi = 19467; // in unit of 1/pb
      

    //
    // set up the looper
    //

    SmurfDYLooper *looper = new SmurfDYLooper(analysis, option, runEra);
    //looper->setGoodRunList("../runlists/runlist_1092.txt");
    looper->setLumiScale(lumi, lumi);

    //
    // set up samples
    //
    int ana = (int)analysis;

    SmurfSample *sample_data = new SmurfSample(option, DATA, kBlack, "Data");
    SmurfSample *sample_dyll = new SmurfSample(option, ZLL, kBlue, "DYLL");
    SmurfSample *sample_vv = new SmurfSample(option, VV, kGreen, "VV");

    // examples of using the unskimed files 
    char *dataDir = "/smurf/jaehyeok/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets_19p5ifb_new/dyskim/";
    sample_data->add(Form("%s/data.root", dataDir));    
    sample_dyll->add(Form("%s/dyll.root", dataDir));
    sample_vv->add(Form("%s/wz.root", dataDir));
    sample_vv->add(Form("%s/zz.root", dataDir));
    
	// 
    // do the looping
    //

    looper->loop(sample_data);
    looper->loop(sample_dyll);
    looper->loop(sample_vv);

    //
    // save histograms  
    //
    std::string outFile = Form("dyhistos_ww_%i.root", int(analysis));
    if ( option == WW_OPT_SMURFXSECSEL ) {
      outFile = Form("dyhistos_ww_%i_wwxsec.root", int(analysis));
    }
    if ( option == HWW_OPT_SMURFCUTSEL ) {
      outFile = Form("dyhistos_ww_%i_cut.root", int(analysis));
    }
    if ( option ==  HWW_OPT_SMURFMVASEL ) {
      outFile = Form("dyhistos_ww_%i_mva.root", int(analysis));
      std::cout << "outputFile Name = " << outFile.c_str() << "\n";
    }
    saveHist(outFile.c_str());  
    deleteHistos();

    
    std::string debugFileName = Form("dyest_mH%i_%.0fpb.txt", int(analysis), lumi);
    if ( option == WW_OPT_SMURFXSECSEL) 
      debugFileName = Form("dyest_mH%i_%.0fpb_wwxsec.txt", int(analysis), lumi);
    if ( option == HWW_OPT_SMURFCUTSEL) 
      debugFileName = Form("dyest_mH%i_%.0fpb_cut.txt", int(analysis), lumi);
    if ( option == HWW_OPT_SMURFMVASEL) 
      debugFileName = Form("dyest_mH%i_%.0fpb_mva.txt", int(analysis), lumi);
    
    FILE *debugtext = fopen(debugFileName.c_str(), "w"); 
    
    // 
    // fill the R 
    // 
    std::string ratioFileName = Form("Routin_mH%i_%.0fpb.root", int(analysis), lumi);
    if ( option == WW_OPT_SMURFXSECSEL) {
      ratioFileName = Form("Routin_mH%i_%.0fpb_wwxsec.root", int(analysis), lumi);
    }
    if ( option == HWW_OPT_SMURFCUTSEL) {
      ratioFileName = Form("Routin_mH%i_%.0fpb_cut.root", int(analysis), lumi);
    }
    if ( option == HWW_OPT_SMURFMVASEL) {
      ratioFileName = Form("Routin_mH%i_%.0fpb_mva.root", int(analysis), lumi);
    }

    fillRoutin(outFile.c_str(), ratioFileName.c_str(), debugtext);
    dyest(analysis, option, outFile.c_str(), ratioFileName.c_str(), debugtext, mHiggs, 
	  DYBkgScaleFactorWWPreselection, DYBkgScaleFactorWWPreselectionKappa, 
	  DYBkgScaleFactorHiggsSelection, DYBkgScaleFactorHiggsSelectionKappa,
	  DYBkgScaleFactorHiggsSelectionMVA, DYBkgScaleFactorHiggsSelectionKappaMVA );
    
    fclose(debugtext);
    
    //
    // clean up
    //

    delete looper;
    delete sample_data;
    delete sample_dyll;
    delete sample_vv;
    
}


