#include "../core/Enums.h"
#include <vector>
#include <iostream>
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TROOT.h"

// main function
void doTopEstimation(RunEra runEra = RUN2012)
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
    gROOT->ProcessLine(".L libSmurfTopLooper.so");


    // Option option = WW_OPT_SMURFXSECSEL;
    // Option option = HWW_OPT_SMURFPRESEL;
    Option option = HWW_OPT_SMURFCUTSEL;

    // 
    // initialize the scale factor and errors
    // 
    double TopBkgScaleFactor[3];
    double TopBkgScaleFactorKappa[3];

    for (int i=0; i<3; i++) {
        TopBkgScaleFactor[i] = 0.0;
        TopBkgScaleFactorKappa[i] = 0.0;
    }

    std::string topestFileName = Form("TopBkgScaleFactors_8TeV.h");
    if ( option == HWW_OPT_SMURFPRESEL) 
        topestFileName = "TopBkgScaleFactors_wwpresel_8TeV.h";
    if ( option == WW_OPT_SMURFXSECSEL)
        topestFileName = "TopBkgScaleFactors_wwxsec_8TeV.h";
    if ( option == HWW_OPT_SMURFCUTSEL)
        topestFileName = ("TopBkgScaleFactors_8TeV.h");
    FILE *topsftext = fopen(topestFileName.c_str(), "w"); 



    // 
    // This is for the Top estimation in the VBF channel 
    //

    const int nHiggs = 1;
    Int_t mHiggs[nHiggs] = {0};
    //const int nHiggs = ;
    //Int_t mHiggs[nHiggs] = {0,110,115,120,125,130,135,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
    Double_t TopVBFBkgScaleFactor[nHiggs];
    Double_t TopVBFBkgScaleFactorKappa[nHiggs];

    for ( int i = 0; i < nHiggs; i++) {
        TopVBFBkgScaleFactor[i] = 1.0;
        TopVBFBkgScaleFactorKappa[i] = 1.0;
    }
    const std::string topestvbfFileName = Form("TopVBFBkgScaleFactors_8TeV.h");
    FILE *topsftextvbf = fopen(topestvbfFileName.c_str(), "w"); 

    // 
    // Configuration
    // 

    bool runWW = 1;
    bool runVBF = 0;

    if ( runWW ) {
        doMassPoint(0., option, runEra, TopBkgScaleFactor, TopBkgScaleFactorKappa);
    }

    if ( runVBF) {
        for ( int i = 0; i < nHiggs ; i++) {
            float higgsMass = float(mHiggs[i]);
            doMassPoint(higgsMass, option, runEra, TopBkgScaleFactor, TopBkgScaleFactorKappa);
            TopVBFBkgScaleFactor[i] = TopBkgScaleFactor[2];
            TopVBFBkgScaleFactorKappa[i] = TopBkgScaleFactorKappa[2];
        }
    }

    //
    // write out files
    // 

    // Top scale factor in the WW level
/*
    fputs("Double_t TopBkgScaleFactorKappa(Int_t jetBin) {\n", topsftext);
    fputs("  assert(jetBin >=0 && jetBin <= 2);\n", topsftext);
    fputs(Form("  Double_t TopBkgScaleFactorKappa[3] = { %.5f, %.5f, %.5f   };\n", 
                TopBkgScaleFactorKappa[0], TopBkgScaleFactorKappa[1], TopBkgScaleFactorKappa[2]), topsftext);
    fputs("  return TopBkgScaleFactorKappa[jetBin];\n", topsftext);
    fputs("}\n", topsftext);
*/
    fputs("Double_t TopBkgScaleFactor(Int_t jetBin) {\n", topsftext);
    fputs("  assert(jetBin >=0 && jetBin <= 2);\n", topsftext);
    fputs(Form("  Double_t TopBkgScaleFactor[3] = { %.5f, %.5f, %.5f   };\n", 
                TopBkgScaleFactor[0], TopBkgScaleFactor[1], TopBkgScaleFactor[2]), topsftext);
    fputs("  return TopBkgScaleFactor[jetBin];\n", topsftext);
    fputs("}\n", topsftext);


    fputs("Double_t TopBkgScaleFactorKappa(Int_t jetBin) {\n", topsftext);
    fputs("  assert(jetBin >=0 && jetBin <= 2);\n", topsftext);
    fputs(Form("  Double_t TopBkgScaleFactorKappa[3] = { %.5f, %.5f, %.5f   };\n", 
                TopBkgScaleFactorKappa[0], TopBkgScaleFactorKappa[1], TopBkgScaleFactorKappa[2]), topsftext);
    fputs("  return TopBkgScaleFactorKappa[jetBin];\n", topsftext);
    fputs("}\n", topsftext);


    // Top Scale factor for the VBF analysis

    fputs("static Double_t TopVBFBkgScaleFactor(Int_t mH) {\n", topsftextvbf);
    fputs("  Int_t mHiggs[33] = {0,110,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300,350,400,450,500,550,600,700,800,900,1000};\n", topsftextvbf);
    fputs("  Double_t TopVBFBkgSF[33] = { \n", topsftextvbf);
    for (int i = 0 ; i < nHiggs-1; i++ )
        fputs(Form("%.5f, ",  TopVBFBkgScaleFactor[i]), topsftextvbf );
    fputs(Form("%.5f};\n",  TopVBFBkgScaleFactor[nHiggs-1]), topsftextvbf );
    fputs("  Int_t massIndex = -1;\n", topsftextvbf);
    fputs("  for (UInt_t m=0; m < 33 ; ++m) {\n", topsftextvbf);
    fputs("    if (mH == mHiggs[m]) massIndex = m;\n", topsftextvbf);
    fputs("  }\n", topsftextvbf);
    fputs("  if (massIndex >= 0) {\n", topsftextvbf);
    fputs("    return TopVBFBkgSF[massIndex];\n", topsftextvbf);
    fputs("  } else { assert(0); }\n", topsftextvbf);
    fputs("  return 1.0;\n", topsftextvbf);
    fputs("}\n", topsftextvbf);


    fputs("static Double_t TopVBFBkgScaleFactorKappa(Int_t mH) {\n", topsftextvbf);
    fputs("  Int_t mHiggs[33] = {0,110,115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200,250,300,350,400,450,500,550,600,700,800,900,1000};\n", topsftextvbf);
    fputs("  Double_t TopVBFBkgKappa[33] = { \n", topsftextvbf);
    for (int i = 0 ; i < nHiggs-1; i++ )
        fputs(Form("%.5f, ",  TopVBFBkgScaleFactorKappa[i]), topsftextvbf );
    fputs(Form("%.5f};\n",  TopVBFBkgScaleFactorKappa[nHiggs-1]), topsftextvbf );
    fputs("  Int_t massIndex = -1;\n", topsftextvbf);
    fputs("  for (UInt_t m=0; m < 33 ; ++m) {\n", topsftextvbf);
    fputs("    if (mH == mHiggs[m]) massIndex = m;\n", topsftextvbf);
    fputs("  }\n", topsftextvbf);
    fputs("  if (massIndex >= 0) {\n", topsftextvbf);
    fputs("    return TopVBFBkgKappa[massIndex];\n", topsftextvbf);
    fputs("  } else { assert(0); }\n", topsftextvbf);
    fputs("  return 1.0;\n", topsftextvbf);
    fputs("}\n", topsftextvbf);



}

void doMassPoint(const float analysis, Option option, RunEra runEra, double TopBkgScaleFactor[3], double TopBkgScaleFactorKappa[3])
{  


    //
    // Loop files and get the histograms
    //

    SmurfSample *sample_data    = new SmurfSample(option, DATA,      kBlack,    "Data",   analysis);
    SmurfSample *sample_ttbar   = new SmurfSample(option, TT,        kMagenta,  "TT",     analysis);
    SmurfSample *sample_tw      = new SmurfSample(option, TW,        kMagenta,  "TW",     analysis);

    // other non-top backgrounds
    SmurfSample *sample_qqww    = new SmurfSample(option, QQWW,      kYellow+2, "qqWW",   analysis);
    SmurfSample *sample_ggww    = new SmurfSample(option, GGWW,      kYellow+2, "ggWW",   analysis);
    SmurfSample *sample_wz      = new SmurfSample(option, WZ,        kGreen,    "WZ",     analysis);
    SmurfSample *sample_zz      = new SmurfSample(option, ZZ,        kGreen,    "ZZ",     analysis);
    SmurfSample *sample_dyll    = new SmurfSample(option, ZLL,       kBlue,     "DYLL",   analysis);
    SmurfSample *sample_wjets   = new SmurfSample(option, WJETSDATA, kCyan,     "Wjets",  analysis);
    SmurfSample *sample_wgamma  = new SmurfSample(option, WGAMMA,    kCyan+2,   "Wgamma", analysis);
    SmurfSample *sample_ztt     = new SmurfSample(option, ZTT,       kBlue+2,   "Ztt",    analysis);

    std::vector<SmurfSample*> samples;
    samples.push_back(sample_data);
    samples.push_back(sample_ttbar);
    samples.push_back(sample_tw);
    samples.push_back(sample_qqww);
    samples.push_back(sample_ggww);
    samples.push_back(sample_wz);
    samples.push_back(sample_zz);
    samples.push_back(sample_dyll);
    samples.push_back(sample_wjets);
    samples.push_back(sample_wgamma);
    samples.push_back(sample_ztt);

    //
    // looper 
    // 

    float lumi = 19467;
    SmurfTopLooper *looper = new SmurfTopLooper(analysis, option);
    looper->setLumiScale(lumi, lumi);
    std::string outFile = Form("tophistos_ww_presel.root");
    if ( option == WW_OPT_SMURFXSECSEL ) 
        outFile = Form("tophistos_ww_wwxsec.root");
    if ( option == HWW_OPT_SMURFCUTSEL)
        outFile = Form("tophistos_ww_%.0f.root", analysis);

    bool skimData = true;

    if ( !skimData) {
        char *dataDir = "/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets/";
        sample_data->add(Form("%s/data.root", dataDir));
        sample_ttbar->add(Form("%s/ttbar_powheg.root", dataDir));
        sample_tw->add(Form("%s/tw.root", dataDir));
        sample_qqww->add(Form("%s/qqww.root", dataDir));
        sample_ggww->add(Form("%s/ggww.root", dataDir));
        sample_dyll->add(Form("%s/dyll.root", dataDir));
        sample_wz->add(Form("%s/wz.root", dataDir));
        sample_zz->add(Form("%s/zz.root", dataDir));
        sample_wgamma->add(Form("%s/zgamma.root", dataDir));
        sample_wgamma->add(Form("%s/wgamma.root", dataDir));
        sample_wgamma->add(Form("%s/wglll.root", dataDir));
        sample_wgamma->add(Form("%s/www.root", dataDir));
        sample_wjets->add(Form("%s/data.root", dataDir));
        sample_wjets->add(Form("%s/qqww.root", dataDir));
        sample_wjets->add(Form("%s/ggww.root", dataDir));
        sample_wjets->add(Form("%s/ttbar_powheg.root", dataDir));
        sample_wjets->add(Form("%s/tw.root", dataDir));
        sample_wjets->add(Form("%s/wz.root", dataDir));
        sample_wjets->add(Form("%s/zz.root", dataDir));
    }

    // Yanyan's skims
    if ( skimData ) {

        char *dataDir = "/smurf/jaehyeok/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets_19p5ifb_new/WW/"; // TAS
        //char *dataDir  = "/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/WW/"; // UAF

        for ( int njet = 0; njet < 3; njet++) {
            sample_data->add(Form("%s/%ij/data.root", dataDir, njet));
            sample_ttbar->add(Form("%s/%ij/ttbar_powheg.root", dataDir, njet));
            sample_tw->add(Form("%s/%ij/tw.root", dataDir, njet));
            sample_qqww->add(Form("%s/%ij/qqww.root", dataDir, njet));
            sample_ggww->add(Form("%s/%ij/ggww.root", dataDir, njet));
            sample_wz->add(Form("%s/%ij/wz.root", dataDir, njet));
            sample_zz->add(Form("%s/%ij/zz.root", dataDir, njet));
            sample_dyll->add(Form("%s/%ij/dyll.root", dataDir, njet));
            sample_wgamma->add(Form("%s/%ij/wgamma.root", dataDir, njet));
            sample_wgamma->add(Form("%s/%ij/zgamma.root", dataDir, njet));
            sample_wgamma->add(Form("%s/%ij/wglll.root", dataDir, njet));
            sample_wgamma->add(Form("%s/%ij/www.root", dataDir, njet));
            sample_wjets->add(Form("%s/%ij/data_PassFail.root", dataDir, njet));
            sample_wjets->add(Form("%s/%ij/qqww_PassFail.root", dataDir, njet));
            sample_wjets->add(Form("%s/%ij/ggww_PassFail.root", dataDir, njet));
            sample_wjets->add(Form("%s/%ij/ttbar_powheg_PassFail.root", dataDir, njet));
            sample_wjets->add(Form("%s/%ij/tw_PassFail.root", dataDir, njet));
            sample_wjets->add(Form("%s/%ij/wz_PassFail.root", dataDir, njet));
            sample_wjets->add(Form("%s/%ij/zz_PassFail.root", dataDir, njet));
            sample_wjets->add(Form("%s/%ij/wgamma_PassFail.root", dataDir, njet));
            sample_wjets->add(Form("%s/%ij/wglll_PassFail.root", dataDir, njet));
            sample_wjets->add(Form("%s/%ij/dyll_PassFail.root", dataDir, njet));
            sample_wjets->add(Form("%s/%ij/www_PassFail.root", dataDir, njet));
            sample_wjets->add(Form("%s/%ij/zgamma_PassFail.root", dataDir, njet));
            
            sample_ztt->add(Form("%s/%ij/data_ztt.root", dataDir, njet));
        }
    }

    looper->loop(sample_data);
    looper->loop(sample_ttbar);
    looper->loop(sample_tw);
    looper->loop(sample_qqww);
    looper->loop(sample_ggww);
    looper->loop(sample_dyll);
    looper->loop(sample_wz);
    looper->loop(sample_zz);
    looper->loop(sample_ztt);
    looper->loop(sample_wgamma);
    looper->loop(sample_wjets);

    // save the histograms for later analysis
    saveHist(outFile.c_str());
    deleteHistos();

    // 
    // Do the top estimation
    // 

    std::string debugFileName = Form("topest_%.0fpb.txt", lumi);
    if ( option == HWW_OPT_SMURFPRESEL ) 
        debugFileName = Form("topest_%.0fpb_wwpresel.txt", lumi);
    if ( option == WW_OPT_SMURFXSECSEL ) 
        debugFileName = Form("topest_%.0fpb_wwxsec.txt", lumi);
    if ( option == HWW_OPT_SMURFCUTSEL ) 
        debugFileName = Form("topest_%.0fpb_mH%.0f.txt", lumi, analysis);

    FILE *debugtext = fopen(debugFileName.c_str(), "w"); 
    topest(analysis, outFile.c_str(), debugtext, TopBkgScaleFactor, TopBkgScaleFactorKappa);
    fclose(debugtext);

    //
    // tidy up
    // 

    delete sample_data;
    delete sample_ttbar;
    delete sample_tw;
    delete sample_qqww;
    delete sample_ggww;
    delete sample_dyll;
    delete sample_wz;
    delete sample_zz;
    delete sample_wgamma;
    delete sample_wjets;
    delete sample_ztt;
}
