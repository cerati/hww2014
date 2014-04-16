
#include "../core/Enums.h"

#include "TROOT.h"

#include <vector>
#include <iostream>

void doWWEstimation(RunEra runEra = RUN2011AB){

    //
    // load the libraries
    //

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gROOT->ProcessLine(".L ../../Smurf/Core/SmurfTree.h+");
    gROOT->ProcessLine(".L ../../Smurf/Core/LeptonScaleLookup.cc+");
    gROOT->ProcessLine(".L ../../Tools/goodrun.cc+");
    gROOT->ProcessLine(".L libSmurfWWLooper.so");

    // 
    // define the  complete mass points considered and dy background file
    // 

    const int nHiggs = 19;
    int mHiggs[nHiggs] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};

    // cut-based scale factor and errors
    double WWBkgScaleFactorCutBased[2][nHiggs]; 
    double WWBkgScaleFactorKappaCutBased[2][nHiggs]; // this is 1 + relative systematic errors
    // initialize the values 
    for (int njet = 0; njet < 2; njet++) {
        for ( int mH = 0; mH < 15; mH++) {
            WWBkgScaleFactorCutBased[njet][mH] = 1.0;
            WWBkgScaleFactorKappaCutBased[njet][mH] = 1.0;
        }
    }

    // cut-based scale factor and errors
    double WWBkgScaleFactorMVA[2][nHiggs]; 
    double WWBkgScaleFactorKappaMVA[2][nHiggs]; // this is 1 + relative systematic errors
    // initialize the values 
    for (int njet = 0; njet < 2; njet++) {
        for ( int mH = 0; mH < nHiggs; mH++) {
            WWBkgScaleFactorMVA[njet][mH] = 1.0;
            WWBkgScaleFactorKappaMVA[njet][mH] = 1.0;
        }
    }


    // 
    // run WW estimation for masses considered above
    //
 
    for ( int i = 0; i < nHiggs ; i++) {
        doMassPoint(mHiggs[i], HWW_OPT_SMURFCUTSEL, runEra, mHiggs, WWBkgScaleFactorCutBased, WWBkgScaleFactorKappaCutBased);
        doMassPoint(mHiggs[i], HWW_OPT_MT2DMLL, runEra, mHiggs, WWBkgScaleFactorMVA, WWBkgScaleFactorKappaMVA);
    }


    //
    // write the data/MC scale factors into the code
    // 

    const std::string wwestFileName = "WWBkgScaleFactors_8TeV.h";
    FILE *wwsftext = fopen(wwestFileName.c_str(), "w"); 

    // scale factors

    fputs("static Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {\n", wwsftext);
    fputs("assert(jetBin >= 0 && jetBin <= 1);\n", wwsftext);
    fputs("  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};\n",wwsftext);
    fputs("  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { \n", wwsftext);
    for ( int njet = 0; njet < 2; njet++) {
        fputs("    {", wwsftext);
        for (int i = 0; i < nHiggs - 1 ; i++) 
            fputs(Form("%.5f,",  WWBkgScaleFactorCutBased[njet][i]), wwsftext);
        if ( njet == 0 )
            fputs(Form("%.5f},\n",  WWBkgScaleFactorCutBased[njet][nHiggs-1]), wwsftext);
        if ( njet == 1 )
            fputs(Form("%.5f} };\n",  WWBkgScaleFactorCutBased[njet][nHiggs-1]), wwsftext);
    }

    fputs("  Int_t massIndex = -1;\n", wwsftext);
    fputs("  for (UInt_t m=0; m < 19 ; ++m) {\n", wwsftext);
    fputs("    if (mH == mHiggs[m]) massIndex = m;\n", wwsftext);
    fputs("  }\n", wwsftext);
    fputs("  if (massIndex >= 0) {\n", wwsftext);
    fputs("    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];\n", wwsftext);
    fputs("  } else {\n", wwsftext);
    fputs("    return 1.0;\n", wwsftext);
    fputs("  }\n", wwsftext);
    fputs("}\n\n", wwsftext);



    fputs("static Double_t WWBkgScaleFactorMVA(Int_t mH, Int_t jetBin) {\n", wwsftext);
    fputs("assert(jetBin >= 0 && jetBin <= 1);\n", wwsftext);
    fputs("  Int_t mHiggs[19] = {115,118,120,122,125,124,126,128,130,135,140,145,150,155,160,170,180,190,200};\n",wwsftext);
    fputs("  Double_t WWBkgScaleFactorHiggsSelection[2][19] = { \n", wwsftext);
    for ( int njet = 0; njet < 2; njet++) {
        fputs("    {", wwsftext);
        for (int i = 0; i < nHiggs-1; i++) 
            fputs(Form("%.5f,",  WWBkgScaleFactorMVA[njet][i]), wwsftext);
            //fputs(Form("%.5f,",  WWBkgScaleFactorMVA[njet][0]), wwsftext);
        if ( njet == 0 )
            fputs(Form("%.5f},\n",  WWBkgScaleFactorMVA[njet][nHiggs-1]), wwsftext);
            //fputs(Form("%.5f},\n",  WWBkgScaleFactorMVA[njet][0]), wwsftext);
        if ( njet == 1 )
            fputs(Form("%.5f} };\n",  WWBkgScaleFactorMVA[njet][nHiggs-1]), wwsftext);
            //fputs(Form("%.5f} };\n",  WWBkgScaleFactorMVA[njet][0]), wwsftext);
    }

    fputs("  Int_t massIndex = -1;\n", wwsftext);
    fputs("  for (UInt_t m=0; m < 19 ; ++m) {\n", wwsftext);
    fputs("    if (mH == mHiggs[m]) massIndex = m;\n", wwsftext);
    fputs("  }\n", wwsftext);
    fputs("  if (massIndex >= 0) {\n", wwsftext);
    fputs("    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];\n", wwsftext);
    fputs("  } else {\n", wwsftext);
    fputs("    return 1.0;\n", wwsftext);
    fputs("  }\n", wwsftext);
    fputs("}\n\n", wwsftext);

    // scale factor errors

    fputs("static Double_t WWBkgScaleFactorKappaCutBased(Int_t mH, Int_t jetBin) {\n", wwsftext);
    fputs("assert(jetBin >= 0 && jetBin <= 1);\n", wwsftext);
    fputs("  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};\n",wwsftext);
    fputs("  Double_t WWBkgScaleFactorKappaHiggsSelection[2][19] = { \n", wwsftext);
    for ( int njet = 0; njet < 2; njet++) {
        fputs("    {", wwsftext);
        for (int i = 0; i < nHiggs -1 ; i++) 
            fputs(Form("%.5f,",  WWBkgScaleFactorKappaCutBased[njet][i]), wwsftext);
        if ( njet == 0 )
            fputs(Form("%.5f},\n",  WWBkgScaleFactorKappaCutBased[njet][nHiggs-1]), wwsftext);
        if ( njet == 1 )
            fputs(Form("%.5f} };\n",  WWBkgScaleFactorKappaCutBased[njet][nHiggs-1]), wwsftext);
    }

    fputs("  Int_t massIndex = -1;\n", wwsftext);
    fputs("  for (UInt_t m=0; m < 19 ; ++m) {\n", wwsftext);
    fputs("    if (mH == mHiggs[m]) massIndex = m;\n", wwsftext);
    fputs("  }\n", wwsftext);
    fputs("  if (massIndex >= 0) {\n", wwsftext);
    fputs("    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];\n", wwsftext);
    fputs("  } else {\n", wwsftext);
    fputs("    return 1.0;\n", wwsftext);
    fputs("  }\n", wwsftext);
    fputs("}\n\n", wwsftext);

    fputs("static Double_t WWBkgScaleFactorKappaMVA(Int_t mH, Int_t jetBin) {\n", wwsftext);
    fputs("assert(jetBin >= 0 && jetBin <= 1);\n", wwsftext);
    fputs("  Int_t mHiggs[19] = {115,118,120,122,124,125,126,128,130,135,140,145,150,155,160,170,180,190,200};\n",wwsftext);
    fputs("  Double_t WWBkgScaleFactorKappaHiggsSelection[2][19] = { \n", wwsftext);
    for ( int njet = 0; njet < 2; njet++) {
        fputs("    {", wwsftext);
        for (int i = 0; i < nHiggs-1; i++) 
            fputs(Form("%.5f,",  WWBkgScaleFactorKappaMVA[njet][i]), wwsftext);
        if ( njet == 0 )
            fputs(Form("%.5f},\n",  WWBkgScaleFactorKappaMVA[njet][nHiggs-1]), wwsftext);
        if ( njet == 1 )
            fputs(Form("%.5f} };\n",  WWBkgScaleFactorKappaMVA[njet][nHiggs-1]), wwsftext);
    }

    fputs("  Int_t massIndex = -1;\n", wwsftext);
    fputs("  for (UInt_t m=0; m < 19 ; ++m) {\n", wwsftext);
    fputs("    if (mH == mHiggs[m]) massIndex = m;\n", wwsftext);
    fputs("  }\n", wwsftext);
    fputs("  if (massIndex >= 0) {\n", wwsftext);
    fputs("    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];\n", wwsftext);
    fputs("  } else {\n", wwsftext);
    fputs("    return 1.0;\n", wwsftext);
    fputs("  }\n", wwsftext);
    fputs("}\n\n", wwsftext);

    fclose(wwsftext);    

}


void doMassPoint(float analysis, Option option, RunEra runEra, int mHiggs[15],  
        double WWBkgScaleFactor[2][15], double WWBkgScaleFactorKappa[2][15])
{

    float lumi = 19467;

    //
    // set up the looper
    //

    SmurfWWLooper *looper = new SmurfWWLooper(analysis, option, runEra);
    //looper->setGoodRunList("../runlists/runlist_1092.txt");
    looper->setLumiScale(lumi, lumi);

    //
    // set up samples
    //
    int ana = (int)analysis;

    SmurfSample *sample_data = new SmurfSample(option, DATA, kBlack, "Data",analysis);
    SmurfSample *sample_ww = new SmurfSample(option, WW , kBlue, "WW",analysis);
    SmurfSample *sample_top = new SmurfSample(option, TOP, kMagenta, "Top",analysis);
    SmurfSample *sample_wjets = new SmurfSample(option, WJETSDATA, kRed, "Wjets",analysis);;
    SmurfSample *sample_vv = new SmurfSample(option, VV, kRed, "VV",analysis);;
    SmurfSample *sample_dyll = new SmurfSample(option, ZLL, kBlue, "DYLL",analysis);
    SmurfSample *sample_ztt = new SmurfSample(option, ZTT, kBlue+2, "ZTT",analysis);
    SmurfSample *sample_wgamma = new SmurfSample(option, WGAMMA, kCyan, "Wgamma",analysis);

    // examples of using the skimed files 
    bool skimData = true;


    if ( skimData) {
        char *dataDir = "/smurf/cerati/skims/Run2012_Summer12_SmurfV9_53X-wwxsecfull8tev/WW/"; // TAS
        //char *dataDir  = "/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/WW/"; // UAF

        for (int jetbin = 0; jetbin < 2; jetbin++) {
            sample_data->add(Form("%s/%ij/data.root", dataDir, jetbin));    

            sample_ww->add(Form("%s/%ij/qqww.root", dataDir, jetbin));
            sample_ww->add(Form("%s/%ij/ggww.root", dataDir, jetbin));

            sample_top->add(Form("%s/%ij/ttbar_powheg.root", dataDir, jetbin));
            sample_top->add(Form("%s/%ij/tw.root", dataDir, jetbin));

            sample_vv->add(Form("%s/%ij/wz.root", dataDir, jetbin));
            sample_vv->add(Form("%s/%ij/zz.root", dataDir, jetbin));

            sample_dyll->add(Form("%s/%ij/dyll.root", dataDir, jetbin));

            sample_wgamma->add(Form("%s/%ij/wgamma.root", dataDir, jetbin));
            sample_wgamma->add(Form("%s/%ij/zgamma.root", dataDir, jetbin));
            sample_wgamma->add(Form("%s/%ij/wglll.root", dataDir, jetbin));
            sample_wgamma->add(Form("%s/%ij/www.root", dataDir, jetbin));

            sample_wjets->add(Form("%s/%ij/data_PassFail.root", dataDir, jetbin));
            sample_wjets->add(Form("%s/%ij/qqww_PassFail.root", dataDir, jetbin));
            sample_wjets->add(Form("%s/%ij/ggww_PassFail.root", dataDir, jetbin));
            sample_wjets->add(Form("%s/%ij/ttbar_powheg_PassFail.root", dataDir, jetbin));
            sample_wjets->add(Form("%s/%ij/tw_PassFail.root", dataDir, jetbin));
            sample_wjets->add(Form("%s/%ij/wz_PassFail.root", dataDir, jetbin));
            sample_wjets->add(Form("%s/%ij/zz_PassFail.root", dataDir, jetbin));
            sample_wjets->add(Form("%s/%ij/wgamma_PassFail.root", dataDir, jetbin));
            sample_wjets->add(Form("%s/%ij/zgamma_PassFail.root", dataDir, jetbin));
            sample_wjets->add(Form("%s/%ij/www_PassFail.root", dataDir, jetbin));
            sample_wjets->add(Form("%s/%ij/wglll_PassFail.root", dataDir, jetbin));
            sample_wjets->add(Form("%s/%ij/dyll_PassFail.root", dataDir, jetbin));
            
            sample_ztt->add(Form("%s/%ij/data_ztt.root", dataDir, jetbin));
        }
    } else {

        char *dataDir = "/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/";
        // data
        sample_data->add(Form("%s/data.root", dataDir));    

        // bkgs
        sample_ww->add(Form("%s/qqww.root", dataDir));
        sample_ww->add(Form("%s/ggww.root", dataDir));
        sample_top->add(Form("%s/ttbar_powheg.root", dataDir));
        sample_top->add(Form("%s/tw.root", dataDir));

        sample_vv->add(Form("%s/wz.root", dataDir));
        sample_vv->add(Form("%s/zz.root", dataDir));

        sample_dyll->add(Form("%s/dyll.root", dataDir));
        sample_wgamma->add(Form("%s/wgamma.root", dataDir));
        sample_wgamma->add(Form("%s/wglll.root", dataDir));

        sample_wjets->add(Form("%s/data.root", dataDir));    
        sample_wjets->add(Form("%s/qqww.root", dataDir));    
        sample_wjets->add(Form("%s/ggww.root", dataDir));    
        sample_wjets->add(Form("%s/ttbar.root", dataDir));    
        sample_wjets->add(Form("%s/tw.root", dataDir));    
        sample_wjets->add(Form("%s/wz.root", dataDir));    
        sample_wjets->add(Form("%s/zz.root", dataDir));    
        sample_wjets->add(Form("%s/wgamma.root", dataDir));    
        sample_wjets->add(Form("%s/wglll.root", dataDir));    
        sample_wjets->add(Form("%s/dyll.root", dataDir));    



    }

    // do the looping
    //

    looper->loop(sample_data);
    looper->loop(sample_ww);
    looper->loop(sample_top);
    looper->loop(sample_wjets);
    looper->loop(sample_vv);
    looper->loop(sample_dyll);
    looper->loop(sample_ztt);
    looper->loop(sample_wgamma);

    //
    // save histograms  
    //

    char *analysisname;
    if ( option == HWW_OPT_SMURFCUTSEL ) analysisname = "cutbased";
    if ( option == HWW_OPT_MT2DMLL )    analysisname  = "2d";

    const std::string outFile = Form("wwhistos_ww_%i_%s.root", int(analysis), analysisname);
    saveHist(outFile.c_str());  
    deleteHistos();

    // 
    // do WW background estimation
    // 

    const std::string debugFileName = Form("wwest_mH%i_%.0fpb_%s.txt", int(analysis), lumi, analysisname);
    FILE *debugtext = fopen(debugFileName.c_str(), "w"); 
    wwest(analysis, outFile.c_str(), debugtext, mHiggs,  WWBkgScaleFactor, WWBkgScaleFactorKappa);
    fclose(debugtext);

    //
    // clean up
    //

    delete looper;
    delete sample_data;
    delete sample_ww;
    delete sample_top;
    delete sample_wjets;
    delete sample_vv;
    delete sample_dyll;
    delete sample_ztt;
    delete sample_wgamma;
}


