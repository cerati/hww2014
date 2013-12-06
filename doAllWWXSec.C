
#include "TROOT.h"
#include "core/Enums.h"
#include <vector>
#include <iostream>
#include <string>

//
// function prototypes
//

void doAnalysis(const Option option, const RunEra runEra, bool doPlots = false, float analysis);

//
// main function
//

void doAllWWXSec(RunEra runEra = RUN2012)
{

    //
    // load the libraries
    //

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gROOT->ProcessLine(".L ../../../../Smurf/Core/SmurfTree.h+");
    gROOT->ProcessLine(".L ../../../../Smurf/Core/LeptonScaleLookup.cc+");
    gROOT->ProcessLine(".L libSmurfLooper.so");

    //
    // choose which analyses to run
    //

    bool runWWXSec      = 1;
    bool runHWWPresel   = 0;
    bool runHWWSSCTL    = 0;

    //
    // run the analyses
    //

    runEra = RUN2012HCP;
    if (runWWXSec)      doAnalysis(WW_OPT_SMURFXSECSEL, runEra, true, 0.0);
    if (runHWWPresel)   doAnalysis(HWW_OPT_SMURFPRESEL, runEra, false, 0.0);
    if (runHWWSSCTL)    doAnalysis(HWW_OPT_SSCTL, runEra, false, 0.0);

}

//
// do analysis
//

void doAnalysis(const Option option, const RunEra runEra, bool doPlots, float analysis)
{

    //
    // set up the looper
    //

    float lumi = 19467;
    SmurfLooper *looper = new SmurfLooper(analysis, option, runEra);
    looper->setLumiScale(lumi, lumi);

    //
    // set up the samples
    //

    // data
    SmurfSample *sample_data = new SmurfSample(option, DATA, kBlack, "Data");

    // ww signals
    SmurfSample *sample_qqww = new SmurfSample(option, QQWW, kYellow+2, "qqWW");
    SmurfSample *sample_ggww = new SmurfSample(option, GGWW, kYellow+2, "ggWW");
    SmurfSample *sample_dpsww = new SmurfSample(option, DPSWW, kYellow+2, "dpsWW");

    // higgs samples
    SmurfSample *sample_gghww = new SmurfSample(option, GGHWW, kCyan, "ggH");
    SmurfSample *sample_qqhww = new SmurfSample(option, QQHWW, kCyan, "qqH");
    SmurfSample *sample_zhww = new SmurfSample(option, ZHWW, kCyan, "ZH");
    SmurfSample *sample_whww = new SmurfSample(option, WHWW, kCyan, "WH");

    // backgrounds
    SmurfSample *sample_wz = new SmurfSample(option, WZ, kGreen, "WZ");
    SmurfSample *sample_zz = new SmurfSample(option, ZZ, kGreen, "ZZ");
    SmurfSample *sample_top = new SmurfSample(option, TOP, kMagenta, "Top");
    SmurfSample *sample_dyll = new SmurfSample(option, ZLL, kBlue, "DYLL");
    SmurfSample *sample_dyll_loosemet = new SmurfSample(option, ZLLLOOSEMET, kBlue, "DYLLLooseMET");
    SmurfSample *sample_zjets = new SmurfSample(option, ZJETS, kBlue, "Zjets");
    SmurfSample *sample_wjetsEle    = new SmurfSample(option, WJETSELEDATA, kCyan, "WjetsE",    analysis);
    SmurfSample *sample_wjetsMu     = new SmurfSample(option, WJETSMUDATA,  kCyan, "WjetsM",    analysis);
    SmurfSample *sample_wjets     = new SmurfSample(option, WJETSDATA,  kCyan, "Wjets",    analysis);
    SmurfSample *sample_wgamma = new SmurfSample(option, WGAMMANORM, kCyan, "Wgamma"); // use WGAMMANORM flag
    SmurfSample *sample_ztt = new SmurfSample(option, ZTT, kBlue+2, "Ztt");
    sample_dpsww->add("/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/ww_dps.root");

    // configure which data to use
    bool skimss = false;
    if(option == HWW_OPT_SSCTL) skimss = true;

    // add the data files
    for ( int jetbin = 0; jetbin <= 2 ; jetbin++) {

        std::string dataPrefix  = Form("/smurf/jaehyeok/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets_19p5ifb_new/WW/%ij/", jetbin);
        std::string dataSuffix  = ".root";
        if (skimss) {
            dataSuffix  = "_SS.root";
        }

        //
        // add the files
        //

        if ((analysis > 0. || (option == HWW_OPT_SMURFPRESEL)) && !skimss) {
            sample_gghww->add(Form("%s/hww%i%s", dataPrefix.c_str(), 125, dataSuffix.c_str()));
            sample_qqhww->add(Form("%s/hww%i%s", dataPrefix.c_str(), 125, dataSuffix.c_str()));
            sample_whww->add(Form("%s/hww%i%s", dataPrefix.c_str(), 125, dataSuffix.c_str()));
            sample_zhww->add(Form("%s/hww%i%s", dataPrefix.c_str(), 125, dataSuffix.c_str()));
        }

        sample_data->add(Form("%s/data%s", dataPrefix.c_str(), dataSuffix.c_str()));
        sample_qqww->add(Form("%s/qqww%s", dataPrefix.c_str(), dataSuffix.c_str()));
        sample_ggww->add(Form("%s/ggww%s", dataPrefix.c_str(), dataSuffix.c_str()));
        sample_dyll->add(Form("%s/dyll%s", dataPrefix.c_str(), dataSuffix.c_str()));
        sample_dyll_loosemet->add(Form("%s/dyll%s", dataPrefix.c_str(), dataSuffix.c_str()));
        sample_zjets->add(Form("%s/dyll%s", dataPrefix.c_str(), dataSuffix.c_str()));
        sample_top->add(Form("%s/ttbar_powheg%s", dataPrefix.c_str(), dataSuffix.c_str()));
        sample_top->add(Form("%s/tw%s", dataPrefix.c_str(), dataSuffix.c_str()));
        sample_wz->add(Form("%s/wz%s", dataPrefix.c_str(), dataSuffix.c_str()));
        sample_zz->add(Form("%s/zz%s", dataPrefix.c_str(), dataSuffix.c_str()));
        sample_wgamma->add(Form("%s/wgamma%s", dataPrefix.c_str(), dataSuffix.c_str()));
        sample_wgamma->add(Form("%s/wglll%s", dataPrefix.c_str(), dataSuffix.c_str()));

        // use directly the unmerged files
        std::string dataSuffixFR  = "_PassFail.root";
        if (skimss) {
            dataSuffixFR  = "_SS.root";
        }
        sample_wjetsEle->add(Form("%s/data%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsEle->add(Form("%s/qqww%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsEle->add(Form("%s/ggww%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsEle->add(Form("%s/ttbar_powheg%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsEle->add(Form("%s/tw%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsEle->add(Form("%s/wz%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsEle->add(Form("%s/zz%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsEle->add(Form("%s/wgamma%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsEle->add(Form("%s/zgamma%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsEle->add(Form("%s/wglll%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsEle->add(Form("%s/dyll%s", dataPrefix.c_str(), dataSuffixFR.c_str()));

        sample_wjetsMu->add(Form("%s/data%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsMu->add(Form("%s/qqww%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsMu->add(Form("%s/ggww%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsMu->add(Form("%s/ttbar_powheg%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsMu->add(Form("%s/tw%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsMu->add(Form("%s/wz%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsMu->add(Form("%s/zz%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsMu->add(Form("%s/wgamma%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsMu->add(Form("%s/zgamma%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsMu->add(Form("%s/wglll%s", dataPrefix.c_str(), dataSuffixFR.c_str()));
        sample_wjetsMu->add(Form("%s/dyll%s", dataPrefix.c_str(), dataSuffixFR.c_str()));

    }

    //
    // do the looping
    //

    looper->loop(sample_data);
    looper->loop(sample_qqww);
    looper->loop(sample_ggww);
    //looper->loop(sample_dpsww);
    looper->loop(sample_wz);
    looper->loop(sample_zz);
    looper->loop(sample_top);
    looper->loop(sample_dyll);
    looper->loop(sample_dyll_loosemet);
    looper->loop(sample_zjets);
    looper->loop(sample_wjetsEle);
    looper->loop(sample_wjetsMu);
    //looper->loop(sample_wjets);

    looper->loop(sample_wgamma);
    looper->loop(sample_ztt);

    if (analysis > 0. || (option == HWW_OPT_SMURFPRESEL)) {
        looper->loop(sample_gghww);
        looper->loop(sample_qqhww);
        looper->loop(sample_whww);
        looper->loop(sample_zhww);
    }

    //
    // analyse cross section results
    //

    std::vector<SmurfSample*> signalSamples;
    if (analysis > 0. || (option == HWW_OPT_SMURFPRESEL)) {
        signalSamples.push_back(sample_gghww);
        signalSamples.push_back(sample_qqhww);
        signalSamples.push_back(sample_whww);
        signalSamples.push_back(sample_zhww);
    }

    std::vector<SmurfSample*> samplesToTabulate;
    if (analysis > 0. || (option == HWW_OPT_SMURFPRESEL)) {
        samplesToTabulate.push_back(sample_gghww);
        samplesToTabulate.push_back(sample_qqhww);
        samplesToTabulate.push_back(sample_whww);
        samplesToTabulate.push_back(sample_zhww);
    }
    samplesToTabulate.push_back(sample_data);
    samplesToTabulate.push_back(sample_qqww);
    samplesToTabulate.push_back(sample_ggww);
    //samplesToTabulate.push_back(sample_dpsww);
    samplesToTabulate.push_back(sample_wz);
    samplesToTabulate.push_back(sample_zz);
    samplesToTabulate.push_back(sample_top);
    samplesToTabulate.push_back(sample_dyll);
    samplesToTabulate.push_back(sample_zjets);
    samplesToTabulate.push_back(sample_wjetsEle);
    samplesToTabulate.push_back(sample_wjetsMu);
    //samplesToTabulate.push_back(sample_wjets);
    samplesToTabulate.push_back(sample_wgamma);
    samplesToTabulate.push_back(sample_ztt);

    gSystem->mkdir("../wwresults", true);

    Tabulate(option, samplesToTabulate, Form("ww_analysis%i_%i", option, int(analysis)));



    //
    // save histograms  
    //

    const std::string outFile = Form("histos_ww_analysis%i_%i.root", option, int(analysis));
    saveHist(outFile.c_str());  
    deleteHistos();

    //
    // make the plots
    //

    const unsigned int kLeptonTypes = 7;
    const unsigned int kJetBins = 3;
    const char*          types[7]                = {"mm","me","em","ee", "sf", "of", "incl"};
    const char*          jetbin_names[3]         = { "0j", "1j", "2j"};

    gROOT->ProcessLine(".L tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    gROOT->ForceStyle();

    if (doPlots) {

        gSystem->mkdir("../wwplots", true);
        TFile *f = new TFile(Form("histos_ww_analysis%i_%i.root", option, int(analysis)), "READ");
        gROOT->cd();


        // scale to measured cross section
        const float wwSF = 69.23 / 57.25;  

        for (unsigned int j = 0; j < kJetBins; ++j) {
            for (unsigned int i = 0; i < kLeptonTypes; ++i) {

                if (i == fOF || i == fSF) continue;

                TCanvas *c_pt1 = makeWWXSecStack(option, analysis, sample_data, sample_qqww, sample_ggww, sample_wz, sample_zz, sample_top, sample_wjetsEle, sample_wjetsMu, sample_wgamma, sample_dyll, sample_zjets, signalSamples,
                        f, i, j, "ww_pt1", "p_{T} (leading lepton) [GeV]", lumi, false, 1, wwSF);
                TCanvas *c_pt2 = makeWWXSecStack(option, analysis, sample_data, sample_qqww, sample_ggww, sample_wz, sample_zz, sample_top, sample_wjetsEle, sample_wjetsMu,  sample_wgamma, sample_dyll, sample_zjets, signalSamples,
                        f, i, j, "ww_pt2", "p_{T} (trailing lepton) [GeV]", lumi, false, 1, wwSF);
                TCanvas *c_mll = makeWWXSecStack(option, analysis, sample_data, sample_qqww, sample_ggww, sample_wz, sample_zz, sample_top, sample_wjetsEle, sample_wjetsMu,  sample_wgamma, sample_dyll, sample_zjets, signalSamples,
                        f, i, j, "ww_mll", "M_{l,l} [GeV/c^{2}]", lumi, false, 1, wwSF);
                TCanvas *c_ptll = makeWWXSecStack(option, analysis, sample_data, sample_qqww, sample_ggww, sample_wz, sample_zz, sample_top, sample_wjetsEle, sample_wjetsMu,  sample_wgamma, sample_dyll, sample_zjets, signalSamples,
                        f, i, j, "ww_ptll", "P_{T} (ll) [GeV]", lumi, false, 1, wwSF);
                TCanvas *c_mt = makeWWXSecStack(option, analysis, sample_data, sample_qqww, sample_ggww, sample_wz, sample_zz, sample_top, sample_wjetsEle, sample_wjetsMu,  sample_wgamma, sample_dyll, sample_zjets, signalSamples,
                        f, i, j, "ww_mt", "Transverse Higgs Mass [GeV]", lumi, false, 1, wwSF);
                TCanvas *c_met = makeWWXSecStack(option, analysis, sample_data, sample_qqww, sample_ggww, sample_wz, sample_zz, sample_top, sample_wjetsEle, sample_wjetsMu,  sample_wgamma, sample_dyll, sample_zjets, signalSamples,
                        f, i, j, "ww_met", "Missing Energy Variable [GeV]", lumi, false, 1, wwSF);
                TCanvas *c_dphi = makeWWXSecStack(option, analysis, sample_data, sample_qqww, sample_ggww, sample_wz, sample_zz, sample_top, sample_wjetsEle, sample_wjetsMu,  sample_wgamma, sample_dyll, sample_zjets, signalSamples,
                        f, i, j, "ww_dphi", "#Delta#phi_{l,l} [rad]", lumi, false, 1, wwSF);

                TCanvas *c1_big = new TCanvas("c1_big", "c1_big", 800, 1000);
                c1_big->Divide(2, 2);
                c1_big->cd(1);
                (TPad*)c_pt1->GetListOfPrimitives()->FindObject("p_main")->Draw();
                (TPad*)c_pt1->GetListOfPrimitives()->FindObject("p_pull")->Draw();
                c1_big->cd(2);
                (TPad*)c_pt2->GetListOfPrimitives()->FindObject("p_main")->Draw();
                (TPad*)c_pt2->GetListOfPrimitives()->FindObject("p_pull")->Draw();
                c1_big->cd(3);
                (TPad*)c_ptll->GetListOfPrimitives()->FindObject("p_main")->Draw();
                (TPad*)c_ptll->GetListOfPrimitives()->FindObject("p_pull")->Draw();
                c1_big->cd(4);
                (TPad*)c_mll->GetListOfPrimitives()->FindObject("p_main")->Draw();
                (TPad*)c_mll->GetListOfPrimitives()->FindObject("p_pull")->Draw();
                c1_big->SaveAs(Form("../wwplots/ww_analysis%i_%i_ALL_%s_%s.eps", option, int(analysis), types[i], jetbin_names[j]));
                c1_big->SaveAs(Form("../wwplots/ww_analysis%i_%i_ALL_%s_%s.png", option, int(analysis), types[i], jetbin_names[j]));

            }
        }

        for (unsigned int i = 0; i < kLeptonTypes; ++i) {

            TCanvas *c_detajj = makeWWXSecStack(option, analysis, sample_data, sample_qqww, sample_ggww, sample_wz, sample_zz, sample_top, sample_wjetsEle, sample_wjetsMu,  sample_wgamma, sample_dyll, sample_zjets, signalSamples,
                    f, i, 2, "ww_detajj", "#Delta#eta_{jet1, jet2}", lumi, false, 1, wwSF);
            TCanvas *c_mjj = makeWWXSecStack(option, analysis, sample_data, sample_qqww, sample_ggww, sample_wz, sample_zz, sample_top, sample_wjetsEle, sample_wjetsMu,  sample_wgamma, sample_dyll, sample_zjets, signalSamples,
                    f, i, 2, "ww_mjj", "M (jet1, jet2) [GeV]", lumi, false, 1, wwSF);

            c_detajj->SaveAs(Form("../wwplots/ww_analysis%i_detajj_%s_2j.png", option, types[i]));
            c_mjj->SaveAs(Form("../wwplots/ww_analysis%i_mjj_%s_2j.png", option, types[i]));

        }


        delete f;

    }

    //
    // clean up
    //

    delete looper;
    delete sample_data;
    delete sample_qqww;
    delete sample_ggww;
    delete sample_wz;
    delete sample_zz;
    delete sample_top;
    delete sample_dyll;
    delete sample_dyll_loosemet;
    delete sample_ztt;
    delete sample_zjets;
    delete sample_wjetsEle;
    delete sample_wjetsMu;
    delete sample_wjets;
    delete sample_wgamma;
    delete sample_dyll_loosemet;
}

