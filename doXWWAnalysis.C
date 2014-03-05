
#include "TROOT.h"
#include "core/Enums.h"
#include <vector>
#include <iostream>


void doXWWAnalysis(RunEra runEra = RUN2012)
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
    gROOT->ProcessLine(".L ../../../NtupleMacros/Tools/goodrun.cc+");
    gROOT->ProcessLine(".L libSmurfLooper.so");

    //
    // configuration  
    //

    doMassPoint(125.0, HWW_OPT_MT2DMLL_JCP, runEra);
    doMassPoint(125.0, XWW_OPT_MT2DMLL_JCP, runEra);
}

void doMassPoint(float analysis, Option option, RunEra runEra)
{

    bool blind = false;

    //
    // set up the looper
    //

    SmurfLooper *looper = new SmurfLooper(analysis, option, runEra);
    //looper->setGoodRunList("../runlists/runlist_1092.txt");
    const float lumi = 19467;
    looper->setLumiScale(lumi, lumi); 


    //
    // set up samples
    //
    int ana = (int)analysis;

    SmurfSample *sample_data = new SmurfSample(option, DATA, kBlack, "Data", analysis);

    // higgs samples
    SmurfSample *sample_gghww = new SmurfSample(option, GGHWW, kCyan, "ggH", analysis);
    SmurfSample *sample_qqhww = new SmurfSample(option, QQHWW, kCyan, "qqH", analysis);
    SmurfSample *sample_zhww = new SmurfSample(option, ZHWW, kCyan, "ZH", analysis);
    SmurfSample *sample_whww = new SmurfSample(option, WHWW, kCyan, "WH", analysis);

    // backgrounds
    SmurfSample *sample_qqww = new SmurfSample(option, QQWW, kYellow+2, "qqWW", analysis);
    SmurfSample *sample_ggww = new SmurfSample(option, GGWW, kYellow+2, "ggWW", analysis);

    SmurfSample *sample_vv = new SmurfSample(option, VV, kGreen, "VV", analysis);
    SmurfSample *sample_top = new SmurfSample(option, TOP, kMagenta, "Top", analysis);
    // note this considers only the Drell Yan part
    SmurfSample *sample_dyll = new SmurfSample(option, ZLL, kBlue, "DYLL", analysis);
    // NOTE. this zjets is a historial nameing..
    // this corresponds to the WZ/ZZ where both leptons are from the same Z
    //SmurfSample *sample_zjets       = new SmurfSample(option, ZJETS, kBlue, "Zjets", analysis);
    SmurfSample *sample_wjets       = new SmurfSample(option, WJETSDATA,    kCyan, "Wjets",     analysis);
    SmurfSample *sample_wjetsEle    = new SmurfSample(option, WJETSELEDATA, kCyan, "WjetsE",    analysis);
    SmurfSample *sample_wjetsMu     = new SmurfSample(option, WJETSMUDATA,  kCyan, "WjetsM",    analysis);
    SmurfSample *sample_wgamma      = new SmurfSample(option, WGAMMA,       kCyan, "Wgamma",        analysis);
    SmurfSample *sample_wgammanorm  = new SmurfSample(option, WGAMMANORM,   kCyan, "Wgammanorm",    analysis);
    SmurfSample *sample_wg3l        = new SmurfSample(option, WG3L,         kCyan, "Wg3l",          analysis);
    SmurfSample *sample_ztt         = new SmurfSample(option, ZTT,          kBlue+2, "Ztt",         analysis);

    // Below are needed for the central shape
    //SmurfSample *sample_dyll_loosemet = new SmurfSample(option, ZLLLOOSEMET, kBlue, "DYLLLooseMET", analysis);

    // Below are needed for the shape variation studies
    SmurfSample *sample_top_var         = new SmurfSample(option, TOPALTER,     kMagenta,   "TopVar", analysis);
    SmurfSample *sample_wwmcnlo         = new SmurfSample(option, WWMCNLO,      kBlue+2,    "WWMCNLO", analysis);
    SmurfSample *sample_wwmcnlo_up      = new SmurfSample(option, WWMCNLOUP,    kBlue+2,    "WWMCNLOUp", analysis);
    SmurfSample *sample_wwmcnlo_down    = new SmurfSample(option, WWMCNLODOWN,  kBlue+2,    "WWMCNLODown", analysis);
    //SmurfSample *sample_wjets_mc        = new SmurfSample(option, WJETS,        kCyan,      "WjetsMC", analysis);
    //SmurfSample *sample_wjets_mc_loose  = new SmurfSample(option, WJETSMCLOOSE, kCyan,      "WjetsMCLoose", analysis);    
    //SmurfSample *sample_dyll_data       = new SmurfSample(option, ZLLDATA,      kBlue,      "DYLLDATA", analysis);
    SmurfSample *sample_gghww_ref       = new SmurfSample(option, GGHWWREF,     kCyan,      "ggHRef", analysis);
    SmurfSample *sample_gghww_jhu       = new SmurfSample(option, GGHWWJHU,     kCyan,      "ggHJHU", analysis);


    bool skimwithmva = true;

    // Below is the setup for running the smurf ntuples at the ww-preselection skim
    if ( skimwithmva) {
        char *dataDir = "/smurf/jaehyeok/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets_19p5ifb_new/WW/"; // TAS
        //char *dataDir  = "/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/WW/"; // UAF

        // set up the mva histogram ranges
        unsigned int n = 20;
        float min = -1.0;
        float max = 1.0;
        if ( option == HWW_OPT_SMURFMESEL) min = 0.0;

        if ( analysis > 0. ) {
            // higgs
            for (int njet=0; njet < 2; njet++) {
                if ( option == XWW_OPT_MT2DMLL_JCP )  
                    sample_gghww->add(Form("%s/mva/%i/xww2p%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
                else
                    sample_gghww->add(Form("%s/mva/%i/hww%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
            }
            sample_gghww->addShapeVariation2D(LEPEFFVAR,    "histo_ggH_CMS_hww_MVALepEffBounding", false);
            sample_gghww->addShapeVariation2D(LEPRESVAR,    "histo_ggH_CMS_hww_MVALepResBounding", true);
            sample_gghww->addShapeVariation2D(METVAR,       "histo_ggH_CMS_hww_MVAMETResBounding", true);
            sample_gghww->addShapeVariation2D(JETRESVAR,    "histo_ggH_CMS_hww_MVAJESBounding", false);

            for (int njet=0; njet < 2; njet++) 
                sample_qqhww->add(Form("%s/mva/%i/hww%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
            sample_qqhww->addShapeVariation2D(LEPEFFVAR, "histo_qqH_CMS_hww_MVALepEffBounding", false);
            sample_qqhww->addShapeVariation2D(LEPRESVAR, "histo_qqH_CMS_hww_MVALepResBounding", true);
            sample_qqhww->addShapeVariation2D(METVAR, "histo_qqH_CMS_hww_MVAMETResBounding", true);
            sample_qqhww->addShapeVariation2D(JETRESVAR, "histo_qqH_CMS_hww_MVAJESBounding", false);

            for (int njet=0; njet < 2; njet++) 
                sample_zhww->add(Form("%s/mva/%i/hww%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
            sample_zhww->addShapeVariation2D(LEPEFFVAR, "histo_ZH_CMS_hww_MVALepEffBounding", false);
            sample_zhww->addShapeVariation2D(LEPRESVAR, "histo_ZH_CMS_hww_MVALepResBounding", true);
            sample_zhww->addShapeVariation2D(METVAR, "histo_ZH_CMS_hww_MVAMETResBounding", true);
            sample_zhww->addShapeVariation2D(JETRESVAR, "histo_ZH_CMS_hww_MVAJESBounding", false);

            for (int njet=0; njet < 2; njet++) 
                sample_whww->add(Form("%s/mva/%i/hww%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
            sample_whww->addShapeVariation2D(LEPEFFVAR, "histo_WH_CMS_hww_MVALepEffBounding", false);
            sample_whww->addShapeVariation2D(LEPRESVAR, "histo_WH_CMS_hww_MVALepResBounding", true);
            sample_whww->addShapeVariation2D(METVAR, "histo_WH_CMS_hww_MVAMETResBounding", true);
            sample_whww->addShapeVariation2D(JETRESVAR, "histo_WH_CMS_hww_MVAJESBounding", false);
        }

        // bkgs
        // qqWW
        for (int njet=0; njet < 2; njet++) 
            sample_qqww->add(Form("%s/mva/%i/qqww_%ij.root", dataDir, int(analysis), njet));

        if ( analysis < 200) 
            sample_qqww->addShapeVariation2D(LEPEFFVAR, "histo_qqWW_CMS_hww_MVALepEffBounding", true);
        else 
            sample_qqww->addShapeVariation2D(LEPEFFVAR, "histo_qqWW_CMS_hww_MVALepEffBounding", false);
        sample_qqww->addShapeVariation2D(LEPRESVAR, "histo_qqWW_CMS_hww_MVALepResBounding", true);
        sample_qqww->addShapeVariation2D(METVAR, "histo_qqWW_CMS_hww_MVAMETResBounding", true);
        sample_qqww->addShapeVariation2D(JETRESVAR, "histo_qqWW_CMS_hww_MVAJESBounding", false);

        // ggWW
        for (int njet=0; njet < 2; njet++) 
            sample_ggww->add(Form("%s/mva/%i/ggww_%ij.root", dataDir, int(analysis), njet));

        if ( analysis < 200) 
            sample_ggww->addShapeVariation2D(LEPEFFVAR, "histo_ggWW_CMS_hww_MVALepEffBounding", true);
        else 
            sample_ggww->addShapeVariation2D(LEPEFFVAR, "histo_ggWW_CMS_hww_MVALepEffBounding", false);
        sample_ggww->addShapeVariation2D(LEPRESVAR, "histo_ggWW_CMS_hww_MVALepResBounding", true);
        sample_ggww->addShapeVariation2D(METVAR, "histo_ggWW_CMS_hww_MVAMETResBounding", true);
        sample_ggww->addShapeVariation2D(JETRESVAR, "histo_ggWW_CMS_hww_MVAJESBounding", false);

        // Drell-Yan
        for (int njet=0; njet < 2; njet++) {
            sample_dyll->add(Form("%s/mva/%i/dyll_%ij.root", dataDir, int(analysis), njet));
            //sample_dyll_loosemet->add(Form("%s/mva/%i/dyll_%ij.root", dataDir, int(analysis), njet));
            //ample_zjets->add(Form("%s/mva/%i/dyll_%ij.root", dataDir, int(analysis), njet));
        }

        // Top 
        for (int njet=0; njet < 2; njet++) {
            sample_top->add(Form("%s/mva/%i/ttbar_powheg_%ij.root", dataDir, int(analysis), njet));
            sample_top->add(Form("%s/mva/%i/tw_%ij.root", dataDir, int(analysis), njet));
        }

        sample_top->addShapeVariation2D(LEPRESVAR, "histo_Top_CMS_hww_MVALepResBounding", true);
        sample_top->addShapeVariation2D(METVAR, "histo_Top_CMS_hww_MVAMETResBounding", true);
        sample_top->addShapeVariation2D(JETRESVAR, "histo_Top_CMS_hww_MVAJESBounding", false);

        // VV 
        for (int njet=0; njet < 2; njet++) {
            sample_vv->add(Form("%s/mva/%i/wz_%ij.root", dataDir, int(analysis), njet));
            sample_vv->add(Form("%s/mva/%i/zz_%ij.root", dataDir, int(analysis), njet)); 
            sample_vv->add(Form("%s/mva/%i/www_%ij.root", dataDir, int(analysis), njet)); 
        }

        sample_vv->addShapeVariation2D(LEPEFFVAR, "histo_VV_CMS_hww_MVALepEffBounding", false);
        sample_vv->addShapeVariation2D(LEPRESVAR, "histo_VV_CMS_hww_MVALepResBounding", true);
        sample_vv->addShapeVariation2D(METVAR, "histo_VV_CMS_hww_MVAMETResBounding", true);
        sample_vv->addShapeVariation2D(JETRESVAR, "histo_VV_CMS_hww_MVAJESBounding", false);


        // Wgamma : l + gamma sample 
        for (int njet=0; njet < 2; njet++) {
            sample_wgamma->add(Form("%s/mva/%i/wgammafo_%ij.root", dataDir, int(analysis), njet));
            sample_wgamma->add(Form("%s/mva/%i/zgammafo_%ij.root", dataDir, int(analysis), njet));
        }

        // Wgamma Normalization
        for (int njet=0; njet < 2; njet++) {
            sample_wgammanorm->add(Form("%s/mva/%i/wgamma_%ij.root", dataDir, int(analysis), njet));
            sample_wgammanorm->add(Form("%s/mva/%i/zgamma_%ij.root", dataDir, int(analysis), njet));
        }

        // Wg3l ( Wgamma* )
        for (int njet=0; njet < 2; njet++) {
            sample_wg3l->add(Form("%s/mva/%i/wglll_%ij.root", dataDir, int(analysis), njet));
        } 

        // Ztt
        for (int njet=0; njet < 2; njet++) {
            sample_ztt->add(Form("%s/mva/%i/data_ztt_%ij.root", dataDir, int(analysis), njet));
        }

        // Wjets 
        for (int njet=0; njet < 2; njet++) {
            sample_wjetsEle->add(Form("%s/mva/%i/data_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/qqww_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/ggww_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/ttbar_powheg_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/tw_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/wz_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/zz_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/wgamma_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/zgamma_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/wglll_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/dyll_PassFail_%ij.root", dataDir, int(analysis), njet));

            sample_wjetsMu->add(Form("%s/mva/%i/data_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/qqww_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/ggww_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/ttbar_powheg_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/tw_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/wz_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/zz_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/wgamma_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/zgamma_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/wglll_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/dyll_PassFail_%ij.root", dataDir, int(analysis), njet));
        }
        sample_wjetsEle->addShapeVariation2D(WJETSELESHAPEVAR, "histo_WjetsE_CMS_hww_MVAWEBounding", true);
        sample_wjetsMu->addShapeVariation2D(WJETSMUSHAPEVAR, "histo_WjetsM_CMS_hww_MVAWMBounding", true);


        // data
        if ( !blind) {
            for (int njet=0; njet < 2; njet++) 
                sample_data->add(Form("%s/mva/%i/data_%ij.root", dataDir, int(analysis), njet));
        }

        // for Shape variations
        // wjets
        for (int njet=0; njet < 2; njet++) {
            //sample_wjets_mc->add(Form("%s/mva/%i/wjets_%ij.root", dataDir, int(analysis), njet));
            //sample_wjets_mc_loose->add(Form("%s/mva/%i/wjets_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_top_var->add(Form("%s/mva/%i/ttbar_%ij.root", dataDir, int(analysis), njet));
            sample_top_var->add(Form("%s/mva/%i/tw_%ij.root", dataDir, int(analysis), njet));
            sample_wwmcnlo->add(Form("%s/mva/%i/wwmcnlo_%ij.root", dataDir, int(analysis), njet));
            sample_wwmcnlo_up->add(Form("%s/mva/%i/wwmcnloup_%ij.root", dataDir, int(analysis), njet));
            sample_wwmcnlo_down->add(Form("%s/mva/%i/wwmcnlodown_%ij.root", dataDir, int(analysis), njet));
            //sample_dyll_data->add(Form("%s/mva/%i/data_%ij.root", dataDir, int(analysis), njet));
            //sample_dyll_data->add(Form("%s/mva/%i/wz_%ij.root", dataDir, int(analysis), njet));
            //sample_dyll_data->add(Form("%s/mva/%i/zz_%ij.root", dataDir, int(analysis), njet));

            sample_gghww_ref->add(Form("%s/mva/%i/hww%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
            sample_gghww_jhu->add(Form("%s/mva/%i/xww2p%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
        }
    } 

    //
    // do the looping
    //

    looper->loop(sample_data);
    looper->loop(sample_gghww);
    looper->loop(sample_qqhww);
    looper->loop(sample_whww);
    looper->loop(sample_zhww);
    looper->loop(sample_qqww);
    looper->loop(sample_ggww);
    looper->loop(sample_vv);
    looper->loop(sample_top);
    looper->loop(sample_dyll);
    //looper->loop(sample_dyll_loosemet);
    //looper->loop(sample_zjets);
    looper->loop(sample_wjetsEle);
    looper->loop(sample_wjetsMu);
    looper->loop(sample_wgammanorm);
    looper->loop(sample_wgamma);
    looper->loop(sample_wg3l);
    looper->loop(sample_ztt);

    // for shape syst.
    //looper->loop(sample_wjets_mc_loose);
    //looper->loop(sample_wjets_mc);
    looper->loop(sample_top_var);
    looper->loop(sample_wwmcnlo);
    looper->loop(sample_wwmcnlo_up);
    looper->loop(sample_wwmcnlo_down);
    looper->loop(sample_dyll_data);
    looper->loop(sample_gghww_ref);
    looper->loop(sample_gghww_jhu);


    //
    // make tables
    //

    std::vector<SmurfSample*> samplesToTabulate;

    samplesToTabulate.push_back(sample_data);
    samplesToTabulate.push_back(sample_gghww);
    samplesToTabulate.push_back(sample_qqhww);
    samplesToTabulate.push_back(sample_whww);
    samplesToTabulate.push_back(sample_zhww);
    samplesToTabulate.push_back(sample_qqww);
    samplesToTabulate.push_back(sample_ggww);
    samplesToTabulate.push_back(sample_vv);
    samplesToTabulate.push_back(sample_top);
    samplesToTabulate.push_back(sample_dyll);
    samplesToTabulate.push_back(sample_zjets);
    samplesToTabulate.push_back(sample_wjetsEle);
    samplesToTabulate.push_back(sample_wjetsMu);
    samplesToTabulate.push_back(sample_wgammanorm);
    samplesToTabulate.push_back(sample_wgamma);
    samplesToTabulate.push_back(sample_wg3l);
    samplesToTabulate.push_back(sample_ztt);
    printResultsTable(samplesToTabulate, option, true);

    //
    // for shape variations
    // 
    //samplesToTabulate.push_back(sample_dyll_loosemet);
    //samplesToTabulate.push_back(sample_wjets_mc); 
    //samplesToTabulate.push_back(sample_wjets_mc_loose);
    samplesToTabulate.push_back(sample_top_var);
    samplesToTabulate.push_back(sample_wwmcnlo);
    samplesToTabulate.push_back(sample_wwmcnlo_up);
    samplesToTabulate.push_back(sample_wwmcnlo_down);
    samplesToTabulate.push_back(sample_dyll_data);
    samplesToTabulate.push_back(sample_gghww_ref); 
    samplesToTabulate.push_back(sample_gghww_jhu); 

    //  
    // make cards
    //
    printf("\n\n[doAllHWW::doMassPoint] Writing cards\n");

    const std::string cardDir = "../cards/hwwjcp_19p5fb/";

    ShapeVar_t mva_option = (1ll<<STATVAR) | (1ll<<TOPSHAPEVAR) | (1ll<<WWSHAPEVAR) | (1ll<<DYSHAPEVAR)
        | (1ll<<LEPEFFVAR) | (1ll<<METVAR) | (1ll<<LEPRESVAR) | (1ll<<JETRESVAR) | (1ll<<WJETSELESHAPEVAR) | (1ll<<WJETSMUSHAPEVAR);

    for (int jetbin = 0; jetbin < 2; ++jetbin) { 
        printCard(samplesToTabulate, option, jetbin, analysis, cardDir, ((1<<1)|(1<<2)), mva_option, runEra);
        print2DShapeHistograms(samplesToTabulate, option, jetbin, analysis, cardDir, ((1<<1)|(1<<2)), runEra);
    }


    //
    // save histograms  
    //

    const std::string outFile = Form("histos_hww_analysis%i_%i.root", option, int(analysis));
    saveHist(outFile.c_str());  
    deleteHistos();


    //
    // clean up
    //

    delete looper;
    delete sample_data;
    delete sample_gghww;
    delete sample_qqhww;
    delete sample_zhww;
    delete sample_whww;
    delete sample_qqww;
    delete sample_ggww;
    delete sample_vv;
    delete sample_top;
    delete sample_dyll;
    delete sample_ztt;
    //delete sample_zjets;
    delete sample_wjets;
    delete sample_wjetsEle;
    delete sample_wjetsMu;
    delete sample_wgamma;
    delete sample_wgammanorm;
    delete sample_wg3l;

    //delete sample_dyll_loosemet;
    delete sample_top_var;
    delete sample_wwmcnlo;
    delete sample_wwmcnlo_up;
    delete sample_wwmcnlo_down;
    //delete sample_wjets_mc;
    //delete sample_wjets_mc_loose;
    //delete sample_dyll_data;
    delete sample_gghww_ref;
    delete sample_gghww_jhu;


}


