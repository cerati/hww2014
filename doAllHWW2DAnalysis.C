
#include "TROOT.h"
#include "core/Enums.h"
#include <vector>
#include <iostream>

const bool makePlots_ = false;

void doAllHWW2DAnalysis(RunEra runEra = RUN2012)
{

    //
    // load the libraries
    //

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gROOT->ProcessLine(".L ../Smurf/Core/SmurfTree.h+");
    gROOT->ProcessLine(".L ../Smurf/Core/LeptonScaleLookup.cc+");
    gROOT->ProcessLine(".L ../Tools/goodrun.cc+");
    gROOT->ProcessLine(".L libSmurfLooper.so");

    //
    // configuration  
    //
    
    bool runSS2D = 0;
    
    bool run110 = 0; 
    bool run115 = 0;
    bool run120 = 0;
    bool run125 = 1; 
    bool run130 = 0; 
    bool run135 = 0;
    bool run140 = 0; 
    bool run150 = 0; 
    bool run160 = 0; 
    bool run170 = 0;
    bool run180 = 0; 
    bool run190 = 0;
    bool run200 = 0; 
    bool run250 = 0;
    bool run300 = 0;  
    bool run350 = 0;
    bool run400 = 0;
    bool run450 = 0;
    bool run500 = 0;
    bool run550 = 0;
    bool run600 = 0;

    //
    // now run it
    //
    if ( runSS2D )  doMassPoint(125.0, HWW_OPT_SSCTL2D, runEra);

    if ( run110 )   doMassPoint(110.0, HWW_OPT_MT2DMLL, runEra);
    if ( run115 )   doMassPoint(115.0, HWW_OPT_MT2DMLL, runEra);
    if ( run120 )   doMassPoint(120.0, HWW_OPT_MT2DMLL, runEra);
    if ( run125 )   doMassPoint(125.0, HWW_OPT_MT2DMLL, runEra);
    if ( run130 )   doMassPoint(130.0, HWW_OPT_MT2DMLL, runEra);
    if ( run135 )   doMassPoint(135.0, HWW_OPT_MT2DMLL, runEra);
    if ( run140 )   doMassPoint(140.0, HWW_OPT_MT2DMLL, runEra);
    if ( run150 )   doMassPoint(150.0, HWW_OPT_MT2DMLL, runEra);
    if ( run160 )   doMassPoint(160.0, HWW_OPT_MT2DMLL, runEra);
    if ( run170 )   doMassPoint(170.0, HWW_OPT_MT2DMLL, runEra);
    if ( run180 )   doMassPoint(180.0, HWW_OPT_MT2DMLL, runEra);
    if ( run190 )   doMassPoint(190.0, HWW_OPT_MT2DMLL, runEra);
    if ( run200 )   doMassPoint(200.0, HWW_OPT_MT2DMLL, runEra);
    if ( run250 )   doMassPoint(250.0, HWW_OPT_MT2DMLL, runEra);
    if ( run300 )   doMassPoint(300.0, HWW_OPT_MT2DMLL, runEra);
    if ( run350 )   doMassPoint(350.0, HWW_OPT_MT2DMLL, runEra);
    if ( run400 )   doMassPoint(400.0, HWW_OPT_MT2DMLL, runEra);
    if ( run450 )   doMassPoint(450.0, HWW_OPT_MT2DMLL, runEra);
    if ( run500 )   doMassPoint(500.0, HWW_OPT_MT2DMLL, runEra);
    if ( run550 )   doMassPoint(550.0, HWW_OPT_MT2DMLL, runEra);
    if ( run600 )   doMassPoint(600.0, HWW_OPT_MT2DMLL, runEra);
    

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
    //const float lumi = 12100; // HCP
    //const float lumi = 7367;  // postHCP
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
    //SmurfSample *sample_zjets   = new SmurfSample(option, ZJETS, kBlue, "Zjets", analysis);
    SmurfSample *sample_wjetsEle    = new SmurfSample(option, WJETSELEDATA, kCyan, "WjetsE",    analysis);
    SmurfSample *sample_wjetsMu     = new SmurfSample(option, WJETSMUDATA,  kCyan, "WjetsM",    analysis);
    SmurfSample *sample_wgamma      = new SmurfSample(option, WGAMMA, kCyan, "Wgamma", analysis);
    SmurfSample *sample_wgammanorm  = new SmurfSample(option, WGAMMANORM,   kCyan, "Wgammanorm",    analysis);
    SmurfSample *sample_wg3l        = new SmurfSample(option, WG3L, kCyan, "Wg3l", analysis);
    SmurfSample *sample_ztt         = new SmurfSample(option, ZTT, kBlue+2, "Ztt", analysis);

    // Below are needed for the central shape
    //SmurfSample *sample_dyll_loosemet = new SmurfSample(option, ZLLLOOSEMET, kBlue, "DYLLLooseMET", analysis);

    // Below are needed for the shape variation studies
    SmurfSample *sample_top_var = new SmurfSample(option, TOPALTER, kMagenta, "TopVar", analysis);
    SmurfSample *sample_wwmcnlo = new SmurfSample(option, WWMCNLO, kBlue+2, "WWMCNLO", analysis);
    SmurfSample *sample_wwmcnlo_up = new SmurfSample(option, WWMCNLOUP, kBlue+2, "WWMCNLOUp", analysis);
    SmurfSample *sample_wwmcnlo_down = new SmurfSample(option, WWMCNLODOWN, kBlue+2, "WWMCNLODown", analysis);
    //SmurfSample *sample_wjets_mc = new SmurfSample(option, WJETS, kCyan, "WjetsMC", analysis);
    //SmurfSample *sample_wjets_mc_loose = new SmurfSample(option, WJETSMCLOOSE, kCyan, "WjetsMCLoose", analysis);    
    //SmurfSample *sample_dyll_data = new SmurfSample(option, ZLLDATA, kBlue, "DYLLDATA", analysis);

    // CHOOSE ONLE ONE 
    bool skimwithmva    = true;
    bool skimwithmvass  = false;

    if ( skimwithmva) {
        char *dataDir = "/smurf/cerati/skims/Run2012_Summer12_SmurfV9_53X-wwxsecfull8tev/WW/"; // TAS 
        //char *dataDir = "/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/WW/"; // UAF

        if ( analysis > 0. ) {
            // higgs
            for (int njet=0; njet < 3; njet++) 
                sample_gghww->add(Form("%s/mva/%i/hww%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
            sample_gghww->addShapeVariation2D(QCDSCALEVAR, "histo_ggH_CMS_hww_MVAggHBounding", true);
            sample_gghww->addShapeVariation2D(LEPEFFVAR,    "histo_ggH_CMS_hww_MVALepEffBounding", false);
            sample_gghww->addShapeVariation2D(LEPRESVAR,    "histo_ggH_CMS_hww_MVALepResBounding", true);
            sample_gghww->addShapeVariation2D(METVAR,       "histo_ggH_CMS_hww_MVAMETResBounding", true);
            sample_gghww->addShapeVariation2D(JETRESVAR,    "histo_ggH_CMS_hww_MVAJESBounding", false);

            for (int njet=0; njet < 3; njet++) 
                sample_qqhww->add(Form("%s/mva/%i/hww%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
            sample_qqhww->addShapeVariation2D(LEPEFFVAR, "histo_qqH_CMS_hww_MVALepEffBounding", false);
            sample_qqhww->addShapeVariation2D(LEPRESVAR, "histo_qqH_CMS_hww_MVALepResBounding", true);
            sample_qqhww->addShapeVariation2D(METVAR, "histo_qqH_CMS_hww_MVAMETResBounding", true);
            sample_qqhww->addShapeVariation2D(JETRESVAR, "histo_qqH_CMS_hww_MVAJESBounding", false);

            for (int njet=0; njet < 3; njet++) 
                sample_zhww->add(Form("%s/mva/%i/hww%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
            sample_zhww->addShapeVariation2D(LEPEFFVAR, "histo_ZH_CMS_hww_MVALepEffBounding", false);
            sample_zhww->addShapeVariation2D(LEPRESVAR, "histo_ZH_CMS_hww_MVALepResBounding", true);
            sample_zhww->addShapeVariation2D(METVAR, "histo_ZH_CMS_hww_MVAMETResBounding", true);
            sample_zhww->addShapeVariation2D(JETRESVAR, "histo_ZH_CMS_hww_MVAJESBounding", false);

            for (int njet=0; njet < 3; njet++) 
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
        sample_qqww->add(Form("%s/mva/%i/qqww_%ij.root", dataDir, 125, 2));

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
        sample_ggww->add(Form("%s/mva/%i/ggww_%ij.root", dataDir, 125, 2));

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
            //sample_zjets->add(Form("%s/mva/%i/dyll_%ij.root", dataDir, int(analysis), njet));
        }
        sample_dyll->add(Form("%s/mva/%i/dyll_%ij.root", dataDir, 125, 2));
        //sample_dyll_loosemet->add(Form("%s/mva/%i/dyll_%ij.root", dataDir, 125, 2));
        //sample_zjets->add(Form("%s/mva/%i/dyll_%ij.root", dataDir, 125, 2));

        // Top 
        for (int njet=0; njet < 2; njet++) {
            sample_top->add(Form("%s/mva/%i/ttbar_powheg_%ij.root", dataDir, int(analysis), njet));
            sample_top->add(Form("%s/mva/%i/tw_%ij.root", dataDir, int(analysis), njet));
        }
        sample_top->add(Form("%s/mva/%i/ttbar_powheg_%ij.root", dataDir, 125, 2));
        sample_top->add(Form("%s/mva/%i/tw_%ij.root", dataDir, 125, 2));

        sample_top->addShapeVariation2D(LEPRESVAR, "histo_Top_CMS_hww_MVALepResBounding", true);
        sample_top->addShapeVariation2D(METVAR, "histo_Top_CMS_hww_MVAMETResBounding", true);
        sample_top->addShapeVariation2D(JETRESVAR, "histo_Top_CMS_hww_MVAJESBounding", false);

        // VV 
        for (int njet=0; njet < 2; njet++) {
            sample_vv->add(Form("%s/mva/%i/wz_%ij.root", dataDir, int(analysis), njet));
            sample_vv->add(Form("%s/mva/%i/zz_%ij.root", dataDir, int(analysis), njet)); 
           sample_vv->add(Form("%s/mva/%i/www_%ij.root", dataDir, int(analysis), njet)); 
        }
        sample_vv->add(Form("%s/mva/%i/wz_%ij.root", dataDir, 125, 2));
        sample_vv->add(Form("%s/mva/%i/zz_%ij.root", dataDir, 125, 2)); 
        sample_vv->add(Form("%s/mva/%i/www_%ij.root", dataDir, 125, 2)); 

        sample_vv->addShapeVariation2D(LEPEFFVAR, "histo_VV_CMS_hww_MVALepEffBounding", false);
        sample_vv->addShapeVariation2D(LEPRESVAR, "histo_VV_CMS_hww_MVALepResBounding", true);
        sample_vv->addShapeVariation2D(METVAR, "histo_VV_CMS_hww_MVAMETResBounding", true);
        sample_vv->addShapeVariation2D(JETRESVAR, "histo_VV_CMS_hww_MVAJESBounding", false);

        // Wgamma : l + gamma sample 
        for (int njet=0; njet < 2; njet++) {
            sample_wgamma->add(Form("%s/mva/%i/wgammafo_%ij.root", dataDir, int(analysis), njet));
            sample_wgamma->add(Form("%s/mva/%i/zgammafo_%ij.root", dataDir, int(analysis), njet));
        }
        sample_wgamma->add(Form("%s/mva/%i/wgammafo_%ij.root", dataDir, 125, 2));
        sample_wgamma->add(Form("%s/mva/%i/zgammafo_%ij.root", dataDir, 125, 2));

        // Wgamma Normalization
        for (int njet=0; njet < 2; njet++) {
            sample_wgammanorm->add(Form("%s/mva/%i/wgamma_%ij.root", dataDir, int(analysis), njet));
            sample_wgammanorm->add(Form("%s/mva/%i/zgamma_%ij.root", dataDir, int(analysis), njet));
        }
        sample_wgammanorm->add(Form("%s/mva/%i/wgamma_%ij.root", dataDir, 125, 2));
        sample_wgammanorm->add(Form("%s/mva/%i/zgamma_%ij.root", dataDir, 125, 2));

        // Wg3l ( Wgamma* )
        for (int njet=0; njet < 2; njet++) {
            sample_wg3l->add(Form("%s/mva/%i/wglll_%ij.root", dataDir, int(analysis), njet));
        }
        sample_wg3l->add(Form("%s/mva/%i/wglll_%ij.root", dataDir, 125, 2));

        // Ztt
        for (int njet=0; njet < 2; njet++) {
            sample_ztt->add(Form("%s/mva/%i/data_ztt_%ij.root", dataDir, int(analysis), njet));
        }
        sample_ztt->add(Form("%s/mva/%i/data_ztt_%ij.root", dataDir, 125, 2));
        

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
            sample_wjetsEle->add(Form("%s/mva/%i/data_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/qqww_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/ggww_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/ttbar_powheg_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/tw_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/wz_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/zz_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/wgamma_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/zgamma_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/wglll_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/dyll_PassFail_%ij.root", dataDir, 125, 2));
        
            sample_wjetsMu->add(Form("%s/mva/%i/data_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/qqww_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/ggww_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/ttbar_powheg_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/tw_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/wz_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/zz_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/wgamma_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/zgamma_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/wglll_PassFail_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/dyll_PassFail_%ij.root", dataDir, 125, 2));
        
      sample_wjetsEle->addShapeVariation2D(WJETSELESHAPEVAR, "histo_WjetsE_CMS_hww_MVAWEBounding", true);
      sample_wjetsMu->addShapeVariation2D(WJETSMUSHAPEVAR, "histo_WjetsM_CMS_hww_MVAWMBounding", true);
        
        // data
        if ( !blind) {
            for (int njet=0; njet < 2; njet++) 
                sample_data->add(Form("%s/mva/%i/data_%ij.root", dataDir, int(analysis), njet));
            sample_data->add(Form("%s/mva/%i/data_%ij.root", dataDir, 125, 2));
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
        }
            //sample_wjets_mc->add(Form("%s/mva/%i/wjets_%ij.root", dataDir, 125, 2));
            //sample_wjets_mc_loose->add(Form("%s/mva/%i/wjets_PassFail_%ij.root", dataDir, 125, 2));
            sample_top_var->add(Form("%s/mva/%i/ttbar_%ij.root", dataDir, 125, 2));
            sample_top_var->add(Form("%s/mva/%i/tw_%ij.root", dataDir, 125, 2));
            sample_wwmcnlo->add(Form("%s/mva/%i/wwmcnlo_%ij.root", dataDir, 125, 2));
            sample_wwmcnlo_up->add(Form("%s/mva/%i/wwmcnloup_%ij.root", dataDir, 125, 2));
            sample_wwmcnlo_down->add(Form("%s/mva/%i/wwmcnlodown_%ij.root", dataDir, 125, 2));
            //sample_dyll_data->add(Form("%s/mva/%i/data_%ij.root", dataDir, 125, 2));
            //sample_dyll_data->add(Form("%s/mva/%i/wz_%ij.root", dataDir, 125, 2));
            //sample_dyll_data->add(Form("%s/mva/%i/zz_%ij.root", dataDir, 125, 2));
    
    
    } 

    //  
    // SS with skim + weights 
    //  
    if ( skimwithmvass) {
        char *dataDir = "/smurf/cerati/skims/Run2012_Summer12_SmurfV9_53X-wwxsecfull8tev/WW/"; // TAS 
        //char *dataDir  = "/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/WW/"; // UAF

        if ( analysis > 0. ) {
            // higgs
            for (int njet=0; njet < 3; njet++) 
                sample_gghww->add(Form("%s/mva/%i/hww%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
            sample_gghww->addShapeVariation2D( QCDSCALEVAR,  "histo_ggH_CMS_hww_MVAggHBounding",     true);
            sample_gghww->addShapeVariation2D( LEPEFFVAR,    "histo_ggH_CMS_hww_MVALepEffBounding",  false);
            sample_gghww->addShapeVariation2D( LEPRESVAR,    "histo_ggH_CMS_hww_MVALepResBounding",  true);
            sample_gghww->addShapeVariation2D( METVAR,       "histo_ggH_CMS_hww_MVAMETResBounding",  true);
            sample_gghww->addShapeVariation2D( JETRESVAR,    "histo_ggH_CMS_hww_MVAJESBounding",     false);

            for (int njet=0; njet < 3; njet++) 
                sample_qqhww->add(Form("%s/mva/%i/hww%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
            sample_qqhww->addShapeVariation2D( LEPEFFVAR,   "histo_qqH_CMS_hww_MVALepEffBounding",  false);
            sample_qqhww->addShapeVariation2D( LEPRESVAR,   "histo_qqH_CMS_hww_MVALepResBounding",  true);
            sample_qqhww->addShapeVariation2D( METVAR,      "histo_qqH_CMS_hww_MVAMETResBounding",  true);
            sample_qqhww->addShapeVariation2D( JETRESVAR,   "histo_qqH_CMS_hww_MVAJESBounding",     false);

            for (int njet=0; njet < 3; njet++) 
                sample_zhww->add(Form("%s/mva/%i/hww%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
            sample_zhww->addShapeVariation2D( LEPEFFVAR,    "histo_ZH_CMS_hww_MVALepEffBounding",   false);
            sample_zhww->addShapeVariation2D( LEPRESVAR,    "histo_ZH_CMS_hww_MVALepResBounding",   true);
            sample_zhww->addShapeVariation2D( METVAR,       "histo_ZH_CMS_hww_MVAMETResBounding",   true);
            sample_zhww->addShapeVariation2D( JETRESVAR,    "histo_ZH_CMS_hww_MVAJESBounding",      false);

            for (int njet=0; njet < 3; njet++) 
                sample_whww->add(Form("%s/mva/%i/hww%i_%ij.root", dataDir, int(analysis), int(analysis), njet));
            sample_whww->addShapeVariation2D( LEPEFFVAR,    "histo_WH_CMS_hww_MVALepEffBounding",   false);
            sample_whww->addShapeVariation2D( LEPRESVAR,    "histo_WH_CMS_hww_MVALepResBounding",   true);
            sample_whww->addShapeVariation2D( METVAR,       "histo_WH_CMS_hww_MVAMETResBounding",   true);
            sample_whww->addShapeVariation2D( JETRESVAR,    "histo_WH_CMS_hww_MVAJESBounding",      false);
        }

        // bkgs
        // qqWW
        for (int njet=0; njet < 2; njet++) 
            sample_qqww->add(Form("%s/mva/%i/qqww_SS_%ij.root", dataDir, int(analysis), njet));
        sample_qqww->add(Form("%s/mva/%i/qqww_SS_%ij.root", dataDir, 125, 2));


        if ( analysis < 200) 
            sample_qqww->addShapeVariation2D(LEPEFFVAR, "histo_qqWW_CMS_hww_MVALepEffBounding", true);
        else 
            sample_qqww->addShapeVariation2D(LEPEFFVAR, "histo_qqWW_CMS_hww_MVALepEffBounding", false);
        sample_qqww->addShapeVariation2D(LEPRESVAR, "histo_qqWW_CMS_hww_MVALepResBounding", true);
        sample_qqww->addShapeVariation2D(METVAR, "histo_qqWW_CMS_hww_MVAMETResBounding", true);
        sample_qqww->addShapeVariation2D(JETRESVAR, "histo_qqWW_CMS_hww_MVAJESBounding", false);

        // ggWW
        for (int njet=0; njet < 2; njet++) 
            sample_ggww->add(Form("%s/mva/%i/ggww_SS_%ij.root", dataDir, int(analysis), njet));
        sample_ggww->add(Form("%s/mva/%i/ggww_SS_%ij.root", dataDir, 125, 2));

        if ( analysis < 200) 
            sample_ggww->addShapeVariation2D(LEPEFFVAR, "histo_ggWW_CMS_hww_MVALepEffBounding", true);
        else 
            sample_ggww->addShapeVariation2D(LEPEFFVAR, "histo_ggWW_CMS_hww_MVALepEffBounding", false);
        sample_ggww->addShapeVariation2D(LEPRESVAR, "histo_ggWW_CMS_hww_MVALepResBounding", true);
        sample_ggww->addShapeVariation2D(METVAR, "histo_ggWW_CMS_hww_MVAMETResBounding", true);
        sample_ggww->addShapeVariation2D(JETRESVAR, "histo_ggWW_CMS_hww_MVAJESBounding", false);

        // Drell-Yan
        for (int njet=0; njet < 2; njet++) {
            sample_dyll->add(Form("%s/mva/%i/dyll_SS_%ij.root", dataDir, int(analysis), njet));
            //sample_dyll_loosemet->add(Form("%s/mva/%i/dyll_SS_%ij.root", dataDir, int(analysis), njet));
            //sample_zjets->add(Form("%s/mva/%i/dyll_SS_%ij.root", dataDir, int(analysis), njet));
        }
        sample_dyll->add(Form("%s/mva/%i/dyll_SS_%ij.root", dataDir, 125, 2));
        //sample_dyll_loosemet->add(Form("%s/mva/%i/dyll_SS_%ij.root", dataDir, 125, 2));
        //sample_zjets->add(Form("%s/mva/%i/dyll_SS_%ij.root", dataDir, 125, 2));

        // Top 
        for (int njet=0; njet < 2; njet++) {
            sample_top->add(Form("%s/mva/%i/ttbar_powheg_SS_%ij.root", dataDir, int(analysis), njet));
            sample_top->add(Form("%s/mva/%i/tw_SS_%ij.root", dataDir, int(analysis), njet));
        }
        sample_top->add(Form("%s/mva/%i/ttbar_powheg_SS_%ij.root", dataDir, 125, 2));
        sample_top->add(Form("%s/mva/%i/tw_SS_%ij.root", dataDir, 125, 2));

        sample_top->addShapeVariation2D(LEPRESVAR, "histo_Top_CMS_hww_MVALepResBounding", true);
        sample_top->addShapeVariation2D(METVAR, "histo_Top_CMS_hww_MVAMETResBounding", true);
        sample_top->addShapeVariation2D(JETRESVAR, "histo_Top_CMS_hww_MVAJESBounding", false);

        // VV 
        for (int njet=0; njet < 2; njet++) {
            sample_vv->add(Form("%s/mva/%i/wz_SS_%ij.root", dataDir, int(analysis), njet));
            sample_vv->add(Form("%s/mva/%i/zz_SS_%ij.root", dataDir, int(analysis), njet)); 
            sample_vv->add(Form("%s/mva/%i/www_SS_%ij.root", dataDir, int(analysis), njet)); 
        }
        sample_vv->add(Form("%s/mva/%i/wz_SS_%ij.root", dataDir, 125, 2));
        sample_vv->add(Form("%s/mva/%i/zz_SS_%ij.root", dataDir, 125, 2)); 
        sample_vv->add(Form("%s/mva/%i/www_SS_%ij.root", dataDir, 125, 2)); 

        sample_vv->addShapeVariation2D(LEPEFFVAR, "histo_VV_CMS_hww_MVALepEffBounding", false);
        sample_vv->addShapeVariation2D(LEPRESVAR, "histo_VV_CMS_hww_MVALepResBounding", true);
        sample_vv->addShapeVariation2D(METVAR, "histo_VV_CMS_hww_MVAMETResBounding", true);
        sample_vv->addShapeVariation2D(JETRESVAR, "histo_VV_CMS_hww_MVAJESBounding", false);

        // Wgamma : l + gamma sample 
        for (int njet=0; njet < 2; njet++) {
            sample_wgamma->add(Form("%s/mva/%i/wgammafo_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wgamma->add(Form("%s/mva/%i/zgammafo_SS_%ij.root", dataDir, int(analysis), njet));
        }
        sample_wgamma->add(Form("%s/mva/%i/wgammafo_SS_%ij.root", dataDir, 125, 2));
        sample_wgamma->add(Form("%s/mva/%i/zgammafo_SS_%ij.root", dataDir, 125, 2));

        // Wgamma Normalization
        for (int njet=0; njet < 2; njet++) {
            sample_wgammanorm->add(Form("%s/mva/%i/wgamma_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wgammanorm->add(Form("%s/mva/%i/zgamma_SS_%ij.root", dataDir, int(analysis), njet));
        }
        sample_wgammanorm->add(Form("%s/mva/%i/wgamma_SS_%ij.root", dataDir, 125, 2));
        sample_wgammanorm->add(Form("%s/mva/%i/zgamma_SS_%ij.root", dataDir, 125, 2));

        // Wg3l ( Wgamma* )
        for (int njet=0; njet < 2; njet++) {
            sample_wg3l->add(Form("%s/mva/%i/wglll_SS_%ij.root", dataDir, int(analysis), njet));
        }
        sample_wg3l->add(Form("%s/mva/%i/wglll_SS_%ij.root", dataDir, 125, 2));

        // Ztt
        for (int njet=0; njet < 2; njet++) {
            sample_ztt->add(Form("%s/mva/%i/data_ztt_SS_%ij.root", dataDir, int(analysis), njet));
        }
        sample_ztt->add(Form("%s/mva/%i/data_ztt_SS_%ij.root", dataDir, 125, 2));
        
        // Wjets 
        for (int njet=0; njet < 2; njet++) {
            sample_wjetsEle->add(Form("%s/mva/%i/data_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/qqww_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/ggww_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/ttbar_powheg_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/tw_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/wz_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/zz_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/wgamma_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/zgamma_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/wglll_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsEle->add(Form("%s/mva/%i/dyll_SS_%ij.root", dataDir, int(analysis), njet));
        
            sample_wjetsMu->add(Form("%s/mva/%i/data_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/qqww_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/ggww_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/ttbar_powheg_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/tw_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/wz_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/zz_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/wgamma_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/zgamma_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/wglll_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wjetsMu->add(Form("%s/mva/%i/dyll_SS_%ij.root", dataDir, int(analysis), njet));
    }
            sample_wjetsEle->add(Form("%s/mva/%i/data_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/qqww_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/ggww_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/ttbar_powheg_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/tw_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/wz_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/zz_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/wgamma_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/zgamma_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/wglll_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsEle->add(Form("%s/mva/%i/dyll_SS_%ij.root", dataDir, 125, 2));
        
            sample_wjetsMu->add(Form("%s/mva/%i/data_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/qqww_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/ggww_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/ttbar_powheg_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/tw_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/wz_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/zz_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/wgamma_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/zgamma_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/wglll_SS_%ij.root", dataDir, 125, 2));
            sample_wjetsMu->add(Form("%s/mva/%i/dyll_SS_%ij.root", dataDir, 125, 2));
       
      sample_wjetsEle->addShapeVariation2D(WJETSELESHAPEVAR, "histo_WjetsE_CMS_hww_MVAWEBounding", true);
      sample_wjetsMu->addShapeVariation2D(WJETSMUSHAPEVAR, "histo_WjetsM_CMS_hww_MVAWMBounding", true);

        // data
        if ( !blind) {
            for (int njet=0; njet < 2; njet++) 
                sample_data->add(Form("%s/mva/%i/data_SS_%ij.root", dataDir, int(analysis), njet));
            sample_data->add(Form("%s/mva/%i/data_SS_%ij.root", dataDir, 125, 2));
        }

        // for Shape variations
        // wjets
        for (int njet=0; njet < 2; njet++) {
            //sample_wjets_mc->add(Form("%s/mva/%i/wjets_%ij.root", dataDir, int(analysis), njet));
            //sample_wjets_mc_loose->add(Form("%s/mva/%i/wjets_PassFail_%ij.root", dataDir, int(analysis), njet));
            sample_top_var->add(Form("%s/mva/%i/ttbar_SS_%ij.root", dataDir, int(analysis), njet));
            sample_top_var->add(Form("%s/mva/%i/tw_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wwmcnlo->add(Form("%s/mva/%i/wwmcnlo_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wwmcnlo_up->add(Form("%s/mva/%i/wwmcnloup_SS_%ij.root", dataDir, int(analysis), njet));
            sample_wwmcnlo_down->add(Form("%s/mva/%i/wwmcnlodown_SS_%ij.root", dataDir, int(analysis), njet));
            //sample_dyll_data->add(Form("%s/mva/%i/data_SS_%ij.root", dataDir, int(analysis), njet));
            //sample_dyll_data->add(Form("%s/mva/%i/wz_SS_%ij.root", dataDir, int(analysis), njet));
            //sample_dyll_data->add(Form("%s/mva/%i/zz_SS_%ij.root", dataDir, int(analysis), njet));
        }
            //sample_wjets_mc->add(Form("%s/mva/%i/wjets_%ij.root", dataDir, 125, 2));
            //sample_wjets_mc_loose->add(Form("%s/mva/%i/wjets_PassFail_%ij.root", dataDir, 125, 2));
            sample_top_var->     add(Form("%s/mva/%i/ttbar_SS_%ij.root",       dataDir, 125, 2));
            sample_top_var->     add(Form("%s/mva/%i/tw_SS_%ij.root",          dataDir, 125, 2));
            sample_wwmcnlo->     add(Form("%s/mva/%i/wwmcnlo_SS_%ij.root",     dataDir, 125, 2));
            sample_wwmcnlo_up->  add(Form("%s/mva/%i/wwmcnloup_SS_%ij.root",   dataDir, 125, 2));
            sample_wwmcnlo_down->add(Form("%s/mva/%i/wwmcnlodown_SS_%ij.root", dataDir, 125, 2));
            //sample_dyll_data->   add(Form("%s/mva/%i/data_SS_%ij.root",        dataDir, 125, 2));
            //sample_dyll_data->   add(Form("%s/mva/%i/wz_SS_%ij.root",          dataDir, 125, 2));
            //sample_dyll_data->   add(Form("%s/mva/%i/zz_SS_%ij.root",          dataDir, 125, 2));

    }

    // 
    // do the looping
    //

    looper->loop(sample_data);
    if ( analysis > 0. ) {
        looper->loop(sample_gghww);
        looper->loop(sample_qqhww);
        looper->loop(sample_whww);
        looper->loop(sample_zhww);
    }
    looper->loop(sample_qqww);
    looper->loop(sample_ggww);
    looper->loop(sample_vv);
    looper->loop(sample_top);
    looper->loop(sample_dyll);
    //looper->loop(sample_dyll_loosemet);
    //looper->loop(sample_zjets);
    looper->loop(sample_wjetsEle);
    looper->loop(sample_wjetsMu);
    looper->loop(sample_wgamma);
    looper->loop(sample_wgammanorm);
    looper->loop(sample_wg3l);
    looper->loop(sample_ztt);

    // for shape syst.
    //looper->loop(sample_wjets_mc_loose);
    //looper->loop(sample_wjets_mc);
    looper->loop(sample_top_var);
    looper->loop(sample_wwmcnlo);
    looper->loop(sample_wwmcnlo_up);
    looper->loop(sample_wwmcnlo_down);
    //looper->loop(sample_dyll_data);


    //
    // make tables
    //

    std::vector<SmurfSample*> samplesToTabulate;

    samplesToTabulate.push_back(sample_data);

    if ( analysis > 0. ) {
        samplesToTabulate.push_back(sample_gghww);
        samplesToTabulate.push_back(sample_qqhww);
        samplesToTabulate.push_back(sample_whww);
        samplesToTabulate.push_back(sample_zhww);
    }

    samplesToTabulate.push_back(sample_qqww);
    samplesToTabulate.push_back(sample_ggww);
    samplesToTabulate.push_back(sample_vv);
    samplesToTabulate.push_back(sample_top);
    samplesToTabulate.push_back(sample_dyll);
    //samplesToTabulate.push_back(sample_zjets);
    samplesToTabulate.push_back(sample_wjetsEle);
    samplesToTabulate.push_back(sample_wjetsMu);
    samplesToTabulate.push_back(sample_wgamma);
    samplesToTabulate.push_back(sample_wgammanorm);
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
    //samplesToTabulate.push_back(sample_dyll_data);

    //  
    // make cards
    //
    printf("\n\n[doAllHWW::doMassPoint] Writing ards\n");

    const std::string cardDir = "../cards/test/";

    // full list of shape variations for BDT analysis 
        ShapeVar_t mva_option = (1ll<<STATVAR) | (1ll<<TOPSHAPEVAR) | (1ll<<WWSHAPEVAR) | (1ll<<DYSHAPEVAR)
            | (1ll<<LEPEFFVAR) | (1ll<<METVAR) | (1ll<<LEPRESVAR) | (1ll<<JETRESVAR) | (1ll<<WJETSELESHAPEVAR) | (1ll<<WJETSMUSHAPEVAR);


    for (int jetbin = 0; jetbin < 2; ++jetbin) { 
        
/*
        // same flavor
        printCard(samplesToTabulate, option, jetbin, analysis, cardDir, ((1<<0)|(1<<3)), mva_option, runEra);
        print2DShapeHistograms(samplesToTabulate, option, jetbin, analysis, cardDir,((1<<0)|(1<<3)), runEra);
*/
        // opposite flavor
        printCard(samplesToTabulate, option, jetbin, analysis, cardDir, ((1<<1)|(1<<2)), mva_option, runEra);
        print2DShapeHistograms(samplesToTabulate, option, jetbin, analysis, cardDir, ((1<<1)|(1<<2)), runEra);

/*
       // me 
        printCard(samplesToTabulate, option, jetbin, analysis, cardDir, (1<<1), mva_option, runEra);
        print2DShapeHistograms(samplesToTabulate, option, jetbin, analysis, cardDir, (1<<1), runEra);
        // em 
        printCard(samplesToTabulate, option, jetbin, analysis, cardDir, (1<<2), mva_option, runEra);
        print2DShapeHistograms(samplesToTabulate, option, jetbin, analysis, cardDir, (1<<2), runEra);
        // mm 
        printCard(samplesToTabulate, option, jetbin, analysis, cardDir, (1<<0), mva_option, runEra);
        print2DShapeHistograms(samplesToTabulate, option, jetbin, analysis, cardDir, (1<<0), runEra);
        // ee 
        printCard(samplesToTabulate, option, jetbin, analysis, cardDir, (1<<3), mva_option, runEra);
        print2DShapeHistograms(samplesToTabulate, option, jetbin, analysis, cardDir, (1<<3), runEra);
*/
    }


    //
    // save histograms  
    //

    const std::string outFile = Form("histos_hww_analysis%i_%i.root", option, int(analysis));
    saveHist(outFile.c_str());  
    deleteHistos();

    //
    // make plots
    //


    if (makePlots_) {

        const unsigned int kLeptonTypes = 7;
        const unsigned int kJetBins = 2;
        const char*          types[7]                = {"mm","me","em","ee", "sf", "of", "incl"};
        const char*          jetbin_names[3]         = { "0j", "1j", "2j"};

        gROOT->ProcessLine(".L tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");
        gROOT->ForceStyle();

        const char *dir = "../hwwplots";
        gSystem->mkdir(plotDir, true);
        TFile *f = new TFile(Form("histos_hww_analysis%i_%i.root", option, int(analysis)), "READ");
        gROOT->cd();

        for (unsigned int j = 0; j < kJetBins; ++j) {
            for (unsigned int i = 0; i < kLeptonTypes; ++i) {

                //if (! (i == fTOTAL || i == fEE || i == fEM || i == fME || i == fMM)) continue;
                if (i < 4 && i!=0) continue;

                //
                // all bdt output
                //

                TCanvas *c_pt1 = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, plotDir, "ww_pt1", "p_{T} (leading lepton) [GeV/c]", lumi);
                TCanvas *c_pt2 = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, plotDir, "ww_pt2", "p_{T} (trailing lepton) [GeV/c]", lumi);
                TCanvas *c_eta1 = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, plotDir, "ww_eta1", "|#eta| (leading lepton) [GeV/c]", lumi);
                TCanvas *c_eta2 = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, plotDir, "ww_eta2", "|#eta| (trailing lepton) [GeV/c]", lumi);
                TCanvas *c_met = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, plotDir, "ww_met", "MET [GeV]", lumi);
                TCanvas *c_mt = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, plotDir, "ww_mt", "MT [GeV]", lumi);
                TCanvas *c_mll = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, plotDir, "ww_mll", "M_{l, l} [GeV/c^{2}]", lumi);
                TCanvas *c_dphi = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, plotDir, "ww_dphi", "#Delta#phi (dilepton)", lumi);
                
                c_pt1->SaveAs(Form("%s/individual/hww_analysis%i_%i_ALL_%s_%s_pt1.pdf", plotDir, option, int(analysis), types[i], jetbin_names[j]));
                c_pt2->SaveAs(Form("%s/individual/hww_analysis%i_%i_ALL_%s_%s_pt2.pdf", plotDir, option, int(analysis), types[i], jetbin_names[j]));
                c_eta1->SaveAs(Form("%s/individual/hww_analysis%i_%i_ALL_%s_%s_eta1.pdf", plotDir, option, int(analysis), types[i], jetbin_names[j]));
                c_eta2->SaveAs(Form("%s/individual/hww_analysis%i_%i_ALL_%s_%s_eta2.pdf", plotDir, option, int(analysis), types[i], jetbin_names[j]));
                c_dphi->SaveAs(Form("%s/individual/hww_analysis%i_%i_ALL_%s_%s_dphi.pdf", plotDir, option, int(analysis), types[i], jetbin_names[j]));
                c_mll->SaveAs(Form("%s/individual/hww_analysis%i_%i_ALL_%s_%s_mll.pdf", plotDir, option, int(analysis), types[i], jetbin_names[j]));
                c_mt->SaveAs(Form("%s/individual/hww_analysis%i_%i_ALL_%s_%s_mt.pdf", plotDir, option, int(analysis), types[i], jetbin_names[j]));
                c_met->SaveAs(Form("%s/individual/hww_analysis%i_%i_ALL_%s_%s_met.pdf", plotDir, option, int(analysis), types[i], jetbin_names[j]));

                TCanvas *c1_big = new TCanvas("c1_big", "c1_big", 1600, 800);
                c1_big->Divide(4, 2);
                c1_big->cd(1);
                (TPad*)c_pt1->GetListOfPrimitives()->FindObject("p_main")->Draw();
                c1_big->cd(2);
                (TPad*)c_pt2->GetListOfPrimitives()->FindObject("p_main")->Draw(); 
                c1_big->cd(3);
                (TPad*)c_eta1->GetListOfPrimitives()->FindObject("p_main")->Draw();
                c1_big->cd(4);
                (TPad*)c_eta2->GetListOfPrimitives()->FindObject("p_main")->Draw();
                c1_big->cd(5);
                (TPad*)c_dphi->GetListOfPrimitives()->FindObject("p_main")->Draw();
                c1_big->cd(6);
                (TPad*)c_mll->GetListOfPrimitives()->FindObject("p_main")->Draw();
                c1_big->cd(7);
                (TPad*)c_mt->GetListOfPrimitives()->FindObject("p_main")->Draw();
                c1_big->cd(8);
                (TPad*)c_met->GetListOfPrimitives()->FindObject("p_main")->Draw();
                //c1_big->SaveAs(Form("%s/hww_analysis%i_%i_ALL_%s_%s.eps", plotDir, option, int(analysis), types[i], jetbin_names[j]));
                //c1_big->SaveAs(Form("%s/hww_analysis%i_%i_ALL_%s_%s.png", plotDir, option, int(analysis), types[i], jetbin_names[j]));
                c1_big->SaveAs(Form("%s/hww_analysis%i_%i_ALL_%s_%s.pdf", plotDir, option, int(analysis), types[i], jetbin_names[j]));

                //delete c_bdt; 
                delete c_pt1; delete c_pt2; delete c_eta1; delete c_eta2; 
                delete c_mll; delete c_mt; delete c_met; delete c_dphi; 
                delete c1_big;
/*
                //
                // LOW SCORE
                //

                TH1F *h1_dy_bdt_tmp = GetHistogram(f, sample_dyll_loosemet, Form("ww_bdt_%s_%s", jetbin_names[j], types[i]));
                float fLowDY = h1_dy_bdt_tmp->Integral(0, h1_dy_bdt_tmp->FindBin(-0.4))
                    / h1_dy_bdt_tmp->Integral(0, 999);
                delete h1_dy_bdt_tmp;

                TCanvas *c_bdt_lo = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLLLOOSEMET,
                        f, i, j, "../hwwplots/", "ww_bdt", "BDT Output", lumi);
                TCanvas *c_pt1_lo = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, "../hwwplots/", "hww_bdtlo_pt1", "p_{T} (leading lepton) [GeV/c]", lumi, fLowDY);
                TCanvas *c_pt2_lo = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, "../hwwplots/", "hww_bdtlo_pt2", "p_{T} (trailing lepton) [GeV/c]", lumi, fLowDY);
                TCanvas *c_mt_lo = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, "../hwwplots/", "hww_bdtlo_mt", "MT [GeV]", lumi, fLowDY);
                TCanvas *c_mll_lo = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, "../hwwplots/", "hww_bdtlo_mll", "M_{l, l} [GeV/c^{2}]", lumi, fLowDY);
                TCanvas *c_dphi_lo = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, "../hwwplots/", "hww_bdtlo_dphi", "#Delta#phi (dilepton)", lumi, fLowDY);
                TCanvas *cbig_bdt_lo = new TCanvas("cbig_bdt_lo", "cbig_bdt_lo", 800, 1200);
                cbig_bdt_lo->Divide(2, 3);
                cbig_bdt_lo->cd(1);
                (TPad*)c_pt1_lo->GetListOfPrimitives()->FindObject("p_main")->Draw();
                cbig_bdt_lo->cd(2);
                (TPad*)c_pt2_lo->GetListOfPrimitives()->FindObject("p_main")->Draw();
                cbig_bdt_lo->cd(3);
                (TPad*)c_dphi_lo->GetListOfPrimitives()->FindObject("p_main")->Draw();
                cbig_bdt_lo->cd(4);
                (TPad*)c_mll_lo->GetListOfPrimitives()->FindObject("p_main")->Draw();
                cbig_bdt_lo->cd(5);
                (TPad*)c_mt_lo->GetListOfPrimitives()->FindObject("p_main")->Draw();
                cbig_bdt_lo->cd(6);
                (TPad*)c_bdt_lo->GetListOfPrimitives()->FindObject("p_main")->Draw();
                cbig_bdt_lo->SaveAs(Form("../hwwplots/hww_bdtlo_analysis%i_%i_ALL_%s_%s.png", option, int(analysis), types[i], jetbin_names[j]));
                cbig_bdt_lo->SaveAs(Form("../hwwplots/hww_bdtlo_analysis%i_%i_ALL_%s_%s.eps", option, int(analysis), types[i], jetbin_names[j]));

                delete c_bdt_lo; delete c_pt1_lo; delete c_pt2_lo; delete c_mt_lo;
                delete c_mll_lo; delete c_dphi_lo; delete cbig_bdt_lo;

                //
                // HIGH SCORE
                //

                TCanvas *c_bdt_hi = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLLLOOSEMET,
                        f, i, j, "../hwwplots/", "ww_bdt", "BDT Output", lumi);
                TCanvas *c_pt1_hi = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, "../hwwplots/", "hww_bdthi_pt1", "p_{T} (leading lepton) [GeV/c]", lumi, (1.0 - fLowDY));
                TCanvas *c_pt2_hi = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, "../hwwplots/", "hww_bdthi_pt2", "p_{T} (trailing lepton) [GeV/c]", lumi, (1.0 - fLowDY));
                TCanvas *c_mt_hi = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, "../hwwplots/", "hww_bdthi_mt", "MT [GeV]", lumi, (1.0 - fLowDY));
                TCanvas *c_mll_hi = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, "../hwwplots/", "hww_bdthi_mll", "M_{l, l} [GeV/c^{2}]", lumi, (1.0 - fLowDY));
                TCanvas *c_dphi_hi = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                        f, i, j, "../hwwplots/", "hww_bdthi_dphi", "#Delta#phi (dilepton)", lumi, (1.0 - fLowDY));
                TCanvas *cbig_bdt_hi = new TCanvas("cbig_bdt_hi", "cbig_bdt_hi", 800, 1200);
                cbig_bdt_hi->Divide(2, 3);
                cbig_bdt_hi->cd(1);
                (TPad*)c_pt1_hi->GetListOfPrimitives()->FindObject("p_main")->Draw();
                cbig_bdt_hi->cd(2);
                (TPad*)c_pt2_hi->GetListOfPrimitives()->FindObject("p_main")->Draw();
                cbig_bdt_hi->cd(3);
                (TPad*)c_dphi_hi->GetListOfPrimitives()->FindObject("p_main")->Draw();
                cbig_bdt_hi->cd(4);
                (TPad*)c_mll_hi->GetListOfPrimitives()->FindObject("p_main")->Draw();
                cbig_bdt_hi->cd(5);
                (TPad*)c_mt_hi->GetListOfPrimitives()->FindObject("p_main")->Draw();
                cbig_bdt_hi->cd(6);
                (TPad*)c_bdt_hi->GetListOfPrimitives()->FindObject("p_main")->Draw();
                cbig_bdt_hi->SaveAs(Form("../hwwplots/hww_bdthi_analysis%i_%i_ALL_%s_%s.png", option, int(analysis), types[i], jetbin_names[j]));
                cbig_bdt_hi->SaveAs(Form("../hwwplots/hww_bdthi_analysis%i_%i_ALL_%s_%s.eps", option, int(analysis), types[i], jetbin_names[j]));

                delete c_bdt_hi; delete c_pt1_hi; delete c_pt2_hi; delete c_mt_hi;
                delete c_mll_hi; delete c_dphi_hi; delete cbig_bdt_hi;
*/
            }
        }
/*
        for (unsigned int i = 0; i < kLeptonTypes; ++i) {
            if (i != fTOTAL) continue;

            TCanvas *c_detajj = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                    f, i, 2, "../hwwplots/", "ww_detajj", "#Delta#eta_{jet1, jet2}", lumi);
            TCanvas *c_mjj = makeHWWAnalysisStack(option, analysis, samplesToTabulate, ZLL,
                    f, i, 2, "../hwwplots/", "ww_mjj", "M_{jet1, jet2}", lumi);

            TCanvas *c1_big = new TCanvas();
            c1_big->Divide(2, 1);
            c1_big->cd(1);
            (TPad*)c_detajj->GetListOfPrimitives()->FindObject("p_main")->Draw();
            c1_big->cd(2);    
            (TPad*)c_mjj->GetListOfPrimitives()->FindObject("p_main")->Draw();

            c1_big->SaveAs(Form("../hwwplots/hww_analysis%i_%i_ALLjj_%s_%s.eps", option, int(analysis), types[i], jetbin_names[2]));
            c1_big->SaveAs(Form("../hwwplots/hww_analysis%i_%i_ALLjj_%s_%s.png", option, int(analysis), types[i], jetbin_names[2]));

            delete c_detajj; delete c_mjj; delete c1_big;

        }
*/

        // tidy up...
        deleteHistos();

    }

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

}


