
#include "SmurfLooper.h"

#include "../Smurf/Core/SmurfTree.h"
#include "../Smurf/Core/LeptonScaleLookup.h"
#include "../Tools/goodrun.h"
#include "../Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"
#include "../Smurf/Analysis/HWWlvlv/factors.h"

#include "core/Selections.h"

#include "core/SmurfPlotUtilities.h"
#include "core/SmurfSample.h"
#include "SmurfScaleFactors.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom.h"

#include <cmath>
#include <iostream>
#include <cassert>
#include <algorithm>

SmurfLooper::SmurfLooper(float analysis, Option option, RunEra runEra) 
{

    std::cout << std::endl;
    std::cout << "[SmurfLooper::SmurfLooper] Doing mass point " << analysis << std::endl;
    std::cout << "[SmurfLooper::SmurfLooper] Doing analysis option " << option << std::endl;

    analysis_ = analysis;
    runEra_ = runEra;
    option_ = option;
    loadWeightHistograms();
    debugtext_ = fopen("debug.txt", "w");
    runlistIsSet_ = false;
    gRandom->SetSeed(48151623);

    // set default lumi
    eeLumi_ = 2121;
    mmLumi_ = 2121;

    std::cout << "[SmurfLooper::SmurfLooper] default lumi is " << eeLumi_ << ", " << mmLumi_ << std::endl;
    std::cout << "[SmurfLooper::SmurfLooper] random seed is " << gRandom->GetSeed() << std::endl;

}

SmurfLooper::~SmurfLooper()
{
}

void SmurfLooper::loop(SmurfSample *sample)
{

    // ----------------------------------------------------------------------------------------
    //      Set up histograms
    // ----------------------------------------------------------------------------------------

    gROOT->cd();

    TH1F *h1_ww_mll[kJetBins][kLeptonTypes];
    TH1F *h1_ww_pt1[kJetBins][kLeptonTypes];
    TH1F *h1_ww_pt2[kJetBins][kLeptonTypes];
    TH1F *h1_ww_eta1[kJetBins][kLeptonTypes];
    TH1F *h1_ww_eta2[kJetBins][kLeptonTypes];
    TH1F *h1_ww_met[kJetBins][kLeptonTypes];
    TH1F *h1_ww_ptll[kJetBins][kLeptonTypes];
    TH1F *h1_ww_ptWW[kJetBins][kLeptonTypes];
    TH1F *h1_ww_mt[kJetBins][kLeptonTypes];
    TH1F *h1_ww_dphi[kJetBins][kLeptonTypes];   
    TH1F *h1_ww_nvtx[kJetBins][kLeptonTypes];   

    for (unsigned int i = 0; i < kJetBins; ++i) {

        FormatHist(h1_ww_mll[i],    sample, Form("ww_mll_%s",   jetbin_names[i]),   "mll",  20,  0.0,   300.0);
        FormatHist(h1_ww_pt1[i],    sample, Form("ww_pt1_%s",   jetbin_names[i]),   "pt1",  20,  0.0,   200.0);
        FormatHist(h1_ww_pt2[i],    sample, Form("ww_pt2_%s",   jetbin_names[i]),   "pt2",  20,  0.0,   160.0);
        FormatHist(h1_ww_eta1[i],   sample, Form("ww_eta1_%s",  jetbin_names[i]),   "eta1",  5,  0.0,   2.5);
        FormatHist(h1_ww_eta2[i],   sample, Form("ww_eta2_%s",  jetbin_names[i]),   "eta2",  5,  0.0,   2.5);
        FormatHist(h1_ww_met[i],    sample, Form("ww_met_%s",   jetbin_names[i]),   "met",  20,  0.0,   200.0);
        FormatHist(h1_ww_ptll[i],   sample, Form("ww_ptll_%s",  jetbin_names[i]),   "ptll", 20,  0.0,   200.0);
        FormatHist(h1_ww_ptWW[i],   sample, Form("ww_ptWW_%s",  jetbin_names[i]),   "ptWW", 20,  0.0,   200.0);
        FormatHist(h1_ww_mt[i],     sample, Form("ww_mt_%s",    jetbin_names[i]),   "mt",   20,  60.0,  280.0);
        FormatHist(h1_ww_dphi[i],   sample, Form("ww_dphi_%s",  jetbin_names[i]),   "dphi", 20,  0.0,   3.2);
        FormatHist(h1_ww_nvtx[i],   sample, Form("ww_nvtx_%s",  jetbin_names[i]),   "nvtx", 40,  0.5,   40.5); 

    } 

    TH1F *h1_ww_detajj[kLeptonTypes];
    TH1F *h1_ww_mjj[kLeptonTypes];
    FormatHist(h1_ww_detajj, sample, "ww_detajj_2j",    "detajj",   20, 0.0, 8);
    FormatHist(h1_ww_mjj,    sample, "ww_mjj_2j",       "mjj",      20, 0.0, 1000.0);

    // ----------------------------------------------------------------------------------------
    //      File loop
    // ----------------------------------------------------------------------------------------

    TObjArray *listOfFiles = sample->getChain()->GetListOfFiles();
    TIter fileIter(listOfFiles);
    if (listOfFiles->GetEntries() == 0) {
        std::cout << "[SmurfLooper::loop] " << sample->getName() << " is not defined" << std::endl;
        return;
    }
    else
        std::cout << "[SmurfLooper::loop] " << sample->getName() << std::endl;

    unsigned int nEventsChain=0;
    unsigned int nEvents = sample->getChain()->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;
    int i_permille_old = 0;

    while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

        SmurfTree *tree = new SmurfTree();
        tree->LoadTree(currentFile->GetTitle());
        tree->InitTree(0);
        
        // ----------------------------------------------------------------------------------------
        //      Load extra variables 
        // ----------------------------------------------------------------------------------------

        //  extra variables for the various correction factors
        float sfWeightPU_ = 1.0;
        float sfWeightTrig_ = 1.0;
        float sfWeightEff_ = 1.0;
        float sfWeightFR_ = 1.0;
        float sfWeightHPt_ = 1.0;

        if (tree->tree_->GetBranchStatus("sfWeightPU"))
            tree->tree_->SetBranchAddress("sfWeightPU", &sfWeightPU_);
        if (tree->tree_->GetBranchStatus("sfWeightTrig"))
            tree->tree_->SetBranchAddress("sfWeightTrig", &sfWeightTrig_);
        if (tree->tree_->GetBranchStatus("sfWeightEff"))
            tree->tree_->SetBranchAddress("sfWeightEff", &sfWeightEff_);
        if (tree->tree_->GetBranchStatus("sfWeightFR"))
            tree->tree_->SetBranchAddress("sfWeightFR", &sfWeightFR_);
        if (tree->tree_->GetBranchStatus("sfWeightHPt"))
            tree->tree_->SetBranchAddress("sfWeightHPt", &sfWeightHPt_);

        // Extra variables for 2D alternative shapes
        float mt_lepup_     = 0.0;
        float mt_lepdown_   = 0.0;
        float mll_lepup_    = 0.0;
        float mll_lepdown_  = 0.0;
        float mt_metup_     = 0.0;
        float mll_metup_    = 0.0;
        if (tree->tree_->GetBranchStatus("mt_lepup"))
            tree->tree_->SetBranchAddress("mt_lepup", &mt_lepup_);
        if (tree->tree_->GetBranchStatus("mt_lepdown"))
            tree->tree_->SetBranchAddress("mt_lepdown", &mt_lepdown_);
        if (tree->tree_->GetBranchStatus("mll_lepup"))
            tree->tree_->SetBranchAddress("mll_lepup", &mll_lepup_);
        if (tree->tree_->GetBranchStatus("mll_lepdown"))
            tree->tree_->SetBranchAddress("mll_lepdown", &mll_lepdown_);
        if (tree->tree_->GetBranchStatus("mt_metup"))
            tree->tree_->SetBranchAddress("mt_metup", &mt_metup_);
        if (tree->tree_->GetBranchStatus("mll_metup"))
            tree->tree_->SetBranchAddress("mll_metup", &mll_metup_);

        // DY MVA : the default value is set to pass the dymva if this branch does not exist
        float dymva_ = 1.; 
        if (tree->tree_->GetBranchStatus("dymva"))
            tree->tree_->SetBranchAddress("dymva", &dymva_);

        // ----------------------------------------------------------------------------------------
        //      Event loop
        // ----------------------------------------------------------------------------------------

        ULong64_t nEvents = tree->tree_->GetEntries();
        for(ULong64_t event = 0; event < nEvents; ++event) {
            tree->tree_->GetEntry(event);

            // ----------------------------------------------------------------------------------------
            //      Incrimenet counters
            // ----------------------------------------------------------------------------------------

            ++nEventsTotal;
            int i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
            if (i_permille != i_permille_old) {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1)) {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                            "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
                    fflush(stdout);
                }
                i_permille_old = i_permille;
            }
            
            // ----------------------------------------------------------------------------------------
            //      Preliminary selections 
            // ----------------------------------------------------------------------------------------
            
            // For WW xsec, use (20, 20)
            if (option_ == WW_OPT_SMURFXSECSEL ) {
                if ( tree->lep2_.Pt() < 20) continue;
            }

            // For 2D analysis
            if (option_ == HWW_OPT_MT2DMLL && tree->njets_ > 3) continue;
            if (option_ == HWW_OPT_SSCTL2D && tree->njets_ > 3) continue;
            if (option_ == HWW_OPT_TOPTAG && tree->njets_ > 3) continue;
            // For Spin-parity analysis, use only 0/1 jet events  
            if (option_ == HWW_OPT_MT2DMLL_JCP || option_ == XWW_OPT_MT2DMLL_JCP ) {
              if ( tree->njets_ > 1 ) continue;
            }

            // ----------------------------------------------------------------------------------------
            //      Separate out different signal sources
            // ----------------------------------------------------------------------------------------

            if ((1ll<<sample->getDataType()) & data_gghiggs) {
                if (tree->processId_ != 10010) continue;
            }
            else if ((1ll<<sample->getDataType()) & data_qqhiggs) {
                if (tree->processId_ != 10001) continue;
            }
            else if (sample->getDataType() == WHWW) {
                if (tree->processId_ != 26) continue;
            }
            else if (sample->getDataType() == ZHWW) {
                if (tree->processId_ != 24) continue;
            }

           
            // ----------------------------------------------------------------------------------------
            //  Make sure events get assigned appropriately :  
            //  This function is defined in core/Selections.cc 
            // ----------------------------------------------------------------------------------------
            if(!hww_assign_this_event(tree, sample->getDataType())) continue;

            // ----------------------------------------------------------------------------------------
            //      Set up weights and binning
            // ----------------------------------------------------------------------------------------

            double weight = 1.0;
            double weight_err = 0.0;
            unsigned int type = tree->type_;
            unsigned int njets = tree->njets_;

            // 3-events are categorized as 2-jet events 
            if (tree->njets_ == 3) njets = 2;

            // ----------------------------------------------------------------------------------------
            //      Apply good run list
            // ----------------------------------------------------------------------------------------

            //if (((1ll<<sample->getDataType()) & data_data) && tree->run_ > 203002) continue; // HCP dataset
            //if (((1ll<<sample->getDataType()) & data_data) && tree->run_ <= 203002) continue; // post-HCP dataset

            if (((1ll<<sample->getDataType()) & data_data) && tree->scale1fb_ == 1.0)  {
                if (runlistIsSet_) {
                    if (!goodrun_json(tree->run_, tree->lumi_)) continue;
                }
            } else {
                if (tree->type_ == 0) weight = tree->scale1fb_*(mmLumi_/1000.0);
                if (tree->type_ == 3) weight = tree->scale1fb_*(eeLumi_/1000.0);
                if (tree->type_ == 1 || tree->type_ == 2) weight = tree->scale1fb_*((eeLumi_+mmLumi_)/(2*1000.0));
            }

            // ----------------------------------------------------------------------------------------
            //      Apply lepton scale factor and PU re-weighting
            // ----------------------------------------------------------------------------------------

            if ( ((1ll<<sample->getDataType()) & (data_allmc | (1ll<<WJETSELEDATA) | (1ll<<WJETSMUDATA) | (1ll<<WJETSDATA) | (1ll<<WJETSMCLOOSE) | (1ll<<ZLLDATA))) 
                    && (sample->getDataType() != ZLLGAMMA && sample->getDataType() != ZVVGAMMA 
                        && sample->getDataType() != WJETSGAMMA)) {
                
                /*
                // recalculate efficiency scale factors : use numbers from the efficiency file 
                double offlineSF1 = leptonSF_->GetExpectedLeptonSF(tree->lep1_.Eta(), tree->lep1_.Pt(), tree->lid1_);
                double offlineSF2 = leptonSF_->GetExpectedLeptonSF(tree->lep2_.Eta(), tree->lep2_.Pt(), tree->lid2_);
                double sfWeightEff = offlineSF1 * offlineSF2;
                double sfWeightTrig = leptonSF_->GetExpectedTriggerEfficiency(fabs(tree->lep1_.eta()),tree->lep1_.pt() ,
                        fabs(tree->lep2_.eta()), tree->lep2_.pt(),
                        TMath::Abs( tree->lid1_), TMath::Abs(tree->lid2_));
                double sfWeightPU =  nPUScaleFactor2012(fhDPU_,tree->npu_);
                */

                // apply weights except data
                if ( tree->dstype_ != SmurfTree::data ) {
                    weight = weight * sfWeightTrig_ * sfWeightEff_ * sfWeightPU_;
                }
            }
            

            // ----------------------------------------------------------------------------------------
            //      Apply data-driven background scale factors
            // ----------------------------------------------------------------------------------------
            
            // Top ------------------------------------------------------------------------------------
            if (((1ll<<sample->getDataType()) & (1ll<<TOP))) {
                weight = weight * TopScaleFactor_[njets];
            }

            // DY ------------------------------------------------------------------------------------
            // note that when make cards, yields are taken from  
            // ../../../../Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h
            // where Data-driven estimation is contained 
            // (It's not scale factors, but the final yield)
            // So, these following lines are not important => FIXME double check  
            if (((1ll<<sample->getDataType()) & (1ll<<ZLL))) {
                if ( (tree->type_ == 0 || tree->type_ == 3) && tree->dstype_ != SmurfTree::dytt ){
                    weight = weight * ZScaleFactor_[njets];
                }
            }

            // WW ------------------------------------------------------------------------------------
            if ( ((1ll<<sample->getDataType()) & ( (1ll<<QQWW) | (1ll<<GGWW)))) {
                if (option_ != WW_OPT_SMURFXSECSEL && option_ != HWW_OPT_SMURFPRESEL ) 
                    weight = weight * WWScaleFactor_[njets];
            }
           
            // W+jets --------------------------------------------------------------------------------

            // flag for mc id of leptons 
            bool isRealLepton = false;
            if((TMath::Abs(tree->lep1McId_) == 11 || TMath::Abs(tree->lep1McId_) == 13) &&
               (TMath::Abs(tree->lep2McId_) == 11 || TMath::Abs(tree->lep2McId_) == 13)) isRealLepton = true;
          
            // apply fakerates for the Wjets in data
            if ( ((1ll<<sample->getDataType()) & ( (1ll<<WJETSELEDATA) | (1ll<<WJETSMUDATA) | (1ll<<WJETSDATA) | (1ll<<WJETSMCLOOSE) ) )) {

                // for the MC events, invert the weight to be subtracted
                if ( tree->dstype_ != SmurfTree::data && (sample->getDataType() != WJETSMCLOOSE) ) {
                    weight = -1 * weight;  
                }
                    
                if ( tree->dstype_ == SmurfTree::data || tree->dstype_ == SmurfTree::wgamma || isRealLepton || (sample->getDataType() == WJETSMCLOOSE)  ) {

                    // lep1 good, lep2 fake
                    if ( tree->cuts_ & SmurfTree::Lep1FullSelection) {
                        weight *= sample->fakeRate( tree->lep2_.Pt(), tree->lep2_.Eta(), fhDFRMu_, fhDFREl_,
                                (tree->cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2,
                                (tree->cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4);
                    }

                    // lep2 good, lep1 fake
                    else if ( tree->cuts_ & SmurfTree::Lep2FullSelection) {
                        weight *= sample->fakeRate( tree->lep1_.Pt(), tree->lep1_.Eta(), fhDFRMu_, fhDFREl_,
                                (tree->cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2,
                                (tree->cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4);

                    }

                } else continue;
            } 

            // Wgamma ----------------------------------------------------------------------------------
            // This is for 2D template 
            // Yield comes from WGAMMANORM 
            // Template is re-normalizaed to WGAMMANORM later on
            if ( ( (1ll<<sample->getDataType()) & (1ll<<WGAMMA) ) ) { 
                if(!(TMath::Abs(tree->lep1McId_) == 11 || TMath::Abs(tree->lep1McId_) == 13)) 
                        weight *= ratioPhotonElectron(fhDRatioPhotonElectron_,tree->lep1_.Eta());
                if(!(TMath::Abs(tree->lep2McId_) == 11 || TMath::Abs(tree->lep2McId_) == 13)) 
                        weight *= ratioPhotonElectron(fhDRatioPhotonElectron_,tree->lep2_.Eta()); 
            }
/*
            // get the data driven estimate of the dyll by doing OF subtraction and VZ subtraction
            if ( (1ll<<sample->getDataType()) & (1ll<<ZLLDATA)) {
                if ( tree->dstype_ != SmurfTree::data) 
                    weight = -1. * weight;
                else {
                    if ( tree->type_ == SmurfTree::em || tree->type_ == SmurfTree::me ) {
                        float kee = 0.8; // electron to muon reconstruction efficiency ratio
                        weight = weight * 0.5 * (kee + 1.0/kee);
                    }
                }
            }
*/

            // Wgamma* --------------------------------------------------------------------------------
            // hard coded scale factors : k-factor = 1.5(Wmm)  and 1.0(Wee)
            // The difference is already applied in the scale1fb i.e. scale1fb(Wmm)/scale1fb(Wee) = 1.5   
            // Applying global scale factor 1.5 accounts for the difference
            if ( tree->dstype_ == SmurfTree::wgstar) weight = weight * 1.5;
            if ( tree->dstype_ == SmurfTree::wgstar && tree->type_==3 && tree->met_>130) continue; 
            // ---> remove unphysical events(bug in MC) in ee final states at high M(gamma*)
           
            // Ztt ------------------------------------------------------------------------------------
            // To Ztt embedded sample, lepton trigger and selection efficieny need to be applied.
            // At this stage, sfWeightTrig_ and sfWeightEff_ are already applied to weights.
            // ----> weight = weight * sfWeightTrig_ * sfWeightEff_ * sfWeightPU_ * scale1fb * lumi; 
            // So, only need to remove sfWeightPU_ and put 1 for sfWeightEff_ and sfWeightTrig_.
            // use lumi = 1 in ZttScaleFactor because scale1fb is already applied to weight
            if ( sample->getDataType() == ZTT) {
                weight = weight * ZttScaleFactor(2, 1, 1, 1) / sfWeightPU_;  
            }  



            // ----------------------------------------------------------------------------------------
            //      Jet binning systematics
            // ----------------------------------------------------------------------------------------

            // define two alternative jet bins for the jet energy differences
            unsigned int jetbin_up = 999;
            unsigned int jetbin_down = 999;

            // scale all jets up by 5%, considering 0->1 and 1->2 jet bin migration
            if ( njets == 0 && tree->jet1_.Pt() * 1.05 < 30 && TMath::Abs(tree->jet1_.Eta()) < 5.0 ) jetbin_up = 0;
            if ( njets == 0 && tree->jet1_.Pt() * 1.05 > 30 && TMath::Abs(tree->jet1_.Eta()) < 5.0 && tree->jet2_.Pt() * 1.05 < 30) jetbin_up = 1;
            if ( njets == 1 && tree->jet2_.Pt() * 1.05 < 30 && TMath::Abs(tree->jet2_.Eta()) < 5.0 ) jetbin_up = 1;
            if ( njets == 1 && tree->jet2_.Pt() * 1.05 > 30 && TMath::Abs(tree->jet2_.Eta()) < 5.0 ) jetbin_up = 2;
            if ( njets >= 2 ) jetbin_up = 2;

            // scale down jets by 5%,  considering only 1->0 jet migration
            if ( njets == 0 ) jetbin_down = 0;
            if ( njets == 1 && tree->jet1_.Pt() * 0.95 < 30) jetbin_down = 0;
            if ( njets == 1 && tree->jet1_.Pt() * 0.95 > 30 && TMath::Abs(tree->jet1_.Eta()) < 5.0 ) jetbin_down = 1;
            if ( njets == 2 && tree->jet2_.Pt() * 0.95 < 30) jetbin_down = 1;
            if ( njets == 2 && tree->jet2_.Pt() * 0.95 > 30 && TMath::Abs(tree->jet2_.Eta()) < 5.0 ) jetbin_down = 2;
            if ( njets >= 3 ) jetbin_down = 2;


            // ----------------------------------------------------------------------------------------
            //      Check what cuts passed
            // ----------------------------------------------------------------------------------------

            Cuts_t cuts_passed = testCuts(tree, sample->getDataType(), njets, dymva_);  


            // ----------------------------------------------------------------------------------------
            //      Record results
            // ----------------------------------------------------------------------------------------
            
            // do not care about 0 weight events 
            if(weight == 0.) continue; 

            // preselection --------------------------------------------------------------------------
            const Cuts_t hww_pass_preselection = (1ll<<HWW_PASS_PRESEL);
            if ( (cuts_passed & hww_pass_preselection) == hww_pass_preselection && 
                 (option_ == HWW_OPT_SMURFPRESEL || option_ == WW_OPT_SMURFXSECSEL || option_ == HWW_OPT_SSCTL || option_ == HWW_OPT_TOPTAG) ) {
                sample->fillResults(njets, tree->type_, weight, weight_err);
            }

            // Cut-based analysis  -------------------------------------------------------------------
            const Cuts_t hww_pass_all_cut = (1ll<<HWW_PASS_PRESEL) | (1ll<<HWW_PASS_CUTSEL);
            if ((cuts_passed & hww_pass_all_cut) == hww_pass_all_cut  && option_ == HWW_OPT_SMURFCUTSEL) {
                sample->fillResults(njets, tree->type_, weight, weight_err);
            }

            // 2D analysis ---------------------------------------------------------------------------
            if ( ( option_ == HWW_OPT_MT2DMLL || option_ == HWW_OPT_MT2DMLL_JCP || option_ == XWW_OPT_MT2DMLL_JCP || option_ == HWW_OPT_SSCTL2D ) 
                 && hww_pass_2DSelection(tree, analysis_, njets, option_) && (cuts_passed & (1ll<<HWW_PASS_PRESEL) ) == (1ll<<HWW_PASS_PRESEL) ) {
    
                if(njets<2) {    
                    sample->fillResults(njets, tree->type_, weight, weight_err);
                    sample->fill2DMVAShape(min(sample->getXMax()-0.0001, (double)tree->mt_), min(sample->getYMax()-0.0001,(double)tree->dilep_.M()), njets, tree->type_, weight); 
                

                    fillAlternateQCD(sample, tree, min(sample->getXMax()-0.0001, (double)tree->mt_), min(sample->getYMax()-0.0001,(double)tree->dilep_.M()), type, weight, njets);
                    fillAlternateFR(sample, tree, min(sample->getXMax()-0.0001, (double)tree->mt_), min(sample->getYMax()-0.0001,(double)tree->dilep_.M()), type, weight, njets);
                    fillAlternateLepEff(sample, tree, min(sample->getXMax()-0.0001, (double)tree->mt_), min(sample->getYMax()-0.0001,(double)tree->dilep_.M()), type, weight, njets);
                    fillAlternateLepRes(sample, tree, min(sample->getXMax()-0.0001, (double)mt_lepup_), min(sample->getXMax()-0.0001, (double)mt_lepdown_), 
                                        min(sample->getYMax()-0.0001,(double)mll_lepup_), min(sample->getYMax()-0.0001,(double)mll_lepdown_), type, weight, njets);
                    fillAlternateMet(sample, tree, min(sample->getXMax()-0.0001, (double)mt_metup_) , min(sample->getYMax()-0.0001,(double)mll_metup_),     type, weight, njets); 
                    if(jetbin_up<2 || (jetbin_up==2 && hww_pass_2DSelection(tree, analysis_, njets+1, option_)))  // when jetbin_up ==2, should pass VBF selection  
                            fillAlternateJES(sample, tree, min(sample->getXMax()-0.0001, (double)tree->mt_), min(sample->getYMax()-0.0001,(double)tree->dilep_.M()), type, weight, jetbin_up, jetbin_down ); 
                } else{ 
                    sample->fillResults(njets, tree->type_, weight, weight_err);
                    sample->fill2DMVAShape( min(sample->getVBFXMax()-0.0001, (double)tree->mt_), 
                                            min(sample->getVBFYMax()-0.0001, (double)tree->dilep_.M()), njets, tree->type_, weight); 

                    fillAlternateQCD(sample, tree, min(sample->getVBFXMax()-0.0001,(double)tree->mt_), 
                                        min(sample->getVBFYMax()-0.0001,(double)tree->dilep_.M()), type, weight, njets);
                    fillAlternateFR(sample, tree, min(sample->getVBFXMax()-0.0001,(double)tree->mt_), 
                                        min(sample->getVBFYMax()-0.0001,(double)tree->dilep_.M()), type, weight, njets);
                    fillAlternateLepEff(sample, tree, min(sample->getVBFXMax()-0.0001,(double)tree->mt_), 
                                        min(sample->getVBFYMax()-0.0001,(double)tree->dilep_.M()), type, weight, njets);
                    fillAlternateLepRes(sample, tree, min(sample->getVBFXMax()-0.0001,(double)mt_lepup_), 
                                        min(sample->getVBFXMax()-0.0001,(double)mt_lepdown_), 
                                        min(sample->getVBFYMax()-0.0001,(double)mll_lepup_), 
                                        min(sample->getVBFYMax()-0.0001,(double)mll_lepdown_), type, weight, njets);
                    fillAlternateMet(sample, tree, min(sample->getVBFXMax()-0.0001,(double)mt_metup_),  
                                        min(sample->getVBFYMax()-0.0001,(double)mll_metup_),    type, weight, njets); 
                    fillAlternateJES(sample, tree, min(sample->getVBFXMax()-0.0001,(double)tree->mt_), 
                                        min(sample->getVBFYMax()-0.0001,(double)tree->dilep_.M()), type, weight, jetbin_up, jetbin_down );  
            
                }
            }
           


            // ----------------------------------------------------------------------------------------
            //      Plots
            // ----------------------------------------------------------------------------------------

            bool doPlot = false;
            if ((option_ == WW_OPT_SMURFXSECSEL || option_ == HWW_OPT_SMURFPRESEL || option_ == HWW_OPT_SSCTL || option_ == HWW_OPT_TOPTAG ) 
                    && ((cuts_passed & hww_pass_preselection) == hww_pass_preselection))    doPlot = true;
            if (option_ == HWW_OPT_SMURFCUTSEL
                    && (cuts_passed & hww_pass_all_cut) == hww_pass_all_cut)                doPlot = true; 
            if ( (option_ == HWW_OPT_SSCTL2D || option_ == HWW_OPT_MT2DMLL ) 
                    && hww_pass_2DSelection(tree, analysis_, njets, option_) 
                    && (cuts_passed & (1ll<<HWW_PASS_PRESEL) ) == (1ll<<HWW_PASS_PRESEL) )  doPlot = true;

            if (doPlot) {

                // values
                float metValue = min(tree->pmet_, tree->pTrackMet_);
                if ( njets == 2 ) metValue = tree->met_;

                // fill plots 
                //if( tree->dilep_.M() > 25. && tree->dilep_.M() < 50. && tree->mt_>120 && tree->mt_<140) { 
                // --->  fill for only interesting kinematic range 

                FillHist(h1_ww_mll[njets],     type, tree->dilep_.M(),      weight);
                FillHist(h1_ww_pt1[njets],     type, tree->lep1_.Pt(),      weight);
                FillHist(h1_ww_pt2[njets],     type, tree->lep2_.Pt(),      weight);
                FillHist(h1_ww_eta1[njets],    type, TMath::Abs(tree->lep1_.Eta()),      weight);
                FillHist(h1_ww_eta2[njets],    type, TMath::Abs(tree->lep2_.Eta()),      weight);
                FillHist(h1_ww_met[njets],     type, metValue,              weight);
                FillHist(h1_ww_ptll[njets],    type, tree->dilep_.Pt(),     weight);
                float ptWWx = tree->dilep_.Px() + tree->met_*TMath::Cos(tree->metPhi_);
                float ptWWy = tree->dilep_.Py() + tree->met_*TMath::Sin(tree->metPhi_);
                FillHist(h1_ww_ptWW[njets],    type, TMath::Sqrt(ptWWx*ptWWx + ptWWy*ptWWy),     weight);
                FillHist(h1_ww_mt[njets],      type, tree->mt_,             weight);
                FillHist(h1_ww_dphi[njets],    type, tree->dPhi_,           weight);
                FillHist(h1_ww_nvtx[njets],    type, tree->nvtx_,           weight);

                //} //if( tree->dilep_.M() > 25. && tree->dilep_.M() < 50. && tree->mt_>120 && tree->mt_<140)  

            }

        } // end loop on events in file
    
        delete tree;

    } // end loop on files in chain


    gROOT->cd();

}

Cuts_t SmurfLooper::testCuts(SmurfTree *tree, DataType dataType, const unsigned int jetbin, const float& dymva)
{

    Cuts_t cuts_passed = 0;

    // 
    // SF special cuts 
    // 

    bool passDY = true;
    bool passMET = true;
    bool passMETSB = true; // met side band for determining the DY shape 

    // apply the addtional cuts to the SF events for all types
    // also apply it to the OF type of the ZLLDATA
    // apply DYMVA everywhere for the 0/1 Jets and met for the 2-jet
    if ( tree->type_ == 0 || tree->type_ == 3 ) { // || (dataType == ZLLDATA && tree->dstype_ == SmurfTree::data)) { // FIXME : erase
      if ( tree->njets_ > 1) {
        if ( ! hww_dy_selection(tree) ) passDY = false;
        if ( ! hww_sfmet_selection(tree) )   passMET = false;
        if ( hww_sfmet_selection(tree) )   passMETSB = false;
      } else {
        if ( ! (hww_sfdymva_selection(tree, dymva)) ) passMET = false;
        if ( hww_sfdymva_selection(tree, dymva) ) passMETSB = false;
        if ( dymva < -0.9 ) passMETSB = false;
      }
    }

    //
    // HWW Selections
    //
    if ( (option_ == HWW_OPT_SSCTL ||  option_ == HWW_OPT_SSCTL2D) &&
            (tree->processId_!=10010 && tree->processId_!=10001 && tree->processId_!=24 &&
             tree->processId_!=26 && tree->processId_!=121 && tree->processId_!=122) ) {
        if ( dataType == WJETSELEDATA || dataType == WJETSMUDATA || dataType == WJETSDATA) {
            if (hww_pass_wwSSPassFailSelection(tree, option_) && passDY && passMET )  cuts_passed |= (1ll<<HWW_PASS_PRESEL);
        }
        if ( dataType == WJETSMCLOOSE   && hww_pass_wwSSPassFailSelection(tree, option_)    && passDY && passMET )      cuts_passed |= (1ll<<HWW_PASS_PRESEL);
        //if ( dataType == ZLLLOOSEMET    && hww_pass_wwSSSelection(tree, option_)            && passDY && passMETSB )    cuts_passed |= (1ll<<HWW_PASS_PRESEL);
        //if ( dataType == ZLLDATA        && hww_pass_wwSSSelection(tree, option_)            && passDY && passMETSB )    cuts_passed |= (1ll<<HWW_PASS_PRESEL);
        if ( dataType == WGAMMA         && hww_pass_wgammaSSSelection(tree, option_)        && passDY && passMET )      cuts_passed |= (1ll<<HWW_PASS_PRESEL);
        if ( dataType != WJETSDATA      && dataType != WJETSELEDATA     && dataType != WJETSMUDATA  &&
                dataType != WJETSMCLOOSE   /*&& dataType != ZLLLOOSEMET      && dataType != ZLLDATA*/      && dataType != WGAMMA  &&
                hww_pass_wwSSSelection(tree, option_) && passDY && passMET )  cuts_passed |= (1ll<<HWW_PASS_PRESEL);
    } else {
        if ( dataType == WJETSELEDATA || dataType == WJETSMUDATA || dataType == WJETSDATA) {
            if (hww_pass_wwPassFailSelection(tree, option_) && passDY && passMET )  cuts_passed |= (1ll<<HWW_PASS_PRESEL);
        }
        if ( dataType == WJETSMCLOOSE && hww_pass_wwPassFailSelection(tree, option_) && passDY && passMET ) cuts_passed |= (1ll<<HWW_PASS_PRESEL); 
        //if ( dataType == ZLLLOOSEMET && hww_pass_wwSelection(tree, option_) && passDY && passMETSB )  cuts_passed |= (1ll<<HWW_PASS_PRESEL);
        //if ( dataType == ZLLDATA && hww_pass_wwSelection(tree, option_) && passDY && passMETSB )  cuts_passed |= (1ll<<HWW_PASS_PRESEL);
        if ( dataType == WGAMMA && hww_pass_wgammaSelection(tree, option_) && passDY && passMET )  cuts_passed |= (1ll<<HWW_PASS_PRESEL);
        if ( dataType != WJETSDATA && dataType != WJETSELEDATA && dataType != WJETSMUDATA &&dataType != WJETSMCLOOSE /*&& dataType != ZLLLOOSEMET 
                && dataType != ZLLDATA*/ && dataType != WGAMMA  && hww_pass_wwSelection(tree, option_) && passDY && passMET )  cuts_passed |= (1ll<<HWW_PASS_PRESEL);
    }

    if (hww_pass_cutSelection(tree, analysis_, jetbin))             cuts_passed |= (1ll<<HWW_PASS_CUTSEL);

    return cuts_passed;
}

void SmurfLooper::loadWeightHistograms() 
{

    //
    // higgs pt dependent re-weighting
    //  

    float ana = analysis_;  
    // use the same file for mH(110) with the mH(115)
    if ( analysis_ <= 115) ana = 115;
    if ( analysis_ > 115 && analysis_ < 125) ana = 120;
    if ( analysis_ >= 125 && analysis_ < 140) ana = 130;


    if (analysis_ == 0 || analysis_ == 1 || analysis_ == 2) ana = 250.0;
    // for the intermediate mass points
    if ( analysis_ <= 300 && analysis_ > 250) ana = 300;
    if ( analysis_ <= 350 && analysis_ > 300) ana = 350;
    if ( analysis_ <= 400 && analysis_ > 350) ana = 400;
    if ( analysis_ <= 500 && analysis_ > 400) ana = 500;
    if ( analysis_ <= 600 && analysis_ > 500) ana = 600;

    const char *HiggsPtKFactorFileName = "/smurf/data/Winter11_4700ipb/auxiliar/ggHWW_KFactors_PowhegToHQT_WithAdditionalMassPoints.root"; // TAS
    //const char *HiggsPtKFactorFileName = "/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/aux/ggHWW_KFactors_PowhegToHQT_WithAdditionalMassPoints.root"; // UAF

    TFile *fHiggsPtKFactorFile = TFile::Open(HiggsPtKFactorFileName, "READ");
    std::string kfactorHistName;
    kfactorHistName = Form("KFactor_PowhegToHQT_mH%i", int(ana));
    HiggsPtKFactor_ = (TH1D*)(fHiggsPtKFactorFile->Get(kfactorHistName.c_str()));
    HiggsPtKFactor_QCDscaleSys1_ = (TH1D*)(fHiggsPtKFactorFile->Get(Form("KFactor_PowhegToHQT_mH%i_QCDscaleSys1", int(ana))));
    HiggsPtKFactor_QCDscaleSys2_ = (TH1D*)(fHiggsPtKFactorFile->Get(Form("KFactor_PowhegToHQT_mH%i_QCDscaleSys2", int(ana))));
    HiggsPtKFactor_QCDscaleSys3_ = (TH1D*)(fHiggsPtKFactorFile->Get(Form("KFactor_PowhegToHQT_mH%i_QCDscaleSys3", int(ana))));
    HiggsPtKFactor_QCDscaleSys4_ = (TH1D*)(fHiggsPtKFactorFile->Get(Form("KFactor_PowhegToHQT_mH%i_QCDscaleSys4", int(ana))));
    HiggsPtKFactor_QCDscaleSys5_ = (TH1D*)(fHiggsPtKFactorFile->Get(Form("KFactor_PowhegToHQT_mH%i_QCDscaleSys5", int(ana))));
    HiggsPtKFactor_QCDscaleSys6_ = (TH1D*)(fHiggsPtKFactorFile->Get(Form("KFactor_PowhegToHQT_mH%i_QCDscaleSys6", int(ana))));
    HiggsPtKFactor_QCDscaleSys7_ = (TH1D*)(fHiggsPtKFactorFile->Get(Form("KFactor_PowhegToHQT_mH%i_QCDscaleSys7", int(ana))));
    HiggsPtKFactor_QCDscaleSys8_ = (TH1D*)(fHiggsPtKFactorFile->Get(Form("KFactor_PowhegToHQT_mH%i_QCDscaleSys8", int(ana))));
    if (HiggsPtKFactor_) {
        HiggsPtKFactor_->SetDirectory(0);
        HiggsPtKFactor_QCDscaleSys1_->SetDirectory(0);
        HiggsPtKFactor_QCDscaleSys2_->SetDirectory(0);
        HiggsPtKFactor_QCDscaleSys3_->SetDirectory(0);
        HiggsPtKFactor_QCDscaleSys4_->SetDirectory(0);
        HiggsPtKFactor_QCDscaleSys5_->SetDirectory(0);
        HiggsPtKFactor_QCDscaleSys6_->SetDirectory(0);
        HiggsPtKFactor_QCDscaleSys7_->SetDirectory(0);
        HiggsPtKFactor_QCDscaleSys8_->SetDirectory(0);
    }

    // assert(HiggsPtKFactor_);
    if ( HiggsPtKFactor_ == 0x0) 
        std::cout << "Warning: HiggsPtKFactor files for mH = " << analysis_ << " is not found!\n";
    fHiggsPtKFactorFile->Close();
    delete fHiggsPtKFactorFile; 


    // 
    // Lepton efficiencies
    //
    leptonSF_ = new LeptonScaleLookup("/smurf/dlevans/Efficiencies/V00-02-09/summary_Moriond_V1.root"); // TAS
    //leptonSF_ = new LeptonScaleLookup("/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/aux/summary_Moriond_V1.root"); // UAF
    
    // lepton efficiency uncertainties
    TFile *fLeptonEffError = 0;
    fLeptonEffError = TFile::Open("/smurf/dlevans/Efficiencies/V00-02-09/systematics_Moriond_V1.root"); // TAS
    //fLeptonEffError = TFile::Open("/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/aux/systematics_Moriond_V1.root"); // UAF
    
    fhDMuonEffError_ = (TH2D*)(fLeptonEffError->Get("h2_nm1_syst_muon_selection")); 
    assert(fhDMuonEffError_);
    fhDMuonEffError_->SetDirectory(0);
    
    fhDElectronEffError_ = (TH2D*)(fLeptonEffError->Get("h2_nm1_syst_electron_selection")); 
    assert(fhDElectronEffError_);
    fhDElectronEffError_->SetDirectory(0);
    
    fLeptonEffError->Close();
    delete fLeptonEffError;
    
    //
    // Fakerate for the hww analysis
    //

    TFile *fLeptonFRFileM = 0;
    TFile *fLeptonFRFileE = 0;

    // load electron fake rate histograms
    fLeptonFRFileE = TFile::Open("/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root"); // TAS
    //fLeptonFRFileE = TFile::Open("/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/aux/summary_fakes_Moriond2012.root"); // UAF
    fhDFREl_ = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta")); 
    assert(fhDFREl_);
    fhDFREl_->SetDirectory(0);

    fhDFREl_systvar_ = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold20_PtEta"));
    assert(fhDFREl_systvar_);
    fhDFREl_systvar_->SetDirectory(0);

    fLeptonFRFileE->Close();
    delete fLeptonFRFileE;

    // load muon fake rate histograms
    fLeptonFRFileM = TFile::Open("/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root"); // TAS
    //fLeptonFRFileM = TFile::Open("/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/aux/summary_fakes_Moriond2012.root"); // UAF
    fhDFRMu_ = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold30_PtEta"));  
    assert(fhDFRMu_);
    fhDFRMu_->SetDirectory(0);

    fhDFRMu_systvar_ = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
    assert(fhDFRMu_systvar_);
    fhDFRMu_systvar_->SetDirectory(0);

    fLeptonFRFileM->Close();
    delete fLeptonFRFileM; 

    // photon -> electron conversion ratio 
    TFile *fRatioPhotonElectron = TFile::Open("/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/ratio_photon_electron.root"); // TAS
    //TFile *fRatioPhotonElectron = TFile::Open("/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/aux/ratio_photon_electron.root"); // UAF
    fhDRatioPhotonElectron_ = (TH1D*)(fRatioPhotonElectron->Get("hDRatioPhotonElectron"));
    assert(fhDRatioPhotonElectron_);
    fhDRatioPhotonElectron_->SetDirectory(0);
    fRatioPhotonElectron->Close();
    delete fRatioPhotonElectron;

    // pileup 
    //TFile *fPUFile = TFile::Open("/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_12p1ifb.root"); // HCP
    //TFile *fPUFile = TFile::Open("/smurf/jaehyeok/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets_19p5ifb/auxillar/puWeights_Summer12_53x_True_postHCP.root"); // postHCP 
    //fhDPU_ = (TH1D*)(fPUFile->Get("puWeights"));
    //assert(fhDPU_);              
    //fhDPU_->SetDirectory(0);    
    //delete fPUFile;

    // 
    // scale factors for the hww analysis
    // 

    getTopScaleFactor(TopScaleFactor_, TopScaleFactorError_, option_, analysis_);
    getWWScaleFactor(WWScaleFactor_, WWScaleFactorError_, option_, analysis_);
    getZScaleFactor(ZScaleFactor_, ZScaleFactorError_, option_, analysis_, "sf");

}

void SmurfLooper::setGoodRunList(const char *runlist)
{
    set_goodrun_file(runlist);
    runlistIsSet_ = true;
}

void SmurfLooper::setLumiScale(float eeLumi, float mmLumi)
{
    std::cout << "[SmurfLooper::setLumiScale] re-setting lumi to " << eeLumi << ", " << mmLumi << std::endl;
    eeLumi_ = eeLumi;
    mmLumi_ = mmLumi;
}

float SmurfLooper::getScaleQCD(float mr, float mf, SmurfTree *tree) {

    int bin = HiggsPtKFactor_QCDscaleSys1_->GetXaxis()->FindFixBin(tree->higgsPt_);
    if (bin > HiggsPtKFactor_QCDscaleSys1_->GetNbinsX()) bin = HiggsPtKFactor_QCDscaleSys1_->GetNbinsX();

    if (mr == 0.5 && mf == 0.5) {
        return HiggsPtKFactor_QCDscaleSys1_->GetBinContent(bin);
    } else if (mr == 0.5 && mf == 1.0) {
        return HiggsPtKFactor_QCDscaleSys2_->GetBinContent(bin);
    } else if (mr == 1.0 && mf == 0.5) {
        return HiggsPtKFactor_QCDscaleSys3_->GetBinContent(bin);
    } else if (mr == 1.0 && mf == 2.0) {
        return HiggsPtKFactor_QCDscaleSys4_->GetBinContent(bin);
    } else if (mr == 2.0 && mf == 1.0) {
        return HiggsPtKFactor_QCDscaleSys5_->GetBinContent(bin);
    } else if (mr == 2.0 && mf == 2.0) {
        return HiggsPtKFactor_QCDscaleSys6_->GetBinContent(bin);
    } else if (mr == 0.5 && mf == 2.0) {
        return HiggsPtKFactor_QCDscaleSys7_->GetBinContent(bin);
    } else if (mr == 2.0 && mf == 0.5) {
        return HiggsPtKFactor_QCDscaleSys8_->GetBinContent(bin);
    } else {
        std::cout << "scaleQCD: Invalid variation" << std::endl;
        return 0.0;
    }

}

void SmurfLooper::fillAlternateQCD(SmurfSample *sample, SmurfTree *tree, 
        const float &varx, const float &vary, unsigned int type, const float &weight, unsigned int jetbin) {

    //
    // find out what extra shapes 
    // this sample has associated with it
    //

    ShapeVar_t availableShapeSystematics = sample->getAvailableShapeSystematicsMask();

    // QCD Scale
    if (availableShapeSystematics & (1ll<<QCDSCALEVAR)) {

        // calculate new weight
        int bin = HiggsPtKFactor_->GetXaxis()->FindFixBin(tree->higgsPt_);
        if (bin > HiggsPtKFactor_->GetNbinsX()) bin = HiggsPtKFactor_->GetNbinsX();
        float HiggsPtKFactor = HiggsPtKFactor_->GetBinContent(bin);
        float newWeightUp   = (weight / HiggsPtKFactor) * getScaleQCD(0.5, 2.0, tree);
        float newWeightDown = (weight / HiggsPtKFactor) * getScaleQCD(2.0, 0.5, tree); 

        // fill the alternate shapes
        if ((1ll<<option_) & HWW_MT2DMLL ) {
            sample->fillShapeVariation2D(QCDSCALEVAR, true, varx, vary, jetbin, type, newWeightUp);
            sample->fillShapeVariation2D(QCDSCALEVAR, false, varx, vary, jetbin, type, newWeightDown);
        } else {
            //sample->fillShapeVariation1D(QCDSCALEVAR, true, varx, jetbin, type, newWeightUp);
            //sample->fillShapeVariation1D(QCDSCALEVAR, false, varx, jetbin, type, newWeightDown);     
        }

    }
}

void SmurfLooper::fillAlternateMet(SmurfSample *sample, SmurfTree *tree,   
        const float &varx, const float &vary, unsigned int type, const float &weight, unsigned int jetbin) {

    unsigned int availableShapeSystematics = sample->getAvailableShapeSystematicsMask();

    // fill the alternate shapes
    if (availableShapeSystematics & (1ll<<METVAR)) {
        // note - this is ultimately going to be a mirror of the 
        // central and the up histograms, so the down histogram is reduntant 
        // and not filled.
       if ((1ll<<option_)  & HWW_MT2DMLL ) {
            sample->fillShapeVariation2D(METVAR, true, varx, vary, jetbin, type, weight);
        } else {
            //sample->fillShapeVariation1D(METVAR, true, varx, jetbin, type, weight);
        }
    }
}

void SmurfLooper::fillAlternateLepRes(SmurfSample *sample, SmurfTree *tree,
        const float &var_upx, const float &var_downx, const float &var_upy, const float &var_downy, 
        unsigned int type, const float &weight, unsigned int jetbin) 
{

    unsigned int availableShapeSystematics = sample->getAvailableShapeSystematicsMask();

    // fill the alternate shapes
    if (availableShapeSystematics & (1ll<<LEPRESVAR)) {
        if ( (1ll<<option_)  & HWW_MT2DMLL ) {
            sample->fillShapeVariation2D(LEPRESVAR, true, var_upx, var_upy, jetbin, type, weight);
            sample->fillShapeVariation2D(LEPRESVAR, false, var_downx, var_downy, jetbin, type, weight);
        } else {
            //sample->fillShapeVariation1D(LEPRESVAR, true, var_upx, jetbin, type, weight);
            //sample->fillShapeVariation1D(LEPRESVAR, false, var_downx, jetbin, type, weight);
        }
    }

}

void SmurfLooper::fillAlternateLepEff(SmurfSample *sample, SmurfTree *tree,   
        const float &varx, const float &vary, unsigned int type, const float &weight, unsigned int jetbin)
{
    //
    // find out what extra shapes 
    // this sample has associated with it
    //

    ShapeVar_t availableShapeSystematics = sample->getAvailableShapeSystematicsMask();

    // lepton efficiency variation
    if (availableShapeSystematics & (1ll<<LEPEFFVAR)) {

        // get existing lepton scale
        double offlineSF1 = leptonSF_->GetExpectedLeptonSF(tree->lep1_.Eta(), tree->lep1_.Pt(), tree->lid1_);
        double offlineSF2 = leptonSF_->GetExpectedLeptonSF(tree->lep2_.Eta(), tree->lep2_.Pt(), tree->lid2_);

        double stat1 = leptonSF_->GetExpectedLeptonSFErr(tree->lep1_.Eta(), tree->lep1_.Pt(), tree->lid1_);
        double stat2 = leptonSF_->GetExpectedLeptonSFErr(tree->lep2_.Eta(), tree->lep2_.Pt(), tree->lid2_);
        double syst1 = TMath::Abs(tree->lid1_)==11 ? 0.02 : 0.015;
        double syst2 = TMath::Abs(tree->lid2_)==11 ? 0.02 : 0.015; 
        double syst1_ext = sample->LepEffError( tree->lep1_.Pt(), tree->lep1_.Eta(), fhDMuonEffError_, fhDElectronEffError_,  
                                                TMath::Abs(tree->lid1_)==13, TMath::Abs(tree->lid1_)==11 ); // additional syst from Eff measurement
        double syst2_ext = sample->LepEffError( tree->lep2_.Pt(), tree->lep2_.Eta(), fhDMuonEffError_, fhDElectronEffError_,  
                                                TMath::Abs(tree->lid2_)==13, TMath::Abs(tree->lid2_)==11 ); // additional syst from Eff measurement

        //cout << "pt1: " << tree->lep1_.Pt() << " eta1 : " << tree->lep1_.Eta() << " ::: " << syst1_ext << endl;
        //cout << "pt2: " << tree->lep2_.Pt() << " eta2 : " << tree->lep2_.Eta() << " ::: " << syst2_ext << endl;

        syst1 = TMath::Sqrt(syst1*syst1 + syst1_ext*syst1_ext);
        syst2 = TMath::Sqrt(syst2*syst2 + syst2_ext*syst2_ext);

        // get the lepton scale uncertainties
        //double offlineSF1Err = leptonSF_->GetExpectedLeptonSFErr(tree->lep1_.Eta(), tree->lep1_.Pt(), tree->lid1_)+syst1;
        //double offlineSF2Err = leptonSF_->GetExpectedLeptonSFErr(tree->lep2_.Eta(), tree->lep2_.Pt(), tree->lid2_)+syst2;
        double offlineSF1Err = TMath::Sqrt(stat1*stat1+syst1*syst1);
        double offlineSF2Err = TMath::Sqrt(stat2*stat2+syst2*syst2);

        // calculate new weight
        double newWeightUp (0.0);
        double newWeightDown (0.0);
        if ( offlineSF1 * offlineSF2 > 0) {
            newWeightUp = (weight / (offlineSF1 * offlineSF2)) * (offlineSF1+offlineSF1Err) * (offlineSF2+offlineSF2Err);
            newWeightDown = (weight / (offlineSF1 * offlineSF2)) * (offlineSF1-offlineSF1Err) * (offlineSF2-offlineSF2Err);
        }

        // fill the alternate shapes
        if ((1ll<<option_) & HWW_MT2DMLL) {
            sample->fillShapeVariation2D(LEPEFFVAR, true, varx, vary, jetbin, type, newWeightUp);
            sample->fillShapeVariation2D(LEPEFFVAR, false, varx, vary, jetbin, type, newWeightDown);
        } else {
            //sample->fillShapeVariation1D(LEPEFFVAR, true, varx, jetbin, type, newWeightUp);
            //sample->fillShapeVariation1D(LEPEFFVAR, false, varx, jetbin, type, newWeightDown);
        }

    }

}

void SmurfLooper::fillAlternateFR(SmurfSample *sample, SmurfTree *tree, 
        const float &varx, const float &vary, unsigned int type, const float &weight, unsigned int jetbin) 
{

    //
    // find out what extra shapes 
    // this sample has associated with it
    //

    ShapeVar_t availableShapeSystematics = sample->getAvailableShapeSystematicsMask();
    if ( (availableShapeSystematics & (1ll<<WJETSELESHAPEVAR)) || (availableShapeSystematics & (1ll<<WJETSMUSHAPEVAR))  ) {

        float central_weight = 1.0;
        float up_weight = 1.0;

        if ( tree->dstype_ == SmurfTree::data || tree->dstype_ == SmurfTree::wgamma || (TMath::Abs(tree->lep1McId_*tree->lep2McId_)>0)) {

            // lep1 good, lep2 fake
            if ( tree->cuts_ & SmurfTree::Lep1FullSelection) {
                central_weight = sample->fakeRate( tree->lep2_.Pt(), tree->lep2_.Eta(), fhDFRMu_, fhDFREl_,
                        (tree->cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2,
                        (tree->cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4);
                up_weight = sample->fakeRate( tree->lep2_.Pt(), tree->lep2_.Eta(), fhDFRMu_systvar_, fhDFREl_systvar_,
                        (tree->cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2,
                        (tree->cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4);
            }

            // lep2 good, lep1 fake
            else if ( tree->cuts_ & SmurfTree::Lep2FullSelection) {
                central_weight = sample->fakeRate( tree->lep1_.Pt(), tree->lep1_.Eta(), fhDFRMu_, fhDFREl_,
                        (tree->cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2,
                        (tree->cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4);
                up_weight = sample->fakeRate( tree->lep1_.Pt(), tree->lep1_.Eta(), fhDFRMu_systvar_, fhDFREl_systvar_,
                        (tree->cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2,
                        (tree->cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4);

            }
        }

        // fill the alternate shapes
        if ((1ll<<option_) & HWW_MT2DMLL) {
            sample->fillShapeVariation2D(WJETSELESHAPEVAR, true,  varx, vary, jetbin, type, weight*up_weight/central_weight);
            sample->fillShapeVariation2D(WJETSELESHAPEVAR, false, varx, vary, jetbin, type, weight*up_weight/central_weight);
            sample->fillShapeVariation2D(WJETSMUSHAPEVAR, true,  varx, vary, jetbin, type, weight*up_weight/central_weight);
            sample->fillShapeVariation2D(WJETSMUSHAPEVAR, false, varx, vary, jetbin, type, weight*up_weight/central_weight);
        } else {
            //sample->fillShapeVariation1D(WJETSELESHAPEVAR, true,  varx, jetbin, type, weight*up_weight/central_weight);
            //sample->fillShapeVariation1D(WJETSELESHAPEVAR, false, varx, jetbin, type, weight*up_weight/central_weight);
            //sample->fillShapeVariation1D(WJETSMUSHAPEVAR, true,  varx, jetbin, type, weight*up_weight/central_weight);
            //sample->fillShapeVariation1D(WJETSMUSHAPEVAR, false, varx, jetbin, type, weight*up_weight/central_weight);
        } 

    }
}

void SmurfLooper::fillAlternateJES(SmurfSample *sample, SmurfTree *tree,   
        const float &varx, const float &vary, unsigned int type, const float &weight, unsigned int jetbin_up, unsigned int jetbin_down) 
{

    unsigned int availableShapeSystematics = sample->getAvailableShapeSystematicsMask();

    // fill the alternate shapes
    if (availableShapeSystematics & (1ll<<JETRESVAR)) {
        if ((1ll<<option_) & HWW_MT2DMLL) {
            if ( jetbin_up == 0 || jetbin_up == 1 || jetbin_up == 2)      sample->fillShapeVariation2D(JETRESVAR, true, varx, vary, jetbin_up, type, weight);
            if ( jetbin_down == 0 || jetbin_down == 1 || jetbin_down == 2)  sample->fillShapeVariation2D(JETRESVAR, false, varx, vary, jetbin_down, type, weight);
        } else {
            //if ( jetbin_up == 0 || jetbin_up == 1)      sample->fillShapeVariation1D(JETRESVAR, true, varx, jetbin_up, type, weight);
            //if ( jetbin_down == 0 || jetbin_down == 1)  sample->fillShapeVariation1D(JETRESVAR, false, varx, jetbin_down, type, weight);
        }
    }
}

void SmurfLooper::fillAlternateJESVBF(SmurfSample *sample, SmurfTree *tree,
        const float &var_upx, const float &var_downx, const float &var_upy, const float &var_downy,
        unsigned int type, const float &weight, unsigned int jetbin) 
{

    unsigned int availableShapeSystematics = sample->getAvailableShapeSystematicsMask();

    // fill the alternate shapes
    if (availableShapeSystematics & (1ll<<JETRESVAR)) {
        if ((1ll<<option_) & HWW_MT2DMLL) {
            sample->fillShapeVariation2D(JETRESVAR, true, var_upx, var_upy, jetbin, type, weight);
            sample->fillShapeVariation2D(JETRESVAR, false, var_downx, var_downy, jetbin, type, weight);
        } else {
            //sample->fillShapeVariation1D(JETRESVAR, true, var_upx, jetbin, type, weight);
            //sample->fillShapeVariation1D(JETRESVAR, false, var_downx, jetbin, type, weight);
        }
    }

}

