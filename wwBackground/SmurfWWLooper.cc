#include "SmurfWWLooper.h"

#include "../../../../../Smurf/Core/SmurfTree.h"
#include "../../../NtupleMacros/Tools/goodrun.h"
#include "../core/Selections.h"
#include "../core/SmurfPlotUtilities.h"
#include "../core/SmurfSample.h"
#include "../core/Enums.h"
#include "../SmurfScaleFactors.h"
#include "../../../../../Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"

#include "TChainElement.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom.h"

#include <cmath>
#include <iostream>
#include <cassert>

SmurfWWLooper::SmurfWWLooper(float analysis, Option option, RunEra runEra) 
{

    std::cout << std::endl;
    std::cout << "[SmurfWWLooper::SmurfWWLooper] Doing mass point " << analysis << std::endl;
    std::cout << "[SmurfWWLooper::SmurfWWLooper] Doing analysis option " << option << std::endl;

    analysis_ = analysis;
    runEra_ = runEra;
    option_ = option;
    runlistIsSet_ = false;
    eeLumi_ = 4630;
    mmLumi_ = 4630;;
    std::cout << "[SmurfWWLooper::SmurfWWLooper] default lumi is " << eeLumi_ << ", " << mmLumi_ << std::endl;    
    loadWeightHistograms();
}

void SmurfWWLooper::loop(SmurfSample *sample)
{
    TObjArray *listOfFiles = sample->getChain()->GetListOfFiles();
    TIter fileIter(listOfFiles);
    if (listOfFiles->GetEntries() == 0) {
        std::cout << "[SmurfWWLooper::loop] " << sample->getName() << " is not defined" << std::endl;
        return;
    }
    else
        std::cout << "[SmurfWWLooper::loop] " << sample->getName() << std::endl;

    //
    // set up histograms
    //

    int nbins(20);
    float xmin = 10.;
    float xmax = 200.;

    TH1F *h1_mll_sig[kJetBins][kLeptonTypes];
    TH1F *h1_mll_bkg[kJetBins][kLeptonTypes];

    for (unsigned int j = 0; j < kJetBins; ++j) {
        FormatHist(h1_mll_sig[j], sample, Form("mll_sig_%s", jetbin_names[j]), "mll in Signal Region", nbins, xmin, xmax);
        FormatHist(h1_mll_bkg[j], sample, Form("mll_bkg_%s", jetbin_names[j]), "mll in Control Region", nbins, xmin, xmax);
    }

    gROOT->cd();

    //
    // file loop
    //

    unsigned int nEventsChain =sample->getChain()->GetEntries();
    unsigned int nEventsTotal = 0;
    int i_permille_old = 0;

    unsigned int nf = 0;
    while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
        SmurfTree *tree = new SmurfTree();
        tree->LoadTree(currentFile->GetTitle());
        tree->InitTree(0);

        // extra variables for the various correction factors
        float sfWeightPU_ = 1.0;
        float sfWeightTrig_ = 1.0;
        float sfWeightEff_ = 1.0;
        float dymva_ = 1.0;
        if (tree->tree_->GetBranchStatus("sfWeightPU"))
            tree->tree_->SetBranchAddress("sfWeightPU", &sfWeightPU_);
        if (tree->tree_->GetBranchStatus("sfWeightTrig"))
            tree->tree_->SetBranchAddress("sfWeightTrig", &sfWeightTrig_);
        if (tree->tree_->GetBranchStatus("sfWeightEff"))
            tree->tree_->SetBranchAddress("sfWeightEff", &sfWeightEff_);
        if (tree->tree_->GetBranchStatus("dymva"))
            tree->tree_->SetBranchAddress("dymva", &dymva_);

        //
        //
        // event loop
        //

        ULong64_t nEvents = tree->tree_->GetEntries();
        for(ULong64_t event = 0; event < nEvents; ++event) {
            tree->tree_->GetEntry(event);

            //
            // incrimenet counters
            //

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

            //
            // set up weights and binning
            //

            double weight = 1.0;
            double weight_err = 0.0;

            unsigned int type = tree->type_;
            unsigned int njets = tree->njets_;
            if ( tree->njets_ > 1 ) continue;

            //
            // apply good run list and MC weights
            //

            if ( tree->dstype_ == SmurfTree::data)  {
                if (runlistIsSet_) {
                    if (!goodrun_json(tree->run_, tree->lumi_)) continue;
                }
            } else {
                if (tree->type_ == 0) weight = tree->scale1fb_*(mmLumi_/1000.0);
                if (tree->type_ == 3) weight = tree->scale1fb_*(eeLumi_/1000.0);
                if (tree->type_ == 1 || tree->type_ == 2) weight = tree->scale1fb_*((eeLumi_+mmLumi_)/(2*1000.0));
                weight = weight * sfWeightTrig_ * sfWeightEff_ * sfWeightPU_;
            }

            // ww-level selections
            if ( sample->getDataType() == WJETSDATA ) {
                if (! hww_pass_wwPassFailSelection(tree, option_)) continue;     
                // for the MC events, invert the weight to be subtracted
                if ( tree->dstype_ != SmurfTree::data ) {
                    weight = -1 * weight;
                }
                if ( tree->dstype_ == SmurfTree::data || tree->dstype_ == SmurfTree::wgamma || (TMath::Abs(tree->lep1McId_*tree->lep2McId_)>0)) {
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
                }
            } else {
                if (! hww_pass_wwSelection(tree, option_) ) continue;
            }

            //
            // SF selections
            // 
            if ( tree->type_ == 0 || tree->type_ == 3 ) {
                if ( ! (hww_sfdymva_selection(tree, dymva_)) ) continue;
            }
            //
            // Apply data/MC scale factors
            //
            if (((1ll<<sample->getDataType()) & (1ll<<TOP))) {
                /*
                // use scale factor by hand
                if ( njets == 0 ) 
                weight = weight * 1.03;
                if ( njets == 1 ) 
                weight = weight * 1.09;
                 */
                weight = weight * TopScaleFactor_[njets];  
                // std::cout << "njets = " << njets << ": TopScaleFactor " << TopScaleFactor_[njets] << "\n";
            }       

            if (((1ll<<sample->getDataType()) & (1ll<<ZLL))) {
                if ( tree->type_ == 0 || tree->type_ == 3 ) {
                    /*
                    // use scale factor by hand
                    if ( njets == 0 ) 
                    weight = weight * 4.26;
                    if ( njets == 1 ) 
                    weight = weight * 3.81;
                     */
                    weight = weight * ZScaleFactor_[njets]; 
                    // std::cout << "njets = " << njets << ": ZScaleFactor " << ZScaleFactor_[njets] << "\n";
                }
            }       

            if ( sample->getDataType() == ZTT) {                                                                                      
                //weight = ZttScaleFactor(2, 1) * (eeLumi_ + mmLumi_) /2/1000.;                      
                weight = weight * ZttScaleFactor(2, 1, 1, 1) / sfWeightPU_ * (eeLumi_ + mmLumi_) /2/1000.;  
            } 

            if ( tree->dstype_ == SmurfTree::wgstar) weight = weight * 1.5;                                                           

            // apply all higgs mass dependent cuts EXCEPT mll/mT/dphi
            float lep1ptCut(20.), lep2ptCut(10.), mllCut(9999.), dPhiCut(9999.), mtLowCut(0.), mtHighCut(9999.),  mllLooseCut(9999.);
            evaluate_hww_cuts(analysis_, lep1ptCut, lep2ptCut, mllCut, dPhiCut, mtLowCut, mtHighCut, mllLooseCut);
            if ( option_ == HWW_OPT_SMURFCUTSEL ) {
                if ( tree->lep1_.Pt() < lep1ptCut) continue;
                if ( tree->lep2_.Pt() < lep2ptCut) continue;
            }

            //
            // fill the control region histograms
            // 

            double mllcrcut = 100.;
            if ( option_ == HWW_OPT_SMURFMVASEL ) 
                //mllcrcut = max(100., double(mllLooseCut));
                mllcrcut = -999.;

            if ( tree->dilep_.mass() > mllcrcut ) {
                FillHist(h1_mll_bkg[njets], type, tree->dilep_.mass(), weight);
            }

            if ( option_ == HWW_OPT_SMURFCUTSEL && hww_pass_cutSelection(tree, analysis_, njets) ) {
                FillHist(h1_mll_sig[njets], type, tree->dilep_.mass(), weight);
                sample->fillResults(njets, tree->type_, weight, weight_err);
            }

            if ( option_ == HWW_OPT_SMURFMVASEL && hww_pass_mvaSelection(tree, analysis_, njets) ) {
                FillHist(h1_mll_sig[njets], type, tree->dilep_.mass(), weight);
                sample->fillResults(njets, tree->type_, weight, weight_err);
            }

        } // end loop on events in file
        nf++;
        delete tree;
    } // end loop on files in chain
    gROOT->cd();
}



void SmurfWWLooper::setGoodRunList(const char *runlist)
{
    set_goodrun_file(runlist);
    runlistIsSet_ = true;
}

void SmurfWWLooper::setLumiScale(float eeLumi, float mmLumi)
{
    std::cout << "[SmurfWWLooper::setLumiScale] re-setting lumi to " << eeLumi << ", " << mmLumi << std::endl;
    eeLumi_ = eeLumi;
    mmLumi_ = mmLumi;

}

void SmurfWWLooper::loadWeightHistograms() 
{

    //
    // Fakerate
    //

    TFile *fLeptonFRFileM = 0;
    TFile *fLeptonFRFileE = 0;

    // load electron fake rate histograms
    fLeptonFRFileE = TFile::Open("/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root");
    fhDFREl_ = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
    assert(fhDFREl_);
    fhDFREl_->SetDirectory(0);
    fLeptonFRFileE->Close();
    delete fLeptonFRFileE;

    // load muon fake rate histograms
    fLeptonFRFileM = TFile::Open("/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root");
    fhDFRMu_ = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold30_PtEta"));
    assert(fhDFRMu_);
    fhDFRMu_->SetDirectory(0);
    fLeptonFRFileM->Close();
    delete fLeptonFRFileM; 


    // 
    // Top scale factors
    // 
    getTopScaleFactor(TopScaleFactor_, TopScaleFactorError_, option_, 0);
    getZScaleFactor(ZScaleFactor_, ZScaleFactorError_, option_, 0, "sf");
}
