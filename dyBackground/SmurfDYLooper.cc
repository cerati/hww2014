#include "SmurfDYLooper.h"
#include "../../Smurf/Core/SmurfTree.h"
#include "../../Tools/goodrun.h"
#include "../core/Selections.h"
#include "../core/SmurfPlotUtilities.h"
#include "../core/SmurfSample.h"

#include "TChainElement.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom.h"

#include <cmath>
#include <iostream>
#include <cassert>

SmurfDYLooper::SmurfDYLooper(float analysis, Option option, RunEra runEra) 
{
 
    std::cout << std::endl;
    std::cout << "[SmurfDYLooper::SmurfDYLooper] Doing mass point " << analysis << std::endl;
    std::cout << "[SmurfDYLooper::SmurfDYLooper] Doing analysis option " << option << std::endl;
    
    analysis_ = analysis;
    runEra_ = runEra;
    option_ = option;
    runlistIsSet_ = false;
    eeLumi_ = 4630;
    mmLumi_ = 4630;;
    std::cout << "[SmurfDYLooper::SmurfDYLooper] default lumi is " << eeLumi_ << ", " << mmLumi_ << std::endl;
}

void SmurfDYLooper::loop(SmurfSample *sample)
{
    TObjArray *listOfFiles = sample->getChain()->GetListOfFiles();
    TIter fileIter(listOfFiles);
    if (listOfFiles->GetEntries() == 0) {
        std::cout << "[SmurfDYLooper::loop] " << sample->getName() << " is not defined" << std::endl;
        return;
    }
    else
        std::cout << "[SmurfDYLooper::loop] " << sample->getName() << std::endl;

    //
    // set up histograms
    //

    // for the DY estimations

    TH1F *h1_met_in_ww[kJetBins][kLeptonTypes];;
    TH1F *h1_met_in[kJetBins][kLeptonTypes];;
    TH1F *h1_met_out[kJetBins][kLeptonTypes];;
    TH1F *h1_met_in_sig[kJetBins][kLeptonTypes];;
    TH1F *h1_met_out_sig[kJetBins][kLeptonTypes];;
    

    for (unsigned int j = 0; j < kJetBins; ++j) {
      
      int nbins = 4;  
      float bins[5] = {20, 25, 30, 45, 50};
      if ( j < 2 ) {
	  bins[0] = -0.9;
	  bins[1] = -0.85;
	  bins[2] = -0.6;
	  bins[3] = 0.88;
	  bins[4] = 1.0;
	  if ( j == 1 ) bins[3] = 0.84;
      }  
      FormatHist(h1_met_in_ww[j], sample, Form("met_in_ww_%s", jetbin_names[j]), "MET for in-window for WW selection", nbins, bins);
      FormatHist(h1_met_in[j], sample, Form("met_in_%s", jetbin_names[j]), "MET for in-window", nbins, bins);
      FormatHist(h1_met_out[j], sample, Form("met_out_%s", jetbin_names[j]), "MET for out-window", nbins, bins);
      FormatHist(h1_met_in_sig[j], sample, Form("met_in_sig_%s", jetbin_names[j]), "MET in-Z with Signal Region", nbins, bins);
      FormatHist(h1_met_out_sig[j], sample, Form("met_out_sig_%s", jetbin_names[j]), "MET in Signal Region", nbins, bins);
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

	if (tree->tree_->GetBranchStatus("sfWeightPU"))
	  tree->tree_->SetBranchAddress("sfWeightPU", &sfWeightPU_);
	if (tree->tree_->GetBranchStatus("sfWeightTrig"))
	  tree->tree_->SetBranchAddress("sfWeightTrig", &sfWeightTrig_);
	if (tree->tree_->GetBranchStatus("sfWeightEff"))
	  tree->tree_->SetBranchAddress("sfWeightEff", &sfWeightEff_);
	// DY MVA
	// the default value is set to pass the dymva if this branch does not exist
	float dymva_ = 1.; 
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
	  if ( tree->njets_ == 3 )     njets = 2;
	  
	  float metValue;
	  if ( njets == 2 ) {
	    metValue = tree->met_;
	  }  else {
	    metValue = dymva_;
	  }
	  
	  //
	  // apply good run list and MC weights
	  //
	  
	  if (((1ll<<sample->getDataType()) & data_data) && tree->dstype_ == SmurfTree::data)  {
	    if (runlistIsSet_) {
	      if (!goodrun_json(tree->run_, tree->lumi_)) continue;
	    }
	  } else {
	    if (tree->type_ == 0) weight = tree->scale1fb_*(mmLumi_/1000.0);
	    if (tree->type_ == 3) weight = tree->scale1fb_*(eeLumi_/1000.0);
	    if (tree->type_ == 1 || tree->type_ == 2) weight = tree->scale1fb_*((eeLumi_+mmLumi_)/(2*1000.0));
	    weight = weight * sfWeightTrig_ * sfWeightEff_ * sfWeightPU_;
	  }
	  
	  // 
	  // apply all SF selections except for mT, mll and met
	  // 
	  
	  const unsigned int basic_selection =
	    SmurfTree::BaseLine |
	    SmurfTree::ChargeMatch |
	    SmurfTree::Lep1FullSelection |
	    SmurfTree::Lep2FullSelection |
	    SmurfTree::TopVeto |
	    SmurfTree::ExtraLeptonVeto;   
	  if ( (tree->cuts_ & basic_selection) != basic_selection ) continue;
	  if ( ! hww_pass_wwBaseline(tree,option_ ) ) continue; 
	  if ( TMath::Min(tree->pmet_, tree->pTrackMet_) < 20. ) continue;

	  // fill the histograms with loose selections for efficiency ratio evaluation k
	  if ( TMath::Abs(tree->dilep_.mass() - 91.1876) < 7.5) 
	    FillHist(h1_met_in_ww[njets], type, metValue, weight);

	  bool passDY = true; // toplogoical cuts based on dphi
	  bool passMET = true;
	  
	  if (  njets > 1 ) {//gc  || (option_ ==  WW_OPT_SMURFXSECSEL)
	    if ( ! hww_dy_selection(tree) ) passDY = false;
	    if ( ! hww_sfmet_selection(tree) )   passMET = false;
	  } else {
	    if ( ! (hww_sfdymva_selection(tree, dymva_)) ) passMET = false;
	  } 
	  if ( !passDY ) continue;
	  
	  //
	  // MC matching for the VV types
	  //

	  if ( tree->dstype_ == SmurfTree::wz || tree->dstype_ == SmurfTree::zz) {      
	    if ( tree->lep1MotherMcId_ != 23 || tree->lep2MotherMcId_ != 23) continue;
	  }
	  
	  // 
	  // apply all higgs mass dependent cuts EXCEPT mll/mT/dphi
	  // 
	  float lep1ptCut(20.), lep2ptCut(10.), mllCut(9999.), dPhiCut(9999.), mtLowCut(0.), mtHighCut(9999.),  mllLooseCut(9999.);
	  evaluate_hww_cuts(analysis_, lep1ptCut, lep2ptCut, mllCut, dPhiCut, mtLowCut, mtHighCut, mllLooseCut);
	  
	  // setting the bit for passVBF
	  bool passVBF = true;
	  if ( njets == 2) {
		  
		  if ( !hww_vbf_selection(tree) ) passVBF = false;
		  LorentzVector dijet = tree->jet1_ + tree->jet2_;
		  if (dijet.M() <= 500.0) passVBF = false;
		  if (fabs( tree->jet1_.Eta() - tree->jet2_.Eta() ) <= 3.5) passVBF = false;
	  }
/*	  
	  if ( (option_ == HWW_OPT_SMURFMVASEL) && analysis_ > 0. ) {
	    mtLowCut = 80.;
	    mtHighCut = analysis_;
	    mllCut = mllLooseCut; 
	  }
*/	  
	  if ( option_ == HWW_OPT_SMURFCUTSEL && analysis_ > 0. && njets == 2 ) {
	    mtLowCut = 30.;
	    mtHighCut = analysis_;
	    if ( analysis_ <= 200. ) 
	      mllCut = TMath::Min(float(100.), mllLooseCut);
	    else 
	      mllCut = mllLooseCut;
	  }

	  // for cut-based analysis
	  if ( option_ ==  HWW_OPT_SMURFCUTSEL ) {
	      
		  if ( tree->lep1_.Pt() < lep1ptCut) continue;
	      if ( tree->lep2_.Pt() < lep2ptCut) continue;
	      if ( (tree->dPhi_ * 180. / TMath::Pi() ) > dPhiCut) continue;
		 
		  if ( njets == 2 && analysis_ > 0.) {
	   	 	  if ( ! passVBF ) continue;
	   	  }
	  }
	  
	  // for the WW xec measurement
	  if ( analysis_ == 0  && tree->lep2_.Pt() < 20. && (option_ ==  WW_OPT_SMURFXSECSEL) ) continue;

      if ( option_ == HWW_OPT_SMURFPRESEL || (option_ ==  WW_OPT_SMURFXSECSEL) ) {//gc add SMURFXSECSEL
          if( tree->dilep_.Pt() < 45.) continue;
      }  

	  // 
	  // fill the in/out/sig met histograms
	  // 
	  
	  /*
	  // check the cuts
	  if ( njets == 2  ) 
	    std::cout << "mtLowCut " << mtLowCut << "\t mtHighCut = " << mtHighCut << "\t mll Cut = " << mllCut << "\n";
	  */

	  // fill the in-Z peak plots
	  if ( TMath::Abs(tree->dilep_.mass() - 91.1876) < 7.5) {
	    // for the Rout/in
	    FillHist(h1_met_in[njets], type, metValue, weight);
	    // for the estimation of Nin applying all signal selections
	    if ( tree->mt_ > mtLowCut && tree->mt_ < mtHighCut && passDY && passMET) {
		FillHist(h1_met_in_sig[njets], type, metValue, weight);	
	    }
	  }
	  
	  // fill the out region plots, defined as the signal selections backing out the mt selections
	  
	  if ( TMath::Abs(tree->dilep_.mass() - 91.1876) > 15  && tree->dilep_.mass() < mllCut) {
	    // for the Nout defintions for Rout/in
	    FillHist(h1_met_out[njets], type, metValue, weight);
	    
	    // for the Nout in the signal selections
	    // 0/1-Jet Selections depend on the analysis type
	    if ( njets < 2) {

	      // MVA selections, using DY MVA
          /*
	      if ( option_ == HWW_OPT_SMURFMVASEL ) {
		if ( hww_pass_mvaSelection(tree, analysis_, njets) && hww_pass_wwSelection(tree, option_) && passDY && passMET ) {
		  FillHist(h1_met_out_sig[njets], type, metValue, weight);
		  sample->fillResults(njets, tree->type_, weight, weight_err);
		}
	      } 
	     */ 
	      // Cut-based selections
	      if ( option_ == HWW_OPT_SMURFCUTSEL) {
		if ( hww_pass_cutSelection(tree, analysis_, njets) && hww_pass_wwSelection(tree, option_) && passDY && passMET ) {
		  FillHist(h1_met_out_sig[njets], type, metValue, weight);
		  sample->fillResults(njets, tree->type_, weight, weight_err);
		}
	      }

	      if ( option_ == WW_OPT_SMURFXSECSEL || option_ == HWW_OPT_SMURFPRESEL) {
		if ( hww_pass_wwSelection(tree, option_) && passDY && passMET ) {
		  FillHist(h1_met_out_sig[njets], type, metValue, weight);
		  sample->fillResults(njets, tree->type_, weight, weight_err);
		}
	      }
	    }

		// 2-Jet bin selections
		if ( njets == 2 ) {
			if ( option_ == HWW_OPT_SMURFCUTSEL && analysis_ > 0.) {
				if (  hww_pass_cutSelection(tree, analysis_, njets) && hww_pass_wwSelection(tree, option_) && passDY && passMET) {
					FillHist(h1_met_out_sig[njets], type, metValue, weight);
					sample->fillResults(njets, tree->type_, weight, weight_err);
				}
			}
			if ( analysis_ == 0.) {
				if (  hww_pass_wwSelection(tree, option_) && passDY && passMET) {
					FillHist(h1_met_out_sig[njets], type, metValue, weight);
					sample->fillResults(njets, tree->type_, weight, weight_err);
				}
			}
		}
	  } 
	}// end loop on events in file
	nf++;
        delete tree;
    } // end loop on files in chain
    gROOT->cd();
}



void SmurfDYLooper::setGoodRunList(const char *runlist)
{
    set_goodrun_file(runlist);
    runlistIsSet_ = true;
}

void SmurfDYLooper::setLumiScale(float eeLumi, float mmLumi)
{
    std::cout << "[SmurfDYLooper::setLumiScale] re-setting lumi to " << eeLumi << ", " << mmLumi << std::endl;
    eeLumi_ = eeLumi;
    mmLumi_ = mmLumi;

}
