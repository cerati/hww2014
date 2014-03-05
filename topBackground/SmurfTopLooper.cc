#include "SmurfTopLooper.h"

#include "../../../../../Smurf/Core/SmurfTree.h"
#include "../../../../../Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"
#include "../../../NtupleMacros/Tools/goodrun.h"
#include "../core/Selections.h"
#include "../core/SmurfPlotUtilities.h"
#include "../core/SmurfSample.h"
#include "../core/Enums.h"
#include "../SmurfScaleFactors.h"

#include "TChainElement.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom.h"

#include <cmath>
#include <iostream>
#include <cassert>

SmurfTopLooper::SmurfTopLooper(float analysis, Option option)
{
    std::cout << std::endl;
    std::cout << "[SmurfTopLooper::SmurfTopLooper] Doing analysis option " << option << std::endl;
    
    runlistIsSet_ = false;
    eeLumi_ = 4630;
    mmLumi_ = 4630;;
    option_ = option;
    analysis_ = analysis;
    std::cout << "[SmurfTopLooper::SmurfTopLooper] default lumi is " << eeLumi_ << ", " << mmLumi_ << std::endl;   
    loadWeightHistograms() ; 
}

void SmurfTopLooper::loop(SmurfSample *sample)
{
    TObjArray *listOfFiles = sample->getChain()->GetListOfFiles();
    TIter fileIter(listOfFiles);
    if (listOfFiles->GetEntries() == 0) {
        std::cout << "[SmurfTopLooper::loop] " << sample->getName() << " is not defined" << std::endl;
        return;
    }
    else
        std::cout << "[SmurfTopLooper::loop] " << sample->getName() << std::endl;

    //
    // set up histograms
    //
    
    int nbins(5);
    float xmin = 0.;
    float xmax = 200.;
    
    TH1F *h1_mll_sig[kJetBins][kLeptonTypes];
    TH1F *h1_mll_control[kJetBins][kLeptonTypes];
    TH1F *h1_mll_cali_denum[kJetBins][kLeptonTypes];
    TH1F *h1_mll_cali_num[kJetBins][kLeptonTypes];

    // central jet eta (for 2-jet only) 
    TH1F *h1_eta_cjet_sig[kLeptonTypes];
    TH1F *h1_eta_cjet_control[kLeptonTypes];
    TH1F *h1_eta_cjet_cali_denum[kLeptonTypes];
    TH1F *h1_eta_cjet_cali_num[kLeptonTypes];


    for (unsigned int j = 0; j < 3; ++j) {
      FormatHist(h1_mll_sig[j], sample, Form("mll_sig_%s_mH%.0f", jetbin_names[j], analysis_), "mll in Signal Region", nbins, xmin, xmax);
      FormatHist(h1_mll_control[j], sample, Form("mll_control_%s_mH%.0f", jetbin_names[j], analysis_), "mll in Control Region", nbins, xmin, xmax);
      FormatHist(h1_mll_cali_denum[j], sample, Form("mll_cali_denum_%s_mH%.0f", jetbin_names[j], analysis_), "mll in top calibration Region", nbins, xmin, xmax);
      FormatHist(h1_mll_cali_num[j], sample, Form("mll_cali_num_%s_mH%.0f", jetbin_names[j], analysis_), "mll in top calibration Region", nbins, xmin, xmax);
    }
    
    int etabins(5);
    float etaMin(0.), etaMax(2.5);
    FormatHist(h1_eta_cjet_sig, sample, Form("eta_cjet_sig_2j_mH%.0f", analysis_), "central jet #eta in Signal Region", etabins, etaMin, etaMax);
    FormatHist(h1_eta_cjet_control, sample, Form("eta_cjet_control_2j_mH%.0f", analysis_), "central jet #eta in Control Region", etabins, etaMin, etaMax);
    FormatHist(h1_eta_cjet_cali_denum, sample, Form("eta_cjet_cali_denum_2j_mH%.0f", analysis_), "central jet #eta in Calibration Denum Region", etabins, etaMin, etaMax);
    FormatHist(h1_eta_cjet_cali_num, sample, Form("eta_cjet_cali_num_2j_mH%.0f", analysis_), "central jet #eta in Calibration Num Region", etabins, etaMin, etaMax);
    
    
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
	float dymva_ = 1.;                                                                   
	if (tree->tree_->GetBranchStatus("sfWeightPU"))
	  tree->tree_->SetBranchAddress("sfWeightPU", &sfWeightPU_);
	if (tree->tree_->GetBranchStatus("sfWeightTrig"))
	  tree->tree_->SetBranchAddress("sfWeightTrig", &sfWeightTrig_);
	if (tree->tree_->GetBranchStatus("sfWeightEff"))
	  tree->tree_->SetBranchAddress("sfWeightEff", &sfWeightEff_);
        if (tree->tree_->GetBranchStatus("dymva"))
          tree->tree_->SetBranchAddress("dymva", &dymva_);

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
	  
	  if ( tree->njets_ == 3 ) njets = 2;
	  
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
	  
	  // 
	  // apply all ww selections except for the top veto
	  // 
	  
	  const unsigned int basic_selection =
	    SmurfTree::BaseLine |
	    SmurfTree::ChargeMatch |
	    SmurfTree::ZVeto |
	    SmurfTree::ExtraLeptonVeto;   
	  if ((tree->cuts_ & basic_selection) != basic_selection) continue;
	  if ( ! hww_pass_wwBaseline(tree, option_) ) continue;
	  if ( tree->type_ == 0 || tree->type_ == 3 ) {
	    if ( njets == 0 && dymva_ < 0.88 ) continue;   
	    if ( njets == 1 && dymva_ < 0.84 ) continue;    
	    if ( njets == 2 ) {
		    if (! hww_dy_selection(tree) ) continue;
		    if (! hww_sfmet_selection(tree) ) continue;
		}	
	  }
	  
	  // apply special ww selection value
	  if (  option_ == WW_OPT_SMURFXSECSEL &&  tree->lep2_.Pt() < 20. ) continue;
	  
	  // 
	  // apply rest of the selections
	  //
	  // flag for mc id of leptons : need this for Wjets
	  bool isRealLepton = false;
	  if((TMath::Abs(tree->lep1McId_) == 11 || TMath::Abs(tree->lep1McId_) == 13) &&
			  (TMath::Abs(tree->lep2McId_) == 11 || TMath::Abs(tree->lep2McId_) == 13)) isRealLepton = true;
	  
	  if ( sample->getDataType() == WJETSDATA ) {
	    
	    // apply pass fail selections
	    // skip events with no lepton pass the full selection
	    if ( ( (tree->cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)  &&  
		 ( (tree->cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ) continue;
	    
	    // skip events with both leptons hat pass the final selection
	    if ( ((tree->cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
		 && ((tree->cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection)) continue;
	    
	    // if lep1 pass full selection, but none of the lep2 pass the FO definition, skip the event
	    if ( tree->cuts_ & SmurfTree::Lep1FullSelection ) {
	      if ( ! ( (tree->cuts_ & SmurfTree::Lep2LooseEleV4) || (tree->cuts_ & SmurfTree::Lep2LooseMuV2) ) ) continue;
	    }
	    
	    // if lep2 pass full selection, but none of the lep1 pass the FO definition, skip the event
	    if ( tree->cuts_ & SmurfTree::Lep2FullSelection ) {
	      if ( ! ( (tree->cuts_ & SmurfTree::Lep1LooseEleV4) || (tree->cuts_ & SmurfTree::Lep1LooseMuV2) ) ) continue;
	    }
	    
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

	    } // end of wjets fake rate reading
	  } else {
	      const unsigned int lepton_selection = SmurfTree::Lep1FullSelection | SmurfTree::Lep2FullSelection;
	      if ((tree->cuts_ & lepton_selection) != lepton_selection) continue;
	  }


	  //
	  // apply necessary scale factors
	  // 
	  // apply the Z scale factors only to the ee/mm final states of the DY samples
	  if (((1ll<<sample->getDataType()) & (1ll<<ZLL)) && ( tree->type_ == 0 || tree->type_ == 3 ) ) {
	    // input by hand the scale factors
		  if ( option_ == WW_OPT_SMURFXSECSEL ) {
			  if ( njets == 0 ) 
				  weight = weight * 3.93;
			  if ( njets == 1 ) 
				  weight = weight * 3.59;
			  if ( njets == 2 ) 
				  weight = weight * 1.63;
		  } else 
	      weight = weight * ZScaleFactor_[njets]; 
	  }
	  
	  if ( tree->dstype_ == SmurfTree::wgstar) weight = weight * 1.5;                                                           
      if ( tree->dstype_ == SmurfTree::wgstar && tree->type_==3 && tree->met_>130) continue; 
	  
      if ( sample->getDataType() == ZTT) weight = weight * ZttScaleFactor(2, 1, 1, 1) / sfWeightPU_;  
	  
      //
	  // b-tag variables on top of the SmurfTree
	  // 
	  float btag_centraljet;
	  float btag_fwdjet;
	  float eta_centraljet;
	  
	  float btag_jet1 = tree->jet1Btag_;
	  float btag_jet2 = tree->jet2Btag_;
	  float btag_jet3 = tree->jet3Btag_;
	  const float btagCut = 2.1;
	  
	  //	  
	  LorentzVector dijet = tree->jet1_ + tree->jet2_; 
	  
	  if ( njets == 2) {
	    int index_centraljet = TMath::Abs(tree->jet1_.Eta()) < TMath::Abs(tree->jet2_.Eta()) ? 1 : 2;
	    if ( index_centraljet == 1) {
	      btag_centraljet = btag_jet1;
	      btag_fwdjet = btag_jet2;
	      eta_centraljet = TMath::Min(TMath::Abs(tree->jet1_.Eta()), 2.499);
	    } else {
	      btag_centraljet = btag_jet2;
	      btag_fwdjet = btag_jet1;
	      eta_centraljet = TMath::Min(TMath::Abs(tree->jet2_.Eta()), 2.499);
	    }
	  }

	  // 
	  // define signal selections
	  // (0, 2) jets the same as in the analysis
	  // 1-Jet: veto the events tagged without using the leading jet
	  //        veto events with leading jet tagged with b-tagging only
	  // 
	  
	  float lep1ptCut(20.), lep2ptCut(10.), mllCut(9999.), dPhiCut(9999.), mtLowCut(0.), mtHighCut(9999.),  mllLooseCut(9999.);
	  evaluate_hww_cuts(analysis_, lep1ptCut, lep2ptCut, mllCut, dPhiCut, mtLowCut, mtHighCut, mllLooseCut);

	  bool top_sig_selection = true;

	  if ( njets == 0 ) {
		  if ( ! hww_topveto(tree) ) top_sig_selection = false;
	  } else if ( njets == 1 ) {
		  if ( tree->cuts_ & SmurfTree::TopTagNotInJets) top_sig_selection = false;
		  if ( btag_jet1 > btagCut)   top_sig_selection = false;
	  } else if ( njets == 2 ) {
		  if ( option_ == HWW_OPT_SMURFCUTSEL) {
			  if ( ! hww_vbf_selection(tree)) top_sig_selection = false;
			  if (dijet.M() <= 500.0) top_sig_selection = false;
			  if (fabs( tree->jet1_.Eta() - tree->jet2_.Eta() ) <= 3.5) top_sig_selection = false;
			  if ( analysis_ > 0. ) {
				  if ( analysis_ <= 200. && tree->dilep_.mass() >= 100.) top_sig_selection = false;
				  if ( tree->dilep_.mass() > mllLooseCut) top_sig_selection = false;
				  if ( tree->mt_ < 30.) top_sig_selection = false;
				  if ( tree->mt_ > analysis_) top_sig_selection = false;
			  }
		  }
		  if ( btag_jet1 > btagCut  || btag_jet2 > btagCut ) top_sig_selection = false;
		  if ( tree->njets_  == 3 && btag_jet3 > btagCut ) top_sig_selection = false; 
		  if ( tree->cuts_ & SmurfTree::TopTagNotInJets) top_sig_selection = false;
		  if ( tree->nSoftMuons_ > 0 ) top_sig_selection = false; 
	  }
	  
	  // 
	  // define control region
	  // 
	  
	  bool top_control_selection = true;
	  if ( njets == 0 ) {
		  if ( tree->cuts_ & SmurfTree::TopVeto)  top_control_selection = false;
	  } else if ( njets == 1 ) {
		  if ( tree->cuts_ & SmurfTree::TopTagNotInJets)  top_control_selection = false;
		  if ( btag_jet1 < btagCut) top_control_selection = false; 
	  }  else if ( njets == 2 ) { // same as signal except that central jet being tagged
		  if ( option_ == HWW_OPT_SMURFCUTSEL) {
			  if ( ! hww_vbf_selection(tree)) top_control_selection = false;
			  if (dijet.M() <= 500.0) top_control_selection =  false;
			  if (fabs( tree->jet1_.Eta() - tree->jet2_.Eta() ) <= 3.5) top_control_selection = false;
			  if ( analysis_ > 0. ) {
				  if ( analysis_ <= 200. && tree->dilep_.mass() >= 100.) top_control_selection = false;
				  if ( tree->dilep_.mass() > mllLooseCut) top_control_selection = false;
				  if ( tree->mt_ < 30.) top_control_selection = false;
				  if ( tree->mt_ > analysis_) top_control_selection = false;
			  }
		  }
		  if ( btag_fwdjet > btagCut ) top_control_selection = false;
		  if ( btag_centraljet < btagCut) top_control_selection = false;   
		  if ( tree->cuts_ & SmurfTree::TopTagNotInJets) top_control_selection = false;
		  if ( tree->njets_  == 3 && btag_jet3 > btagCut ) top_control_selection = false; 
		  if ( tree->nSoftMuons_ > 0 )  top_control_selection = false; 
	  }
	  
	  // 
	  // define calibration selections to get the top-tagging efficiency
	  // denom vs num
	  // 

	  // Measure 0-Jet top-tag efficiency from 1-jet events
	  bool top_cali_denum_0j(false), top_cali_num_0j(false);
	  if ( tree->njets_ == 1 ) {
	    if ( tree->cuts_ & SmurfTree::OneBJet )     {
	      top_cali_denum_0j = true;
	      if ( tree->cuts_ & SmurfTree::TopTagNotInJets ) 
		top_cali_num_0j = true;
	    }
	  }
	  
	  // Measure 1-Jet top-tag efficiency from 2-jet events
	  // make use of the *highest jet pT tagging*
	  bool top_cali_denum_1j(false), top_cali_num_1j(false);
	  if ( tree->njets_ == 2 ) {
		  if ( ! (tree->cuts_ & SmurfTree::TopTagNotInJets) && (btag_jet2 > btagCut)) {
			  top_cali_denum_1j = true;
			  if ( btag_jet1 > btagCut) top_cali_num_1j = true;
		  }
	  }
	  
	  // Measure 2-Jet top-tag efficiency from 2-jet events
	  // Maasure the central jet tagging efficiency
	  // Get 2-jet events where the forward jet is tagged after vbf cuts
	  
	  bool top_cali_denum_2j(false), top_cali_num_2j(false);
	  if ( njets == 2 &&  ( tree->njets_ == 2 || btag_jet3 < btagCut ) && tree->nSoftMuons_ == 0 ) {
	    if ( !( (tree->cuts_) & SmurfTree::TopTagNotInJets) && btag_fwdjet < btagCut) {
	      top_cali_denum_2j = true;
	      if ( btag_centraljet > btagCut ) 
		top_cali_num_2j = true;
	    }
	  }
	  
	  //
	  // fill the histograms
	  // 
	  
	  if ( top_sig_selection) {
	    FillHist(h1_mll_sig[njets], type, tree->dilep_.mass(), weight);
	    if ( njets == 2) 
	      FillHist(h1_eta_cjet_sig, type, eta_centraljet, weight);
	    sample->fillResults(njets, tree->type_, weight, weight_err);
	  }
	  
	  if ( top_control_selection) {
	    FillHist(h1_mll_control[njets], type, tree->dilep_.mass(), weight);
	    if ( njets == 2 ) {
	      FillHist(h1_eta_cjet_control, type, eta_centraljet, weight);
	      // if ( sample->getDataType() == DATA ) 
	      // std::cout << tree->event_ << " " << eta_centraljet << "\n";
	    }
	  }
	  
	  if ( top_cali_denum_0j) 
	    FillHist(h1_mll_cali_denum[0], type, tree->dilep_.mass(), weight);
	  if ( top_cali_num_0j) 
	    FillHist(h1_mll_cali_num[0], type, tree->dilep_.mass(), weight);

	  if ( top_cali_denum_1j) 
	    FillHist(h1_mll_cali_denum[1], type, tree->dilep_.mass(), weight);
	  if ( top_cali_num_1j) 
	    FillHist(h1_mll_cali_num[1], type, tree->dilep_.mass(), weight);
	  
	  if ( top_cali_denum_2j) {
	    FillHist(h1_eta_cjet_cali_denum, type, eta_centraljet, weight);
	  }
	  if ( top_cali_num_2j) {
	    FillHist(h1_eta_cjet_cali_num, type, eta_centraljet, weight);
	  }


	} // end loop on events in file
	nf++;
        delete tree;
    } // end loop on files in chain
    gROOT->cd();
}



void SmurfTopLooper::setGoodRunList(const char *runlist)
{
    set_goodrun_file(runlist);
    runlistIsSet_ = true;
}

void SmurfTopLooper::setLumiScale(float eeLumi, float mmLumi)
{
    std::cout << "[SmurfTopLooper::setLumiScale] re-setting lumi to " << eeLumi << ", " << mmLumi << std::endl;
    eeLumi_ = eeLumi;
    mmLumi_ = mmLumi;

}

void SmurfTopLooper::loadWeightHistograms() 
{
     
    //
    // Fakerate
    //
  
    TFile *fLeptonFRFileM = 0;
    TFile *fLeptonFRFileE = 0;
    
    // load electron fake rate histograms
    //fLeptonFRFileE = TFile::Open("/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root");
    fLeptonFRFileE = TFile::Open("/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/aux/summary_fakes_Moriond2012.root");
    fhDFREl_ = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
    assert(fhDFREl_);
    fhDFREl_->SetDirectory(0);
    fLeptonFRFileE->Close();
    delete fLeptonFRFileE;

    // load muon fake rate histograms
    //fLeptonFRFileM = TFile::Open("/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root");
    fLeptonFRFileM = TFile::Open("/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/aux/summary_fakes_Moriond2012.root");
    fhDFRMu_ = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold30_PtEta"));
    assert(fhDFRMu_);
    fhDFRMu_->SetDirectory(0);
    fLeptonFRFileM->Close();
    delete fLeptonFRFileM; 
    
    getZScaleFactor(ZScaleFactor_, ZScaleFactorError_, HWW_OPT_SMURFCUTSEL, 0, "sf");

}
