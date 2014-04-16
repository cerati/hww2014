#include "Selections.h"
#include "../../Smurf/Core/SmurfTree.h"

#include "TMath.h"
#include <iostream> 
using namespace std; 

//
// Top veto 
//
bool hww_topveto(SmurfTree *tree) 
{
    if ( (tree->cuts_ & SmurfTree::TopVeto) != SmurfTree::TopVeto ) return false; 
    return true;
}

//
// Selection for SF events 
//
bool hww_dy_selection(SmurfTree *tree) 
{
    if ( tree->njets_ < 2 && tree->jet1_.Pt() > 15.0 && tree->dPhiDiLepJet1_ > (165./180.)*TMath::Pi())  return false;
    if ( tree->njets_ >= 2 && ( fabs(ROOT::Math::VectorUtil::DeltaPhi((tree->jet1_+tree->jet2_),tree->dilep_))*180.0/TMath::Pi() > 165.)) return false; 
    return true;
}

//
// MET selection for SF events : used for 2jet
//
bool hww_sfmet_selection(SmurfTree *tree) 
{
    if ( tree->njets_ < 2 ) {
        if ( std::min(tree->pmet_, tree->pTrackMet_)  <= 45. )  return false;
    } else {
        if ( tree->met_ < 45. ) return false;
    }
    return true;
}

//
// DY MVA selection for SF events
//
bool hww_sfdymva_selection(SmurfTree *tree, const float& dymva)
{

// Relaxed MET requirement for SS mm test : threshold to kill DY completely
// DY template has very large weights with very little entries 
// So, it's better to kill this process completely
//  if ( std::min(tree->pmet_, tree->pTrackMet_)  <= 35. )  return false; 

    if ( tree->njets_ == 0  && dymva < 0.88 ) return false;  
    if ( tree->njets_ == 1  && dymva < 0.84 ) return false;

    return true;
}

bool hww_vbf_selection(SmurfTree *tree)
{
    if ( tree->njets_ == 2 || tree->njets_ == 3 ) {
        // centrality cuts
        float largestEta = tree->jet1_.Eta();
        float smallestEta = tree->jet2_.Eta();
        if (tree->jet2_.Eta() > tree->jet1_.Eta()) {
            largestEta = tree->jet2_.Eta();
            smallestEta = tree->jet1_.Eta();
        }
        if (! (tree->lep1_.Eta() > smallestEta && tree->lep1_.Eta() < largestEta)) return false;
        if (! (tree->lep2_.Eta() > smallestEta && tree->lep2_.Eta() < largestEta)) return false;
    }
    return true;
}

bool hww_pass_wwBaseline(SmurfTree *tree, const Option option) 
{
    //if (tree->lq1_ * tree->lq2_ > 0 ) return false; // comment out for SS selections
    if (tree->lep1_.Pt() < 20. || tree->lep2_.Pt() < 10.) return false;
    if ( tree->dilep_.Pt() <=  30)  return false;   
    if ( tree->type_ == 0 || tree->type_ == 3 ) { 
        if ( tree->dilep_.Pt() <=  45)  return false;   
    }
    if ( option != XWW_OPT_MT2DMLL_JCP && option != HWW_OPT_MT2DMLL_JCP && option != HWW_OPT_MT2DMLL
	 && option !=HWW_OPT_SSCTL2D && option != HWW_OPT_SMURFPRESEL && option !=WW_OPT_SMURFXSECSEL)//gc add SMURFXSECSEL
        if ( tree->dilep_.Pt() <=  45)  return false;  
    if (tree->dilep_.mass() <= 12.0) return false; 
    if (std::min(tree->pmet_, tree->pTrackMet_) < 20.) return false;

    if (tree->njets_ > 3) return false;
    if (tree->dstype_ == SmurfTree::data) {
        if ( (tree->cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger ) return false;
    }

    if ( tree->njets_ == 3 ) {
        // 3-jet veto
        float largestEta = tree->jet1_.Eta();
        float smallestEta = tree->jet2_.Eta();
        if (tree->jet2_.Eta() > tree->jet1_.Eta()) {
            largestEta = tree->jet2_.Eta();
            smallestEta = tree->jet1_.Eta();
        }
        if ( tree->jet3_.Eta() > smallestEta && tree->jet3_.Eta() < largestEta ) return false; 
    }
    return true;
}

// Full ID plus the preliminary cuts
bool hww_pass_wwSelection(SmurfTree *tree, const Option option)  
{   
    const unsigned int basic_selection =  
        SmurfTree::BaseLine |
        SmurfTree::ChargeMatch |
        SmurfTree::Lep1FullSelection |
        SmurfTree::Lep2FullSelection |
        SmurfTree::ZVeto |
        //        SmurfTree::TopVeto |
        SmurfTree::ExtraLeptonVeto;   
    if ((tree->cuts_ & basic_selection) != basic_selection) return false;
    if ( ! hww_pass_wwBaseline(tree, option) ) return false; 
    if ( option==HWW_OPT_TOPTAG ){ 
        if( hww_topveto(tree) ) return false; 
    } else {
        if ( ! hww_topveto(tree) ) return false;
    }
    return true;
}

// Same sign with full ID plus the preliminary cuts
bool hww_pass_wwSSSelection(SmurfTree *tree, const Option option)  
{   
    const unsigned int basic_selection =  
        SmurfTree::BaseLine |
        //SmurfTree::ChargeMatch |
        SmurfTree::Lep1FullSelection |
        SmurfTree::Lep2FullSelection |
        SmurfTree::ZVeto |              // Remove Zveto for SS mm test 
        //        SmurfTree::TopVeto |
        SmurfTree::ExtraLeptonVeto;   

    if ( (tree->cuts_ & basic_selection) != basic_selection ) return false;
    if ( tree->lq1_ * tree->lq2_ < 0 ) return false;
    if ( ! hww_pass_wwBaseline(tree, option) ) return false;
    if ( ! hww_topveto(tree) ) return false;
    return true;
}

// Wgamma FO + the preliminary cuts for opposite sign events 
bool hww_pass_wgammaSelection(SmurfTree *tree, const Option option)  
{   
    const unsigned int basic_selection =   
        SmurfTree::BaseLine |
        SmurfTree::ChargeMatch |
        SmurfTree::ZVeto |
        //        SmurfTree::TopVeto |
        SmurfTree::ExtraLeptonVeto;   
    if ((tree->cuts_ & basic_selection) != basic_selection) return false;   
    // require Tight + FO
    if ( ! ((tree->cuts_ & SmurfTree::Lep1LooseEleV2) && (tree->cuts_ & SmurfTree::Lep2FullSelection)) && 
         ! ((tree->cuts_ & SmurfTree::Lep2LooseEleV2) && (tree->cuts_ & SmurfTree::Lep1FullSelection)) ) return false;
    if ( ! hww_pass_wwBaseline(tree, option) ) return false;
    if ( option==HWW_OPT_TOPTAG ){ 
        if( hww_topveto(tree) ) return false; 
    } else {
        if ( ! hww_topveto(tree) ) return false;
    }
    return true;
}

// Wgamma FO + the preliminary cuts for same sign events
bool hww_pass_wgammaSSSelection(SmurfTree *tree, const Option option)  
{   
    const unsigned int basic_selection =   
        SmurfTree::BaseLine |
        //SmurfTree::ChargeMatch |
        SmurfTree::ZVeto |              // Remove Zveto for SS mm test  
        //        SmurfTree::TopVeto |
        SmurfTree::ExtraLeptonVeto;   
    if ((tree->cuts_ & basic_selection) != basic_selection) return false;   
    if ( tree->lq1_ * tree->lq2_ < 0 ) return false;
    // require Tight + FO
    if ( ! ((tree->cuts_ & SmurfTree::Lep1LooseEleV2) && (tree->cuts_ & SmurfTree::Lep2FullSelection)) && 
         ! ((tree->cuts_ & SmurfTree::Lep2LooseEleV2) && (tree->cuts_ & SmurfTree::Lep1FullSelection)) ) return false;
    if ( ! hww_pass_wwBaseline(tree, option) ) return false;
    if ( ! hww_topveto(tree) ) return false;
    return true;
}

// Pass+Fail Lepton selections plus the full MET
bool hww_pass_wwPassFailSelection(SmurfTree *tree, const Option option)
{
    const unsigned int basic_selection =
        SmurfTree::BaseLine |
        SmurfTree::ChargeMatch |
        SmurfTree::ZVeto |
        // SmurfTree::TopVeto |
        SmurfTree::ExtraLeptonVeto;   
    if ((tree->cuts_ & basic_selection) != basic_selection) return false;
    if ( ! hww_pass_wwBaseline(tree, option) ) return false;
    if ( option==HWW_OPT_TOPTAG ){
        if( hww_topveto(tree) ) return false;
    } else {
        if ( ! hww_topveto(tree) ) return false;
    }
    // Pass fail selections
    // skip events with no lepton pass the full selection
    if ( ( (tree->cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)  &&  
            ( (tree->cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ) return false;

    // skip events with both leptons passing the final selection
    if ( ((tree->cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
            && ((tree->cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))  return false;

    // if lep1 passes the full selection, but lep2 does not pass none of the FO definitions, skip the event
    if ( tree->cuts_ & SmurfTree::Lep1FullSelection ) {
        if ( ! ( (tree->cuts_ & SmurfTree::Lep2LooseEleV4) || (tree->cuts_ & SmurfTree::Lep2LooseMuV2) ) ) return false;
    }

    // if lep2 passes the full selection, but lep1 does not pass none of the FO definitions, skip the event
    if ( tree->cuts_ & SmurfTree::Lep2FullSelection ) {
        if ( ! ( (tree->cuts_ & SmurfTree::Lep1LooseEleV4) || (tree->cuts_ & SmurfTree::Lep1LooseMuV2) ) ) return false;
    }

    return true;
}

// Same Sign with Pass+Fail Lepton selections plus the full MET
bool hww_pass_wwSSPassFailSelection(SmurfTree *tree, const Option option)
{
    const unsigned int basic_selection =
        SmurfTree::BaseLine |
        //SmurfTree::ChargeMatch |
        SmurfTree::ZVeto |              // Remove Zveto for SS mm test  
        // SmurfTree::TopVeto |
        SmurfTree::ExtraLeptonVeto;    
    if ( (tree->cuts_ & basic_selection) != basic_selection ) return false;
    if ( tree->lq1_ * tree->lq2_ < 0 ) return false;
    if ( ! hww_pass_wwBaseline(tree, option) ) return false;
    if ( ! hww_topveto(tree) ) return false;

    // Pass fail selections
    // skip events with no lepton pass the full selection
    if ( ( (tree->cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)  &&  
            ( (tree->cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ) return false;

    // skip events with both leptons hat pass the final selection
    if ( ((tree->cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection)
            && ((tree->cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))  return false;

    // if lep1 pass full selection, but none of the lep2 pass the FO definition, skip the event
    if ( tree->cuts_ & SmurfTree::Lep1FullSelection ) {
        if ( ! ( (tree->cuts_ & SmurfTree::Lep2LooseEleV4) || (tree->cuts_ & SmurfTree::Lep2LooseMuV2) ) ) return false;
    }

    // if lep2 pass full selection, but none of the lep1 pass the FO definition, skip the event
    if ( tree->cuts_ & SmurfTree::Lep2FullSelection ) {
        if ( ! ( (tree->cuts_ & SmurfTree::Lep1LooseEleV4) || (tree->cuts_ & SmurfTree::Lep1LooseMuV2) ) ) return false;
    }

    return true;
}

void evaluate_hww_cuts(const float &analysis, float &lep1ptCut, float &lep2ptCut, float &mllCut, float &dPhiCut, float &mtLowCut, float &mtHighCut, float &mllLooseCut) 
{

    if ( analysis == 0.) return;
    // default cuts 
    mtLowCut = 80.;
    mtHighCut = analysis;
    mllLooseCut = analysis;

    // differences to the 3 default cuts are specified below

    if ( analysis <= 115.) {
        lep1ptCut = 20.;
        lep2ptCut = 10.;
        mllCut = 40;
        dPhiCut = 115;
        mtHighCut = 110;
        mllLooseCut = 70;
        return;
    }

    if ( analysis == 118.) {
        lep1ptCut = 20.;
        lep2ptCut = 10.;
        mllCut = 40;
        dPhiCut = 115;
        mtHighCut = 115;
        mllLooseCut = 70;
        return;
    }

    if ( analysis == 120.) {
        lep1ptCut = 20.;
        lep2ptCut = 10.;
        mllCut = 40;
        dPhiCut = 115;
        mtHighCut = 120;
        mllLooseCut = 70;
        return;
    }

    if ( analysis == 122.) {
        lep1ptCut = 21.;
        lep2ptCut = 10.;
        mllCut = 41;
        dPhiCut = 110;
        mtHighCut = 121;
        mllLooseCut = 70;
        return;
    }


    if ( analysis == 124. ) {
        lep1ptCut = 22.;
        lep2ptCut = 10.;
        mllCut = 42;
        dPhiCut = 105;
        mtHighCut = 122;
        mllLooseCut = 70;
        return;
    }

    if ( analysis == 126. || analysis == 125) {
        lep1ptCut = 23.;
        lep2ptCut = 10.;
        mllCut = 43;
        dPhiCut = 100;
        mtHighCut = 123;  
        mllLooseCut = 80;
        return;
    }

    if ( analysis == 128.) {
        lep1ptCut = 24.;
        lep2ptCut = 10.;
        mllCut = 44;
        dPhiCut = 95;
        mtHighCut = 124;
        mllLooseCut = 80;
        return;
    }

    if ( analysis == 130.) {
        lep1ptCut = 25.;
        lep2ptCut = 10.;
        mllCut = 45;
        dPhiCut = 90;
        mtHighCut = 125;
        mllLooseCut = 80;
        return;
    }

    if ( analysis == 135.) {
        lep1ptCut = 25.;
        lep2ptCut = 12.;
        mllCut = 45;
        dPhiCut = 90;
        mtHighCut = 128;
        mllLooseCut = 90;
        return;
    }

    if ( analysis == 140.) {
        lep1ptCut = 25.;
        lep2ptCut = 15.;
        mllCut = 45;
        dPhiCut = 90;
        mtHighCut = 130;
        mllLooseCut = 90;
        return;
    }

    if ( analysis == 145.) {
        lep1ptCut = 25.;
        lep2ptCut = 15.;
        mllCut = 45;
        dPhiCut = 90;
        mtHighCut = 130;
        mllLooseCut = 100;
        return;
    }

    if ( analysis == 150. || analysis == 155.) {
        lep1ptCut = 27.;
        lep2ptCut = 25.;
        mllCut = 50;
        dPhiCut = 90;
        mllLooseCut = 100;
        return;
    }

    if ( analysis == 160.) {
        lep1ptCut = 30.;
        lep2ptCut = 25.;
        mllCut = 50;
        dPhiCut = 60;
        mtLowCut = 90.;
        mllLooseCut = 100;
        return;
    }

    if ( analysis == 170.) {
        lep1ptCut = 34.;
        lep2ptCut = 25.;
        mllCut = 50;
        dPhiCut = 60;
        mtLowCut = 110.;
        mllLooseCut = 100;
        return;
    }

    if ( analysis == 180.) {
        lep1ptCut = 36.;
        lep2ptCut = 25.;
        mllCut = 60;
        dPhiCut = 70;
        mtLowCut = 120.;
        mllLooseCut = 110;
        return;
    }

    if ( analysis == 190.) {
        lep1ptCut = 38.;
        lep2ptCut = 25.;
        mllCut = 80;
        dPhiCut = 90;
        mtLowCut = 120.;
        mllLooseCut = 120;
        return;
    }

    if ( analysis == 200.) {
        lep1ptCut = 40.;
        lep2ptCut = 25.;
        mllCut = 90;
        dPhiCut = 100;
        mtLowCut = 120.;
        mllLooseCut = 130;
        return;
    }

    if ( analysis == 250.) {
        lep1ptCut = 55.;
        lep2ptCut = 25.;
        mllCut = 150;
        dPhiCut = 140;
        mtLowCut = 120.;
        return;
    }

    if ( analysis == 300.) {
        lep1ptCut = 70.;
        lep2ptCut = 25.;
        mllCut = 200;
        dPhiCut = 175;
        mtLowCut = 120.;
        return;
    }

    if ( analysis == 350.) {
        lep1ptCut = 80.;
        lep2ptCut = 25.;
        mllCut = 250;
        dPhiCut = 175;
        mtLowCut = 120.;
        return;
    }

    if ( analysis == 400.) {
        lep1ptCut = 90.;
        lep2ptCut = 25.;
        mllCut = 300;
        dPhiCut = 175;
        mtLowCut = 120.;
        return;
    }

    if ( analysis == 450.) {
        lep1ptCut = 110.;
        lep2ptCut = 25.;
        mllCut = 350;
        dPhiCut = 175;
        mtLowCut = 120.;
        return;
    }

    if ( analysis == 500.) {
        lep1ptCut = 120.;
        lep2ptCut = 25.;
        mllCut = 400;
        dPhiCut = 175;
        mtLowCut = 120.;
        return;
    }

    if ( analysis == 550.) {
        lep1ptCut = 130.;
        lep2ptCut = 25.;
        mllCut = 450;
        dPhiCut = 175;
        mtLowCut = 120.;
        return;
    }

    if ( analysis == 600.) {
        lep1ptCut = 140.;
        lep2ptCut = 25.;
        mllCut = 500;
        dPhiCut = 175;
        mtLowCut = 120.;
        return;
    }
}

// only the higgs related cuts, to be used with 
// ww, ww_PassFail, ww_LooseMET together to complete the cuts
// WARNING, no higgs level SF cuts defined here, they need to applied at the looper level

bool hww_pass_cutSelection(SmurfTree *tree, const float &analysis, const unsigned int jetbin)
{  
    float lep1ptCut(20.), lep2ptCut(10.), mllCut(9999.), dPhiCut(9999.), mtLowCut(0.), mtHighCut(9999.),  mllLooseCut(9999.);
    evaluate_hww_cuts(analysis, lep1ptCut, lep2ptCut, mllCut, dPhiCut, mtLowCut, mtHighCut, mllLooseCut);
 
    // cuts applied for 0/1 Jets
    if (jetbin < 2) {
        if ( tree->lep1_.Pt() < lep1ptCut) return false;
        if ( tree->lep2_.Pt() < lep2ptCut) return false;
        if ( tree->dilep_.mass() > mllCut) return false;
        if ( (tree->dPhi_ * 180. / TMath::Pi() ) > dPhiCut) return false;
        if ( tree->mt_ < mtLowCut) return false;
        if ( tree->mt_ > mtHighCut) return false;
    }

    if ( jetbin == 2 ) { 
        // 0/1jet cuts for VBF   
        if ( tree->lep1_.Pt() < lep1ptCut) return false;
        if ( tree->lep2_.Pt() < lep2ptCut) return false;
        if ( tree->dilep_.mass() > mllCut) return false;
        if ( (tree->dPhi_ * 180. / TMath::Pi() ) > dPhiCut) return false;
        if ( tree->mt_ > mtHighCut) return false; 

        if ( ! hww_vbf_selection(tree)) return false;
        LorentzVector dijet = tree->jet1_ + tree->jet2_;
        if (dijet.M() <= 500.0) return false;
        if (fabs( tree->jet1_.Eta() - tree->jet2_.Eta() ) <= 3.5) return false;
        if ( analysis <= 200. && tree->dilep_.mass() >= 100.) return false;
        if ( tree->dilep_.mass() > mllLooseCut) return false;
        if ( tree->mt_ < 30.) return false;
        if ( tree->mt_ > analysis) return false;
    }
    
    return true;
}

bool hww_pass_mvaSelection(SmurfTree *tree, const float &analysis, const unsigned int jetbin)
{ 
    float lep1ptCut(20.), lep2ptCut(10.), mllCut(9999.), dPhiCut(9999.), mtLowCut(0.), mtHighCut(9999.),  mllLooseCut(9999.);
    evaluate_hww_cuts(analysis, lep1ptCut, lep2ptCut, mllCut, dPhiCut, mtLowCut, mtHighCut, mllLooseCut);
    if ( jetbin < 2) {  

        // new BDT : same [mT, mll] range as 2D
        if(analysis < 300) {
            if ( tree->mt_ < 80) return false;   
            if ( tree->mt_ > 280) return false;     
            if ( tree->dilep_.mass() > 200) return false;   
        } else {
            if ( tree->mt_ < 80) return false;   
            if ( tree->mt_ > 600) return false;
            if ( tree->dilep_.mass() > 600) return false;
            if ( tree->lep1_.Pt() < 50) return false;
        }

        // ICHEP BDT
        //if ( tree->mt_ < 80) return false;
        //if ( tree->mt_ > analysis) return false;   
        //if ( tree->dilep_.mass() > mllLooseCut) return false;  
    }

    if ( jetbin == 2 ) {
        if ( ! hww_vbf_selection(tree)) return false;
        if ( analysis <= 200. && tree->dilep_.mass() >= 100.) return false;
        if ( tree->dilep_.mass() > mllLooseCut) return false;
        if ( tree->mt_ < 30.) return false;
        if ( tree->mt_ > analysis) return false;
    }
    return true;                        
}

bool hww_pass_2DSelection(SmurfTree *tree, const float &analysis, const unsigned int jetbin,  const Option option)
{ 
    float lep1ptCut(20.), lep2ptCut(10.), mllCut(9999.), dPhiCut(9999.), mtLowCut(0.), mtHighCut(9999.),  mllLooseCut(9999.);
    evaluate_hww_cuts(analysis, lep1ptCut, lep2ptCut, mllCut, dPhiCut, mtLowCut, mtHighCut, mllLooseCut);

    //  for the nominal analysis 
    if ( option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D ) {
      if ( jetbin<2 ) {
        if ( tree->mt_ < 60) return false; 
      } else {
        if ( tree->mt_ < 30) return false;   
      }
      
      // different cuts for mH < 300 and mH >= 300
      if(analysis < 300) {
        if ( tree->mt_ > 280) return false;     
        if ( tree->dilep_.mass() > 200) return false;   
        //if ( tree->mt_ < 80) return false;            // for HCP2012 analysis
        //if ( tree->dilep_.Pt() < 45) return false;    // for HCP2012 analysis
      } else {
        if ( tree->mt_ > 600) return false;
        if ( tree->dilep_.mass() > 600) return false;
        if ( tree->lep1_.Pt() < 50) return false;
        //if ( tree->dilep_.Pt() < 45) return false; // to sync with Guillelmo
      } 
      
      if ( jetbin == 2 ) {
        if ( ! hww_vbf_selection(tree)) return false; 
        LorentzVector dijet = tree->jet1_ + tree->jet2_;
        if (dijet.M() <= 500.0) return false;
        if (fabs( tree->jet1_.Eta() - tree->jet2_.Eta() ) <= 3.5) return false;
      }
    }

    //  for the JCP studies
    if ( option == HWW_OPT_MT2DMLL_JCP || option == XWW_OPT_MT2DMLL_JCP ) {
      if ( jetbin > 1 ) return false; 
      if ( tree->mt_ < 60. ) return false;
      if ( tree->mt_ > 280. ) return false;
      if ( tree->dilep_.mass() > 200 ) return false;
      //if ( tree->mt_ > 200. ) return false;
      //if ( tree->dilep_.mass() > 150. ) return false;
    }
    
    return true;                        
}

bool hww_assign_this_event(SmurfTree *tree, DataType dataType)
{
    // ZLL collects only the SF of the dyee and dymm
    if (dataType == ZLL) {
        if ( ! (tree->dstype_ == SmurfTree::dyee || tree->dstype_ == SmurfTree::dymm)) return false;
        if ( tree->type_ == 1 || tree->type_ == 2 ) return false;
    }

    // Zjets includes the OF of Z->ee/mm and the peaking component of the WZ/ZZ 
    if (dataType == ZJETS) {
        return false; 
        /*
        if ( tree->dstype_ == SmurfTree::dyee || tree->dstype_ == SmurfTree::dymm) {
            if ( tree->type_ == 0 || tree->type_ == 3 ) return false;
        }
        if ( (tree->dstype_ == SmurfTree::wz) || (tree->dstype_ == SmurfTree::zz) ){
            if (tree->lep1MotherMcId_ != 23 || tree->lep2MotherMcId_ !=  23)   return false;
        }
        */
    }

    // VV means both leptons must not 
    // come from the same Z
    if (dataType == VV) {
        //if (tree->lep1MotherMcId_ == 23 && tree->lep2MotherMcId_ == 23) return false;
    }
    
    //
    // Wjets data
    //
    if (dataType == WJETSELEDATA || dataType == WJETSMUDATA || dataType == WJETSDATA  ) {
        if ( tree->dstype_ != SmurfTree::data && tree->dstype_ != SmurfTree::wgamma && tree->lep1McId_ * tree->lep2McId_ == 0) return false;
    }

    if (dataType == WJETSELEDATA) {
        if( !(tree->cuts_ & SmurfTree::Lep1LooseEleV4 || tree->cuts_ & SmurfTree::Lep2LooseEleV4)  )  return false;
    }

    if (dataType == WJETSMUDATA) {
        if( !(tree->cuts_ & SmurfTree::Lep1LooseMuV2 || tree->cuts_ & SmurfTree::Lep2LooseMuV2)  )  return false;
    }

    return true;
}

