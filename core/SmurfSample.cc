
#include "TROOT.h"

#include "SmurfSample.h"
#include "../../Smurf/Core/SmurfTree.h"
#include "core/Selections.h"

SmurfSample::SmurfSample()
{
}

SmurfSample::~SmurfSample()
{
    for (unsigned int j = 0; j < kJetBins; ++j) {
        for (unsigned int i = 0; i < kLeptonTypes; ++i) {
//            if (h1k_shape_[j][i]) delete h1k_shape_[j][i];
        }
    }
}

SmurfSample::SmurfSample(Option option, DataType dataType, Color_t colour, std::string name, float analysis)
{

    dataType_ = dataType;
    name_ = name;
    colour_ = colour;
    chain_ = new TChain("tree");

    // expanded non-uniform Binning
    arrayXBinning_ = 0;
    arrayYBinning_ = 0;
    arrayXBinningVBF_ = 0;
    arrayYBinningVBF_ = 0;
    
    //
    // 2D shape analysis
    //
    if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D || option == HWW_OPT_MT2DMLL_JCP || option == XWW_OPT_MT2DMLL_JCP ) {
        float lep1ptCut(20.), lep2ptCut(10.), mllCut(9999.), dPhiCut(9999.), mtLowCut(0.), mtHighCut(9999.),  mllLooseCut(9999.);
        evaluate_hww_cuts(analysis, lep1ptCut, lep2ptCut, mllCut, dPhiCut, mtLowCut, mtHighCut, mllLooseCut);
        

        // X: mT / Y: mll
        if (analysis < 300) { 

            // for 0/1 jet bin
            nBins2DShapeX_          = 14;
            nBins2DShapeY_          = 9;
            arrayXBinning_ = new double[15]();
            arrayYBinning_ = new double[10]();
            double  mTBinning[15]    = {60,70,80,90,100,110,120,140,160,180,200,220,240,260,280};
            double  mllBinning[10]   = {12,30,45,60,75,100,125,150,175,200};
            setBinningArray(mTBinning, nBins2DShapeX_+1, arrayXBinning_, mllBinning, nBins2DShapeY_+1, arrayYBinning_);

            // for VBF
            nBins2DVBFShapeX_ = 4;
            nBins2DVBFShapeY_ = 4;  
            arrayXBinningVBF_ = new double[5]();
            arrayYBinningVBF_ = new double[5]();
            double   mTBinningVBF[5]    = {30, 92.5, 155, 217.5, 280};
            double   mllBinningVBF[5]   = {0, 50, 100, 150, 200};
            setBinningArray(mTBinningVBF, nBins2DVBFShapeX_+1, arrayXBinningVBF_, mllBinningVBF, nBins2DVBFShapeY_+1, arrayYBinningVBF_);

        } 
        else {
            // for 0/1 jet bin
            nBins2DShapeX_ = 10;
            nBins2DShapeY_ = 8;  
            arrayXBinning_ = new double[11]();
            arrayYBinning_ = new double[9]();
            double   mTBinning[11]   = {80, 110, 140, 170, 200, 230, 260, 290, 320, 350, 380};
            double   mllBinning[9]   = {0, 56.25, 112.5, 168.75, 225, 281.25, 337.5, 393.75, 450};
            setBinningArray(mTBinning, nBins2DShapeX_+1, arrayXBinning_, mllBinning, nBins2DShapeY_+1, arrayYBinning_);
    
            // for VBF
            nBins2DVBFShapeX_ = 2;
            nBins2DVBFShapeY_ = 3; 
            arrayXBinningVBF_ = new double[3]();
            arrayYBinningVBF_ = new double[4]();
            double   mTBinningVBF[3]    = {30, 180, 330};
            double   mllBinningVBF[4]   = {0, 150, 300, 450};
            setBinningArray(mTBinningVBF, nBins2DVBFShapeX_+1, arrayXBinningVBF_, mllBinningVBF, nBins2DVBFShapeY_+1, arrayYBinningVBF_);

        }
    }

/*
    // 
    // 2D shape analysis for the JCP studies
    // 
    if ( option == HWW_OPT_MT2DMLL_JCP || option == XWW_OPT_MT2DMLL_JCP ) {

            // 0/1-jet
            nBins2DShapeX_ = 14;
            nBins2DShapeY_ = 9; 
            arrayXBinning_ = new double[15]();
            arrayYBinning_ = new double[10]();
            double  mTBinning[15]    = {60,70,80,90,100,110,120,140,160,180,200,220,240,260,280};    
            double  mllBinning[10]   = {12,30,45,60,75,100,125,150,175,200}; 
            setBinningArray(mTBinning, nBins2DShapeX_+1, arrayXBinning_, mllBinning, nBins2DShapeY_+1, arrayYBinning_);

            // VBF - note: defining same as normal low mass analysis
            //             as need binning to be defined even if not used
            nBins2DVBFShapeX_ = 4;
            nBins2DVBFShapeY_ = 4;
            arrayXBinningVBF_ = new double[5]();
            arrayYBinningVBF_ = new double[5]();
            double   mTBinningVBF[5]    = {30, 92.5, 155, 217.5, 280};
            double   mllBinningVBF[5]   = {0, 50, 100, 150, 200};
            setBinningArray(mTBinningVBF, nBins2DVBFShapeX_+1, arrayXBinningVBF_, mllBinningVBF, nBins2DVBFShapeY_+1, arrayYBinningVBF_);

    }

    gROOT->cd();
*/
    // 
    // 2D shape results
    //

    if ( option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D || option ==HWW_OPT_MT2DMLL_JCP || option == XWW_OPT_MT2DMLL_JCP ) {
        for (unsigned int j = 0; j < kJetBins; ++j) {
            for (unsigned int i = 0; i < kLeptonTypes; ++i) {
                // 2d shape histogram 
                if(j<2) {
                   h2_shape_[j][i] = new TH2F(Form("%s_histo2_shape_%s_%s", getName().c_str(), jetbin_names[j], types[i]), 
                        "h2_shape", nBins2DShapeX_, arrayXBinning_, nBins2DShapeY_, arrayYBinning_);

                } else {
                   h2_shape_[j][i] = new TH2F(Form("%s_histo2_shape_%s_%s", getName().c_str(), jetbin_names[j], types[i]), 
                        "h2_shape", nBins2DVBFShapeX_, arrayXBinningVBF_, nBins2DVBFShapeY_, arrayYBinningVBF_);
                }
                h2_shape_[j][i]->Sumw2();
            }
        }
    }
 /*   
    // 
    // JCP shape results
    //
    } else if ( ( option ==HWW_OPT_MT2DMLL_JCP ) || ( option == XWW_OPT_MT2DMLL_JCP) ) {

        for (unsigned int j = 0; j < kJetBins; ++j) {
            for (unsigned int i = 0; i < kLeptonTypes; ++i) { 

                if (j<2) { 
                    h2_shape_[j][i] = new TH2F(Form("%s_histo2_shape_%s_%s", getName().c_str(), jetbin_names[j], types[i]), 
                        "h2_shape", nBins2DShapeX_, arrayXBinning_, nBins2DShapeY_, arrayYBinning_);
                } else {
                   h2_shape_[j][i] = new TH2F(Form("%s_histo2_shape_%s_%s", getName().c_str(), jetbin_names[j], types[i]), 
                        "h2_shape", nBins2DVBFShapeX_, arrayXBinningVBF_, nBins2DVBFShapeY_, arrayYBinningVBF_); 
                }
                h2_shape_[j][i]->Sumw2();

            }

        }
    
    }
*/
    //
    // cut results
    //

    for (unsigned int i = 0; i < kLeptonTypes; ++i) {

        // record yields
        h1_results_cutsel_[i]    = new TH1F(Form("%s_h1_results_%s", getName().c_str(), types[i]), "h1_results", 3, -0.5, 2.5);

        // now keep track of the weights directly
        // this is needed for gamma uncertainty
        h1_sumw_cutsel_[i]    = new TH1F(Form("%s_h1_sumw_%s", getName().c_str(), types[i]), "h1_sumw", 3, -0.5, 2.5);

        h1_results_cutsel_[i]->Sumw2();
        h1_sumw_cutsel_[i]->Sumw2();

    }

}

void SmurfSample::add(std::string fileName) 
{
    chain_->Add(fileName.c_str());
}

TChain *SmurfSample::getChain()
{
    return chain_;
}

DataType SmurfSample::getDataType() 
{
    return dataType_;
}

std::string SmurfSample::getName()
{
    return name_;
}

Color_t SmurfSample::getColour() 
{
    return colour_;
}

void SmurfSample::fillWeighted(TH1F *hist, TH1F *sum, unsigned int bin, double weight, double weight_err)
{
    int binidx = hist->FindBin(bin);
    double thiscontent = hist->GetBinContent(binidx);
    double thiserr = hist->GetBinError(binidx);
    double newerr = sqrt(thiserr*thiserr + weight_err*weight_err); 
    hist->SetBinContent(binidx, thiscontent + weight);
    hist->SetBinError(binidx, newerr);
    sum->Fill(bin, weight);
}

void SmurfSample::getWeighted(TH1F *hist, TH1F *sum, unsigned int binMin, unsigned int binMax, double &yield, double &err)
{
    yield = 0.0;
    double sum_err = 0.0;
    double wt_err = 0.0;
    yield = hist->IntegralAndError(binMin, binMax, wt_err);
    yield = sum->IntegralAndError(binMin, binMax, sum_err); 
    err = sqrt (wt_err*wt_err + sum_err*sum_err);
}

void SmurfSample::fillResults(unsigned int bin, 
        unsigned int type, double weight, double weight_err)
{

    //
    // store the results
    //

    if (dataType_ == GAMMA || dataType_ == WJETSELEDATA || dataType_ == WJETSMUDATA) { 
        fillWeighted(h1_results_cutsel_[type], h1_sumw_cutsel_[type], bin, weight, weight_err);
        fillWeighted(h1_results_cutsel_[6], h1_sumw_cutsel_[6], bin, weight, weight_err);                               // inclusive
        if (type == 0 || type == 3) fillWeighted(h1_results_cutsel_[4], h1_sumw_cutsel_[4], bin, weight, weight_err);   // sf
        if (type == 1 || type == 2) fillWeighted(h1_results_cutsel_[5], h1_sumw_cutsel_[5], bin, weight, weight_err);   // of
    } else {
        h1_results_cutsel_[type]->Fill(bin, weight);
        h1_sumw_cutsel_[type]->Fill(bin);
        h1_results_cutsel_[6]->Fill(bin, weight);
        h1_sumw_cutsel_[6]->Fill(bin);
        if (type == 0 || type == 3) {
            h1_results_cutsel_[4]->Fill(bin, weight);
            h1_sumw_cutsel_[4]->Fill(bin);
        }
        if (type == 1 || type == 2) {
            h1_results_cutsel_[5]->Fill(bin, weight);
            h1_sumw_cutsel_[5]->Fill(bin);
        }
    }

}

double SmurfSample::getNWeights(unsigned int binMin, unsigned int binMax, unsigned int type)
{
    double nWeights = h1_sumw_cutsel_[type]->Integral(binMin, binMax);
    return nWeights;
}

void SmurfSample::getResults(unsigned int binMin, unsigned int binMax, 
        unsigned int type, double &yield, double &err)
{

    //
    // look up the results
    //

    if (dataType_ == GAMMA || dataType_ == WJETSELEDATA || dataType_ == WJETSMUDATA)
        getWeighted(h1_results_cutsel_[type], h1_sumw_cutsel_[type], binMin, binMax, yield, err);
    else
        yield = h1_results_cutsel_[type]->IntegralAndError(binMin, binMax, err);

}

//
// 2D shape analysis
//

void SmurfSample::get2DResults(unsigned int binMin, unsigned int binMax, 
        unsigned int type, double &yield, double &err)
{
    // binMin = binMax = jetbin+1 
    unsigned int jetbin = binMin-1;
    unsigned int nbinsX = h2_shape_[jetbin][type]->GetXaxis()->GetNbins();
    unsigned int nbinsY = h2_shape_[jetbin][type]->GetYaxis()->GetNbins();

    yield = h2_shape_[jetbin][type]->IntegralAndError(1, nbinsX, 1, nbinsY, err); 
}

void SmurfSample::fill2DMVAShape(double x, double y, unsigned int jetbin, unsigned int type, double weight)
{

    // fill the normal shape
    h2_shape_[jetbin][type]->Fill(x, y, weight);
    h2_shape_[jetbin][6]->Fill(x, y, weight);
    // fill the normal shape
    if (type == 0 || type == 3) {
        h2_shape_[jetbin][4]->Fill(x, y, weight);
    }
    if (type == 1 || type == 2) {
        h2_shape_[jetbin][5]->Fill(x, y, weight); 
    }
}

TH2F *SmurfSample::get2DMVAShape(unsigned int jetbin, unsigned int type)
{

    TH2F *temp = 0;

    // same flavor
    if ((type & ((1<<0)|(1<<3))) == ((1<<0)|(1<<3))) {
        temp = (TH2F*)h2_shape_[jetbin][0]->Clone();
        temp->Add(h2_shape_[jetbin][3]);
        temp->SetName(Form("histo2_%s", getName().c_str()));
    }
    // opposite flavor
    else if ((type & ((1<<1)|(1<<2))) == ((1<<1)|(1<<2))) {
        temp = (TH2F*)h2_shape_[jetbin][1]->Clone();
        temp->Add(h2_shape_[jetbin][2]);
        temp->SetName(Form("histo2_%s", getName().c_str()));
    }
    // mm
    else if ((type & (1<<0)) == (1<<0)) {
        temp = (TH2F*)h2_shape_[jetbin][0]->Clone();
        temp->SetName(Form("histo2_%s", getName().c_str()));
    }
    // me
    else if ((type & (1<<1)) == (1<<1)) {
        temp = (TH2F*)h2_shape_[jetbin][1]->Clone();
        temp->SetName(Form("histo2_%s", getName().c_str()));
    }
    // em 
    else if ((type & (1<<2)) == (1<<2)) {
        temp = (TH2F*)h2_shape_[jetbin][2]->Clone();
        temp->SetName(Form("histo2_%s", getName().c_str()));
    }
    // ee 
    else if ((type & (1<<3)) == (1<<3)) {
        temp = (TH2F*)h2_shape_[jetbin][3]->Clone();
        temp->SetName(Form("histo2_%s", getName().c_str()));
    }
    // undefined
    else {
        std::cout << "[SmurfSample::get2DMVAShape] Undefined flavor combination" << std::endl;
        return temp;
    }
/*
    // need to check any negative bins and if negative set that bin to be 0
    for (unsigned int x = 1; x <= temp->GetXaxis()->GetNbins(); ++x) {
        for (unsigned int y = 1; y <= temp->GetYaxis()->GetNbins(); ++y) {
            if( temp->GetBinContent(x, y) < 0) {
                temp->SetBinContent(x, y, 0.0);
                //temp->SetBinError(x, y, 0.0);
            }
        }
    }
*/
    
    return temp;
}

//
// utility functions
//

double SmurfSample::fakeRate(double pt, double eta, TH2D *&fhDFRMu, TH2D *&fhDFREl, int fm, int fe){
    // fm == apply muon fake rate, fe == apply electron fake rate
    if(fm == 0 && fe == 0) return 1.0;
    double mypt   = TMath::Min(pt,34.999);
    double myeta  = TMath::Min(fabs(eta),2.4999);
    double prob = 1.0;
    if     (fm == 1){
        Int_t ptbin = fhDFRMu->GetXaxis()->FindBin(mypt);
        Int_t etabin = fhDFRMu->GetYaxis()->FindBin(myeta);  
        prob = fhDFRMu->GetBinContent(ptbin,etabin);
    }
    else if(fe == 1){
        Int_t ptbin = fhDFREl->GetXaxis()->FindBin(mypt);
        Int_t etabin = fhDFREl->GetYaxis()->FindBin(myeta);  
        prob = fhDFREl->GetBinContent(ptbin,etabin);
    }
    return prob/(1-prob);
}

double SmurfSample::fakeRateError(double pt, double eta, TH2D *&fhDFRMu, TH2D *&fhDFREl, int fm, int fe){
    // fm == apply muon fake rate, fe == apply electron fake rate
    if(fm == 0 && fe == 0) return 0.0;
    double mypt   = TMath::Min(pt,34.999);
    double myeta  = TMath::Min(fabs(eta),2.4999);
    double error = 1.0;
    double prob = 1.0;
    if     (fm == 1){
        Int_t ptbin = fhDFRMu->GetXaxis()->FindBin(mypt);
        Int_t etabin = fhDFRMu->GetYaxis()->FindBin(myeta);  
        error = fhDFRMu->GetBinError(ptbin,etabin);
        prob = fhDFRMu->GetBinContent(ptbin,etabin);
    }
    else if(fe == 1){
        Int_t ptbin = fhDFREl->GetXaxis()->FindBin(mypt);
        Int_t etabin = fhDFREl->GetYaxis()->FindBin(myeta);  
        error = fhDFREl->GetBinError(ptbin,etabin);
        prob = fhDFREl->GetBinContent(ptbin,etabin);
    }

    // fractional uncert : sigma_X / X   
    return TMath::Abs(1/(1-prob)/prob*error); 
}

double SmurfSample::LepEffError(double pt, double eta, TH2D *&fhDEffMu, TH2D *&fhDEffEl, int fm, int fe){
    // fm == apply muon fake rate, fe == apply electron fake rate
    if(fm == 0 && fe == 0) return 0.0;
    double mypt   = TMath::Min(pt,99.999);
    double myeta  = fm ? TMath::Min(fabs(eta),2.3999): TMath::Min(fabs(eta),2.4999);
    double error = 1.0;
    if     (fm == 1){
        Int_t ptbin = fhDEffMu->GetXaxis()->FindBin(mypt);
        Int_t etabin = fhDEffMu->GetYaxis()->FindBin(myeta);  
        error = fhDEffMu->GetBinContent(ptbin,etabin);
    }
    else if(fe == 1){
        Int_t ptbin = fhDEffEl->GetXaxis()->FindBin(mypt);
        Int_t etabin = fhDEffEl->GetYaxis()->FindBin(myeta);  
        error = fhDEffEl->GetBinContent(ptbin,etabin);
    }

    return error; 
}

/*
void SmurfSample::addShapeVariation1D(ShapeSyst syst, std::string name, bool normaliseToCentral)
{

    // flag this shape as enabled for this sample
    alternateShapesVector_.insert(syst);

    // create the histograms
    gROOT->cd();

    // qcd scale
    if (((1ll<<syst) & (1ll<<QCDSCALEVAR)) == (1ll<<QCDSCALEVAR)) {
        normalise_qcdscale_ = normaliseToCentral;
        initAlternateShape1D(h1_shape_qcdscale_up_, h1_shape_qcdscale_down_, 
                name.c_str(), name.c_str(), nBinsShape_, minRangeShape_, maxRangeShape_);

        // lepton efficiency
    } else if (((1ll<<syst) & (1ll<<LEPEFFVAR)) == (1ll<<LEPEFFVAR)) {
        normalise_lepeff_ = normaliseToCentral;
        initAlternateShape1D(h1_shape_lepeff_up_, h1_shape_lepeff_down_, 
                name.c_str(), name.c_str(), nBinsShape_, minRangeShape_, maxRangeShape_);

        // jet resulution
    } else if (((1ll<<syst) & (1ll<<JETRESVAR)) == (1ll<<JETRESVAR)) {
        normalise_jetres_ = normaliseToCentral;
        initAlternateShape1D(h1_shape_jetres_up_, h1_shape_jetres_down_, 
                name.c_str(), name.c_str(), nBinsShape_, minRangeShape_, maxRangeShape_);

        // wjets variation for electron
    } else if (((1ll<<syst) & (1ll<<WJETSELESHAPEVAR)) == (1ll<<WJETSELESHAPEVAR)) {
        normalise_wjets_ = normaliseToCentral;
        initAlternateShape1D(h1_shape_wjets_up_, h1_shape_wjets_down_, 
                name.c_str(), name.c_str(), nBinsShape_, minRangeShape_, maxRangeShape_);

        // wjets variation for muon
    } else if (((1ll<<syst) & (1ll<<WJETSMUSHAPEVAR)) == (1ll<<WJETSMUSHAPEVAR)) {
        normalise_wjets_ = normaliseToCentral;
        initAlternateShape1D(h1_shape_wjets_up_, h1_shape_wjets_down_, 
                name.c_str(), name.c_str(), nBinsShape_, minRangeShape_, maxRangeShape_);

        // lepton resolution
    } else if (((1ll<<syst) & (1ll<<LEPRESVAR)) == (1ll<<LEPRESVAR)) {
        normalise_lepres_ = normaliseToCentral;
        initAlternateShape1D(h1_shape_lepres_up_, h1_shape_lepres_down_, 
                name.c_str(), name.c_str(), nBinsShape_, minRangeShape_, maxRangeShape_);

        // met resolution
    } else if (((1ll<<syst) & (1ll<<METVAR)) == (1ll<<METVAR)) {
        normalise_met_ = normaliseToCentral;
        initAlternateShape1D(h1_shape_met_up_, h1_shape_met_down_, 
                name.c_str(), name.c_str(), nBinsShape_, minRangeShape_, maxRangeShape_);

        // invalid option
    } else {
        std::cout << "[SmurfSample::addShapeVariation] Shape variation " << syst << " is not implemented" << std::endl;
    }

}
*/

void SmurfSample::addShapeVariation2D(ShapeSyst syst, std::string name, bool normaliseToCentral)
{

    // flag this shape as enabled for this sample
    alternateShapesVector_.insert(syst);

    // create the histograms
    gROOT->cd();
    
    // qcd scale
    if (((1ll<<syst) & (1ll<<QCDSCALEVAR)) == (1ll<<QCDSCALEVAR)) { 
        normalise_qcdscale_ = normaliseToCentral;
        initAlternateShape2D(h2_shape_qcdscale_up_, h2_shape_qcdscale_down_, name.c_str(), name.c_str(),
            nBins2DShapeX_, arrayXBinning_, nBins2DShapeY_, arrayYBinning_);

    // lepton efficiency
    } else if (((1ll<<syst) & (1ll<<LEPEFFVAR)) == (1ll<<LEPEFFVAR)) {
        normalise_lepeff_ = normaliseToCentral;
        initAlternateShape2D(h2_shape_lepeff_up_, h2_shape_lepeff_down_, name.c_str(), name.c_str(),
            nBins2DShapeX_, arrayXBinning_, nBins2DShapeY_, arrayYBinning_);

    // jet resulution
    } else if (((1ll<<syst) & (1ll<<JETRESVAR)) == (1ll<<JETRESVAR)) {
        normalise_jetres_ = normaliseToCentral;
        initAlternateShape2D(h2_shape_jetres_up_, h2_shape_jetres_down_, name.c_str(), name.c_str(),
            nBins2DShapeX_, arrayXBinning_, nBins2DShapeY_, arrayYBinning_);

    // wjets variation for electron
    } else if (((1ll<<syst) & (1ll<<WJETSELESHAPEVAR)) == (1ll<<WJETSELESHAPEVAR)) {
        normalise_wjets_ = normaliseToCentral;
        initAlternateShape2D(h2_shape_wjets_up_, h2_shape_wjets_down_, name.c_str(), name.c_str(),
            nBins2DShapeX_, arrayXBinning_, nBins2DShapeY_, arrayYBinning_);

    // wjets variation for muon
    } else if (((1ll<<syst) & (1ll<<WJETSMUSHAPEVAR)) == (1ll<<WJETSMUSHAPEVAR)) {
        normalise_wjets_ = normaliseToCentral;
        initAlternateShape2D(h2_shape_wjets_up_, h2_shape_wjets_down_, name.c_str(), name.c_str(), 
            nBins2DShapeX_, arrayXBinning_, nBins2DShapeY_, arrayYBinning_);
  
   // lepton resolution
   } else if (((1ll<<syst) & (1ll<<LEPRESVAR)) == (1ll<<LEPRESVAR)) {
        normalise_lepres_ = normaliseToCentral;
        initAlternateShape2D(h2_shape_lepres_up_, h2_shape_lepres_down_, name.c_str(), name.c_str(),
            nBins2DShapeX_, arrayXBinning_, nBins2DShapeY_, arrayYBinning_);

   // met resolution
   } else if (((1ll<<syst) & (1ll<<METVAR)) == (1ll<<METVAR)) {
       normalise_met_ = normaliseToCentral;
       initAlternateShape2D(h2_shape_met_up_, h2_shape_met_down_, name.c_str(), name.c_str(),
               nBins2DShapeX_, arrayXBinning_, nBins2DShapeY_, arrayYBinning_);

    // invalid option
    } else {
        std::cout << "[SmurfSample::addShapeVariation2D] Shape variation " << syst << " is not implemented" << std::endl;
    }   
    
}

/*
void SmurfSample::fillShapeVariation1D(ShapeSyst syst, bool up, double val, unsigned int jetbin, unsigned int type, double weight)
{

    TH1F **h1_shape;

    // qcd scale
    if (((1ll<<syst) & (1ll<<QCDSCALEVAR)) == (1ll<<QCDSCALEVAR)) {
        if (up) h1_shape = h1_shape_qcdscale_up_[jetbin];
        else    h1_shape = h1_shape_qcdscale_down_[jetbin];

        // lepton efficiency
    } else if (((1ll<<syst) & (1ll<<LEPEFFVAR)) == (1ll<<LEPEFFVAR)) {
        if (up) h1_shape = h1_shape_lepeff_up_[jetbin];
        else    h1_shape = h1_shape_lepeff_down_[jetbin];

        // jet resolution
    } else if (((1ll<<syst) & (1ll<<JETRESVAR)) == (1ll<<JETRESVAR)) {
        if (up) h1_shape = h1_shape_jetres_up_[jetbin];
        else    h1_shape = h1_shape_jetres_down_[jetbin];

        //wjets variation  for electron
    } else if (((1ll<<syst) & (1ll<<WJETSELESHAPEVAR)) == (1ll<<WJETSELESHAPEVAR)) {
        if (up) h1_shape = h1_shape_wjets_up_[jetbin];
        else    h1_shape = h1_shape_wjets_down_[jetbin];

        //wjets variation  for muon
    } else if (((1ll<<syst) & (1ll<<WJETSMUSHAPEVAR)) == (1ll<<WJETSMUSHAPEVAR)) {
        if (up) h1_shape = h1_shape_wjets_up_[jetbin];
        else    h1_shape = h1_shape_wjets_down_[jetbin];

        // lepton resolution
    } else if (((1ll<<syst) & (1ll<<LEPRESVAR)) == (1ll<<LEPRESVAR)) {
        if (up) h1_shape = h1_shape_lepres_up_[jetbin];
        else    h1_shape = h1_shape_lepres_down_[jetbin];

        // met resolution
    } else if (((1ll<<syst) & (1ll<<METVAR)) == (1ll<<METVAR)) {
        if (up) h1_shape = h1_shape_met_up_[jetbin];
        else    h1_shape = h1_shape_met_down_[jetbin];

        // invalid option
    } else {
        std::cout << "[SmurfSample::fillShapeVariation1D] Shape variation " << syst << " is not implemented" << std::endl;
    }

    // fill this alternate shape
    h1_shape[type]->Fill(val, weight);
    h1_shape[4]->Fill(val, weight);

}

TH1F *SmurfSample::getShapeVariation1D(ShapeSyst syst, bool up, unsigned int jetbin, unsigned int type)
{

    TH1F *temp = 0;
    TH1F **h1_shape;
    bool normalise = false;

    // qcd scale
    if (((1ll<<syst) & (1ll<<QCDSCALEVAR)) == (1ll<<QCDSCALEVAR)) {
        if (up) h1_shape = h1_shape_qcdscale_up_[jetbin];
        else    h1_shape = h1_shape_qcdscale_down_[jetbin];
        normalise = normalise_qcdscale_;

        // lepton efficiency
    } else if (((1ll<<syst) & (1ll<<LEPEFFVAR)) == (1ll<<LEPEFFVAR)) {
        if (up) h1_shape = h1_shape_lepeff_up_[jetbin];
        else    h1_shape = h1_shape_lepeff_down_[jetbin];
        normalise = normalise_lepeff_;

        // jet resolution
    } else if (((1ll<<syst) & (1ll<<JETRESVAR)) == (1ll<<JETRESVAR)) {
        if (up) h1_shape = h1_shape_jetres_up_[jetbin];
        else    h1_shape = h1_shape_jetres_down_[jetbin];
        normalise = normalise_jetres_;

        // wjets variation for electron
    } else if (((1ll<<syst) & (1ll<<WJETSELESHAPEVAR)) == (1ll<<WJETSELESHAPEVAR)) {
        if (up) h1_shape = h1_shape_wjets_up_[jetbin];
        else    h1_shape = h1_shape_wjets_down_[jetbin];
        normalise = normalise_wjets_;

        // wjets variation for muon
    } else if (((1ll<<syst) & (1ll<<WJETSMUSHAPEVAR)) == (1ll<<WJETSMUSHAPEVAR)) {
        if (up) h1_shape = h1_shape_wjets_up_[jetbin];
        else    h1_shape = h1_shape_wjets_down_[jetbin];
        normalise = normalise_wjets_;

        // lepton resolution
    } else if (((1ll<<syst) & (1ll<<LEPRESVAR)) == (1ll<<LEPRESVAR)) {
        if (up) h1_shape = h1_shape_lepres_up_[jetbin];
        else    h1_shape = h1_shape_lepres_down_[jetbin];
        normalise = normalise_lepres_;

        // met resolution
    } else if (((1ll<<syst) & (1ll<<METVAR)) == (1ll<<METVAR)) {
        if (up) h1_shape = h1_shape_met_up_[jetbin];
        else    h1_shape = h1_shape_met_down_[jetbin];
        normalise = normalise_met_;

        // invalid option
    } else {
        std::cout << "[SmurfSample::getShapeVariation1D] Shape variation " << syst << " is not implemented" << std::endl;
    }

    // same flavor
    if ((type & ((1<<0)|(1<<3))) == ((1<<0)|(1<<3))) {
        temp = (TH1F*)h1_shape[0]->Clone();
        temp->Add(h1_shape[3]);
    }
    // opposite flavor
    else if ((type & ((1<<1)|(1<<2))) == ((1<<1)|(1<<2))) {
        temp = (TH1F*)h1_shape[1]->Clone();
        temp->Add(h1_shape[2]);
    }
    // mm 
    else if ((type & (1<<0)) == (1<<0)) {
        temp = (TH1F*)h1_shape[0]->Clone();
    }
    // me
    else if ((type & (1<<1)) == (1<<1)) {
        temp = (TH1F*)h1_shape[1]->Clone();
    }
    // em 
    else if ((type & (1<<2)) == (1<<2)) {
        temp = (TH1F*)h1_shape[2]->Clone();
    }
    // ee 
    else if ((type & (1<<3)) == (1<<3)) {
        temp = (TH1F*)h1_shape[3]->Clone();
    }
    // undefined
    else {
        std::cout << "[SmurfSample::getShapeVariation1D] Undefined flavor combination" << std::endl;
        return temp;
    }

    // if required, normalise to the central yield
    if (normalise) {
        TH1F *central = (TH1F*)getMVAShape(jetbin, type)->Clone("temp_central");
        if (temp->Integral(0, temp->GetNbinsX()+1) != 0.0)
            temp->Scale(central->Integral(0, central->GetNbinsX()+1)/temp->Integral(0, temp->GetNbinsX()+1));
        delete central;
    }

    return temp;

}
*/
void SmurfSample::fillShapeVariation2D(ShapeSyst syst, bool up, double valx, double valy, unsigned int jetbin, unsigned int type, double weight)
{

    TH2F **h2_shape;
    
    // qcd scale
    if (((1ll<<syst) & (1ll<<QCDSCALEVAR)) == (1ll<<QCDSCALEVAR)) {
        if (up) h2_shape = h2_shape_qcdscale_up_[jetbin];
        else    h2_shape = h2_shape_qcdscale_down_[jetbin];
        
        // lepton efficiency
    } else if (((1ll<<syst) & (1ll<<LEPEFFVAR)) == (1ll<<LEPEFFVAR)) {
        if (up) h2_shape = h2_shape_lepeff_up_[jetbin];
        else    h2_shape = h2_shape_lepeff_down_[jetbin];

        // jet resolution
    } else if (((1ll<<syst) & (1ll<<JETRESVAR)) == (1ll<<JETRESVAR)) {
        if (up) h2_shape = h2_shape_jetres_up_[jetbin];
        else    h2_shape = h2_shape_jetres_down_[jetbin];

        //wjets variation  for electron
    } else if (((1ll<<syst) & (1ll<<WJETSELESHAPEVAR)) == (1ll<<WJETSELESHAPEVAR)) {
        if (up) h2_shape = h2_shape_wjets_up_[jetbin];
        else    h2_shape = h2_shape_wjets_down_[jetbin];

        //wjets variation  for electron
    } else if (((1ll<<syst) & (1ll<<WJETSMUSHAPEVAR)) == (1ll<<WJETSMUSHAPEVAR)) {
        if (up) h2_shape = h2_shape_wjets_up_[jetbin];
        else    h2_shape = h2_shape_wjets_down_[jetbin];

        // lepton resolution
    } else if (((1ll<<syst) & (1ll<<LEPRESVAR)) == (1ll<<LEPRESVAR)) {
        if (up) h2_shape = h2_shape_lepres_up_[jetbin];
        else    h2_shape = h2_shape_lepres_down_[jetbin];

        // met resolution
    } else if (((1ll<<syst) & (1ll<<METVAR)) == (1ll<<METVAR)) {
        if (up) h2_shape = h2_shape_met_up_[jetbin];
        else    h2_shape = h2_shape_met_down_[jetbin];
        
        // invalid option
    } else {
        std::cout << "[SmurfSample::fillShapeVariation2D] Shape variation " << syst << " is not implemented" << std::endl;
    }

    // fill this alternate shape
    h2_shape[type]->Fill(valx, valy, weight);
    h2_shape[4]->Fill(valx, valy, weight);

}

TH2F *SmurfSample::getShapeVariation2D(ShapeSyst syst, bool up, unsigned int jetbin, unsigned int type)
{

    TH2F *temp = 0;
    TH2F **h2_shape;
    bool normalise = false;
    
    // qcd scale
    if (((1ll<<syst) & (1ll<<QCDSCALEVAR)) == (1ll<<QCDSCALEVAR)) {
        if (up) h2_shape = h2_shape_qcdscale_up_[jetbin];
        else    h2_shape = h2_shape_qcdscale_down_[jetbin];
        normalise = normalise_qcdscale_;
        
        // lepton efficiency
    } else if (((1ll<<syst) & (1ll<<LEPEFFVAR)) == (1ll<<LEPEFFVAR)) {
        if (up) h2_shape = h2_shape_lepeff_up_[jetbin];
        else    h2_shape = h2_shape_lepeff_down_[jetbin];
        normalise = normalise_lepeff_;

        // jet resolution
    } else if (((1ll<<syst) & (1ll<<JETRESVAR)) == (1ll<<JETRESVAR)) {
        if (up) h2_shape = h2_shape_jetres_up_[jetbin];
        else    h2_shape = h2_shape_jetres_down_[jetbin];
        normalise = normalise_jetres_;

        // wjets variation for electron
    } else if (((1ll<<syst) & (1ll<<WJETSELESHAPEVAR)) == (1ll<<WJETSELESHAPEVAR)) {
        if (up) h2_shape = h2_shape_wjets_up_[jetbin];
        else    h2_shape = h2_shape_wjets_down_[jetbin];
        normalise = normalise_wjets_;

        // wjets variation for muon
    } else if (((1ll<<syst) & (1ll<<WJETSMUSHAPEVAR)) == (1ll<<WJETSMUSHAPEVAR)) {
        if (up) h2_shape = h2_shape_wjets_up_[jetbin];
        else    h2_shape = h2_shape_wjets_down_[jetbin];
        normalise = normalise_wjets_;

        // lepton resolution
    } else if (((1ll<<syst) & (1ll<<LEPRESVAR)) == (1ll<<LEPRESVAR)) {
        if (up) h2_shape = h2_shape_lepres_up_[jetbin];
        else    h2_shape = h2_shape_lepres_down_[jetbin];
        normalise = normalise_lepres_;

        // met resolution
    } else if (((1ll<<syst) & (1ll<<METVAR)) == (1ll<<METVAR)) {
        if (up) h2_shape = h2_shape_met_up_[jetbin];
        else    h2_shape = h2_shape_met_down_[jetbin];
        normalise = normalise_met_;
        
        // invalid option
    } else {
        std::cout << "[SmurfSample::getShapeVariation2D] Shape variation " << syst << " is not implemented" << std::endl;
    }

    // same flavor
    if ((type & ((1<<0)|(1<<3))) == ((1<<0)|(1<<3))) {
        temp = (TH2F*)h2_shape[0]->Clone();
        temp->Add(h2_shape[3]);
    }
    // opposite flavor
    else if ((type & ((1<<1)|(1<<2))) == ((1<<1)|(1<<2))) {
        temp = (TH2F*)h2_shape[1]->Clone();
        temp->Add(h2_shape[2]);
    }
    // mm 
    else if ((type & (1<<0)) == (1<<0)) {
        temp = (TH2F*)h2_shape[0]->Clone();
    }
    // me 
    else if ((type & (1<<1)) == (1<<1)) {
        temp = (TH2F*)h2_shape[1]->Clone();
    }
    // em 
    else if ((type & (1<<2)) == (1<<2)) {
        temp = (TH2F*)h2_shape[2]->Clone();
    }
    // mm
    else if ((type & (1<<3)) == (1<<3)) {
        temp = (TH2F*)h2_shape[3]->Clone();
    }
    // undefined
    else {
        std::cout << "[SmurfSample::getShapeVariation1D] Undefined flavor combination" << std::endl;
        return temp;
    }

    // if required, normalise to the central yield
    if (normalise) {   
        TH2F *central = (TH2F*)get2DMVAShape(jetbin, type)->Clone("temp_central"); 
        if (temp->Integral() != 0.0) temp->Scale( central->Integral() / temp->Integral() ); 
        delete central;
    }

    return temp;
}


std::set<ShapeSyst> SmurfSample::getAvailableShapeSystematics()
{
    return alternateShapesVector_;
}

ShapeVar_t SmurfSample::getAvailableShapeSystematicsMask()
{
    ShapeVar_t alternateShapesMask = 0;
    std::set<ShapeSyst>::const_iterator var;
    for (var = alternateShapesVector_.begin(); var != alternateShapesVector_.end(); ++var) {
        alternateShapesMask |= (1ll<<(*var));
    }
    return alternateShapesMask;
}
/*
void SmurfSample::initAlternateShape1D(TH1F *histUp[kJetBins][kLeptonTypes], 
        TH1F *histDown[kJetBins][kLeptonTypes], 
        const char *name, const char *title, int nbins, float min, float max)
{

    for (unsigned int j = 0; j < kJetBins; ++j) {
        for (unsigned int i = 0; i < kLeptonTypes; ++i) {
            histUp[j][i] = new TH1F(Form("%s_histoUP_%s_%s_%s", getName().c_str(), jetbin_names[j], types[i], name), 
                    name, nbins, min, max);
            histDown[j][i] = new TH1F(Form("%s_histoDOWN_%s_%s_%s", getName().c_str(), jetbin_names[j], types[i], name),  
                    name, nbins, min, max);
            histUp[j][i]->Sumw2();
            histDown[j][i]->Sumw2();
        }
    }

}
*/

void SmurfSample::initAlternateShape2D(TH2F *histUp[kJetBins][kLeptonTypes],
        TH2F *histDown[kJetBins][kLeptonTypes],
        const char *name, const char *title, int nbinsx, double *binsx,
        int nbinsy, double *binsy)
{
    for (unsigned int j = 0; j < kJetBins; ++j) {
        for (unsigned int i = 0; i < kLeptonTypes; ++i) { 
            if(j<2) {
                histUp[j][i] = new TH2F(Form("%s_histoUP_%s_%s_%s", getName().c_str(), jetbin_names[j], types[i], name),  
                    name, nbinsx, binsx, nbinsy, binsy);
                histDown[j][i] = new TH2F(Form("%s_histoDOWN_%s_%s_%s", getName().c_str(), jetbin_names[j], types[i], name),
                    name, nbinsx, binsx, nbinsy, binsy);
            } else {
                histUp[j][i] = new TH2F(Form("%s_histoUP_%s_%s_%s", getName().c_str(), jetbin_names[j], types[i], name),  
                    name, nbinsx, binsx, nbinsy, binsy);
                histDown[j][i] = new TH2F(Form("%s_histoDown_%s_%s_%s", getName().c_str(), jetbin_names[j], types[i], name), 
                    name, nbinsx, binsx, nbinsy, binsy); 
            }
            histUp[j][i]->Sumw2();
            histDown[j][i]->Sumw2();  
        }
    }

}

//
// 2D binning 0/1-jet
//

unsigned int SmurfSample::getNBins2DShapeX()     
{ 
    return nBins2DShapeX_; 
}

unsigned int SmurfSample::getNBins2DShapeY()     
{ 
    return nBins2DShapeY_; 
}

double* SmurfSample::getXBinning()   
{   
    return arrayXBinning_;  
}

double* SmurfSample::getYBinning()   
{   
    return arrayYBinning_;  
}

double SmurfSample::getXMax()     
{
    if (arrayXBinning_ !=0) return arrayXBinning_[nBins2DShapeX_];
    else return 0.0;
}

double SmurfSample::getYMax()     
{
    if (arrayYBinning_ !=0) return arrayYBinning_[nBins2DShapeY_];
    else return 0.0; 
}

//
// 2D binning VBF
//
        
unsigned int SmurfSample::getNBins2DVBFShapeX()     
{ 
    return nBins2DVBFShapeX_; 
}

unsigned int SmurfSample::getNBins2DVBFShapeY()     
{ 
    return nBins2DVBFShapeY_; 
}

double* SmurfSample::getVBFXBinning()   
{   
    return arrayXBinningVBF_;  
}

double* SmurfSample::getVBFYBinning()   
{   
    return arrayYBinningVBF_;  
}
double SmurfSample::getVBFXMax()     
{   
    if (arrayXBinningVBF_ !=0) return arrayXBinningVBF_[nBins2DVBFShapeX_];
    else return 0.0;
}

double SmurfSample::getVBFYMax()    
{   
    if (arrayYBinningVBF_ !=0) return arrayYBinningVBF_[nBins2DVBFShapeY_];
    else return 0.0;
}

//
// binning utilities
//

void SmurfSample::setBinningArray(const double* x, const unsigned int nx, double* arrx,
                                  const double* y, const unsigned int ny, double* arry)
{
    for (unsigned int i = 0; i < nx; ++i)
           arrx[i] = x[i];
    for (unsigned int i = 0; i < ny; ++i) 
           arry[i] = y[i];
}

void SmurfSample::printBinning(double* arr, unsigned int n)
{
    for (unsigned int i = 0; i < n; ++i) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
}

