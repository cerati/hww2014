
//
// Dave "the one but not the only" Evans 
//

#include "MyScanChain.h"
#include "core/Selections.h"
#include "../SmurfScaleFactors.h"

// ROOT includes
#include "TROOT.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <cmath>

// LHAPDF
#include "/tas/dlevans/HWW2012/CMSSW_5_2_3/src/LHAPDF-5.8.92b/include/LHAPDF/LHAPDF.h"

MyScanChain::MyScanChain(float analysis, Option option, std::string pdfName, unsigned int pdfSubset) 
{
    analysis_ = analysis;
    option_ = option;
    LHAPDF::setPDFPath("/tas/dlevans/HWW2012/CMSSW_5_2_3/src/LHAPDF-5.8.92b/PDFSets/");
    LHAPDF::initPDFSetM(genset_, pdfName); 
    LHAPDF::initPDFM(genset_, pdfSubset);
}

//
// Main function
//

int MyScanChain::ScanChain(SmurfSample *sample, std::string pdfName) {

    TObjArray *listOfFiles = sample->getChain()->GetListOfFiles();
    if (listOfFiles->GetEntries() == 0) {
        std::cout << "[SmurfLooper::loop] " << sample->getName() << " is not defined" << std::endl;
        return 1;
    }
    else
        std::cout << "[SmurfLooper::loop] " << sample->getName() << std::endl;

    //
    // setup pdf stuff
    //

    if (pdfName != "cteq6ll")   LHAPDF::initPDFSetM(set_, pdfName + ".LHgrid");
    else                        LHAPDF::initPDFSetM(set_, pdfName + ".LHpdf");

    unsigned int nsets = 1 + LHAPDF::numberPDFM(set_);
    if (pdfName == "NNPDF20_as_0116_100" || pdfName == "NNPDF20_as_0122_100") nsets = 5;
    if (pdfName == "NNPDF20_as_0117_100" || pdfName == "NNPDF20_as_0121_100") nsets = 27;
    if (pdfName == "NNPDF20_as_0118_100" || pdfName == "NNPDF20_as_0120_100") nsets = 72;
    if (pdfName == "cteq6mE") nsets = 1;
    if (pdfName == "cteq6ll") nsets = 1;

    //
    // setup histograms
    //

    std::vector <TH2F*> histArr_all_1bin;
    std::vector <TH2F*> histArr_DF_0j;
    std::vector <TH2F*> histArr_DF_1j;
    std::vector <TH2F*> histArr_DF_VBF;

    double *tmpx = sample->getXBinning();
    double *tmpy = sample->getYBinning();
    double *tmpxvbf = sample->getVBFXBinning();
    double *tmpyvbf = sample->getVBFYBinning();

	float mtbins[15]    =   {60,70,80,90,100,110,120,140,160,180,200,220,240,260,280};
	float mllbins[10]   =   {12,30,45,60,75,100,125,150,175,200};

	for (unsigned int i = 0; i < nsets; ++i) {


        histArr_all_1bin.push_back(new TH2F(Form("%s_all_1bin_%s_%i", sample->getName().c_str(), pdfName.c_str(), i),
                    Form("%s_all_1bin_%s_%i", sample->getName().c_str(), pdfName.c_str(), i),
                    1, 0.0, 100.0, 1, 0.0, 100.0));

        histArr_DF_0j.push_back(new TH2F(Form("%s_DF_0j_%s_%i", sample->getName().c_str(), pdfName.c_str(), i),
                    Form("%s_DF_0j_%s_%i", sample->getName().c_str(), pdfName.c_str(), i), 
                    sample->getNBins2DShapeX(), tmpx, sample->getNBins2DShapeY(), tmpy));

        histArr_DF_1j.push_back(new TH2F(Form("%s_DF_1j_%s_%i", sample->getName().c_str(), pdfName.c_str(), i),
                    Form("%s_DF_1j_%s_%i", sample->getName().c_str(), pdfName.c_str(), i), 
                    sample->getNBins2DShapeX(), tmpx, sample->getNBins2DShapeY(), tmpy));

        histArr_DF_VBF.push_back(new TH2F(Form("%s_DF_VBF_%s_%i", sample->getName().c_str(), pdfName.c_str(), i),
                    Form("%s_DF_VBF_%s_%i", sample->getName().c_str(), pdfName.c_str(), i), 
                    sample->getNBins2DVBFShapeX(), tmpxvbf, sample->getNBins2DVBFShapeY(), tmpyvbf));

    }

    //
    // get other scale factors
    //

    double WWScaleFactor_[3];
    double WWScaleFactorError_[3];
    getWWScaleFactor(WWScaleFactor_, WWScaleFactorError_, option_, analysis_);

    //
    // loop over pdf subsets
    //

    for (unsigned int subset = 0; subset < nsets; ++subset)
    {

        std::cout << "doing set, subset: " << set_ << ", " << subset << std::endl;
        LHAPDF::initPDFM(set_, subset);

        //
        // loop over content of sample
        //

        TIter fileIter(listOfFiles);
        while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {

            SmurfTree *t = new SmurfTree();
            t->LoadTree(currentFile->GetTitle());
            t->InitTree(0);

            // get extra variables
            float dymva = 1.;
            if (t->tree_->GetBranchStatus("dymva"))     t->tree_->SetBranchAddress("dymva", &dymva);

            //
            // loop over events in file
            //

            ULong64_t nEvents = t->tree_->GetEntries();
            //nEvents = 1000;
            for(ULong64_t event = 0; event < nEvents; ++event) {

                t->tree_->GetEntry(event);
                unsigned int njets = t->njets_;
                if (t->njets_ == 3) njets = 2;
                if (t->type_ == 0 || t->type_ == 3) continue;
                if (!cuts(t, sample->getDataType(), njets, dymva)) continue;

                // generated pdf values
                double fx1Q0gen = LHAPDF::xfxM(genset_, t->x1_, t->Q_, int(t->id1_)) / t->x1_;
                double fx2Q0gen = LHAPDF::xfxM(genset_, t->x2_, t->Q_, int(t->id2_)) / t->x2_;

                // subset pdf values
                double fx1Qi = LHAPDF::xfxM(set_, t->x1_, t->Q_, int(t->id1_)) / t->x1_;
                double fx2Qi = LHAPDF::xfxM(set_, t->x2_, t->Q_, int(t->id2_)) / t->x2_;

                // calculate weight and fill histogram
                double pdf_weight = ((fx1Qi*fx2Qi)/(fx1Q0gen*fx2Q0gen));
                double experimental_weight = t->scale1fb_ * t->sfWeightTrig_ * t->sfWeightEff_ * t->sfWeightPU_;
                if ( ((1ll<<sample->getDataType()) & ( (1ll<<QQWW) | (1ll<<GGWW)))) {
                    if (option_ != WW_OPT_SMURFXSECSEL && option_ != HWW_OPT_SMURFPRESEL )
                    experimental_weight = experimental_weight * WWScaleFactor_[njets];
                }

                // inclusive uncertainty, one fixed bin
                histArr_all_1bin[subset]->Fill(50.0, 50.0, pdf_weight * experimental_weight);

                // differential uncertainty, in jet bins...
                if (njets == 0) {
                    histArr_DF_0j[subset]->Fill(min(sample->getXMax()-0.0001, (double)t->mt_), min(sample->getYMax()-0.0001,(double)t->dilep_.M()), 
                            pdf_weight * experimental_weight);
                } 
                else if (njets == 1) {
                    histArr_DF_1j[subset]->Fill(min(sample->getXMax()-0.0001, (double)t->mt_), min(sample->getYMax()-0.0001,(double)t->dilep_.M()),
                            pdf_weight * experimental_weight);
                } 
                else if (njets == 2) {
                    histArr_DF_VBF[subset]->Fill(min(sample->getVBFXMax()-0.0001, (double)t->mt_), min(sample->getVBFYMax()-0.0001, (double)t->dilep_.M()),
                            pdf_weight * experimental_weight);
                }

            } // end loop on events

            delete t;

        } // end loop on files in chain

    } // end loop on subsets

    return 0;

}

bool MyScanChain::cuts(SmurfTree *tree, DataType dataType, const unsigned int jetbin, const float& dymva) 
{

    // 
    // SF special cuts 
    // 

    bool passDY = true;
    bool passMET = true;

    // apply the addtional cuts to the SF events for all types
    // also apply it to the OF type of the ZLLDATA
    // apply DYMVA everywhere for the 0/1 Jets and met for the 2-jet
    if ( tree->type_ == 0 || tree->type_ == 3 || (dataType == ZLLDATA && tree->dstype_ == SmurfTree::data)) {
        if ( tree->njets_ > 1) {
            if ( ! hww_dy_selection(tree) ) passDY = false;
            if ( ! hww_sfmet_selection(tree) )   passMET = false;
        } else {
            if ( ! (hww_sfdymva_selection(tree, dymva)) ) passMET = false;
        }
    }

    if (!passDY)    return false;
    if (!passMET)   return false;
    if (!hww_pass_2DSelection(tree, analysis_, jetbin, option_))                 return false;
    if (!hww_pass_wwSelection(tree, option_))    return false;
    return true;
}

