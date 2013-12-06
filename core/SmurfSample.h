#ifndef SMURFSAMPLE_H
#define SMURFSAMPLE_H

#include "TChain.h"
#include "TH1F.h"
#include "TH2.h"

#include "Enums.h"
#include "TH1Keys.h"

#include <set>
#include <iostream>

//
// global constants
//

static const unsigned int   kLeptonTypes            = 7;
static const unsigned int   kJetBins                = 3;
static const unsigned int   kNDataTypes             = 50;
//static const unsigned int   kNarrayXBinning          = 14;
//static const unsigned int   kNarrayYBinning          = 9;
static const char*          types[7]                = {"mm","me","em","ee", "sf", "of", "incl"};
static const char*          jetbin_names[3]         = { "0j", "1j", "2j"};

//
// a class to represent the results
// for a particular physics process
//

class SmurfSample {

    public:

        SmurfSample();
        SmurfSample(Option option, DataType dataType, Color_t colour, std::string name, float analysis);
        SmurfSample(Option option, DataType dataType, Color_t colour, std::string name);
        ~SmurfSample();

        // Add a smurf tree to the chain
        void add(std::string fileName);

        // Get parameters of this sample
        TChain *getChain();
        std::string getName();
	 	int getTreeType();	 
        Color_t getColour();
        DataType getDataType();

        // Add and get results for this sample
        void fillResults(unsigned int bin, unsigned int type, double weight, double weight_err = 0.0);
        void getResults(unsigned int binMin, unsigned int binMax, unsigned int type, double &yield, double &err);
        double getNWeights(unsigned int binMin, unsigned int binMax, unsigned int type);

        // Add and get mva shape for this sample
        void fillMVAShape(double val, unsigned int jetbin, unsigned int type, double weight);
        TH1F *getMVAShape(unsigned int jetbin, unsigned int type);
        TH1F *getKeysMVAShape(unsigned int jetbin, unsigned int type);
        void fillShapeVariation1D(MVAShapeSyst syst, bool up, double val, unsigned int jetbin, unsigned int type, double weight);
        TH1F *getShapeVariation1D(MVAShapeSyst syst, bool up, unsigned int jetbin, unsigned int type);
        void addShapeVariation1D(MVAShapeSyst syst, std::string name, bool normaliseToCentral);

        // Add and get 2D shape variations
        void get2DResults(unsigned int binMin, unsigned int binMax, unsigned int type, double &yield, double &err);
        void fill2DMVAShape(double x, double y, unsigned int jetbin, unsigned int type, double weight);
        TH2F *get2DMVAShape(unsigned int jetbin, unsigned int type);
        void fillShapeVariation2D(MVAShapeSyst syst, bool up, double x, double y, unsigned int jetbin, unsigned int type, double weight);
        TH2F *getShapeVariation2D(MVAShapeSyst syst, bool up, unsigned int jetbin, unsigned int type);
        void addShapeVariation2D(MVAShapeSyst syst, std::string name, bool normaliseToCentral);

        // Get available valid shape variations associated with this sample
        std::set<MVAShapeSyst> getAvailableShapeSystematics();
        ShapeVar_t getAvailableShapeSystematicsMask();

        // apply fakerate
        double fakeRate(double pt, double eta, TH2D *& fhDFRMu, TH2D *& fhDFREl, int fm, int fe);
        double fakeRateError(double pt, double eta, TH2D *& fhDFRMu, TH2D *& fhDFREl, int fm, int fe);
        double LepEffError(double pt, double eta, TH2D *&fhDEFfMu, TH2D *&fhDEffEl, int fm, int fe);
	
        //
        // 2D binning 0/1-jet
        //
	
        unsigned int getNBins2DShapeX();     
        unsigned int getNBins2DShapeY();     
		double* getXBinning();	
		double* getYBinning();	
        double getXMax();     
        double getYMax();    

        //
        // 2D binning VBF
        //

        unsigned int getNBins2DVBFShapeX();
        unsigned int getNBins2DVBFShapeY();
        double* getVBFXBinning();
        double* getVBFYBinning();
        double getVBFXMax();
        double getVBFYMax();

    private:

        //
        // Details of the sample
        //

        DataType dataType_;
        std::string name_;
        TChain *chain_;
		int treeType_;
        Color_t colour_;

        //
        // cut based analysis
        //

        TH1F *h1_results_cutsel_[kLeptonTypes];        // yields
        TH1F *h1_sumw_cutsel_[kLeptonTypes];           // sum of weights

        //
        // shape based analysis
        //

        // binning parameters
        unsigned int nBinsShape_;
        float minRangeShape_;
        float maxRangeShape_;

        bool normalise_met_;
        bool normalise_lepres_;
        bool normalise_wjets_;
        bool normalise_qcdscale_;
        bool normalise_jetres_;
        bool normalise_lepeff_;

        TH1F    *h1_shape_[kJetBins][kLeptonTypes];           // central shape
        TH1Keys *h1k_shape_[kJetBins][kLeptonTypes];          // smoothed central shape

        // lepton efficiency alternate shapes
        TH1F    *h1_shape_lepeff_up_[kJetBins][kLeptonTypes];
        TH1F    *h1_shape_lepeff_down_[kJetBins][kLeptonTypes];

        // jet resulution alternate shapes
        TH1F    *h1_shape_jetres_up_[kJetBins][kLeptonTypes];
        TH1F    *h1_shape_jetres_down_[kJetBins][kLeptonTypes];

        // qcd scale alternate shapes
        TH1F    *h1_shape_qcdscale_up_[kJetBins][kLeptonTypes];
        TH1F    *h1_shape_qcdscale_down_[kJetBins][kLeptonTypes];

        // wjets alternate shape
        TH1F    *h1_shape_wjets_up_[kJetBins][kLeptonTypes];
        TH1F    *h1_shape_wjets_down_[kJetBins][kLeptonTypes];

        // lepton resolution
        TH1F *h1_shape_lepres_up_[kJetBins][kLeptonTypes];
        TH1F *h1_shape_lepres_down_[kJetBins][kLeptonTypes];

        // met resolution
        TH1F *h1_shape_met_up_[kJetBins][kLeptonTypes];
        TH1F *h1_shape_met_down_[kJetBins][kLeptonTypes];
        
        //
        // 2D shape analysis
        //

        // binning parameters 
		// 0/1 jet bins
        unsigned int nBins2DShapeX_;
        unsigned int nBins2DShapeY_;
        double* arrayXBinning_;
        double* arrayYBinning_;

		// VBF bins
        unsigned int nBins2DVBFShapeX_;
        unsigned int nBins2DVBFShapeY_;
        double* arrayXBinningVBF_;
        double* arrayYBinningVBF_;


        TH2F    *h2_shape_[kJetBins][kLeptonTypes];           // central shape

        // lepton efficiency alternate shapes
        TH2F    *h2_shape_lepeff_up_[kJetBins][kLeptonTypes];
        TH2F    *h2_shape_lepeff_down_[kJetBins][kLeptonTypes];
       
        // jet resulution alternate shapes
        TH2F    *h2_shape_jetres_up_[kJetBins][kLeptonTypes];
        TH2F    *h2_shape_jetres_down_[kJetBins][kLeptonTypes];
 
        // qcd scale alternate shapes
        TH2F    *h2_shape_qcdscale_up_[kJetBins][kLeptonTypes];
        TH2F    *h2_shape_qcdscale_down_[kJetBins][kLeptonTypes];
        
        // wjets alternate shape
        TH2F    *h2_shape_wjets_up_[kJetBins][kLeptonTypes];
        TH2F    *h2_shape_wjets_down_[kJetBins][kLeptonTypes];
        
        // lepton resolution
        TH2F *h2_shape_lepres_up_[kJetBins][kLeptonTypes];
        TH2F *h2_shape_lepres_down_[kJetBins][kLeptonTypes];

        // met resolution
        TH2F *h2_shape_met_up_[kJetBins][kLeptonTypes];
        TH2F *h2_shape_met_down_[kJetBins][kLeptonTypes];
		
        // which alternate shapes are used
        std::set<MVAShapeSyst> alternateShapesVector_;

        //
        // Utility functions
        //

        void fillWeighted(TH1F *hist, TH1F *sum, unsigned int bin, double weight, double weight_err);
        void getWeighted(TH1F *hist, TH1F *sum, unsigned int binMin, unsigned int binMax, double &yield, double &err);
        void initAlternateShape1D(TH1F *histUp[kJetBins][kLeptonTypes], 
                TH1F *histDown[kJetBins][kLeptonTypes], 
                const char *name, const char *title, int nbins, float min, float max);
        void initAlternateShape2D(TH2F *histUp[kJetBins][kLeptonTypes],
                TH2F *histDown[kJetBins][kLeptonTypes],
                const char *name, const char *title, int nbinsx, double *binsx,
                int nbinsy, double *binsy);
        void setBinningArray(const double* x, const unsigned int nx, double* arrx,
                const double* y, const unsigned int ny, double* arry);
        void printBinning(double* arr, unsigned int n);

};

#endif

