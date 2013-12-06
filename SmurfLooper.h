#ifndef SMURFLOOPER_H
#define SMURFLOOPER_H

#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"

#include "core/Enums.h"
#include "core/TypeDefs.h"
#include <iostream>

class LeptonScaleLookup;
class SmurfSample;
class SmurfTree;

class SmurfLooper {

    public:
        SmurfLooper(float analysis, Option option, RunEra runEra);
        ~SmurfLooper();

        void setGoodRunList(const char *runlist); 
        void setLumiScale(float eeLumi, float mmLumi);
        void loop(SmurfSample *sample);

    private:

        // analysis parameters
        float analysis_;
        Option option_;
        RunEra runEra_;
        double eeLumi_;
        double mmLumi_;

        bool runlistIsSet_;
        
        LeptonScaleLookup *leptonSF_;

        TH1D *HiggsPtKFactor_;
        TH1D *HiggsPtKFactor_QCDscaleSys1_; // --
        TH1D *HiggsPtKFactor_QCDscaleSys2_; // -c
        TH1D *HiggsPtKFactor_QCDscaleSys3_; // c-
        TH1D *HiggsPtKFactor_QCDscaleSys4_; // c+
        TH1D *HiggsPtKFactor_QCDscaleSys5_; // +c
        TH1D *HiggsPtKFactor_QCDscaleSys6_; // ++
        TH1D *HiggsPtKFactor_QCDscaleSys7_; // -+
        TH1D *HiggsPtKFactor_QCDscaleSys8_; // +-

        TH1D *fhDPUS4_;

        // hww analysis
        double ZScaleFactor_[3];
        double ZScaleFactorError_[3];
        double TopScaleFactor_[3];
        double TopScaleFactorError_[3];
        double WWScaleFactor_[3];
        double WWScaleFactorError_[3];
        TH2D *fhDMuonEffError_;
        TH2D *fhDElectronEffError_;
        TH2D *fhDFRMu_;
        TH2D *fhDFREl_;
        TH2D *fhDFRMu_systvar_;
        TH2D *fhDFREl_systvar_;
        TH1D *fhDRatioPhotonElectron_;
        TH1D *fhDPU_;

        // cuts
        Cuts_t testCuts(SmurfTree *tree, DataType dataType, const unsigned int jetbin, const float & dymva);

        // utilities
        void loadWeightHistograms();
	    int getMEProc(const int mH, Option option);

        // for debug purpose  
        FILE *debugtext_;

        //
        // for systematics studies
        //

        // function to call to fill different types of
        // alternate shapes, if associated with the sample
        void fillAlternateQCD(SmurfSample *sample, SmurfTree *tree,
			      const float &varx, const float &vary, unsigned int type, const float &weight, unsigned int jetbin);
        void fillAlternateMet(SmurfSample *sample, SmurfTree *tree,
                const float &varx, const float &vary, unsigned int type, const float &weight, unsigned int jetbin);
        void fillAlternateLepRes(SmurfSample *sample, SmurfTree *tree,
                const float &var_upx, const float &var_downx, const float &var_upy, const float &var_downy, 
                unsigned int type, const float &weight, unsigned int jetbin);
        void fillAlternateLepEff(SmurfSample *sample, SmurfTree *tree,
				 const float &varx, const float &vary, unsigned int type, const float &weight, unsigned int jetbin);
        void fillAlternateFR(SmurfSample *sample, SmurfTree *tree,
			     const float &varx, const float &vary, unsigned int type, const float &weight, unsigned int jetbin);
        void fillAlternateJES(SmurfSample *sample, SmurfTree *tree,
			      const float &varx, const float &vary, unsigned int type, const float &weight, unsigned int jetbin_up, unsigned int jetbin_down);
        void fillAlternateJESVBF(SmurfSample *sample, SmurfTree *tree,
                const float &var_upx, const float &var_downx, const float &var_upy, const float &var_downy, 
                unsigned int type, const float &weight, unsigned int jetbin);

        // get alternate QCD scale for higgs
        float getScaleQCD(float mr, float mf, SmurfTree *tree);

};

#endif

