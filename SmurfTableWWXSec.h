#ifndef SMURFTABLEWWXSEC_H
#define SMURFTABLEWWXSEC_H

#include <vector>

#include "core/SmurfSample.h"
#include "core/Enums.h"

#include "TCanvas.h"
#include "TFile.h"

#include <fstream>

struct CrossSectionResults {

    public:
        CrossSectionResults() {
            xs_=0.0;
            stat_err_pb_=0.0;
            syst_err_pt_=0.0;
            lumi_err_pb_=0.0;
            total_err_pb_=0.0;
            eff_ww_=0.0;
            eff_ww_err_=0.0;
        }

        float xs_;
        float stat_err_pb_;
        float syst_err_pt_;
        float lumi_err_pb_;
        float total_err_pb_;
        float eff_ww_;
        float eff_ww_err_;
};

// get scale factors
float GetJetVetoEffScaleFactor();
float GetLuminosity();

// get systematics
float GetBackgroundEstimationSystematic(Option option, DataType dataType, unsigned int jetbin, std::string flavor);

// fill yield and uncertainty arrays
void SetYieldsAndUncertainties(Option option, double yield[][kNDataTypes], double nwght[][kNDataTypes],
        double stat2[][kNDataTypes], double syst2[][kNDataTypes], unsigned int jetbin, std::vector<SmurfSample*> samples);

// print tables etc.
void Tabulate(Option option, std::vector<SmurfSample*> samples, std::string filename);

// plots for PAS
TCanvas *makeWWXSecStack(Option option, float analysis, SmurfSample* data, SmurfSample* qqww, SmurfSample* ggww,  SmurfSample *wz, SmurfSample *zz,
        SmurfSample* top, SmurfSample* wjetsEle, SmurfSample *wjetsMu,  SmurfSample* wgamma, SmurfSample* dy, SmurfSample *zjets, std::vector<SmurfSample*> signals,
        TFile *file, const unsigned int flav, const unsigned int njet, const char *name, const char *title, float lumi, bool log, unsigned int rebin, float wwSF);

// tables for PAS

void PrintWWYieldTableVersion1(Option option, unsigned int jetbin,
        std::vector<SmurfSample*> samples, FILE *fout);

void PrintWWYieldTableVersion2(Option option, unsigned int jetbin,
        std::vector<SmurfSample*> samples, FILE *fout);

void PrintSameSignClosureTest(Option option, unsigned int jetbin, unsigned int flavors,
        std::vector<SmurfSample*> samples, FILE *fout);

void PrintWWYieldTable(Option option, unsigned int jetbin, unsigned int flavors,
        std::vector<SmurfSample*> samples, FILE *fout);

void CalculateCrossSection(Option option, std::vector<SmurfSample*> samples, FlavorType flavorType, unsigned int jetbin,
        CrossSectionResults &results);

void DoCrossSection(FlavorType flavorType, double yield[kNDataTypes], double nwght[kNDataTypes],
        double stat2[kNDataTypes], double syst2[kNDataTypes], CrossSectionResults &results);


#endif // SMURFANALYSISWWXSEC_H
