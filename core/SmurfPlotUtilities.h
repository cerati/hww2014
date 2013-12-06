#ifndef SMURFPLOTUTILITIES_H
#define SMURFPLOTUTILITIES_H

#include <vector>
#include "TGraph.h"
#include "TArrow.h"
#include "TH1F.h"

class THStack;
class SmurfSample;
class TCanvas;
class TFile;

typedef TH1F H;

// making histograms

void FillHist(TH1F **hist, const unsigned int type, const float &val, const float &weight);
void FormatHist(TH1F **hist, SmurfSample *sample, const char *title, 
        const char *name, int n, float min, float max);
void FormatHist(TH1F **hist, SmurfSample *sample, const char *name,
        const char *title, int n, float bins[]);


// plotting

TGraph *GetEffRej(TFile *file, std::vector<SmurfSample*> bgSamples,
        std::vector<SmurfSample*> signalSamples,
        const char *name, const char *title, bool increasing);

TCanvas *GetStack(TFile *file, std::vector<SmurfSample*> bgSamples, 
        std::vector<SmurfSample*> signalSamples,
        SmurfSample *dataSample, const char *name, const char *title,
        float lumi, float xmin, float xmax, bool log, int rebin = 1);

TH1F *GetHistogram(TFile *file, std::vector<SmurfSample*> samples, const char *name, int rebin = 1);
TH1F *GetHistogram(TFile *file, SmurfSample* sample, const char *name, int rebin = 1);

// for interactive stuff
// comes from old "histtools.C"

void saveHist(const char* filename, const char* pat="*");
void deleteHistos();
H cumulate (const H &in, bool increasing);
TGraph eff_rej (const H &signal, H &background, bool normalize, bool increasing);

// random utilities
TArrow *GetGraphArrow(TGraph *gr, int bin, Color_t colour);

// test gamma+jet prediction
// to data ratio
TH1F *GetGammaJetExtrap(TFile *file, std::vector<SmurfSample*> bgSamples,
        SmurfSample* gammaSample,
        SmurfSample* dataSample,
        const char *name, const char *title);

// draw two histograms
TCanvas *ComparePlots(TFile *f, const char *hist1, const char *hist2, const char *hist3, const char *label);
TCanvas *ComparePlots(TFile *f1, TFile *f2, float lumi1, float lumi2, const char *hist, unsigned int rebin);
TCanvas *ComparePlots(TFile *f, const char *hist1, const char *hist2, const char *label1, const char *label2, unsigned int rebin);

// rebin but into the same number of bins...
TH1F* SmurfRebin(const TH1F *old, const unsigned int rebin);

#endif

