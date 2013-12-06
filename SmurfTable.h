#ifndef SMURFTABLE_H
#define SMURFTABLE_H

#include <vector>
#include "core/Enums.h"
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "core/SmurfSample.h"

// results tables
void printResultsTable(std::vector<SmurfSample*> samples, Option option, bool doJetBins = false);

// card for statistical tools
void printCard(std::vector<SmurfSample*> samples, Option option, 
	       unsigned int jetbin, float analysis, std::string cdir, unsigned int fcode, unsigned int mva_option, unsigned int runEra);

// root file containing 2d shape hists
void print2DShapeHistograms(std::vector<SmurfSample*> samples, Option option,
              unsigned int jetbin, float analysis, std::string cdir, unsigned int fcode, unsigned int runEra);

// root file containing shape hists
void printShapeHistograms(std::vector<SmurfSample*> samples, Option option,
			  unsigned int jetbin, float analysis, std::string cdir, unsigned int fcode, unsigned int runEra);
void print2DShapeHistograms(std::vector<SmurfSample*> samples, Option option,
			  unsigned int jetbin, float analysis, std::string cdir, unsigned int fcode, unsigned int runEra);

// make a plot for HWW
TCanvas *makeHWWAnalysisStack(Option option, float analysis, std::vector<SmurfSample *> samples, DataType dyType,
    TFile *file, const unsigned int flav, const unsigned int njet, const char *dir, const char *name, const char *title, float lumi, float dyScale = 1.0);

#endif

