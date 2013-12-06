#ifndef SMURFDYESTIMATION_H
#define SMURFDYESTIMATION_H

#include <vector>
#include "../core/Enums.h"
#include <string>
#include "TH1F.h"

class SmurfSample;
class TFile;
class TH1F;

// main functions
void fillRoutin(const char* inputFileName, const char* routinFileName, FILE *debugtext);
void dyest(const float analysis, Option option, const char* inputFileName, const char* rouinFileName, FILE *debugtext, int mHiggs[20],  
	   double DYBkgScaleFactorWWPreselection[3], double DYBkgScaleFactorWWPreselectionKappa[3], 
	   double DYBkgScaleFactorHiggsSelection[3][20], double DYBkgScaleFactorHiggsSelectionKappa[3][20],
	   double DYBkgScaleFactorHiggsSelectionMVA[3][20], double DYBkgScaleFactorHiggsSelectionKappaMVA[3][20]);

// utility functions
void calcR(double Nout, double NoutE, double Nin, double NinE, double & R, double & RE);
void ofsubtraction (double Nll, double NllE, double Nem, double NemE, double kll, double kllE, double & Nll_subt, double & Nll_subtE);
void ofsubt_single(const double ndata[4], double k_ee, double k_eeE, double & nee_subt, double &nee_subtE, 
		   double &nmm_subt, double &nmm_subtE, double &ntot_subt, double &ntot_subtE, 
		   FILE *debugtext, bool forRatio = false, bool verbose=false);
void combll(double nee, double neeE, double nmm, double nmmE, double & n, double & nE); 
void lookupR(int njet, const char* fileName, const char *suffix , double & R_ee, double & R_eeE, double & R_eeE_syst, 
	     double & R_mm, double & R_mmE, double & R_mmE_syst, double & R, double & RE, double & RE_syst);
void ratio_syst(TH1F* & ratio_vs_met, double & R, double & RE_stat, double & RE_syst);
void getEstimates(double Di_subt, double  Di_subtE, double R, double RE, double RE_syst,
		  double & pred, double & predE, double & predE_syst);

#endif

