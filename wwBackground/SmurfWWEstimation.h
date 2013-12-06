#ifndef SMURFWWESTIMATION_H
#define SMURFWWESTIMATION_H

#include <vector>
#include "../core/Enums.h"
#include <string>
#include "TH1F.h"

class SmurfSample;
class TFile;
class TH1F;

// main functions
void wwest(const float analysis, const char* inputFileName, FILE *debugtext, int mHiggs[19],  
	   double WWBkgScaleFactor[2][19], double WWBkgScaleFactorKappa[2][19]); 

// utility functions
void calcalpha(TH1F* &h1_sr, TH1F *& h1_cr, double &alpha, double &alphaerr); 
void getndataincr(int njet, const char* inputFileName, const char *flavor,  double & nWWCR, double &nWWCRErr, FILE *& debugtext);
void calcwwbkg(const double nCR, const double nCRErr, const double alpha, const double alphaErr, double & nSR, double & nSRErr);

#endif

