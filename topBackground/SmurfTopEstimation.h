#ifndef SMURFTOPESTIMATION_H
#define SMURFTOPESTIMATION_H

#include <vector>
#include "../core/Enums.h"
#include <string>
#include "TH1F.h"

class SmurfSample;
class TFile;
class TH1F;

// main functions

void topest( const float mH, const char* inputFileName, FILE *debugtext,  
	     double TopBkgScaleFactor[3], double TopBkgScaleFactorKappa[3]);
 

// utility functions

void getTopData( const float analysis, const char* inputFileName, FILE *debugtext, int njet, const char *histLabel, 
		 const char *flavor, int bin_min, int bin_max, double & ndata, double & err_ndata, double & nbkg, double & err_nbkg, bool isTWBkg, bool verbose);
void calcEff(const double& num, const double& numerr, const double& denum, const double& denumerr, double& eff, double& err_eff);
void calctopeff0j(const double& f_ttbar, const double& err_f_ttbar, const double& eff_1leg, const double& err_eff_1leg,
		  double & eff, double & err_eff);
void calctopbkg(const double& ncr, const double& err_ncr, const double& eff, const double& err_eff, double &yield, double &err_yield);
void getIntegralAndError(TH1F *hist, const int bin_min, const int bin_max, double& n, double& nerr);
#endif

