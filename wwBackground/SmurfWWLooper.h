#ifndef SMURFWWLOOPER_H
#define SMURFWWLOOPER_H

#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"

#include "../core/Enums.h"
#include "../core/TypeDefs.h"
#include <iostream>

class SmurfSample;
class SmurfTree;

class SmurfWWLooper {
 
 public:
  
  SmurfWWLooper(){};
  SmurfWWLooper(float analysis, Option option, RunEra runEra);
  ~SmurfWWLooper(){};
  
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
  
  // fake rate histograms
  TH2D *fhDFRMu_;
  TH2D *fhDFREl_;
  
  // top scale factors
  double TopScaleFactor_[3];
  double TopScaleFactorError_[3];
  double ZScaleFactor_[3];
  double ZScaleFactorError_[3];

  // utility histograms
  void loadWeightHistograms();
  
};

#endif
