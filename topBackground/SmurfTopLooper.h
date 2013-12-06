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

class SmurfTopLooper {
 
 public:
  
  SmurfTopLooper();
  ~SmurfTopLooper(){};
  SmurfTopLooper(float analysis, Option option);
  
  void setGoodRunList(const char *runlist); 
  void setLumiScale(float eeLumi, float mmLumi);
  void loop(SmurfSample *sample);
  
 private:
  
  // analysis parameters
  double eeLumi_;
  double mmLumi_;
  bool runlistIsSet_;
  Option option_;
  float analysis_;

  // fake rate histograms
  TH2D *fhDFRMu_;
  TH2D *fhDFREl_;

  // scale factors
  double ZScaleFactor_[3];
  double ZScaleFactorError_[3];

  // utility histograms
  void loadWeightHistograms();


};

#endif
