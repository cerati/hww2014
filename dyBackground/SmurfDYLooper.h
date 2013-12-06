#ifndef SMURFDYLOOPER_H
#define SMURFDYLOOPER_H

#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"

#include "../core/Enums.h"
#include "../core/TypeDefs.h"
#include "../core/SmurfSample.h"
#include <iostream>

class SmurfDYLooper {
 
 public:
  
  SmurfDYLooper(){};
  SmurfDYLooper(float analysis, Option option, RunEra runEra);
  ~SmurfDYLooper(){};
  
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
  
};

#endif
