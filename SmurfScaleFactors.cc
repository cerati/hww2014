#include "SmurfScaleFactors.h"
#include "TString.h"
#include "../Smurf/Analysis/HWWlvlv/TopBkgScaleFactors_8TeV.h"
#include "../Smurf/Analysis/HWWlvlv/DYBkgScaleFactors_8TeV.h"
#include "../Smurf/Analysis/HWWlvlv/WWBkgScaleFactors_8TeV.h"
#include "../Smurf/Analysis/HWWlvlv/TopVBFBkgScaleFactors_8TeV.h"

#include <cstdlib>
#include <algorithm>

void getLumiScaleFactor(double &sf, double &err, Option option)
{
  sf = 1.000;
  err = 0.044;
}

void getWjetsScaleFactor(double *sf, double *err, Option option)
{
  sf[0] = 1.00;
  err[0] = 0.36;
  sf[1] = 1.00;
  err[1] = 0.36;
  sf[2] = 1.00;
  err[2] = 0.36;
  
}


void getTopScaleFactor(double *sf, double *err, Option option, float mass)
{ 
  
  for ( int i = 0; i < 2 ; i++) {
    sf[i] = TopBkgScaleFactor(i);
    err[i] = TopBkgScaleFactorKappa(i) -1;
  }
  sf[2] = TopVBFBkgScaleFactor(0); // use index=0 for VBF 
  err[2] = TopVBFBkgScaleFactorKappa(0) -1; 
 
  if ( option == WW_OPT_SMURFXSECSEL) {
    //gc fixme
    //these are obtained from data driven estimation, need to make those not hardcoded
    sf[0] = 1.10959;
    err[0] = 1.13549-1.;
    sf[1] = 1.08424;
    err[1] = 1.02863-1.;
    sf[2] = 1.28938;
    err[2] = 1.02916-1.;
    //     sf[0] = 1.182;
    //     err[0] = 0.205;
    //     sf[1] = 1.086;
    //     err[1] = 0.055;
    //     sf[2] = 1.106;
    //     err[2] = 0.057;
  }

  // Apply the ww level scale factors for the 2-jet bin 
  // for ww preselection or SS closure tests.
  if ( option == HWW_OPT_SMURFPRESEL || option == HWW_OPT_SSCTL ) {
    sf[2] = TopBkgScaleFactor(2);
    err[2] = TopBkgScaleFactorKappa(2) -1;
  }
}


void getZScaleFactor(double *sf, double *err, Option option, double mass, std::string flavor)
{
	if ( mass < 115 && mass > 0 ) mass = 115;
	// these values are for the same flavor
	for ( int i = 0 ; i < 3 ; i++) {

		if (flavor == "sf") {

			sf[i] = DYBkgScaleFactor(0, i);
			err[i] = DYBkgScaleFactorKappa(0, i) -1;

			if (option == HWW_OPT_SMURFCUTSEL) {
				sf[i] = DYBkgScaleFactor(int(mass), i);
				err[i] = DYBkgScaleFactorKappa(int(mass), i) -1;
			} 

			if ( (1ll<<option) & HWW_SHAPE ) {
				sf[i] = DYBkgScaleFactorBDT(int(mass), i);
				err[i] = DYBkgScaleFactorBDTKappa(int(mass), i) -1;
			} 
		}
	}

	// redefine the scale factors if it is of
	if (flavor == "of") {
		sf[0] = 1.000;
		err[0] = 0.100;
		sf[1] = 1.000;
		err[1] = 0.100;
		sf[2] = 1.000;
		err[2] = 0.100;
	}

	if ( option == WW_OPT_SMURFXSECSEL && flavor == "sf") {
	  //gc fixme
	  //these are obtained from data driven estimation, need to make those not hardcoded
	  sf[0] =  5.01;//3.62;
	  err[0] = 41.7/139.1;//17.47/59.46; // only the systematic part
	  sf[1] = 3.78;//3.66;
	  err[1] = 19.2/64.1;//11.38/87.38; // only the systematic part
	  //sf[2] = 2.03;
	  err[2] = 24.16/178.71; // only the systematic part
	}

}

void getWWScaleFactor(double *sf, double *err, Option option, double mass) 
{
  if ( mass < 115) mass = 115; 
  if (option == HWW_OPT_SMURFPRESEL || option == HWW_OPT_SSCTL ) {  
    //sf[0] = 1.0;    sf[1] = 1.0;    sf[2] = 1.0;
    //err[0] = 0.0;   err[1] = 0.0;   err[2] = 0.0; 
    sf[0] = WWBkgScaleFactorMVA(125, 0);
    sf[1] = WWBkgScaleFactorMVA(125, 1);
    sf[2] = WWBkgScaleFactorMVA(125, 1);
    err[0] = WWBkgScaleFactorKappaMVA(125, 0) - 1.;
    err[1] = WWBkgScaleFactorKappaMVA(125, 1) - 1.;
    err[2] = WWBkgScaleFactorKappaMVA(125, 1) - 1.;
    return;
  }

  // fill 0/1-j scale factors from Guillelmo/si code
  for ( int i = 0; i < 2; i++) {  
    if (option == HWW_OPT_SMURFCUTSEL) {
      sf[i] = WWBkgScaleFactorCutBased(std::min(int(mass),int(200)), i);
      err[i] = WWBkgScaleFactorKappaCutBased(std::min(int(mass),int(200)), i) -1; 
    }
    
    if ( (1ll<<option) & HWW_SHAPE ) {
      sf[i] = WWBkgScaleFactorMVA(std::min(int(mass),int(200)), i);
      err[i] = WWBkgScaleFactorKappaMVA(std::min(int(mass),int(200)), i) -1;
    }
/*
    // special treatment of mH > 200 cases
    if ( mass > 200) {
      sf[0] = 1.0;
      err[0] = 0.0;
      sf[1] = 1.0;
      err[1] = 0.0;
    }
*/
  }
  
  sf[2] = sf[1]*2.0;
  err[2] = 0.5; 

}

//
// for scaling leptons and met
//

LorentzVector scaleLepton(const LorentzVector &normal, const int &id, bool up)
{

    // scale
    float scale = 0.0;
    if (abs(id) == 13) {
        scale = up ? 0.01 : -0.01;
    } else if (abs(id) == 11) {
        scale = up ? 0.02 : -0.02;
    } else {
        std::cout << "scaleLepton: Invalid flavor" << std::endl;
        return LorentzVector(0.0, 0.0, 0.0, 0.0);
    }
    
    // compute new px and py and pz
    float newPx = normal.Px() + scale*normal.Px();
    float newPy = normal.Py() + scale*normal.Py();
    float newPz = normal.Pz() + scale*normal.Pz();
    float newMag = sqrt(newPx*newPx + newPy*newPy + newPz*newPz);

    // return scaled vector
    LorentzVector scaled(newPx, newPy, newPz, newMag);
    return scaled;

}

void scaleMet(const LorentzVector &ll, const float &met, const float &metPhi, 
        float &newMet, float &newMetPhi, bool up)
{

    // scale
    float delta = 0.05;
    float scale = up ? (1+delta) : (1-delta);

    // get the hadronic part of MET
    float methx = met*cos(metPhi) + ll.Px();
    float methy = met*sin(metPhi) + ll.Py();

    // scale the met components without the leptons
    // and put the leptons back in to get the new scaled met
    float newMetX = (methx*scale - ll.Px());
    float newMetY = (methy*scale - ll.Py());
    
    // return the modified values
    newMet = sqrt(newMetX*newMetX + newMetY*newMetY);

}

