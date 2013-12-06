#ifndef SMURFSCALEFACTORS_H
#define SMURFSCALEFACTORS_H

#include "core/Enums.h"
#include "core/TypeDefs.h"

void getLumiScaleFactor(double &sf, double &err, Option option);
void getWjetsScaleFactor(double *sf, double *err, Option option);
void getZScaleFactor(double *sf, double *err, Option option, double mass, std::string flavor);
void getWWScaleFactor(double *sf, double *err, Option option, double mass);
void getTopScaleFactor(double *sf, double *err, Option option, float mass);

//
// smearing factors 
// for lepton and met syst etc.
//

LorentzVector scaleLepton(const LorentzVector &normal, const int &id, bool up);
void scaleMet(const LorentzVector &ll, const float &met, const float &metPhi,
        float &newMet, float &newMetPhi, bool up);


#endif

