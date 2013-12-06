#include "SmurfWWEstimation.h"
#include "../core/SmurfSample.h"
#include <cmath>
#include <set>
#include <iostream>
#include <stdio.h>
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"



void wwest(const float analysis, const char* inputFileName, FILE *debugtext, int mHiggs[19],  
	   double WWBkgScaleFactor[2][19], double WWBkgScaleFactorKappa[2][19]) {

  TFile *inputFile = TFile::Open(inputFileName, "READ");
  assert(inputFile);
  std::cout << "[SmurfWWEstimation::wwest]: opening file " << inputFileName  << std::endl;
  gROOT->cd();
  
  
  for (int njet = 0; njet < 2; njet++) {
    fputs("-------------------------------------------------------\n", debugtext);
    fputs(Form("*** WW Background Estimation for %i-Jet\n", njet), debugtext);
    fputs("-------------------------------------------------------\n", debugtext);   

    // 
    // 1- get SR/CR ratio "alpha"
    // 

    fputs("------------Calculate alpha from WW  MC------------------------\n", debugtext);
    // OF 
    double alpha_of (0.), alphaErr_of(0.);
    TH1F *h1_sr_of_wwmc = (TH1F*) inputFile->Get(Form("WW_mll_sig_%ij_of", njet));
    TH1F *h1_cr_of_wwmc = (TH1F*) inputFile->Get(Form("WW_mll_bkg_%ij_of", njet));
    calcalpha(h1_sr_of_wwmc, h1_cr_of_wwmc, alpha_of, alphaErr_of); 
    fputs(Form("alpha (OF) =  %.3f +/- %.3f\n", alpha_of, alphaErr_of), debugtext);

    // OF 
    double alpha_sf (0.), alphaErr_sf(0.);
    TH1F *h1_sr_sf_wwmc = (TH1F*) inputFile->Get(Form("WW_mll_sig_%ij_sf", njet));
    TH1F *h1_cr_sf_wwmc = (TH1F*) inputFile->Get(Form("WW_mll_bkg_%ij_sf", njet));
    calcalpha(h1_sr_sf_wwmc, h1_cr_sf_wwmc, alpha_sf, alphaErr_sf); 
    fputs(Form("alpha (SF) =  %.3f +/- %.3f\n", alpha_sf, alphaErr_sf), debugtext);

    // SF + OF
    double alpha (0.), alphaErr(0.);
    TH1F *h1_sr_wwmc = (TH1F*) h1_sr_of_wwmc->Clone("h1_sr_wwmc");
    h1_sr_wwmc->Add(h1_sr_sf_wwmc);
    TH1F *h1_cr_wwmc = (TH1F*) h1_cr_of_wwmc->Clone("h1_cr_wwmc");
    h1_cr_wwmc->Add(h1_cr_sf_wwmc);
    calcalpha(h1_sr_wwmc, h1_cr_wwmc, alpha, alphaErr); 
    fputs(Form("alpha (SF+OF) =  %.3f +/- %.3f\n", alpha, alphaErr), debugtext);
    
    // 
    // 2 - count the number of events in the control region
    // 
  
    fputs("------------nEvents in Control region in data------------------\n", debugtext);

    // OF 
    double nWWCR_of (0.), nWWCRErr_of (0.);
    getndataincr(njet, inputFileName, "of",  nWWCR_of, nWWCRErr_of, debugtext);
    fputs(Form("WW contribution in control region (OF) %i-Jet = %.2f +/- %.2f\n", njet, nWWCR_of, nWWCRErr_of), debugtext);
    fputs("--------------------------------------------------------------\n", debugtext);
    // SF 
    double nWWCR_sf (0.), nWWCRErr_sf (0.);
    getndataincr(njet, inputFileName, "sf",  nWWCR_sf, nWWCRErr_sf, debugtext);
    fputs(Form("WW contribution in control region (SF) %i-Jet = %.2f +/- %.2f\n", njet, nWWCR_sf, nWWCRErr_sf), debugtext);
    // Total
    fputs("--------------------------------------------------------------\n", debugtext);
    double nWWCR = nWWCR_of + nWWCR_sf;
    double nWWCRErr = sqrt(pow(nWWCRErr_of,2) + pow(nWWCRErr_sf,2));
    fputs(Form("WW contribution in control region (SF+OF) %i-Jet = %.2f +/- %.2f\n", njet, nWWCR, nWWCRErr), debugtext);
    
        
    // 
    // 3 - calculate WW background in signal region
    // 
    
    fputs("------------WW background in the signal region in data--------------\n", debugtext);
    double nWWSR_of(0.), nWWSRErr_of(0.);
    calcwwbkg(nWWCR_of, nWWCRErr_of, alpha_of, alphaErr_of, nWWSR_of, nWWSRErr_of);
    fputs(Form("WW background in signal region in Data (OF) %i-Jet = %.1f +/- %.1f\n", njet, nWWSR_of, nWWSRErr_of), debugtext);
    double nWWSR_sf(0.), nWWSRErr_sf(0.);
    calcwwbkg(nWWCR_sf, nWWCRErr_sf, alpha_sf, alphaErr_sf, nWWSR_sf, nWWSRErr_sf);
    fputs(Form("WW background in signal region in Data (SF) %i-Jet = %.1f +/- %.1f\n", njet, nWWSR_sf, nWWSRErr_sf), debugtext);
    double nWWSR(0.), nWWSRErr(0.);
    calcwwbkg(nWWCR, nWWCRErr, alpha, alphaErr, nWWSR, nWWSRErr);
    fputs(Form("WW background in signal region in Data (SF+OF) %i-Jet = %.1f +/- %.1f\n", njet, nWWSR, nWWSRErr), debugtext);
    fputs("--------------------------------------------------------------\n", debugtext);

    // 
    // 4 - data/MC scalefactors
    // 
    
    double nWWSRMC_of(0.), nWWSRMCErr_of(0.);
    nWWSRMC_of = ((TH1F*)inputFile->Get(Form("WW_mll_sig_%ij_of", njet)))->IntegralAndError(0, 1000, nWWSRMCErr_of);
    double SF_of = nWWSR_of / nWWSRMC_of;
    double SFErr_of = SF_of * nWWSRErr_of / nWWSRMC_of;
    

    double nWWSRMC_sf(0.), nWWSRMCErr_sf(0.);
    nWWSRMC_sf = ((TH1F*)inputFile->Get(Form("WW_mll_sig_%ij_sf", njet)))->IntegralAndError(0, 1000, nWWSRMCErr_sf);
    double SF_sf = nWWSR_sf / nWWSRMC_sf;
    double SFErr_sf = SF_sf * nWWSRErr_sf / nWWSRMC_sf;
    
    
    double nWWSRMC = nWWSRMC_sf + nWWSRMC_of;
    double nWWSRMCErr = sqrt(pow(nWWSRMCErr_sf,2) + pow(nWWSRMCErr_of,2));
    
    double SF = nWWSR / nWWSRMC;
    double SFErr = SF * nWWSRErr / nWWSRMC;
    
    fputs(Form("WW background in signal region in MC (OF) %i-Jet = %.1f +/- %.1f\n", njet, nWWSRMC_of, nWWSRMCErr_of), debugtext);
    fputs(Form("WW background in signal region in MC (SF) %i-Jet = %.1f +/- %.1f\n", njet, nWWSRMC_sf, nWWSRMCErr_sf), debugtext);
    fputs(Form("WW background in signal region in MC (SF+OF) %i-Jet = %.1f +/- %.1f\n", njet, nWWSRMC, nWWSRMCErr), debugtext);
    
    fputs("--------------------------------------------------------------\n", debugtext);
    
    fputs(Form("data/MC scale factor (OF) %i-Jet = %.2f +/- %.2f\n", njet, SF_of, SFErr_of), debugtext);
    fputs(Form("data/MC scale factor (SF) %i-Jet = %.2f +/- %.2f\n", njet, SF_sf, SFErr_sf), debugtext);
    fputs(Form("data/MC scale factor (SF+OF) %i-Jet = %.2f +/- %.2f\n", njet, SF, SFErr), debugtext);


    // 5 - log in the scale factors
    for (int i = 0; i < 19; i++) {
      if ( analysis == mHiggs[i]) {
	WWBkgScaleFactor[njet][i] = SF;
	WWBkgScaleFactorKappa[njet][i] = 1 + SFErr / SF; 
      }
    }
  }
  fputs("--------------------------------------------------------------\n", debugtext);
  
  inputFile->Close();
}

void calcalpha(TH1F* &h1_sr, TH1F *& h1_cr, double &alpha, double &alphaErr) 
{
  double Nsr(0.), Ncr(0.), NsrErr(0.), NcrErr(0.);
  int bin_start(0), bin_end(1000);
  Nsr = h1_sr->IntegralAndError(bin_start, bin_end, NsrErr);
  Ncr = h1_cr->IntegralAndError(bin_start, bin_end, NcrErr);
  
  // std::cout << Form("Nsr = %.2f +/- %.2f\t Ncr = %.2f +/- %.2f \n", Nsr, NsrErr, Ncr, NcrErr);
  alpha =  Ncr > 0.0 ? Nsr / Ncr : 0.0;
  alphaErr = Ncr > 0.0 ? alpha * sqrt(pow(NsrErr/Nsr,2) + pow(NcrErr/Ncr,2)) : 0.0;
  // add the systematic error
  // alphaErr = sqrt(alphaErr*alphaErr+0.08*0.08);
}


void getndataincr(int njet, const char* inputFileName, const char *flavor,  double & nWWCR, double &nWWCRErr, FILE *& debugtext) 
{
  TFile *inputFile = TFile::Open(inputFileName, "READ");
  assert(inputFile);
  // std::cout << "[SmurfWWEstimation::getndataincr]: file " << inputFileName  << std::endl;
  gROOT->cd();  
  
  TH1F *h1_cr_data = (TH1F*) inputFile->Get(Form("Data_mll_bkg_%ij_%s", njet, flavor));
  double nDataCR = h1_cr_data->Integral(0,1000);    
  double nDataCRErr = sqrt(nDataCR);
  fputs(Form("nEvents in control region in data (%s) = %.0f\n", flavor, nDataCR), debugtext);
  
  TH1F *h1_cr_wjets = (TH1F*) inputFile->Get(Form("Wjets_mll_bkg_%ij_%s", njet, flavor));
  double nWjetsCR (0.), nWjetsCRErr(0.);
  if ( h1_cr_wjets != 0) 
    nWjetsCR = h1_cr_wjets->IntegralAndError(0,1000, nWjetsCRErr);
  fputs(Form("Wjets contribution in control region (%s) = %.2f +/- %.2f\n", flavor, nWjetsCR, nWjetsCRErr), debugtext);
  
  TH1F *h1_cr_top = (TH1F*) inputFile->Get(Form("Top_mll_bkg_%ij_%s", njet, flavor));
  double nTopCR (0.), nTopCRErr(0.);
  if ( h1_cr_top != 0 ) 
    nTopCR = h1_cr_top->IntegralAndError(0,1000, nTopCRErr);
  fputs(Form("Top contribution in control region (%s) = %.2f +/- %.2f\n", flavor, nTopCR, nTopCRErr), debugtext);

  TH1F *h1_cr_vv = (TH1F*) inputFile->Get(Form("VV_mll_bkg_%ij_%s", njet, flavor));
  double nVVCR (0.), nVVCRErr(0.);
  if ( h1_cr_vv != 0 ) 
    nVVCR = h1_cr_vv->IntegralAndError(0,1000, nVVCRErr);
  fputs(Form("VV contribution in control region (%s) = %.2f +/- %.2f\n", flavor, nVVCR, nVVCRErr), debugtext);

  TH1F *h1_cr_dyll = (TH1F*) inputFile->Get(Form("DYLL_mll_bkg_%ij_%s", njet, flavor));
  double nDYLLCR (0.), nDYLLCRErr(0.);
  if ( h1_cr_dyll !=  0 ) 
    nDYLLCR = h1_cr_dyll->IntegralAndError(0,1000, nDYLLCRErr);
  fputs(Form("DYLL contribution in control region (%s) = %.2f +/- %.2f\n", flavor, nDYLLCR, nDYLLCRErr), debugtext);
 
  TH1F *h1_cr_ztt = (TH1F*) inputFile->Get(Form("ZTT_mll_bkg_%ij_%s", njet, flavor));
  double nZTTCR (0.), nZTTCRErr(0.);
  if ( h1_cr_ztt != 0 ) 
    nZTTCR = h1_cr_ztt->IntegralAndError(0,1000, nZTTCRErr);
  fputs(Form("ZTT contribution in control region (%s) = %.2f +/- %.2f\n", flavor, nZTTCR, nZTTCRErr), debugtext);


  TH1F *h1_cr_Wgamma = (TH1F*) inputFile->Get(Form("Wgamma_mll_bkg_%ij_%s", njet, flavor));
  double nWgammaCR (0.), nWgammaCRErr(0.);
  if ( h1_cr_Wgamma != 0 ) 
    nWgammaCR = h1_cr_Wgamma->IntegralAndError(0,1000, nWgammaCRErr);
  
  fputs(Form("WGAMMA contribution in control region (%s) = %.2f +/- %.2f\n", flavor, nWgammaCR, nWgammaCRErr), debugtext);
 
  // subtract non-ww backgrounds
  nWWCR = nDataCR - nWjetsCR - nTopCR - nVVCR - nDYLLCR - nZTTCR - nWgammaCR;
  nWWCRErr = sqrt( pow(nDataCRErr,2) + pow(nWjetsCRErr, 2) + pow(nTopCRErr, 2) + pow(nVVCRErr, 2) + pow(nDYLLCRErr,2) + pow(nZTTCRErr,2) + pow(nWgammaCRErr, 2));
  
  inputFile->Close();
}


void calcwwbkg(const double nCR, const double nCRErr, const double alpha, const double alphaErr, double & nSR, double & nSRErr )
{
  nSR = nCR > 0. ? nCR * alpha : 0;
  nSRErr = nSR > 0. ? nSR * sqrt(pow(nCRErr/nCR, 2) + pow(alphaErr/alpha, 2)) : 0.;
}
