//
// FIXME: current code is very repeatitive...
// 

#include "SmurfTopEstimation.h"
#include "../core/SmurfSample.h"
#include <cmath>
#include <set>
#include <iostream>
#include <stdio.h>
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"

void topest(const float mH, const char* inputFileName, FILE *debugtext,  
	    double TopBkgScaleFactor[3], double TopBkgScaleFactorKappa[3]) {
  
  TFile *inputFile = TFile::Open(inputFileName, "READ");
  assert(inputFile);
  std::cout << "[SmurfTopEstimation::topest]: opening file " << inputFileName  << std::endl;
  gROOT->cd();

  
  int bin_min = 0;
  int bin_max = 1000;

  double ntop_sig_mc[3];
  double err_ntop_sig_mc[3];
  double toptageff[3];
  double err_toptageff[3];

  double ntop_control_data[3];
  double err_ntop_control_data[3];
  double nbkg_control[3];
  double err_nbkg_control[3];
  
  double ntop_sig_data[3];
  double err_ntop_sig_data[3];

  bool verbose = false;
  
  for ( int njet = 0; njet < 3; njet++) {
    
    fputs("-------------------------------------------------------\n", debugtext);
    fputs(Form("*** Top Background Estimation for %i-Jet\n", njet), debugtext);
    fputs("-------------------------------------------------------\n", debugtext);   

    // 
    // 1- Get the Top events in MC in signal region
    //
    double ntt_sig_mc_sf(0.), err_ntt_sig_mc_sf(0.);
    getIntegralAndError( ((TH1F*)inputFile->Get(Form("TT_mll_sig_%ij_mH%.0f_sf", njet, mH))), bin_min, bin_max, ntt_sig_mc_sf, err_ntt_sig_mc_sf);
        
    double ntt_sig_mc_of(0.), err_ntt_sig_mc_of(0.);
    getIntegralAndError( ((TH1F*)inputFile->Get(Form("TT_mll_sig_%ij_mH%.0f_of", njet, mH))), bin_min, bin_max, ntt_sig_mc_of, err_ntt_sig_mc_of);
    
    double ntw_sig_mc_sf(0.), err_ntw_sig_mc_sf(0.);
    getIntegralAndError( ((TH1F*)inputFile->Get(Form("TW_mll_sig_%ij_mH%.0f_sf", njet, mH))), bin_min, bin_max, ntw_sig_mc_sf, err_ntw_sig_mc_sf);

    double ntw_sig_mc_of(0.), err_ntw_sig_mc_of(0.);
    getIntegralAndError( ((TH1F*)inputFile->Get(Form("TW_mll_sig_%ij_mH%.0f_of", njet, mH))), bin_min, bin_max, ntw_sig_mc_of, err_ntw_sig_mc_of);
    
    double ntop_sig_mc_sf = ntt_sig_mc_sf  + ntw_sig_mc_sf;
    double err_ntop_sig_mc_sf = sqrt(pow(err_ntt_sig_mc_sf,2) + pow(err_ntw_sig_mc_sf, 2));
    
    double ntop_sig_mc_of = ntt_sig_mc_of + ntw_sig_mc_of;
    double err_ntop_sig_mc_of = sqrt(pow(err_ntt_sig_mc_of,2) + pow(err_ntw_sig_mc_of, 2));
    
    
    ntop_sig_mc[njet] = ntop_sig_mc_sf + ntop_sig_mc_of;
    err_ntop_sig_mc[njet] = sqrt(pow(err_ntop_sig_mc_of, 2) + pow(err_ntop_sig_mc_sf, 2));

    fputs(Form("Estimated top events in simulation (SF+OF) %i-Jet %.1f +/- %.1f\n", njet, ntop_sig_mc[njet], err_ntop_sig_mc[njet]), debugtext);
    fputs(Form("Estimated top events in simulation (OF) %i-Jet %.1f +/- %.1f\n", njet, ntop_sig_mc_of, err_ntop_sig_mc_of), debugtext);
    fputs(Form("Estimated top events in simulation (SF) %i-Jet %.1f +/- %.1f\n", njet, ntop_sig_mc_sf, err_ntop_sig_mc_sf), debugtext);
    fputs("-------------------------------------------------------\n", debugtext);


    // 
    // 2 - Get the top-tagged events in data (control region)
    // 
    double ntop_control_data_sf(0.), ntop_control_data_of(0.), err_ntop_control_data_sf(0.), err_ntop_control_data_of(0.);
    double nbkg_control_sf(0.), nbkg_control_of(0.), err_nbkg_control_sf(0.), err_nbkg_control_of(0.);
    // background subtracted numbers
    double ntop_control_of(0.), err_ntop_control_of(0.), ntop_control_sf(0.), err_ntop_control_sf(0.), ntop_control(0.), err_ntop_control(0.);
    
    if ( njet < 2 ) {

      ntop_control_data_sf = ((TH1F*) inputFile->Get(Form("Data_mll_control_%ij_mH%.0f_sf", njet, mH)))->GetEntries();
      ntop_control_data_of = ((TH1F*) inputFile->Get(Form("Data_mll_control_%ij_mH%.0f_of", njet, mH)))->GetEntries();
      ntop_control_data[njet] = ntop_control_data_sf + ntop_control_data_of;
      err_ntop_control_data[njet] = sqrt(ntop_control_data[njet]);

      fputs(Form("top-tagged events in data (SF+OF) %i-Jet %.0f\n", njet, ntop_control_data[njet]), debugtext);
      fputs(Form("top-tagged events in data (OF) %i-Jet %.0f\n", njet, ntop_control_data_of), debugtext);
      fputs(Form("top-tagged events in data (SF) %i-Jet %.0f\n", njet, ntop_control_data_sf), debugtext);

      getTopData(mH, inputFileName, debugtext, njet, "mll_control", "of", bin_min, bin_max, ntop_control_data_of, err_ntop_control_data_of, nbkg_control_of, err_nbkg_control_of, false, verbose);
      ntop_control_of = ntop_control_data_of - nbkg_control_of;
      err_ntop_control_of = sqrt(ntop_control_data_of + pow(err_nbkg_control_of, 2));
      
      getTopData(mH, inputFileName, debugtext, njet, "mll_control", "sf", bin_min, bin_max, ntop_control_data_sf, err_ntop_control_data_sf, nbkg_control_sf, err_nbkg_control_sf, false, verbose);
      ntop_control_sf = ntop_control_data_sf - nbkg_control_sf;
      err_ntop_control_sf = sqrt(ntop_control_data_sf + pow(err_nbkg_control_sf, 2));

      ntop_control = ntop_control_sf + ntop_control_of;
      err_ntop_control = sqrt(pow(err_ntop_control_sf, 2) + pow(err_ntop_control_of,2)); 
      nbkg_control[njet] = nbkg_control_of + nbkg_control_sf;
      err_nbkg_control[njet] = sqrt(pow( err_nbkg_control_of, 2) + pow(err_nbkg_control_sf, 2));

      fputs(Form("top-tagged events in data (bkg subtracted) (SF+OF) %i-Jet %.1f +/- %.1f\n", njet, ntop_control, err_ntop_control), debugtext);
      fputs(Form("top-tagged events in data (bkg subtracted) (SF) %i-Jet %.1f +/- %.1f\n", njet, ntop_control_sf, err_ntop_control_sf), debugtext);
      fputs(Form("top-tagged events in data (bkg subtracted) (OF) %i-Jet %.1f +/- %.1f\n", njet, ntop_control_of, err_ntop_control_of), debugtext);

      
      
    }

    //
    // For two-jet bins do not seperate the sf/of
    // 

    const int etabins =  ((TH1F*) inputFile->Get(Form("Data_eta_cjet_sig_2j_mH%.0f_sf", mH)))->GetNbinsX();
    double ntop_control_data_2j[etabins];
    double err_ntop_control_data_2j[etabins];
    double nbkg_control_2j[etabins];
    double err_nbkg_control_2j[etabins];
    double ntop_control_2j[etabins];
    double err_ntop_control_2j[etabins];
    
    if ( njet == 2 ) {
      fputs("top-tagged events in data in 5 bins of eta: \n", debugtext);
      for (int bin = 0; bin < etabins; bin++) {
	double ntop_control_data_of_2j(0.), err_ntop_control_data_of_2j(0.), ntop_control_data_sf_2j(0.), err_ntop_control_data_sf_2j(0.);
	double nbkg_control_of_2j(0.), err_nbkg_control_of_2j(0.), nbkg_control_sf_2j(0.), err_nbkg_control_sf_2j(0.);
	// background subtracted yields
	double ntop_control_of_2j(0.), err_ntop_control_of_2j(0.), ntop_control_sf_2j(0.), err_ntop_control_sf_2j(0.);
	
	getTopData(mH, inputFileName, debugtext, njet, "eta_cjet_control", "of", bin+1, bin+1, ntop_control_data_of_2j, err_ntop_control_data_of_2j, nbkg_control_of_2j, err_nbkg_control_of_2j, false, verbose);
	ntop_control_of_2j = ntop_control_data_of_2j - nbkg_control_of_2j;
	err_ntop_control_of_2j = sqrt(ntop_control_data_of_2j + pow(err_nbkg_control_of_2j, 2));
	
	getTopData(mH, inputFileName, debugtext, njet, "eta_cjet_control", "sf", bin+1, bin+1, ntop_control_data_sf_2j, err_ntop_control_data_sf_2j, nbkg_control_sf_2j, err_nbkg_control_sf_2j, false, verbose);
	ntop_control_sf_2j = ntop_control_data_sf_2j - nbkg_control_sf_2j;
	err_ntop_control_sf_2j = sqrt(ntop_control_data_sf_2j + pow(err_nbkg_control_sf_2j, 2));
	
	ntop_control_data_2j[bin] = ntop_control_data_sf_2j + ntop_control_data_of_2j;
	err_ntop_control_data_2j[bin] = sqrt(ntop_control_data_sf_2j + ntop_control_data_of_2j);

	nbkg_control_2j[bin] = nbkg_control_sf_2j + nbkg_control_of_2j;
	err_nbkg_control_2j[bin] = sqrt(pow(err_nbkg_control_sf_2j,2) + pow(err_nbkg_control_of_2j, 2));

	ntop_control_2j[bin] = ntop_control_sf_2j + ntop_control_of_2j;
	if ( ntop_control_2j[bin] < 0. ) 	  ntop_control_2j[bin] = 0.;
	err_ntop_control_2j[bin] = sqrt(pow(err_ntop_control_sf_2j, 2) + pow(err_ntop_control_of_2j,2)); 
	
	fputs(Form(" %.1f+/-%.1f,  ", ntop_control_data_2j[bin], err_ntop_control_data_2j[bin]), debugtext);
      }
      fputs("\nbackground events in control region in 5 bins of eta: \n", debugtext);
      for ( int bin = 0; bin < etabins; bin++) 
	fputs(Form(" %.1f+/-%.1f,  ", nbkg_control_2j[bin], err_nbkg_control_2j[bin]), debugtext);
      fputs("\n-------------------------------------------------------\n", debugtext);
    }

    // 
    // 3 - Calculate the top tagging efficiency
    // 
    
    double toptageff_sf(0.), err_toptageff_sf(0.), toptageff_of(0.), err_toptageff_of(0.0);
    
    if ( njet < 2 ) {
      
      // get the numerator and denumerator counts

      // OF 
      double ntop_cali_data_num_of(0.), err_ntop_cali_data_num_of(0.);
      double nbkg_cali_num_of(0.), err_nbkg_cali_num_of(0.);
      double ntop_cali_num_of(0.), err_ntop_cali_num_of(0.); // bkg subtracted
      getTopData(mH, inputFileName, debugtext, njet, "mll_cali_num", "of", bin_min, bin_max, ntop_cali_data_num_of, err_ntop_cali_data_num_of, nbkg_cali_num_of, err_nbkg_cali_num_of, true, verbose);
      ntop_cali_num_of = ntop_cali_data_num_of - nbkg_cali_num_of;
      err_ntop_cali_num_of = sqrt( ntop_cali_data_num_of + pow( err_nbkg_cali_num_of, 2) );

      double ntop_cali_data_denum_of(0.), err_ntop_cali_data_denum_of(0.);
      double nbkg_cali_denum_of(0.), err_nbkg_cali_denum_of(0.);
      double ntop_cali_denum_of(0.), err_ntop_cali_denum_of(0.); // bkg subtracted
      getTopData(mH, inputFileName, debugtext, njet, "mll_cali_denum", "of", bin_min, bin_max, ntop_cali_data_denum_of, err_ntop_cali_data_denum_of, nbkg_cali_denum_of, err_nbkg_cali_denum_of, true, verbose);
      ntop_cali_denum_of = ntop_cali_data_denum_of - nbkg_cali_denum_of;
      err_ntop_cali_denum_of = sqrt( ntop_cali_data_denum_of + pow( err_nbkg_cali_denum_of, 2) );
      
      double eff_1leg_of(0.), err_eff_1leg_of(0.);
      calcEff(ntop_cali_num_of, err_ntop_cali_num_of, ntop_cali_denum_of, err_ntop_cali_denum_of, eff_1leg_of, err_eff_1leg_of);
      fputs(Form("one-leg top-tagging efficiency for OF is %.2f +/- %.2f\n", eff_1leg_of, err_eff_1leg_of), debugtext);

      // SF 
      double ntop_cali_data_num_sf(0.), err_ntop_cali_data_num_sf(0.);
      double nbkg_cali_num_sf(0.), err_nbkg_cali_num_sf(0.);
      double ntop_cali_num_sf(0.), err_ntop_cali_num_sf(0.); // bkg subtracted
      getTopData(mH, inputFileName, debugtext, njet, "mll_cali_num", "sf", bin_min, bin_max, ntop_cali_data_num_sf, err_ntop_cali_data_num_sf, nbkg_cali_num_sf, err_nbkg_cali_num_sf, true, verbose);
      ntop_cali_num_sf = ntop_cali_data_num_sf - nbkg_cali_num_sf;
      err_ntop_cali_num_sf = sqrt( ntop_cali_data_num_sf + pow( err_nbkg_cali_num_sf, 2) );

      double ntop_cali_data_denum_sf(0.), err_ntop_cali_data_denum_sf(0.);
      double nbkg_cali_denum_sf(0.), err_nbkg_cali_denum_sf(0.);
      double ntop_cali_denum_sf(0.), err_ntop_cali_denum_sf(0.); // bkg subtracted
      getTopData(mH, inputFileName, debugtext, njet, "mll_cali_denum", "sf", bin_min, bin_max, ntop_cali_data_denum_sf, err_ntop_cali_data_denum_sf, nbkg_cali_denum_sf, err_nbkg_cali_denum_sf, true, verbose);
      ntop_cali_denum_sf = ntop_cali_data_denum_sf - nbkg_cali_denum_sf;
      err_ntop_cali_denum_sf = sqrt( ntop_cali_data_denum_sf + pow( err_nbkg_cali_denum_sf, 2) );

      // SF + OF 
      double ntop_cali_num = ntop_cali_num_of + ntop_cali_num_sf;
      double err_ntop_cali_num = sqrt(pow(err_ntop_cali_num_of,2) + pow(err_ntop_cali_num_sf,2));
      double ntop_cali_denum = ntop_cali_denum_of + ntop_cali_denum_sf;
      double err_ntop_cali_denum = sqrt(pow(err_ntop_cali_denum_of,2) + pow(err_ntop_cali_denum_sf,2));
      
      // 0-jet top-tag efficiency 
      
      if ( njet == 0 ) {
	// 1) get 1-leg efficiency 
	double eff_1leg_sf(0.), err_eff_1leg_sf(0.);
	calcEff(ntop_cali_num_sf, err_ntop_cali_num_sf, ntop_cali_denum_sf, err_ntop_cali_denum_sf, eff_1leg_sf, err_eff_1leg_sf);
	fputs(Form("one-leg top-tagging efficiency for SF is %.2f +/- %.2f\n", eff_1leg_sf, err_eff_1leg_sf), debugtext);
	
	double eff_1leg(0.), err_eff_1leg(0.);
	calcEff(ntop_cali_num, err_ntop_cali_num, ntop_cali_denum, err_ntop_cali_denum, eff_1leg, err_eff_1leg);
	fputs(Form("one-leg top-tagging efficiency for SF+OF is %.2f +/- %.2f\n", eff_1leg, err_eff_1leg), debugtext);
	
	// 2) Get the fraction of ttbar events in simulation (SF+OF);
	double ntt_control_mc_sf(0.), err_ntt_control_mc_sf(0.);
	getIntegralAndError(((TH1F*)inputFile->Get(Form("TT_mll_control_0j_mH%.0f_sf", mH))), bin_min, bin_max, ntt_control_mc_sf, err_ntt_control_mc_sf);
	
	double ntt_control_mc_of(0.), err_ntt_control_mc_of(0.);
	getIntegralAndError(((TH1F*)inputFile->Get(Form("TT_mll_control_0j_mH%.0f_of", mH))), bin_min, bin_max, ntt_control_mc_of, err_ntt_control_mc_of);
	
	double ntt_control_mc = ntt_control_mc_sf + ntt_control_mc_of;
	double err_ntt_control_mc = sqrt(pow(err_ntt_control_mc_sf,2) + pow(err_ntt_control_mc_of, 2));
	
	double ntw_control_mc_sf(0.), err_ntw_control_mc_sf(0.);
	getIntegralAndError(((TH1F*)inputFile->Get(Form("TW_mll_control_0j_mH%.0f_sf", mH))), bin_min, bin_max, ntw_control_mc_sf, err_ntw_control_mc_sf);
	
	double ntw_control_mc_of(0.), err_ntw_control_mc_of(0.);
	getIntegralAndError(((TH1F*)inputFile->Get(Form("TW_mll_control_0j_mH%.0f_of", mH))), bin_min, bin_max, ntw_control_mc_of, err_ntw_control_mc_of);
	
	double ntw_control_mc = ntw_control_mc_sf + ntw_control_mc_of;
	double err_ntw_control_mc = sqrt(pow(err_ntw_control_mc_sf,2) + pow(err_ntw_control_mc_of, 2));
	
	// hard coded number as the fraction of two b final states for the tw
	double x = 0.153; // FIXME
	double f_ttbar = (ntt_control_mc + x * ntw_control_mc) / (ntt_control_mc + ntw_control_mc);
	double err_f_ttbar = ((1 - x) * 0.17) * f_ttbar;
	fputs(Form("f_ttbar = %.2f +/- %.2f\n", f_ttbar, err_f_ttbar), debugtext);
      
	// 3) Calculate the top-tagging efficiency
	calctopeff0j(f_ttbar, err_f_ttbar, eff_1leg_sf, err_eff_1leg_sf, toptageff_sf, err_toptageff_sf);
	calctopeff0j(f_ttbar, err_f_ttbar, eff_1leg_of, err_eff_1leg_of, toptageff_of, err_toptageff_of);
	calctopeff0j(f_ttbar, err_f_ttbar, eff_1leg, err_eff_1leg, toptageff[njet], err_toptageff[njet]);
      }
      
      // 1-jet top tagging efficiency
      if ( njet == 1) {
	calcEff(ntop_cali_num_of, err_ntop_cali_num_of, ntop_cali_denum_of, err_ntop_cali_denum_of, toptageff_of, err_toptageff_of);
	calcEff(ntop_cali_num_sf, err_ntop_cali_num_sf, ntop_cali_denum_sf, err_ntop_cali_denum_sf, toptageff_sf, err_toptageff_sf);
	calcEff(ntop_cali_num, err_ntop_cali_num, ntop_cali_denum, err_ntop_cali_denum, toptageff[njet], err_toptageff[njet]);
      }
    
      fputs(Form("top-tag efficiency (SF+OF) %i-Jet %.3f +/- %.3f\n", njet, toptageff[njet], err_toptageff[njet]), debugtext);
      fputs(Form("top-tag efficiency (SF) %i-Jet %.3f +/- %.3f\n", njet, toptageff_sf, err_toptageff_sf), debugtext);
      fputs(Form("top-tag efficiency (OF) %i-Jet %.3f +/- %.3f\n", njet, toptageff_of, err_toptageff_of), debugtext);
      fputs("-------------------------------------------------------\n", debugtext);
    }
    
    // special for the 2-jet bin cases where numbers are logged in 5 bins of eta
    double toptageff_2j[etabins], err_toptageff_2j[etabins];
    if ( njet == 2 ) {
      fputs("top-tag efficiency in 5 bins of eta: \n", debugtext);
      for ( int bin = 0; bin < etabins; bin++) {

	// numerator
	double ntop_cali_data_num_of(0.), err_ntop_cali_data_num_of(0.), ntop_cali_data_num_sf(0.), err_ntop_cali_data_num_sf(0.);
	double nbkg_cali_num_of(0.), err_nbkg_cali_num_of(0.), nbkg_cali_num_sf(0.), err_nbkg_cali_num_sf(0.);
	double ntop_cali_num_of(0.), err_ntop_cali_num_of(0.), ntop_cali_num_sf(0.), err_ntop_cali_num_sf(0.); // bkg subtracted
	
	getTopData(mH, inputFileName, debugtext, njet, "eta_cjet_cali_num", "of", bin+1, bin+1, ntop_cali_data_num_of, err_ntop_cali_data_num_of, nbkg_cali_num_of, err_nbkg_cali_num_of, true, verbose); 
	ntop_cali_num_of = ntop_cali_data_num_of - nbkg_cali_num_of;
	err_ntop_cali_num_of = sqrt(ntop_cali_num_of + pow(err_nbkg_cali_num_of, 2) );

	getTopData(mH, inputFileName, debugtext, njet, "eta_cjet_cali_num", "sf", bin+1, bin+1, ntop_cali_data_num_sf, err_ntop_cali_data_num_sf, nbkg_cali_num_sf, err_nbkg_cali_num_sf, true, verbose); 
	ntop_cali_num_sf = ntop_cali_data_num_sf - nbkg_cali_num_sf;
	err_ntop_cali_num_sf = sqrt(ntop_cali_num_sf + pow(err_nbkg_cali_num_sf, 2) );

	double ntop_cali_num = ntop_cali_num_of + ntop_cali_num_sf;
	double err_ntop_cali_num = sqrt(pow(err_ntop_cali_num_of,2) + pow(err_ntop_cali_num_sf,2));
	
	// denumerator

	double ntop_cali_data_denum_of(0.), err_ntop_cali_data_denum_of(0.), ntop_cali_data_denum_sf(0.), err_ntop_cali_data_denum_sf(0.);
	double nbkg_cali_denum_of(0.), err_nbkg_cali_denum_of(0.), nbkg_cali_denum_sf(0.), err_nbkg_cali_denum_sf(0.);
	double ntop_cali_denum_of(0.), err_ntop_cali_denum_of(0.), ntop_cali_denum_sf(0.), err_ntop_cali_denum_sf(0.); // bkg subtracted
	
	getTopData(mH, inputFileName, debugtext, njet, "eta_cjet_cali_denum", "of", bin+1, bin+1, ntop_cali_data_denum_of, err_ntop_cali_data_denum_of, nbkg_cali_denum_of, err_nbkg_cali_denum_of, true, verbose); 
	ntop_cali_denum_of = ntop_cali_data_denum_of - nbkg_cali_denum_of;
	err_ntop_cali_denum_of = sqrt(ntop_cali_denum_of + pow(err_nbkg_cali_denum_of, 2) );

	getTopData(mH, inputFileName, debugtext, njet, "eta_cjet_cali_denum", "sf", bin+1, bin+1, ntop_cali_data_denum_sf, err_ntop_cali_data_denum_sf, nbkg_cali_denum_sf, err_nbkg_cali_denum_sf, true, verbose); 
	ntop_cali_denum_sf = ntop_cali_data_denum_sf - nbkg_cali_denum_sf;
	err_ntop_cali_denum_sf = sqrt(ntop_cali_denum_sf + pow(err_nbkg_cali_denum_sf, 2) );

	double ntop_cali_denum = ntop_cali_denum_of + ntop_cali_denum_sf;
	double err_ntop_cali_denum = sqrt(pow(err_ntop_cali_denum_of,2) + pow(err_ntop_cali_denum_sf,2));
	/*
	std::cout << "eta bin " << bin << ": denumerator in data " <<  ntop_cali_data_denum_of  +  ntop_cali_data_denum_sf  << "\n";
	std::cout << "eta bin " << bin << ": numerator in data " <<  ntop_cali_data_num_of  +  ntop_cali_data_num_sf  << "\n";
	std::cout << "eta bin " << bin << ": denumerator in data (bkg) " <<  nbkg_cali_denum_of  +  nbkg_cali_denum_sf  << "\n";
	std::cout << "eta bin " << bin << ": numerator in data (bkg) " <<  nbkg_cali_num_of  +  nbkg_cali_num_sf  << "\n";
	*/

	// efficiency 
	calcEff(ntop_cali_num, err_ntop_cali_num, ntop_cali_denum, err_ntop_cali_denum, toptageff_2j[bin], err_toptageff_2j[bin]);
	fputs(Form("%.2f+/-%.2f, ", toptageff_2j[bin], err_toptageff_2j[bin]), debugtext);
      }
      fputs("\n-------------------------------------------------------\n", debugtext);
    }
    
    // 
    // 4 - Calculate the data-driven estimate
    // 
    double ntop_sig_data_sf(0.), err_ntop_sig_data_sf(0.), ntop_sig_data_of(0.), err_ntop_sig_data_of(0.);
    if ( njet < 2 ) {
      calctopbkg(ntop_control, err_ntop_control, toptageff[njet], err_toptageff[njet], ntop_sig_data[njet], err_ntop_sig_data[njet]);
      calctopbkg(ntop_control_sf, err_ntop_control_sf, toptageff_sf, err_toptageff_sf, ntop_sig_data_sf, err_ntop_sig_data_sf);
      calctopbkg(ntop_control_of, err_ntop_control_of, toptageff_of, err_toptageff_of, ntop_sig_data_of, err_ntop_sig_data_of);
    }

    if ( njet == 2 ) {
      for ( int bin = 0; bin < etabins; bin++ ) {
	double ntop_sig_data_2j(0.), err_ntop_sig_data_2j(0.); 
	calctopbkg(ntop_control_2j[bin], err_ntop_control_2j[bin], toptageff_2j[bin], err_toptageff_2j[bin], ntop_sig_data_2j, err_ntop_sig_data_2j);
	ntop_sig_data[njet] += ntop_sig_data_2j;
	err_ntop_sig_data[njet] = sqrt(pow(err_ntop_sig_data[njet], 2) + pow(err_ntop_sig_data_2j, 2));
      }
    }
    
    fputs(Form("Estimated top events in data (SF+OF) %i-Jet %.1f +/- %.1f\n", njet, ntop_sig_data[njet], err_ntop_sig_data[njet]), debugtext);
    if ( njet < 2 ) {
      fputs(Form("Estimated top events in data (OF) %i-Jet %.1f +/- %.1f\n", njet, ntop_sig_data_of, err_ntop_sig_data_of), debugtext);
      fputs(Form("Estimated top events in data (SF) %i-Jet %.1f +/- %.1f\n", njet, ntop_sig_data_sf, err_ntop_sig_data_sf), debugtext);
    }
    fputs("-------------------------------------------------------\n", debugtext);
     
    // 
    // 5 - Calculate data/MC scale factors
    // 

    double scalefactor(0.), err_scalefactor(0.);
    scalefactor = ntop_sig_mc[njet] > 0. ? ntop_sig_data[njet]/ntop_sig_mc[njet] : 0.;
    // do not double count the MC uncertainty
    err_scalefactor = ntop_sig_mc[njet]*ntop_sig_data[njet] > 0. ? scalefactor*err_ntop_sig_data[njet]/ntop_sig_data[njet] : 0.; 
    
    double scalefactor_sf(0.), err_scalefactor_sf(0.);
    scalefactor_sf = ntop_sig_mc_sf > 0. ? ntop_sig_data_sf/ntop_sig_mc_sf : 0.;
    err_scalefactor_sf = ntop_sig_mc_sf*ntop_sig_data_sf> 0. ? scalefactor_sf * err_ntop_sig_data_sf/ntop_sig_data_sf : 0.;

    double scalefactor_of(0.), err_scalefactor_of(0.);
    scalefactor_of = ntop_sig_mc_of > 0. ? ntop_sig_data_of/ntop_sig_mc_of : 0.;
    err_scalefactor_of = ntop_sig_mc_of * ntop_sig_data_of> 0. ? scalefactor_of * err_ntop_sig_data_of/ntop_sig_data_of : 0.;
    
    fputs(Form("data/MC scale factor (SF+OF) %i-Jet %.2f +/- %.2f\n", njet, scalefactor, err_scalefactor), debugtext);
    if ( njet < 2 ) {
      fputs(Form("data/MC scale factor (SF) %i-Jet %.2f +/- %.2f\n", njet, scalefactor_sf, err_scalefactor_sf), debugtext);
      fputs(Form("data/MC scale factor (OF) %i-Jet %.2f +/- %.2f\n", njet, scalefactor_of, err_scalefactor_of), debugtext);
    }
  
    TopBkgScaleFactor[njet] = scalefactor;
    TopBkgScaleFactorKappa[njet] = 1 + err_scalefactor/scalefactor;
  }

  // 
  // write out easy output
  // 
  
  bool doLatex = false;
  if (!doLatex) {
    std::cout << "--------------------------------------------------------------------------------------------------\n";
    std::cout << Form("| %40s | %-15s | %-15s | %-15s |","Sample","0-jet","1-jet", "2-jet") << endl;
    std::cout << "--------------------------------------------------------------------------------------------------\n";
  } else {
    std::cout << "\\begin{table}[ht!]\n\\begin{center}\n\\begin{tabular}{l c c c}\n";
    std::cout << "\\hline" << endl;
    std::cout << Form(" %40s & %-15s & %-15s & %-15s \\\\","Sample","0-jet","1-jet", "2-jet") << endl;
    std::cout << "\\hline" << endl;
  }
  
  TString formstr = "| %40s | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f |";
  if (doLatex) formstr = " %40s & %5.1f $\\pm$ %-5.1f & %5.1f $\\pm$ %-5.1f &  %5.1f $\\pm$ %-5.1f \\\\";

  TString formstr_short = "| %40s | %5.1f +/- %-5.1f | %5.1f +/- %-5.1f |        -        |";
  if (doLatex) formstr_short = " %40s & %5.1f $\\pm$ %-5.1f & %5.1f $\\pm$ %-5.1f &  - \\\\";
  
  
  std::cout << Form(formstr,
		    "Estimated top events in simulation", 
		    round(10.*ntop_sig_mc[0])/10., round(10.*err_ntop_sig_mc[0])/10., 
		    round(10.*ntop_sig_mc[1])/10., round(10.*err_ntop_sig_mc[1])/10., 
		    round(10.*ntop_sig_mc[2])/10., round(10.*err_ntop_sig_mc[2])/10.) << "\n" ; 
  
  std::cout << Form(formstr_short, 		    
		    "tagging efficiency (\\%)",
		    round(1000.*toptageff[0])/10.,round(1000.*err_toptageff[0])/10.,
		    round(1000.*toptageff[1])/10.,round(1000.*err_toptageff[1])/10.) << "\n" ;

  std::cout << Form(formstr_short,
		    "top-tagged events in data",
		    round(10.*ntop_control_data[0])/10.,round(10.*err_ntop_control_data[0])/10.,
		    round(10.*ntop_control_data[1])/10.,round(10.*err_ntop_control_data[1])/10.) << "\n";

  std::cout << Form(formstr_short,
		    "background events in control region",
		    round(10.*nbkg_control[0])/10.,round(10.*err_nbkg_control[0])/10.,
		    round(10.*nbkg_control[1])/10.,round(10.*err_nbkg_control[1])/10.) << "\n";

  std::cout << Form(formstr,
		    "Data-driven top background estimate",
		    round(10.*ntop_sig_data[0])/10., round(10.*err_ntop_sig_data[0])/10., 
		    round(10.*ntop_sig_data[1])/10., round(10.*err_ntop_sig_data[1])/10., 
		    round(10.*ntop_sig_data[2])/10., round(10.*err_ntop_sig_data[2])/10.) << "\n" ; 
  
  TString formstr_sf = "| %40s | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f | %5.2f +/- %-5.2f |";
  if (doLatex) formstr_sf = " %40s & %5.2f $\\pm$ %-4.2f & %5.2f $\\pm$ %-4.2f & %5.2f $\\pm$ %-4.2f \\\\";

  
  std::cout << Form(formstr_sf,
		    "Scale factors",
		      round(100.*TopBkgScaleFactor[0])/100., round(100.*(TopBkgScaleFactorKappa[0]-1.)*TopBkgScaleFactor[0])/100.,
		      round(100.*TopBkgScaleFactor[1])/100., round(100.*(TopBkgScaleFactorKappa[1]-1.)*TopBkgScaleFactor[1])/100.,
		      round(100.*TopBkgScaleFactor[2])/100., round(100.*(TopBkgScaleFactorKappa[2]-1.)*TopBkgScaleFactor[2])/100.) << "\n";
  
  if (!doLatex)     std::cout << "--------------------------------------------------------------------------------------------------\n"; 
  if ( doLatex ) {
    std::cout << "\\end{tabular}\n\\caption{Top Estimation Table}\n\\end{center}\n\\end{table}\n";
    std::cout << "\\hline" << std::endl;
  }
  
  
  inputFile->Close();
}


// 
// This function get the top estimation in data, subtracting non-top backgrounds
// 
void getTopData( const float mH, const char* inputFileName, FILE *debugtext, int njet,  const char *histLabel, const char *flavor, int bin_min, int bin_max, 
		 double & ndata, double & err_ndata, double & nbkg, double & err_nbkg, bool isTWBkg, bool verbose) {

  if ( verbose ) 
    std::cout << "analyzing histograms with " << histLabel << " in the " << flavor <<  "\n";
  
  TFile *inputFile = TFile::Open(inputFileName, "READ");
  assert(inputFile);
  gROOT->cd();

  TH1F *hdata  = (TH1F*) inputFile->Get(Form("Data_%s_%ij_mH%.0f_%s", histLabel, njet, mH, flavor));
  TH1F *htw = (TH1F*) inputFile->Get(Form("TW_%s_%ij_mH%.0f_%s", histLabel, njet, mH, flavor));
  TH1F *hqqww = (TH1F*) inputFile->Get(Form("qqWW_%s_%ij_mH%.0f_%s", histLabel, njet, mH, flavor));
  TH1F *hggww = (TH1F*) inputFile->Get(Form("ggWW_%s_%ij_mH%.0f_%s", histLabel, njet, mH, flavor));
  TH1F *hdyll = (TH1F*) inputFile->Get(Form("DYLL_%s_%ij_mH%.0f_%s", histLabel, njet, mH, flavor));
  TH1F *hwz = (TH1F*) inputFile->Get(Form("WZ_%s_%ij_mH%.0f_%s", histLabel, njet, mH, flavor));
  TH1F *hzz = (TH1F*) inputFile->Get(Form("ZZ_%s_%ij_mH%.0f_%s", histLabel, njet, mH, flavor));
  TH1F *hztt = (TH1F*) inputFile->Get(Form("Ztt_%s_%ij_mH%.0f_%s", histLabel, njet, mH, flavor));
  TH1F *hwjets = (TH1F*) inputFile->Get(Form("Wjets_%s_%ij_mH%.0f_%s", histLabel, njet, mH, flavor));
  TH1F *hwgamma = (TH1F*) inputFile->Get(Form("Wgamma_%s_%ij_mH%.0f_%s", histLabel, njet, mH, flavor));

  ndata = hdata->Integral(bin_min, bin_max);
  err_ndata = sqrt(ndata);
  
  double ntw(0.), err_ntw(0.);
  getIntegralAndError(htw, bin_min, bin_max, ntw, err_ntw);
  err_ntw = 0.5 * err_ntw;
  
  double nqqww(0.), err_nqqww(0.);
  getIntegralAndError(hqqww, bin_min, bin_max, nqqww, err_nqqww);
  err_nqqww = 0.5 * nqqww;
  
  double nggww(0.), err_nggww(0.);
  getIntegralAndError(hggww, bin_min, bin_max, nggww, err_nggww);
  err_nggww = 0.5 * nggww;
  
  double ndyll(0.), err_ndyll(0.);
  getIntegralAndError(hdyll, bin_min, bin_max, ndyll, err_ndyll);
  err_ndyll = 0.5 * ndyll;

  double nwz(0.), err_nwz(0.);
  getIntegralAndError(hwz, bin_min, bin_max, nwz, err_nwz);
  err_nwz = 0.5 * nwz;

  double nzz(0.), err_nzz(0.);
  getIntegralAndError(hzz, bin_min, bin_max, nzz, err_nzz);
  err_nzz = 0.5 * nzz;
  
  double nztt(0.), err_nztt(0.);
  getIntegralAndError(hztt, bin_min, bin_max, nztt, err_nztt);
  err_nztt = 0.5 * nztt;
  
  double nwjets(0.), err_nwjets(0.);
  getIntegralAndError(hwjets, bin_min, bin_max, nwjets, err_nwjets);
  err_nwjets = sqrt( pow(err_nwjets, 2) + pow(nwjets*0.35, 2));
  
  double nwgamma(0.), err_nwgamma(0.);
  getIntegralAndError(hwgamma, bin_min, bin_max, nwgamma, err_nwgamma);
  err_nwgamma = nwgamma * 0.5;

  nbkg = nqqww + nggww + ndyll + nwz + nzz + nztt + nwjets + nwgamma;
  // err_nbkg =  sqrt(pow((nqqww + nggww + ndyll + nwz + nzz + nztt + nwgamma)*0.5, 2) + pow(err_nwjets, 2));
  err_nbkg = sqrt( pow(err_nqqww,2) + pow(err_nggww,2) + pow(err_ndyll, 2) + pow(err_nwz, 2) 
		   + pow(err_nzz, 2) + pow(err_nztt, 2) + pow(err_nwjets, 2) + pow(err_nwgamma, 2) );

  if ( isTWBkg ) {
    nbkg = nbkg +  ntw;
    err_nbkg = sqrt(pow(err_nbkg, 2) + pow(err_ntw,2));
  }

  if ( verbose ) {
    std::cout << "bin_min = " << bin_min << "\t bin_max = " << bin_max << "\n";
    std::cout << Form("ndata = %.2f +/- %.2f \n", ndata, err_ndata);
    std::cout << Form("ntw = %.2f +/- %.2f \n", ntw, err_ntw);
    std::cout << Form("nqqww = %.2f +/- %.2f \n", nqqww, err_nqqww);
    std::cout << Form("nggww = %.2f +/- %.2f \n", nggww, err_nggww);
    std::cout << Form("ndyll = %.2f +/- %.2f \n", ndyll, err_ndyll);
    std::cout << Form("nwz = %.2f +/- %.2f \n", nwz, err_nwz);
    std::cout << Form("nzz = %.2f +/- %.2f \n", nzz, err_nzz);
    std::cout << Form("nztt = %.2f +/- %.2f \n", nztt, err_nztt);
    std::cout << Form("nwjets = %.2f +/- %.2f \n", nwjets, err_nwjets);
    std::cout << Form("nwgamma = %.2f +/- %.2f \n", nwgamma, err_nwgamma);
    std::cout << Form("total bkg = %.2f +/- %.2f \n", nbkg, err_nbkg);
  }

  inputFile->Close();
}


void calcEff(const double& num, const double& numerr, const double& denum, const double& denumerr, double& eff, double& err_eff) {
  
  if ( denum <= 0.) {
    eff = 0.;
    err_eff = 0.;
    return;
  }
  eff = num / denum;
  err_eff = sqrt(eff*(1-eff)/denum);
  return;
}


void calctopeff0j(const double& f_ttbar, const double& err_f_ttbar, const double& eff_1leg, const double& err_eff_1leg,
		  double & eff, double & err_eff)
{
  eff = f_ttbar * (1 - pow(1 - eff_1leg, 2)) + (1 - f_ttbar) * eff_1leg;
  // add the uncertainty2 due to eff_1leg_data_bgsub
  err_eff = 0.;
  err_eff += pow((err_eff_1leg) * (f_ttbar * ( 1 - 2 * eff_1leg) + 1), 2);
  // add the uncertainty2 due to err_f_ttbar
  err_eff += pow((err_f_ttbar) * (eff_1leg - eff_1leg * eff_1leg), 2);
  // take square root
  err_eff  = sqrt(err_eff);
}


void calctopbkg(const double& ncr, const double& err_ncr, const double& eff, const double& err_eff, double &yield, double &err_yield)
{
    
  yield = (1 - eff) * ncr / eff;
  err_yield = sqrt( pow((1 - eff)/eff, 2) * pow(err_ncr, 2) + pow(ncr, 2)*pow(err_eff, 2)/pow(eff,4) );
}


void getIntegralAndError(TH1F *hist, const int bin_min, const int bin_max, double& n, double& nerr ) {
  
  if ( hist != 0 ) {
    n = hist->IntegralAndError(bin_min, bin_max, nerr);
  }   else {
    n = 0.;
    nerr = 0.;
  }
}
