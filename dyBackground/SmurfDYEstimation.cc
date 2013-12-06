#include "SmurfDYEstimation.h"
#include "../core/SmurfSample.h"
#include <cmath>
#include <set>
#include <iostream>
#include <stdio.h>
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TMath.h"


void fillRoutin(const char* inputFileName, const char* routinFileName, FILE *debugtext)
{
  
  TFile *inputFile = TFile::Open(inputFileName, "READ");
  assert(inputFile);
  std::cout << "[SmurfDYEstimation::fillRoutin]: opening file " << inputFileName  << std::endl;
  gROOT->cd();


  // 
  // declare the histograms needed to do Routin
  // 
  
  TH1F *Ree_mc[3];
  TH1F *Ninee_mc[3];
  TH1F *Noutee_mc[3];
  TH1F *Rmm_mc[3];
  TH1F *Ninmm_mc[3];
  TH1F *Noutmm_mc[3];
  TH1F *R_mc[3];
  
    
  TH1F *Ree_data[3];
  TH1F *Ninee_data[3];
  TH1F *Noutee_data[3];
  TH1F *Rmm_data[3];
  TH1F *Ninmm_data[3];
  TH1F *Noutmm_data[3];
  TH1F *R_data[3];


  for (int njet = 0; njet < 3; njet++) {
    
    Ninee_mc[njet] = (TH1F*)inputFile->Get(Form("DYLL_met_in_%ij_ee", njet));
    Noutee_mc[njet] = (TH1F*)inputFile->Get(Form("DYLL_met_out_%ij_ee", njet));
    Ree_mc[njet] = (TH1F*)Noutee_mc[njet]->Clone(Form("Ree_vs_met_mc_%iJet", njet));
    Ninmm_mc[njet] = (TH1F*)inputFile->Get(Form("DYLL_met_in_%ij_mm", njet));
    Noutmm_mc[njet] = (TH1F*)inputFile->Get(Form("DYLL_met_out_%ij_mm", njet));
    Rmm_mc[njet] = (TH1F*)Noutmm_mc[njet]->Clone(Form("Rmm_vs_met_mc_%iJet", njet));
    R_mc[njet] = (TH1F*)Noutmm_mc[njet]->Clone(Form("R_vs_met_mc_%iJet", njet));
    
    Ninee_data[njet] = (TH1F*)inputFile->Get(Form("Data_met_in_%ij_ee", njet));
    Noutee_data[njet] = (TH1F*)inputFile->Get(Form("Data_met_out_%ij_ee", njet));
    Ree_data[njet] = (TH1F*)Noutee_data[njet]->Clone(Form("Ree_vs_met_data_%iJet", njet));
    Ninmm_data[njet] = (TH1F*)inputFile->Get(Form("Data_met_in_%ij_mm", njet));
    Noutmm_data[njet] = (TH1F*)inputFile->Get(Form("Data_met_out_%ij_mm", njet));
    Rmm_data[njet] = (TH1F*)Noutmm_data[njet]->Clone(Form("Rmm_vs_met_data_%iJet", njet));
    R_data[njet] = (TH1F*)Noutmm_data[njet]->Clone(Form("R_vs_met_data_%iJet", njet));
    
    if ( Ninee_mc[njet] == 0x0 || Noutee_mc[njet] == 0x0 ||  Ninmm_mc[njet] == 0x0 ||  Ninmm_mc[njet] == 0x0 ) {
      std::cout << Form("Error, histograms DYLL_met_in_%ij_ee(mm) or DYLL_met_out_%ij_ee(mm) are not found! \n", njet, njet);
      return;
    }

    if ( Ninee_data[njet] == 0x0 || Noutee_data[njet] == 0x0 ||  Ninmm_data[njet] == 0x0 ||  Ninmm_data[njet] == 0x0 ) {
      std::cout << Form("Error, histograms Data_met_in_%ij_ee(mm) or Data_met_out_%ij_ee(mm) are not found! \n", njet, njet);
      return;
    }
  }
  
  
  // 
  // fill Routin from MC
  // 
  
  for (int njet = 0; njet < 3; njet++) {
  
    fputs("-------------------------------------------------------\n", debugtext);
    fputs(Form("***Calculating R based on MC for %i-Jet\n", njet), debugtext);
    fputs("-------------------------------------------------------\n", debugtext);
    
    for (int i = 1; i < Ninee_mc[njet]->GetNbinsX()+1; i++) {
      float metcut_low = Ninee_mc[njet]->GetBinLowEdge(i);
      float metcut_high = Ninee_mc[njet]->GetBinLowEdge(i) + Ninee_mc[njet]->GetBinWidth(i);
      int bin_start = i;
      int bin_end = i;
      
      if ( i == Ninee_mc[njet]->GetNbinsX()  ) {
	metcut_high = 99999.;
	bin_end =  Ninee_mc[njet]->GetNbinsX()+1;
      }
      
      fputs(Form("***Calculating R for Met Cut (%.2f-%.2f):\n", metcut_low, metcut_high), debugtext);
            
      
      // Calculate R in MM
      double Nin_mm(0.), Nin_mmE(0.0), Nout_mm(0.0), Nout_mmE(0.0), Rmm(0.), RmmE(0.);
      Nin_mm = Ninmm_mc[njet]->IntegralAndError(bin_start,bin_end, Nin_mmE);
      Nout_mm = Noutmm_mc[njet]->IntegralAndError(bin_start,bin_end, Nout_mmE);
      calcR(Nout_mm, Nout_mmE, Nin_mm, Nin_mmE, Rmm, RmmE);
      
      // Calculate R in EE
      double Nin_ee(0.), Nin_eeE(0.0), Nout_ee(0.0), Nout_eeE(0.0), Ree(0.), ReeE(0.);
      Nin_ee = Ninee_mc[njet]->IntegralAndError(bin_start,bin_end, Nin_eeE);
      Nout_ee = Noutee_mc[njet]->IntegralAndError(bin_start,bin_end, Nout_eeE);
      calcR(Nout_ee, Nout_eeE, Nin_ee, Nin_eeE, Ree, ReeE);
      
      // Calculate R in EE+MM
      double Nin(0.), NinE(0.0), Nout(0.0), NoutE(0.0), R(0.), RE(0.0);
      Nin = Nin_ee + Nin_mm;
      NinE = sqrt(Nin_eeE*Nin_eeE + Nin_mmE*Nin_mmE);
      Nout = Nout_ee + Nout_mm;
      NoutE = sqrt(Nout_eeE*Nout_eeE + Nout_mmE*Nout_mmE);
      calcR(Nout, NoutE, Nin, NinE, R, RE);
       
      // text outputs
      fputs(Form("MM: Nin = %.1f +/- %.1f, Nout = %.1f +/- %.1f, R = %.3f +/- %.3f\n", Nin_mm, Nin_mmE, Nout_mm, Nout_mmE, Rmm, RmmE), debugtext);
      fputs(Form("EE: Nin = %.1f +/- %.1f, Nout = %.1f +/- %.1f, R = %.3f +/- %.3f\n", Nin_ee, Nin_eeE, Nout_ee, Nout_eeE, Ree, ReeE), debugtext);
      fputs(Form("EE+MM: Nin = %.1f +/- %.1f, Nout = %.1f +/- %.1f, R = %.3f +/- %.3f\n", Nin, NinE, Nout, NoutE, R, RE), debugtext);
      
      // Fill ratio
      Rmm_mc[njet]->SetBinContent(i, Rmm);
      Rmm_mc[njet]->SetBinError(i, RmmE);
      Ree_mc[njet]->SetBinContent(i, Ree);
      Ree_mc[njet]->SetBinError(i, ReeE);
      R_mc[njet]->SetBinContent(i, R);
      R_mc[njet]->SetBinError(i, RE);
    }
  }

  // 
  // fill Routin from Data, calculate the signal box
  // 
  for ( int njet = 0; njet < 3; njet++) {    
    
    fputs("-------------------------------------------------------\n", debugtext);
    fputs(Form("***Calculating R based on Data for %i-Jet\n", njet), debugtext);
    fputs("-------------------------------------------------------\n", debugtext);

    // get electron/muon efficiency ratio (the same code as in fillRoutin)
    double Ninee(0.), NineeE(0.), Ninmm(0.), NinmmE(0.);
    Ninee = ((TH1F*) inputFile->Get(Form("Data_met_in_ww_%ij_ee", njet)))->GetEntries();
    NineeE = sqrt(Ninee);
    
    Ninmm = ((TH1F*) inputFile->Get(Form("Data_met_in_ww_%ij_mm", njet)))->GetEntries();
    NinmmE = sqrt(Ninmm);
    
    double k_ee = sqrt(Ninee/Ninmm);
    double k_eeE = sqrt(pow(NineeE/Ninmm,2) + pow(NinmmE*Ninee/(Ninmm*Ninmm),2))/2.;

    std::cout << "k_ee = " << k_ee << "\n"; 
    for (int i = 1; i < Ninee_data[njet]->GetNbinsX()+1; i++) {    
      float metcut_low = Ninee_mc[njet]->GetBinLowEdge(i);
      float metcut_high = Ninee_mc[njet]->GetBinLowEdge(i) + Ninee_mc[njet]->GetBinWidth(i);
      int bin_start = i;
      int bin_end = i;
      
      if ( i == Ninee_mc[njet]->GetNbinsX() ) {
	metcut_high = 99999.;
	bin_end =  Ninee_mc[njet]->GetNbinsX()+1;
      }
      
      fputs(Form("***Calculating R for Met Cut (%.2f-%.2f):\n", metcut_low, metcut_high), debugtext);
      double ndata[4];
      
      // get the in-peak of-subtracted yield
      ndata[0] = ((TH1F*)inputFile->Get(Form("Data_met_in_%ij_mm", njet)))->Integral(bin_start,bin_end); 
      ndata[1] = ((TH1F*)inputFile->Get(Form("Data_met_in_%ij_me", njet)))->Integral(bin_start,bin_end); 
      ndata[2] = ((TH1F*)inputFile->Get(Form("Data_met_in_%ij_em", njet)))->Integral(bin_start,bin_end); 
      ndata[3] = ((TH1F*)inputFile->Get(Form("Data_met_in_%ij_ee", njet)))->Integral(bin_start,bin_end); 
      fputs(Form("In-Z peak yield(mm/me/em/ee): %.0f\t %.0f\t %.0f\t %.0f\t\n", ndata[0], ndata[1], ndata[2], ndata[3]), debugtext);
      
      double Nin_mm(0.), Nin_mmE(0.0), Nin_ee(0.), Nin_eeE(0.0), Nin(0.0), NinE(0.0);  
      //ofsubt_single(ndata, k_ee, k_eeE,  Nin_ee, Nin_eeE, Nin_mm, Nin_mmE, Nin, NinE, debugtext, true, true);
      ofsubt_single(ndata, k_ee, k_eeE,  Nin_ee, Nin_eeE, Nin_mm, Nin_mmE, Nin, NinE, debugtext, false, true);
      
      // get the out-peak of-subtracted yield
      ndata[0] = ((TH1F*)inputFile->Get(Form("Data_met_out_%ij_mm", njet)))->Integral(bin_start,bin_end); 
      ndata[1] = ((TH1F*)inputFile->Get(Form("Data_met_out_%ij_me", njet)))->Integral(bin_start,bin_end); 
      ndata[2] = ((TH1F*)inputFile->Get(Form("Data_met_out_%ij_em", njet)))->Integral(bin_start,bin_end); 
      ndata[3] = ((TH1F*)inputFile->Get(Form("Data_met_out_%ij_ee", njet)))->Integral(bin_start,bin_end); 
      fputs(Form("Out-Z peak yield(mm/me/em/ee): %.0f\t %.0f\t %.0f\t %.0f\t\n", ndata[0], ndata[1], ndata[2], ndata[3]), debugtext);
	
      double Nout_mm(0.), Nout_mmE(0.0), Nout_ee(0.), Nout_eeE(0.0), Nout(0.0), NoutE(0.0);
      //ofsubt_single(ndata, k_ee, k_eeE,  Nout_ee, Nout_eeE, Nout_mm, Nout_mmE, Nout, NoutE, debugtext, true, true);
      ofsubt_single(ndata, k_ee, k_eeE,  Nout_ee, Nout_eeE, Nout_mm, Nout_mmE, Nout, NoutE, debugtext, false, true);

      // subtract the vz contribution as well

      double Nin_ee_vz(0.), Nin_ee_vzE(0.), Nin_mm_vz(0.), Nin_mm_vzE(0.), Nin_vz(0.), Nin_vzE(0.);
      Nin_ee_vz = ((TH1F*)inputFile->Get(Form("VV_met_in_%ij_ee", njet)))->IntegralAndError(bin_start,bin_end, Nin_ee_vzE);
      Nin_mm_vz = ((TH1F*)inputFile->Get(Form("VV_met_in_%ij_mm", njet)))->IntegralAndError(bin_start,bin_end, Nin_mm_vzE);
      combll(Nin_mm_vz, Nin_mm_vzE, Nin_ee_vz, Nin_ee_vzE, Nin_vz, Nin_vzE);

      double Nout_ee_vz(0.), Nout_ee_vzE(0.), Nout_mm_vz(0.), Nout_mm_vzE(0.), Nout_vz(0.), Nout_vzE(0.);
      Nout_ee_vz = ((TH1F*)inputFile->Get(Form("VV_met_out_%ij_ee", njet)))->IntegralAndError(bin_start,bin_end, Nout_ee_vzE);
      Nout_mm_vz = ((TH1F*)inputFile->Get(Form("VV_met_out_%ij_mm", njet)))->IntegralAndError(bin_start,bin_end, Nout_mm_vzE);
      combll(Nout_mm_vz, Nout_mm_vzE, Nout_ee_vz, Nout_ee_vzE, Nout_vz, Nout_vzE);
      
	  Nout_mm = Nout_mm - Nout_mm_vz;
      Nout_mmE = sqrt( pow(Nout_mmE, 2) + pow(Nout_mm_vzE, 2) + pow(0.1 * Nout_mm_vz, 2) );

      Nout_ee = Nout_ee - Nout_ee_vz;
      Nout_eeE = sqrt( pow(Nout_eeE, 2) + pow(Nout_ee_vzE, 2) + pow(0.1 * Nout_ee_vz, 2) );
      
      Nout = Nout - Nout_vz;
      NoutE = sqrt( pow(NoutE, 2) + pow(Nout_vzE, 2) + pow(0.1 * Nout_vz, 2) );
     
	  Nin_mm = Nin_mm - Nin_mm_vz;
      Nin_mmE = sqrt( pow(Nin_mmE, 2) + pow(Nin_mm_vzE, 2) + pow(0.1 * Nin_mm_vz, 2) );

      Nin_ee = Nin_ee - Nin_ee_vz;
      Nin_eeE = sqrt( pow(Nin_eeE, 2) + pow(Nin_ee_vzE, 2) + pow(0.1 * Nin_ee_vz, 2) );
      
      Nin = Nin - Nin_vz;
      NinE = sqrt( pow(NinE, 2) + pow(Nin_vzE, 2) + pow(0.1 * Nin_vz, 2) );

      // calculated Routin
      double Rmm(0.), RmmE(0.), Ree(0.), ReeE(0.), R(0.), RE(0.);
      calcR(Nout_mm, Nout_mmE, Nin_mm, Nin_mmE, Rmm, RmmE);
      calcR(Nout_ee, Nout_eeE, Nin_ee, Nin_eeE, Ree, ReeE);
      calcR(Nout, NoutE, Nin, NinE, R, RE);

      // text outputs    
      fputs(Form("MM: Nin = %.1f +/- %.1f, Nout = %.1f +/- %.1f, R = %.3f +/- %.3f\n", Nin_mm, Nin_mmE, Nout_mm, Nout_mmE, Rmm, RmmE), debugtext); 
      fputs(Form("EE: Nin = %.1f +/- %.1f, Nout = %.1f +/- %.1f, R = %.3f +/- %.3f\n", Nin_ee, Nin_eeE, Nout_ee, Nout_eeE, Ree, ReeE), debugtext);
      fputs(Form("EE+MM: Nin = %.1f +/- %.1f, Nout = %.1f +/- %.1f, R = %.3f +/- %.3f\n", Nin, NinE, Nout, NoutE, R, RE), debugtext);

      // Fill ratio
      Rmm_data[njet]->SetBinContent(i, Rmm);
      Rmm_data[njet]->SetBinError(i, RmmE);
      Ree_data[njet]->SetBinContent(i, Ree);
      Ree_data[njet]->SetBinError(i, ReeE);
      R_data[njet]->SetBinContent(i, R);
      R_data[njet]->SetBinError(i, RE);
    }
    int nbins = Ninee_data[njet]->GetNbinsX(); 
    Rmm_data[njet]->SetBinContent(nbins, 0);
    Rmm_data[njet]->SetBinError(nbins, 0);
    Ree_data[njet]->SetBinContent(nbins,0);
    Ree_data[njet]->SetBinError(nbins,0);
    R_data[njet]->SetBinContent(nbins, 0);
    R_data[njet]->SetBinError(nbins, 0);


  }
  
  inputFile->Close();  
  
  std::cout << "creating " << routinFileName << "\n";
  TFile *f = new TFile(routinFileName, "RECREATE");
  f->cd();
  for (int njet = 0; njet < 3; njet++) {
    Ree_mc[njet]->Write();
    Rmm_mc[njet]->Write();
    R_mc[njet]->Write();
    Ree_data[njet]->Write();
    Rmm_data[njet]->Write();
    R_data[njet]->Write();
  }
  f->Close();
  
}

void dyest(const float analysis, Option option, const char* inputFileName, const char* routinFileName, FILE *debugtext, int mHiggs[20],  
	   double DYBkgScaleFactorWWPreselection[3], double DYBkgScaleFactorWWPreselectionKappa[3], 
	   double DYBkgScaleFactorHiggsSelection[3][20], double DYBkgScaleFactorHiggsSelectionKappa[3][20],
	   double DYBkgScaleFactorHiggsSelectionMVA[3][20], double DYBkgScaleFactorHiggsSelectionKappaMVA[3][20])
{  
  TFile *inputFile = TFile::Open(inputFileName, "READ");
  assert(inputFile);
  std::cout << "[SmurfDYEstimation::dyest]: opening file " << inputFileName  << std::endl;
  
  for (int njet = 0; njet < 3; njet++) {
    fputs("-------------------------------------------------------\n", debugtext);
    fputs(Form("***Final Drell-Yan Estimation for %i-Jet\n", njet), debugtext);
    fputs("-------------------------------------------------------\n", debugtext);   

    // - 0 get electron/muon efficiency ratio (the same code as in fillRoutin)
    double Ninee(0.), NineeE(0.), Ninmm(0.), NinmmE(0.);
    Ninee = ((TH1F*) inputFile->Get(Form("Data_met_in_ww_%ij_ee", njet)))->GetEntries();
    NineeE = sqrt(Ninee);
    
    Ninmm = ((TH1F*) inputFile->Get(Form("Data_met_in_ww_%ij_mm", njet)))->GetEntries();
    NinmmE = sqrt(Ninmm);
    
    double k_ee = sqrt(Ninee/Ninmm);
    double k_eeE = sqrt(pow(NineeE/Ninmm,2) + pow(NinmmE*Ninee/(Ninmm*Ninmm),2))/2.;

    // - 1 get in-peak raw events after all signal selections but the mll
    double ndata[4];
    int bin_start = 1;
    int bin_end = 1000;
    
    ndata[0] = ((TH1F*)inputFile->Get(Form("Data_met_in_sig_%ij_mm", njet)))->Integral(bin_start,bin_end);
    ndata[1] = ((TH1F*)inputFile->Get(Form("Data_met_in_sig_%ij_me", njet)))->Integral(bin_start,bin_end); 
    ndata[2] = ((TH1F*)inputFile->Get(Form("Data_met_in_sig_%ij_em", njet)))->Integral(bin_start,bin_end); 
    ndata[3] = ((TH1F*)inputFile->Get(Form("Data_met_in_sig_%ij_ee", njet)))->Integral(bin_start,bin_end); 
    fputs(Form("In-Z peak yield (mm/me/em/ee) after all signal selections %i-Jet:\n %.0f\t %.0f\t %.0f\t %.0f\t\n", 
	       ndata[0], ndata[1], ndata[2], ndata[3], njet), debugtext);

    // - 2 do of-subtraction
    double Di_ee_subt(0.), Di_ee_subtE(0.), Di_mm_subt(0.), Di_mm_subtE(0.), Di_subt(0.), Di_subtE(0.);
    ofsubt_single(ndata, k_ee, k_eeE,  Di_ee_subt, Di_ee_subtE, Di_mm_subt, Di_mm_subtE, Di_subt, Di_subtE, debugtext, false, true);
    fputs(Form("OF subtracted yields in Z window %i-Jet:\n %.1f +/- %.1f(MM)\t %.1f +/- %.1f(EE)\t %.1f +/- %.1f (EE+MM) \n", 
	       njet, Di_mm_subt, Di_mm_subtE, Di_ee_subt, Di_ee_subtE, Di_subt, Di_subtE), debugtext);

    // get the in-peak VZ contribution
    double Di_ee_vz(0.), Di_ee_vzE(0.), Di_mm_vz(0.), Di_mm_vzE(0.), Di_vz(0.), Di_vzE(0.);
    Di_ee_vz = ((TH1F*)inputFile->Get(Form("VV_met_in_sig_%ij_ee", njet)))->IntegralAndError(bin_start,bin_end, Di_ee_vzE);
    Di_mm_vz = ((TH1F*)inputFile->Get(Form("VV_met_in_sig_%ij_mm", njet)))->IntegralAndError(bin_start,bin_end, Di_mm_vzE);
    combll(Di_mm_vz, Di_mm_vzE, Di_ee_vz, Di_ee_vzE, Di_vz, Di_vzE);
    double vzsyst = 0.1; // the uncertainty to the MC cross-section on WZ/ZZ
    double Di_ee_vzE_syst = Di_ee_vz * vzsyst;
    double Di_mm_vzE_syst = Di_mm_vz * vzsyst;
    double Di_vzE_syst = Di_vz * vzsyst;
    fputs(Form("WZ/ZZ contribution in Z window %i-Jet: \n %.1f +/- %.1f +/- %.1f (MM), %.1f +/- %.1f +/- %.1f (EE), %.1f +/- %.1f +/- %.1f (EE+MM)\n", 
	       njet, Di_mm_vz, Di_mm_vzE, Di_mm_vzE_syst, Di_ee_vz, Di_ee_vzE, Di_ee_vzE_syst, Di_vz, Di_vzE, Di_vzE_syst), debugtext);
    
    // - 3 get in-peak yield subtracting the VZ contribution
    Di_ee_subt = Di_ee_subt - Di_ee_vz;
    Di_ee_subtE = sqrt( pow(Di_ee_subtE,2) + pow(Di_ee_vzE,2) + pow(Di_ee_vzE_syst,2));
    
    Di_mm_subt = Di_mm_subt - Di_mm_vz;
    Di_mm_subtE = sqrt( pow(Di_mm_subtE,2) + pow(Di_mm_vzE,2) + pow(Di_mm_vzE_syst,2));
    
    Di_subt = Di_subt - Di_vz;
    Di_subtE = sqrt( pow(Di_subtE,2) + pow(Di_vzE,2) + pow(Di_vzE_syst,2));
    
    // if over-subtracting, fluctuate the event count to 1
    if (Di_subt <= 0.0)
      Di_subt = 1.0;
    fputs(Form("OF/VZ subtracted yields in Z window %i-Jet:\n %.1f +/- %.1f(MM)\t %.1f +/- %.1f(EE)\t %.1f +/- %.1f (EE+MM) \n", 
	       njet, Di_mm_subt, Di_mm_subtE, Di_ee_subt, Di_ee_subtE, Di_subt, Di_subtE), debugtext);
    
    // - 4 read Rout/in and do the data-driven estimate 
    double R_ee(0.), R_eeE(0.), R_eeE_syst(0.0),  R_mm(0.), R_mmE(0.), R_mmE_syst(0.), R(0.), RE(0.), RE_syst(0.);
    
//	if ( njet < 2 ) 
//      lookupR(njet,routinFileName, "mc", R_ee, R_eeE, R_eeE_syst, R_mm, R_mmE, R_mmE_syst, R, RE, RE_syst);
//    else 
      lookupR(njet,routinFileName, "data", R_ee, R_eeE, R_eeE_syst, R_mm, R_mmE, R_mmE_syst, R, RE, RE_syst);
    
	fputs("-------------------------------------------------------\n", debugtext);
    fputs(Form("Choose the nominal R as..\n"), debugtext);
    fputs(Form("R(EE) %i-Jet = %.2f +/- %.2f +/- %.2f\n", njet, R_ee, R_eeE, R_eeE_syst), debugtext);
    fputs(Form("R(MM) %i-Jet = %.2f +/- %.2f +/- %.2f\n", njet, R_mm, R_mmE, R_mmE_syst), debugtext);
    fputs(Form("R(EE+MM) %i-Jet = %.2f +/- %.2f +/- %.2f\n", njet, R, RE, RE_syst), debugtext);

    double pred_ee(0.), pred_eeE(0.), pred_eeE_syst(0.);
    getEstimates( Di_ee_subt, Di_ee_subtE, R_ee, R_eeE, R_eeE_syst, pred_ee, pred_eeE, pred_eeE_syst);

    double pred_mm(0.), pred_mmE(0.), pred_mmE_syst(0.);
    getEstimates( Di_mm_subt, Di_mm_subtE, R_mm, R_mmE, R_mmE_syst, pred_mm, pred_mmE, pred_mmE_syst);

    double pred(0.), predE(0.), predE_syst(0.);
    getEstimates( Di_subt, Di_subtE, R, RE, RE_syst, pred, predE, predE_syst);
    
    fputs("-------------------------------------------------------\n", debugtext);
    fputs(Form("data-driven estimate %i-Jet:\n %.1f +/- %.1f +/- %.1f(MM)\t %.1f +/- %.1f +/- %.1f (EE)\t %.1f +/- %.1f +/- %.1f (EE+MM)\n", 
	       njet, pred_mm, pred_mmE, pred_mmE_syst, pred_ee, pred_eeE, pred_eeE_syst, pred, predE, predE_syst), debugtext);
    
    // - 5 calculate the data/MC scale factors
    double nmmMC(0.), nmmMCE(0.0), neeMC(0.), neeMCE(0.0), nMC(0.0), nMCE(0.);
    nmmMC = ((TH1F*)inputFile->Get(Form("DYLL_met_out_sig_%ij_mm", njet)))->IntegralAndError(bin_start,bin_end, nmmMCE);
    neeMC = ((TH1F*)inputFile->Get(Form("DYLL_met_out_sig_%ij_ee", njet)))->IntegralAndError(bin_start,bin_end, neeMCE);
    combll(nmmMC, nmmMCE, neeMC, neeMCE, nMC, nMCE);
    fputs("-------------------------------------------------------\n", debugtext);
    fputs(Form("MC estimation in signal region %i-Jet:\n %.2f +/- %.2f(MM)\t %.2f +/- %.2f(EE)\t  %.2f +/- %.2f (EE+MM)\n", 
	       njet, nmmMC, nmmMCE, neeMC, neeMCE, nMC, nMCE), debugtext);
    double sf_ee = neeMC > 0. ? pred_ee / neeMC : 1.;
    double sf_eeE_syst = neeMC > 0. ? sf_ee * pred_eeE_syst/pred_ee : 0.;
    double sf_eeE = pred_ee > 0. ? sf_ee * pred_eeE/pred_ee : 0.;
    
    double sf_mm = nmmMC > 0. ? pred_mm / nmmMC : 1.;
    double sf_mmE_syst = nmmMC > 0. ? sf_mm * pred_mmE_syst/pred_mm : 0.;
    double sf_mmE = nmmMC > 0. ? sf_mm *pred_mmE/pred_mm : 0.;
    
    double sf = nMC > 0. ? pred / nMC : 1.;
    double sfE_syst = nMC > 0. ? sf * predE_syst/pred : 0.;
    double sfE = nMC > 0. ? sf * predE/pred : 0.;
    fputs(Form("data/MC scale factor from Rout/in method %i-Jet:\n %.2f +/- %.2f (MM)\t  %.2f +/- %.2f (EE)\t %.2f +/- %.2f (EE+MM)\n",
	       njet, sf_mm, sqrt(sf_mmE*sf_mmE + sf_mmE_syst*sf_mmE_syst), sf_ee, sqrt(sf_eeE*sf_eeE + sf_eeE_syst*sf_eeE_syst), 
	       sf, sqrt(sfE*sfE + sfE_syst*sfE_syst)), debugtext);
    fputs("-------------------------------------------------------\n", debugtext);
  

    // - 6 calculate the VZ contribution in the signal region
    double Do_ee_vz(0.), Do_ee_vzE(0.), Do_mm_vz(0.), Do_mm_vzE(0.), Do_vz(0.), Do_vzE(0.);
    Do_ee_vz = ((TH1F*)inputFile->Get(Form("VV_met_out_sig_%ij_ee", njet)))->IntegralAndError(bin_start,bin_end, Do_ee_vzE);
    Do_mm_vz = ((TH1F*)inputFile->Get(Form("VV_met_in_sig_%ij_mm", njet)))->IntegralAndError(bin_start,bin_end, Do_mm_vzE);
    combll(Do_mm_vz, Do_mm_vzE, Do_ee_vz, Do_ee_vzE, Do_vz, Do_vzE);
    fputs(Form("WZ/ZZ contribution in signal region %i-Jet: \n %.2f +/- %.2f (MM), %.2f +/- %.2f (EE), %.2f +/- %.2f (EE+MM)\n", 
	       njet, Do_mm_vz, Do_mm_vzE, Do_ee_vz, Do_ee_vzE, Do_vz, Do_vzE), debugtext);
	  
    // - 7 fill the information needed for the DY estimations
    if ( analysis == 0 ) {
      DYBkgScaleFactorWWPreselection[njet] = sf;
      DYBkgScaleFactorWWPreselectionKappa[njet] = 1 + sqrt(sfE*sfE+sfE_syst*sfE_syst)/sf;
    } else { 
	for (int i = 0; i < 20; i++) {
	  if ( analysis == mHiggs[i]) {
	    if ( option == HWW_OPT_SMURFCUTSEL ) {
	      DYBkgScaleFactorHiggsSelection[njet][i] = pred;
	      DYBkgScaleFactorHiggsSelectionKappa[njet][i] = sqrt(predE*predE+predE_syst*predE_syst) / pred + 1.;
	    } else if ( option == HWW_OPT_SMURFMVASEL ) {
	      DYBkgScaleFactorHiggsSelectionMVA[njet][i] = pred;
	      DYBkgScaleFactorHiggsSelectionKappaMVA[njet][i] = sqrt(predE*predE+predE_syst*predE_syst) / pred + 1.;
	    }
	  }
	}
    } 
  }
  
  gROOT->cd();
  inputFile->Close();
  
}
// calculate the R from out/in
void calcR(double Nout, double NoutE, double Nin, double NinE, double & R, double & RE) {
  R  = Nin > 0.0 ? Nout/Nin : 0.0;
  RE = Nin > 0.0 ? R*sqrt(pow(NoutE/Nout,2) + pow(NinE/Nin,2)) : 0.0;
}

void ofsubtraction (double Nll, double NllE, double Nem, double NemE, double kll, double kllE, double & Nll_subt, double & Nll_subtE) {
  Nll_subt = Nll - 0.5 * Nem * kll; 
  Nll_subtE = sqrt( pow(NllE,2)  + pow( 0.5 * NemE *kll, 2) +pow( 0.5 * Nem * kllE, 2) ) ;
}

void ofsubt_single(const double ndata[4], double k_ee, double k_eeE, double & nee_subt, double &nee_subtE, 
		   double &nmm_subt, double &nmm_subtE, double &ntot_subt, double &ntot_subtE, 
		   FILE *debugtext, bool forRatio, bool verbose)
{  
  // get the event counts with signal selections
  double ndataerr[4]; 
  for (int i=0; i < 4 ; i++) {
    ndataerr[i] = sqrt(ndata[i]);
  }
  
  // do the opposite flavor subtraction
  // fill the needed information
  double nem = ndata[1] + ndata[2];
  double nemE = sqrt(pow(ndataerr[1],2) + pow(ndataerr[2],2));
  double nee = ndata[3];
  double neeE = ndataerr[3];
  double nmm = ndata[0];
  double nmmE = ndataerr[0];
  double k_mm = 1./k_ee;
  double k_mmE = k_mm * k_eeE / k_ee;

  // for ee 
  ofsubtraction (nee, neeE, nem, nemE, k_ee, k_eeE, nee_subt, nee_subtE);
  // for mm
  ofsubtraction (nmm, nmmE, nem, nemE, k_mm, k_mmE, nmm_subt, nmm_subtE);

  // now do the sum assuming 0 uncertainty on k
  if (forRatio) {
    ntot_subt = nee + k_ee*k_ee*nmm - k_ee*nem;
    ntot_subtE = sqrt( neeE*neeE  + pow(k_ee*k_ee*nmmE,2) + pow(k_ee*nemE, 2));
  }
  else {
    ntot_subt = nee + nmm - 0.5 * nem * (k_ee + 1./k_ee);
    ntot_subtE = sqrt( neeE*neeE  + nmmE*nmmE + pow(0.5*(k_ee+1./k_ee)*nemE,2));
  }
}

void combll(double nee, double neeE, double nmm, double nmmE, double & n, double & nE) 
{
  n = nee + nmm;
  nE = sqrt( neeE*neeE + nmmE*nmmE);
}


void lookupR(int njet, const char *fileName, const char *suffix,
	     double & R_ee, double & R_eeE, double & R_eeE_syst,
	     double & R_mm, double & R_mmE, double & R_mmE_syst,
	     double & R, double & RE, double & RE_syst)
{

  TFile *file = TFile::Open(fileName,"READ");
  if (file == 0x0) 
    std::cout << "CANNOT Find " << fileName << "...! You need to fill the Rout/in..\n";
  //assert(file);
  gROOT->cd();
  TH1F *hRee = (TH1F*) file->Get(Form("Ree_vs_met_%s_%iJet", suffix, njet));  
  cout << "Ree " << njet << endl; 
  ratio_syst(hRee, R_ee, R_eeE, R_eeE_syst);
  
  TH1F *hRmm = (TH1F*) file->Get(Form("Rmm_vs_met_%s_%iJet", suffix, njet));
  cout << "Rmm " << njet << endl; 
  ratio_syst(hRmm, R_mm, R_mmE, R_mmE_syst);

  TH1F *hR = (TH1F*) file->Get(Form("R_vs_met_%s_%iJet", suffix, njet));
  cout << "R " << njet << endl; 
  ratio_syst(hR, R, RE, RE_syst);
  
  file->Close();
}

void ratio_syst(TH1F* & ratio_vs_met, double & R, double & RE_stat, double & RE_syst)
{
  int nbinsX = ratio_vs_met->GetNbinsX();

  // take the R from the bin in the signal box
  R = ratio_vs_met->GetBinContent(nbinsX-1);
  RE_stat = ratio_vs_met->GetBinError(nbinsX-1);
  RE_syst = 0.0;
  
  for(int i=1;i<=nbinsX-2;i++) {
    if ( ratio_vs_met->GetBinContent(i) > 0 ) {
      double rel_err =  ratio_vs_met->GetBinError(i) / ratio_vs_met->GetBinContent(i); 
      // protect against any bin that has low statsitics error
      //if ( i == nbinsX && rel_err > 0.4)  continue;
      RE_syst = RE_syst > TMath::Abs(ratio_vs_met->GetBinContent(i) - R) ? RE_syst : TMath::Abs(ratio_vs_met->GetBinContent(i) - R); 
	  std::cout << "Bin " << i << ": R " << R << " R(this bin) " << ratio_vs_met->GetBinContent(i) << " "  
	  			<< "syst " << RE_syst << "\n";
    }
  }
  // add difference btw Rout/in and Zeta (30 %) 
  RE_syst = TMath::Max( RE_syst, R*0.3);  
  std::cout << "Final syst = " << RE_syst << ", relative " << RE_syst/R << endl;
}


void getEstimates(double Di_subt, double  Di_subtE, double R, double RE, double RE_syst,double & pred, double & predE, double & predE_syst)
{
  pred	= R*Di_subt;
  predE = sqrt(pow(R*Di_subtE,2) +pow(RE*Di_subt,2));
  predE_syst = pred * (RE_syst / R);
}

