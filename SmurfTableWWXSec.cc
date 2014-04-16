
#include "SmurfTableWWXSec.h"
#include "core/SmurfPlotUtilities.h"
#include "SmurfScaleFactors.h"

#include "TROOT.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLine.h"

#include <fstream>

float GetLuminosity()
{

    return 19467;
    //return (5098. - 3553.);
     // first 3.5 fb
//    return 3553;
    // fill 5.1 fb
    //return 5098;
    // 2012A
    //return 619;
    // 2012B
    //return 4431;

}

float GetLumiSystematic()
{
    return 0.044;
}

float GetJetVetoEffScaleFactor()
{
    // no jet veto scale factor
    // is to be applied
    return 1.00;
}

float GetBackgroundEstimationSystematic(Option option ,DataType dataType, unsigned int jetbin, std::string flavor)
{

    //
    // for SS closure test
    //
    if (option == HWW_OPT_SSCTL) {
        if (dataType == GGWW || dataType == QQWW || dataType == VV || dataType == WZ || dataType == ZZ)   return 0.1;
        if (dataType == TOP)    return 0.2;
        if (dataType == WGAMMA || dataType == WGAMMANORM || dataType == WG3L)   return 0.3;
        if (dataType == ZLL || dataType == ZJETS)   return 0.3;
    }

    //
    // systematics for MC derived processes
    //

    double lumi         = GetLumiSystematic();
    double ptScale      = 0.015;
    double met          = 0.020;
    double pu           = 0.023;
    double trigger      = 0.015;
    //double jetVeto      = 0.047;
    double jetVeto      = sqrt(0.042*0.042+0.026*0.026+0.035*0.035);//gc where 4.2% from QCDscale_WW, 2.6% from QCDscale_WW1in, 3.5% from UEPS
    if (jetbin==1)  jetVeto      = sqrt(0.076*0.076+0.086*0.086+0.035*0.035);//gc where 7.6% from QCDscale_WW1in, 8.6% from QCDscale_WW2in, 3.5% from UEPS
    if (jetbin==2)  jetVeto      = 0.2;//gc fimxe
    double perElectron  = 0.020;
    double perMuon      = 0.015;
    double lepton = 2 * perMuon;
    if (flavor == "ee") lepton = 2 * perElectron;
    if (flavor == "em" || flavor == "me") lepton = sqrt(pow(perElectron, 2) + pow(perMuon, 2));
    double mcSyst2 = pow(lepton, 2) + pow(trigger, 2) + pow(jetVeto, 2) 
            + pow(met, 2) + pow(pu, 2) + pow(lumi, 2) + pow(ptScale, 2);
    double mcSystNoLumi2 = pow(lepton, 2) + pow(trigger, 2) + pow(jetVeto, 2)
            + pow(met, 2) + pow(pu, 2) + pow(ptScale, 2);

    //
    // Drell-Yan
    //

    if (dataType == ZLL) {
        double ZScaleFactor_[3] = {0.0, 0.0, 0.0};
	double ZScaleFactorError_[3] = {0.0, 0.0, 0.0};
	if ( flavor == "ee" || flavor == "mm" || flavor  == "sf" ) 
	  getZScaleFactor(ZScaleFactor_, ZScaleFactorError_, option, 0.0, "sf");
	else 
	  getZScaleFactor(ZScaleFactor_, ZScaleFactorError_, option, 0.0, "of");
        return ZScaleFactorError_[jetbin];
    }

    //
    // Z+jets part (not Z->ll)...
    //

    if (dataType == ZJETS) {
        return 0.1;
    }

    //
    // Top
    //

    else if (dataType == TOP) {
        double TopScaleFactor_[3] = {0.0, 0.0, 0.0};
        double TopScaleFactorError_[3] = {0.0, 0.0, 0.0};
        getTopScaleFactor(TopScaleFactor_, TopScaleFactorError_, option, 0.);
        return TopScaleFactorError_[jetbin];
    }

    //
    // W+Jets
    //

    else if (dataType == WJETSDATA || dataType == WJETSELEDATA || dataType == WJETSMUDATA) {
        double WjetsScaleFactor_[3] = {0.0, 0.0, 0.0};
        double WjetsScaleFactorError_[3] = {0.0, 0.0, 0.0};
        getWjetsScaleFactor(WjetsScaleFactor_, WjetsScaleFactorError_, option);
        return WjetsScaleFactorError_[jetbin];
    }

    //
    // WZ
    //

    else if (dataType == WZ) {
        return sqrt(mcSyst2 + 0.065*0.062 + 0.042*0.042);
    }

    //
    // ZZ
    //

    else if (dataType == ZZ) {
        return sqrt(mcSyst2 + 0.048*0.048 + 0.018*0.018);
    }

    //
    // WGamma
    //

    else if (dataType == WGAMMA || dataType == WGAMMANORM) {
        return 0.30;
    }

    //
    // WGamma*
    //

    else if (dataType == WGAMMASTAR) {
        //return sqrt(mcSyst2 + 0.30*0.30);
        return sqrt(mcSyst2 + 0.40*0.40);
    }

    //
    // qqWW
    // note - this does not include lumi systematics
    // because systematic on signal enters efficiency systematic
    // and this should not depend on lumi
    //

    else if (dataType == QQWW) {
        //pdf 2.3%, scale 1.5% //gc
        return sqrt(mcSystNoLumi2 + 0.015*0.015 + 0.023*0.023);
    }

    //
    // ggWW
    // note - this does not include lumi systematics
    // because systematic on signal enters efficiency systematic
    // and this should not depend on lumi
    //

    else if (dataType == GGWW) {
        return sqrt(mcSystNoLumi2 + 0.008*0.008 + 0.3*0.3);
    }

    return 0.0;

}

void SetYieldsAndUncertainties(Option option, double yield[][kNDataTypes], double nwght[][kNDataTypes],
        double stat2[][kNDataTypes], double syst2[][kNDataTypes], 
        unsigned int jetbin, std::vector<SmurfSample*> samples)
{

    unsigned int binMinIdx = 1;
    unsigned int binMaxIdx = 1;
    if (jetbin == 1) {
        binMinIdx = 2;
        binMaxIdx = 2;
    }
    if (jetbin == 2) { 
        binMinIdx = 3;
        binMaxIdx = 3;
    }

    // loop on samples
    for (unsigned int s = 0; s < samples.size(); ++s) {

        DataType dataType = samples[s]->getDataType();

        // loop on flavor combinations
        for (unsigned int i = 0; i < 4; ++i) {

            std::string flav = "mm";
            if (i == 1) flav = "me";
            if (i == 2) flav = "em";
            if (i == 3) flav = "ee";

            // get results
            double ytmp = 0.0;
            double stmp = 0.0;
            samples[s]->getResults(binMinIdx, binMaxIdx, i, ytmp, stmp);
            yield[i][dataType] += ytmp;
            stat2[i][dataType] += stmp*stmp;
            syst2[i][dataType] += pow(ytmp * GetBackgroundEstimationSystematic(option, dataType, jetbin, flav), 2);
            nwght[i][dataType] += samples[s]->getNWeights(binMinIdx, binMaxIdx, i);
            
        } // end loop on flavors

        // fill same flavor
        yield[fSF][dataType] = yield[fEE][dataType] + yield[fMM][dataType];
        stat2[fSF][dataType] = stat2[fEE][dataType] + stat2[fMM][dataType];
        syst2[fSF][dataType] = pow(sqrt(syst2[fEE][dataType]) + sqrt(syst2[fMM][dataType]), 2);
        nwght[fSF][dataType] = nwght[fEE][dataType] + nwght[fMM][dataType];

        // fill e-mu
        yield[fOF][dataType] = yield[fEM][dataType] + yield[fME][dataType];
        stat2[fOF][dataType] = stat2[fEM][dataType] + stat2[fME][dataType];
        syst2[fOF][dataType] = pow(sqrt(syst2[fEM][dataType]) + sqrt(syst2[fME][dataType]), 2);
        nwght[fSF][dataType] = nwght[fEM][dataType] + nwght[fME][dataType];

        // fill total
        yield[fTOTAL][dataType] = yield[fEE][dataType] + yield[fEM][dataType] + yield[fME][dataType] + yield[fMM][dataType];
        stat2[fTOTAL][dataType] = stat2[fEE][dataType] + stat2[fEM][dataType] + stat2[fME][dataType] + stat2[fMM][dataType];
        syst2[fTOTAL][dataType] = pow(sqrt(syst2[fEE][dataType]) + sqrt(syst2[fEM][dataType]) + sqrt(syst2[fME][dataType]) + sqrt(syst2[fMM][dataType]), 2);
        nwght[fTOTAL][dataType] = nwght[fEE][dataType] + nwght[fEM][dataType] + nwght[fME][dataType] + nwght[fMM][dataType];

        // total up backgrounds
        if (dataType != GGWW && dataType != QQWW && dataType != DATA
           //&& dataType != GGHWW && dataType != QQHWW && dataType != ZHWW && dataType != WHWW //gc higgs is background
           ) {

            yield[fEE][0] += yield[fEE][dataType];
            yield[fEM][0] += yield[fEM][dataType];
            yield[fME][0] += yield[fME][dataType];
            yield[fMM][0] += yield[fMM][dataType];
            yield[fSF][0] += yield[fSF][dataType];
            yield[fOF][0] += yield[fOF][dataType];
            yield[fTOTAL][0] += yield[fTOTAL][dataType];

            stat2[fEE][0] += stat2[fEE][dataType];
            stat2[fEM][0] += stat2[fEM][dataType];
            stat2[fME][0] += stat2[fME][dataType];
            stat2[fMM][0] += stat2[fMM][dataType];
            stat2[fSF][0] += stat2[fSF][dataType];
            stat2[fOF][0] += stat2[fOF][dataType];
            stat2[fTOTAL][0] += stat2[fTOTAL][dataType];

            syst2[fEE][0] += syst2[fEE][dataType];
            syst2[fEM][0] += syst2[fEM][dataType];
            syst2[fME][0] += syst2[fME][dataType];
            syst2[fMM][0] += syst2[fMM][dataType];
            syst2[fSF][0] += syst2[fSF][dataType];
            syst2[fOF][0] += syst2[fOF][dataType];
            syst2[fTOTAL][0] += syst2[fTOTAL][dataType];

            nwght[fEE][0] += nwght[fEE][dataType];
            nwght[fEM][0] += nwght[fEM][dataType];
            nwght[fME][0] += nwght[fME][dataType];
            nwght[fMM][0] += nwght[fMM][dataType];
            nwght[fSF][0] += nwght[fSF][dataType];
            nwght[fOF][0] += nwght[fOF][dataType];
            nwght[fTOTAL][0] += nwght[fTOTAL][dataType];
        
        }
            
    } // end loop on saples

}

void Tabulate(Option option, std::vector<SmurfSample*> samples, std::string filename)
{

    //
    // set up files
    //

    FILE *fout_debug;
    FILE *fout_pas;
    fout_debug = fopen (Form("../wwresults/%s.tex_debug", filename.c_str()), "w");
    fout_pas = fopen (Form("../wwresults/%s.tex", filename.c_str()), "w");

    //
    // print debug tables
    //

    PrintWWYieldTableVersion1(option, 0, samples, fout_debug);
    PrintWWYieldTableVersion1(option, 1, samples, fout_debug);
    PrintWWYieldTableVersion1(option, 2, samples, fout_debug);

    PrintWWYieldTableVersion2(option, 0, samples, fout_debug);
    PrintWWYieldTableVersion2(option, 1, samples, fout_debug);
    PrintWWYieldTableVersion2(option, 2, samples, fout_debug);

    //
    // print document header
    //

    printf("[SmurfTableWWXSec::Tabulate] Writing summary tables for PAS %s\n", filename.c_str());

    fprintf(fout_pas, "\\documentclass{cmspaper}\n");
    fprintf(fout_pas, "\\usepackage{graphicx}\n");
    fprintf(fout_pas, "\\begin{document}\n");
    fprintf(fout_pas, "\\title{Summary of WW Cross Section Results}\n");
    fprintf(fout_pas, "\\tableofcontents\n");
    fprintf(fout_pas, "\\clearpage\n");

    //
    // print yield tables
    //

    fprintf(fout_pas, "\\section{Yields}\n");
    fprintf(fout_pas, "\\subsection{Zero jet bin}\n\n");
    PrintWWYieldTable(option, 0, ((1<<fMM)|(1<<fME)|(1<<fEM)|(1<<fEE)), samples, fout_pas);    
    PrintWWYieldTable(option, 0, ((1<<fOF)|(1<<fSF)), samples, fout_pas);//gc
    PrintWWYieldTable(option, 0, (1<<fTOTAL), samples, fout_pas);    
    fprintf(fout_pas, "\\clearpage\n");
    fprintf(fout_pas, "\\subsection{One jet bin}\n\n");
    PrintWWYieldTable(option, 1, ((1<<fMM)|(1<<fME)|(1<<fEM)|(1<<fEE)), samples, fout_pas);    
    PrintWWYieldTable(option, 1, ((1<<fOF)|(1<<fSF)), samples, fout_pas);//gc
    PrintWWYieldTable(option, 1, (1<<fTOTAL), samples, fout_pas);                              
    fprintf(fout_pas, "\\clearpage\n");
    fprintf(fout_pas, "\\subsection{Two jet bin}\n\n");
    PrintWWYieldTable(option, 2, ((1<<fMM)|(1<<fME)|(1<<fEM)|(1<<fEE)), samples, fout_pas);    
    PrintWWYieldTable(option, 2, ((1<<fOF)|(1<<fSF)), samples, fout_pas);//gc
    PrintWWYieldTable(option, 2, (1<<fTOTAL), samples, fout_pas);                              
    fprintf(fout_pas, "\\clearpage\n");

    //
    // compute cross sections
    //
   
    CrossSectionResults results_0j[kLeptonTypes];
    CrossSectionResults results_1j[kLeptonTypes];
    CrossSectionResults results_2j[kLeptonTypes];
 
    // zero jet
    CalculateCrossSection(option, samples, fMM, 0, results_0j[fMM]);
    CalculateCrossSection(option, samples, fME, 0, results_0j[fME]);
    CalculateCrossSection(option, samples, fEM, 0, results_0j[fEM]);
    CalculateCrossSection(option, samples, fEE, 0, results_0j[fEE]);
    CalculateCrossSection(option, samples, fTOTAL, 0, results_0j[fTOTAL]);

    // one jet
    CalculateCrossSection(option, samples, fMM, 1, results_1j[fMM]);
    CalculateCrossSection(option, samples, fME, 1, results_1j[fME]);
    CalculateCrossSection(option, samples, fEM, 1, results_1j[fEM]);
    CalculateCrossSection(option, samples, fEE, 1, results_1j[fEE]);
    CalculateCrossSection(option, samples, fTOTAL, 1, results_1j[fTOTAL]);

    // two jet
    CalculateCrossSection(option, samples, fMM, 2, results_2j[fMM]);
    CalculateCrossSection(option, samples, fME, 2, results_2j[fME]);
    CalculateCrossSection(option, samples, fEM, 2, results_2j[fEM]);
    CalculateCrossSection(option, samples, fEE, 2, results_2j[fEE]);
    CalculateCrossSection(option, samples, fTOTAL, 2, results_2j[fTOTAL]);

    //
    // print cross section table
    //

    fprintf(fout_pas, "\\section{Cross Section Summary (All Channels)}\n\n");
    fprintf(fout_pas, "The measured cross section in the zero jet bin is:\n");
    fprintf(fout_pas, "\\begin{equation}\n");
    fprintf(fout_pas, "\\sigma_{WW}  = %4.2f \\pm %4.2f~\\mathrm{(stat.)} \\pm %4.2f~\\mathrm{(syst.)} \\pm %4.2f~\\mathrm{(lumi.)~pb}\n",
        results_0j[fTOTAL].xs_, results_0j[fTOTAL].stat_err_pb_, results_0j[fTOTAL].syst_err_pt_, results_0j[fTOTAL].lumi_err_pb_);
    fprintf(fout_pas, "\\end{equation}\n");

    fprintf(fout_pas, "\\begin{table}[!ht]\n");
    fprintf(fout_pas, "\\begin{center}\n");
    fprintf(fout_pas, "\\begin{tabular}{|l|c|c|c|}\n");
    fprintf(fout_pas, "\\hline\n");
    fprintf(fout_pas, "Channel              & 0-jet & 1-jet & 2-jet \\\\ \\hline \n");
    fprintf(fout_pas, "$\\sigma_{\\mu\\mu}$   &  $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$  & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ \\\\ \n",
        results_0j[fMM].xs_, results_0j[fMM].stat_err_pb_, results_0j[fMM].syst_err_pt_, results_0j[fMM].lumi_err_pb_,
        results_1j[fMM].xs_, results_1j[fMM].stat_err_pb_, results_1j[fMM].syst_err_pt_, results_1j[fMM].lumi_err_pb_,
        results_2j[fMM].xs_, results_2j[fMM].stat_err_pb_, results_2j[fMM].syst_err_pt_, results_2j[fMM].lumi_err_pb_);
    fprintf(fout_pas, "$\\sigma_{\\mu e}$   &  $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$  & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ \\\\ \n",
        results_0j[fME].xs_, results_0j[fME].stat_err_pb_, results_0j[fME].syst_err_pt_, results_0j[fME].lumi_err_pb_,
        results_1j[fME].xs_, results_1j[fME].stat_err_pb_, results_1j[fME].syst_err_pt_, results_1j[fME].lumi_err_pb_,
        results_2j[fME].xs_, results_2j[fME].stat_err_pb_, results_2j[fME].syst_err_pt_, results_2j[fME].lumi_err_pb_);
    fprintf(fout_pas, "$\\sigma_{e \\mu}$   &  $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$  & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ \\\\ \n",
        results_0j[fEM].xs_, results_0j[fEM].stat_err_pb_, results_0j[fEM].syst_err_pt_, results_0j[fEM].lumi_err_pb_,
        results_1j[fEM].xs_, results_1j[fEM].stat_err_pb_, results_1j[fEM].syst_err_pt_, results_1j[fEM].lumi_err_pb_,
        results_2j[fEM].xs_, results_2j[fEM].stat_err_pb_, results_2j[fEM].syst_err_pt_, results_2j[fEM].lumi_err_pb_);
    fprintf(fout_pas, "$\\sigma_{ee}$   &  $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$  & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ \\\\ \n",
        results_0j[fEE].xs_, results_0j[fEE].stat_err_pb_, results_0j[fEE].syst_err_pt_, results_0j[fEE].lumi_err_pb_,
        results_1j[fEE].xs_, results_1j[fEE].stat_err_pb_, results_1j[fEE].syst_err_pt_, results_1j[fEE].lumi_err_pb_,
        results_2j[fEE].xs_, results_2j[fEE].stat_err_pb_, results_2j[fEE].syst_err_pt_, results_2j[fEE].lumi_err_pb_);
    fprintf(fout_pas, "\\hline \\hline\n");
    fprintf(fout_pas, "$\\sigma_{tot.}$   &  $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$  & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ & $%4.2f\\pm%4.2f\\pm%4.2f\\pm%4.2f$ \\\\ \n",
        results_0j[fTOTAL].xs_, results_0j[fTOTAL].stat_err_pb_, results_0j[fTOTAL].syst_err_pt_, results_0j[fTOTAL].lumi_err_pb_,
        results_1j[fTOTAL].xs_, results_1j[fTOTAL].stat_err_pb_, results_1j[fTOTAL].syst_err_pt_, results_1j[fTOTAL].lumi_err_pb_,
        results_2j[fTOTAL].xs_, results_2j[fTOTAL].stat_err_pb_, results_2j[fTOTAL].syst_err_pt_, results_2j[fTOTAL].lumi_err_pb_);
    fprintf(fout_pas, "\\hline\n");
    fprintf(fout_pas, "\\end{tabular}\n");
    fprintf(fout_pas, "\\caption{Summary of cross section results.  Uncertainties are $\\mathrm{(stat.)} \\pm \\mathrm{(syst.)} \\pm\\mathrm{(lumi.)~pb}$.}\n");
    fprintf(fout_pas, "\\label{tab:xs_summary}\n");
    fprintf(fout_pas, "\\end{center}\n");
    fprintf(fout_pas, "\\end{table}\n");

    //
    // make a summary plot
    //

    gROOT->ProcessLine(".L tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    // theory
    TH1F *h1_theory_systdown = new TH1F ("h1_theory_systdown", "theory_systdown", 15, 0, 15);
    TH1F *h1_theory_systup = new TH1F ("h1_theory_systup", "theory_systup;;Cross Section (pb)", 15, 0, 15);
    h1_theory_systup->SetFillColor(kBlack);
    h1_theory_systup->SetLineColor(kWhite);
    h1_theory_systup->SetFillStyle(3005);
    h1_theory_systdown->SetFillColor(10);
    h1_theory_systdown->SetFillStyle(1001);
    h1_theory_systdown->SetLineColor(kWhite);
    for (unsigned int i = 1; i <= 16; ++i) {
        h1_theory_systdown->SetBinContent(i, (47.0 - 2.0)*1.2151);
        h1_theory_systup->SetBinContent(i, (47.0 + 2.0)*1.2151);
    }
    h1_theory_systup->GetXaxis()->SetBinLabel(1, "#mu#mu");
    h1_theory_systup->GetXaxis()->SetBinLabel(2, "#mue");
    h1_theory_systup->GetXaxis()->SetBinLabel(3, "e#mu");
    h1_theory_systup->GetXaxis()->SetBinLabel(4, "ee");
    h1_theory_systup->GetXaxis()->SetBinLabel(5, "tot.");
    h1_theory_systup->GetXaxis()->SetBinLabel(6, "#mu#mu");
    h1_theory_systup->GetXaxis()->SetBinLabel(7, "#mue");
    h1_theory_systup->GetXaxis()->SetBinLabel(8, "e#mu");
    h1_theory_systup->GetXaxis()->SetBinLabel(9, "ee");
    h1_theory_systup->GetXaxis()->SetBinLabel(10, "tot.");
    h1_theory_systup->GetXaxis()->SetBinLabel(11, "#mu#mu");
    h1_theory_systup->GetXaxis()->SetBinLabel(12, "#mue");
    h1_theory_systup->GetXaxis()->SetBinLabel(13, "e#mu");
    h1_theory_systup->GetXaxis()->SetBinLabel(14, "ee");
    h1_theory_systup->GetXaxis()->SetBinLabel(15, "tot.");

    // 0-jet 
    TH1F *h1_0j = new TH1F ("h1_0j", "0-jet", 15, 0, 15);
    h1_0j->SetMarkerStyle(20);
    h1_0j->SetMarkerSize(2);
    h1_0j->SetLineWidth(2);
    h1_0j->SetBinContent(1, results_0j[fMM].xs_);
    h1_0j->SetBinContent(2, results_0j[fME].xs_);
    h1_0j->SetBinContent(3, results_0j[fEM].xs_);
    h1_0j->SetBinContent(4, results_0j[fEE].xs_);
    h1_0j->SetBinError(1, results_0j[fMM].total_err_pb_);
    h1_0j->SetBinError(2, results_0j[fME].total_err_pb_);
    h1_0j->SetBinError(3, results_0j[fEM].total_err_pb_);
    h1_0j->SetBinError(4, results_0j[fEE].total_err_pb_);
    TH1F *h1_0j_incl = new TH1F ("h1_0j_incl", "0-jet incl", 15, 0, 15);
    h1_0j_incl->SetMarkerStyle(20);
    h1_0j_incl->SetMarkerSize(2);
    h1_0j_incl->SetLineWidth(2);
    h1_0j_incl->SetLineColor(kRed);
    h1_0j_incl->SetMarkerColor(kRed);
    h1_0j_incl->SetBinContent(5, results_0j[fTOTAL].xs_);
    h1_0j_incl->SetBinError(5, results_0j[fTOTAL].total_err_pb_);

    // 1-jet
    TH1F *h1_1j = new TH1F ("h1_1j", "1-jet", 15, 0, 15);
    h1_1j->SetMarkerStyle(21);
    h1_1j->SetMarkerSize(2);
    h1_1j->SetLineWidth(2);
    h1_1j->SetBinContent(6, results_1j[fMM].xs_);
    h1_1j->SetBinContent(7, results_1j[fME].xs_);
    h1_1j->SetBinContent(8, results_1j[fEM].xs_);
    h1_1j->SetBinContent(9, results_1j[fEE].xs_);
    h1_1j->SetBinError(6, results_1j[fMM].total_err_pb_);
    h1_1j->SetBinError(7, results_1j[fME].total_err_pb_);
    h1_1j->SetBinError(8, results_1j[fEM].total_err_pb_);
    h1_1j->SetBinError(9, results_1j[fEE].total_err_pb_);
    TH1F *h1_1j_incl = new TH1F ("h1_1j_incl", "1-jet incl", 15, 0, 15);
    h1_1j_incl->SetMarkerStyle(21);
    h1_1j_incl->SetMarkerSize(2);
    h1_1j_incl->SetLineWidth(2);
    h1_1j_incl->SetLineColor(kRed);
    h1_1j_incl->SetMarkerColor(kRed);
    h1_1j_incl->SetBinContent(10, results_1j[fTOTAL].xs_);
    h1_1j_incl->SetBinError(10, results_1j[fTOTAL].total_err_pb_);

    // 2-jet
    TH1F *h1_2j = new TH1F ("h1_2j", "2-jet", 15, 0, 15);
    h1_2j->SetMarkerStyle(22);
    h1_2j->SetMarkerSize(2);
    h1_2j->SetLineWidth(2);
    h1_2j->SetBinContent(11, results_2j[fMM].xs_);
    h1_2j->SetBinContent(12, results_2j[fME].xs_);
    h1_2j->SetBinContent(13, results_2j[fEM].xs_);
    h1_2j->SetBinContent(14, results_2j[fEE].xs_);
    h1_2j->SetBinError(11, results_2j[fMM].total_err_pb_);
    h1_2j->SetBinError(12, results_2j[fME].total_err_pb_);
    h1_2j->SetBinError(13, results_2j[fEM].total_err_pb_);
    h1_2j->SetBinError(14, results_2j[fEE].total_err_pb_);
    TH1F *h1_2j_incl = new TH1F ("h1_2j_incl", "2-jet incl", 15, 0, 15);
    h1_2j_incl->SetMarkerStyle(22);
    h1_2j_incl->SetMarkerSize(2);
    h1_2j_incl->SetLineWidth(2);
    h1_2j_incl->SetLineColor(kRed);
    h1_2j_incl->SetMarkerColor(kRed);
    h1_2j_incl->SetBinContent(15, results_2j[fTOTAL].xs_);
    h1_2j_incl->SetBinError(15, results_2j[fTOTAL].total_err_pb_);

    TLegend *l1 = new TLegend(0.45, 0.85, 0.98, 0.99);
    l1->SetShadowColor(kWhite);
    l1->SetLineColor(kWhite);
    l1->SetFillColor(10);
    l1->AddEntry(h1_theory_systup, "Theory", "f");
    l1->AddEntry(h1_0j, "0-jet channels +/- (stat.#oplus syst.#oplus lumi.)", "lp");
    l1->AddEntry(h1_1j, "1-jet channels +/- (stat.#oplus syst.#oplus lumi.)", "lp");
    l1->AddEntry(h1_2j, "2-jet channels +/- (stat.#oplus syst.#oplus lumi.)", "lp");

    TLatex *tex = new TLatex(0.05,0.89, Form("#sqrt{s}=8 TeV, #int Ldt = %u pb^{-1}", int(GetLuminosity())));
    tex->SetNDC();
    tex->SetTextSize(0.035);
    tex->SetLineWidth(2);
    TLatex *tex2 = new TLatex(0.05,0.95,"CMS Preliminary 2012");
    tex2->SetNDC();
    tex2->SetTextSize(0.035);
    tex2->SetLineWidth(2);

    TLine *line1 = new TLine(5.0, 0.0, 5.0, 120.0);
    line1->SetLineStyle(kDashed);
    TLine *line2 = new TLine(10.0, 0.0, 10.0, 120.0);
    line2->SetLineStyle(kDashed);

    TCanvas *c1 = new TCanvas("CrossSectionResults");
    c1->cd();
    c1->SetTopMargin(0.20);
    h1_theory_systup->Draw();
    h1_theory_systup->GetYaxis()->SetRangeUser(0.0, 120.0);
    h1_theory_systdown->Draw("SAME");
    h1_0j->Draw("SAME E1");
    h1_0j_incl->Draw("SAME E1");
    h1_1j->Draw("SAME E1");
    h1_1j_incl->Draw("SAME E1");
    h1_2j->Draw("SAME E1");
    h1_2j_incl->Draw("SAME E1");
    l1->Draw();
    tex->Draw("SAME");
    tex2->Draw("SAME");
    line1->Draw();
    line2->Draw();
    c1->RedrawAxis();
    c1->SaveAs(Form("../wwresults/%s_summary.pdf", filename.c_str()));
    c1->SaveAs(Form("../wwresults/%s_summary.png", filename.c_str()));

    delete h1_0j;
    delete h1_0j_incl;
    delete h1_1j;
    delete h1_1j_incl;
    delete h1_2j;
    delete h1_2j_incl;
    delete h1_theory_systdown;
    delete h1_theory_systup;
    delete tex;
    delete tex2;
    delete l1;
    delete line1;
    delete line2;
    delete c1;

    fprintf(fout_pas, "\\vspace{30pt}\n");
    fprintf(fout_pas, "\\begin{figure}[!hbtp]\n");
    fprintf(fout_pas, "\\centering\n");
    fprintf(fout_pas, "\\includegraphics[width=.8\\textwidth]{%s_summary.pdf}\n", filename.c_str());
    fprintf(fout_pas, "\\caption{Summary of all channels. Total uncertainty is shown.}\n");
    fprintf(fout_pas, "\\label{fig:xs_summary_figure}\n");
    fprintf(fout_pas, "\\end{figure}\n");

    //
    // print document footer
    //

    fprintf(fout_pas, "\\end{document}\n");

    // close log file
    fclose(fout_debug);
    fclose(fout_pas);

}

void CalculateCrossSection(Option option, std::vector<SmurfSample*> samples, FlavorType flavorType, unsigned int jetbin,
    CrossSectionResults &results)
{

    //
    // containers
    //

    double yield[kLeptonTypes][kNDataTypes];
    double stat2[kLeptonTypes][kNDataTypes];
    double syst2[kLeptonTypes][kNDataTypes];
    double nwght[kLeptonTypes][kNDataTypes];
    for (unsigned int i = 0; i < kNDataTypes; ++i) {
        for (unsigned int j = 0; j < kLeptonTypes; ++j) {
            yield[j][i] = 0.0;
            stat2[j][i] = 0.0;
            syst2[j][i] = 0.0;
            nwght[j][i] = 0.0;
        }
    }

    SetYieldsAndUncertainties(option, yield, nwght, stat2, syst2, jetbin, samples);
    DoCrossSection(flavorType, yield[flavorType], nwght[flavorType], stat2[flavorType], syst2[flavorType], results);

}

void DoCrossSection(FlavorType flavorType, double yield[kNDataTypes], double nwght[kNDataTypes], 
         double stat2[kNDataTypes], double syst2[kNDataTypes], CrossSectionResults &results)
{

    //
    // calculate the acceptance
    //

    float brwwlnln      = 0.108*0.108*9.0;
    float sigma_qqWW    = 47.0*(1-0.03)*1.2151;
    float sigma_ggWW    = 47.0*0.03*1.2151;
    float lumi          = GetLuminosity();
    float eff_qqWW      = yield[QQWW]/(sigma_qqWW*lumi*brwwlnln);
    eff_qqWW            *= GetJetVetoEffScaleFactor();
    float eff_ggWW      = yield[GGWW]/(sigma_ggWW*lumi*brwwlnln);
    eff_ggWW            *= GetJetVetoEffScaleFactor();
    float err_qqWW      = sqrt(eff_qqWW*(1-eff_qqWW)/ (nwght[QQWW]/eff_qqWW) );
    float err_ggWW      = sqrt(eff_ggWW*(1-eff_ggWW)/ (nwght[GGWW]/eff_ggWW) );
    float eff_WW        = (yield[QQWW] + yield[GGWW])/((sigma_qqWW+sigma_ggWW)*lumi*brwwlnln);
    eff_WW              *= GetJetVetoEffScaleFactor();
    float err_WW        = sqrt(pow(err_qqWW, 2) + pow(err_ggWW, 2));

    //
    // calculate the uncertainties
    // and the cross section
    //

    float sigmaEff          = 0.03 * sqrt(syst2[GGWW])/yield[GGWW]
                                    +  0.97 * sqrt(syst2[QQWW])/yield[QQWW];
    float effE              = eff_WW*sigmaEff;
    float sigmaLumi         = GetLumiSystematic();
    float lumiE             = lumi * sigmaLumi;
    float xs                = (yield[DATA] - yield[0]) / (lumi * eff_WW * brwwlnln);
    float stat_err_pb       = sqrt(stat2[DATA]) / (lumi * eff_WW * brwwlnln);
    float systbg_err_pb     = sqrt(stat2[0] + syst2[0]) / (lumi * eff_WW * brwwlnln);
    float lumi_err_pb       = xs  * lumiE / lumi;
    float systeff_err_pb    = xs * effE / eff_WW;
    float xs_err_pb         = sqrt( pow(stat_err_pb, 2) + pow(systbg_err_pb, 2)
                                    + pow(lumi_err_pb, 2) + pow (systeff_err_pb, 2));

    results.xs_             = xs;
    results.stat_err_pb_    = stat_err_pb;
    results.syst_err_pt_    = sqrt(systbg_err_pb*systbg_err_pb + systeff_err_pb*systeff_err_pb);
    results.lumi_err_pb_    = lumi_err_pb;
    results.eff_ww_         = eff_WW;
    results.eff_ww_err_     = effE;
    results.total_err_pb_   = xs_err_pb;
}

TCanvas *makeWWXSecStack(Option option, float analysis, SmurfSample* data, SmurfSample* qqww, SmurfSample *ggww,  SmurfSample *wz, SmurfSample *zz,
        SmurfSample* top, SmurfSample* wjetsEle, SmurfSample *wjetsMu,  SmurfSample* wgamma, SmurfSample* dy, SmurfSample *zjets, std::vector<SmurfSample*> signals,
        TFile *file, const unsigned int flav, const unsigned int njet, const char *name, const char *title, float lumi, bool log, unsigned int rebin, float wwSF)
{

    const char *fullname = Form("%s_%s_%s", name, jetbin_names[njet], types[flav]);

    TH1F *h1_signal = 0;
    if (signals.size() != 0) {
        h1_signal = GetHistogram(file, signals, fullname, rebin);
        h1_signal->SetLineWidth(1);
        h1_signal->SetLineColor(kRed);
    }

    TH1F *h1_data = GetHistogram(file, data, fullname, rebin);
    h1_data->SetMarkerStyle(20);   
    h1_data->SetLineWidth(2);
        
    TH1F *h1_qqww = GetHistogram(file, qqww, fullname, rebin);
    TH1F *h1_ggww = GetHistogram(file, ggww, fullname, rebin);
    // scale WW by a specified scale factor
    h1_qqww->Scale(wwSF);
    h1_ggww->Scale(wwSF);
    TH1F *h1_ww = (TH1F*)h1_qqww->Clone();
    h1_ww->Add(h1_ggww);
    h1_ww->SetFillColor(kAzure-9);
    h1_ww->SetLineColor(kAzure-9);
    
    TH1F *h1_wz = GetHistogram(file, wz, fullname, rebin);
    TH1F *h1_zz = GetHistogram(file, zz, fullname, rebin);
    
    TH1F *h1_vv = (TH1F*)h1_wz->Clone();
    h1_vv->Add(h1_zz);
    h1_vv->SetFillColor(kAzure-2);
    h1_vv->SetLineColor(kAzure-2);

    TH1F *h1_top = GetHistogram(file, top, fullname, rebin);
    h1_top->SetFillColor(kYellow);
    h1_top->SetLineColor(kYellow);
    
    TH1F *h1_wjetsEle = GetHistogram(file, wjetsEle, fullname, rebin);
    TH1F *h1_wjetsMu = GetHistogram(file, wjetsMu, fullname, rebin);
    TH1F *h1_wgamma = GetHistogram(file, wgamma, fullname, rebin);
    TH1F *h1_fakes = (TH1F*)h1_wjetsEle->Clone();
    h1_fakes->Add(h1_wjetsMu);
    h1_fakes->Add(h1_wgamma);
    h1_fakes->SetFillColor(kGray+1);
    h1_fakes->SetLineColor(kGray+1);
        
    // note - assuming the sf uncertainty for drell-yan in the inclusive plot
    // - this is not really correct, but will make a negligible difference to
    // how the uncertainty band on the plot looks.  Really should construct
    // the inclusive plot uncertainties from the sum in quadrature of the uncertainties
    // on the exclusive channels.  But that's a bit fiddly so neglecting it for now.
    std::string flavor = "sf";
    if (flav == fEM || flav == fME || flav == fOF) flavor = "of";
   
    // DY is tricky because right now there can be zero events in the MC template
    // in some bins, so the SF represents the data estimate directly... 
    TH1F *h1_dy = GetHistogram(file, dy, fullname, rebin);
    if (analysis > 0.0 && analysis <= 300.0) {
        double zSF[3] = {0.0, 0.0, 0.0};
        double zSFErr[3] = {0.0, 0.0, 0.0};
        getZScaleFactor(zSF, zSFErr, option, analysis, flavor.c_str());
        double dy_yield = h1_dy->Integral(0, h1_dy->GetNbinsX() + 1);
        if (dy_yield > 0.0) h1_dy->Scale(zSF[njet]/dy_yield);
    }

    TH1F *h1_zjets = GetHistogram(file, zjets, fullname, rebin);
    h1_dy->Add(h1_zjets);
    h1_dy->SetFillColor(kGreen+2);
    h1_dy->SetLineColor(kGreen+2);
    
    //
    // figure out the uncertainties
    //

    TH1F *h1_band = (TH1F*)h1_data->Clone(Form("%s_up", fullname));
    TH1F *h1_band_rel = (TH1F*)h1_data->Clone(Form("%s_up", fullname));
    TH1F *h1_ratio = (TH1F*)h1_data->Clone(Form("%s_up", fullname));
    for (unsigned int i = 1; i <= h1_data->GetNbinsX(); ++i)
    {
        float syst_qqww = GetBackgroundEstimationSystematic(option, QQWW, njet, flavor) * h1_qqww->GetBinContent(i);
        float syst_ggww = GetBackgroundEstimationSystematic(option, GGWW, njet, flavor) * h1_ggww->GetBinContent(i);
        float syst_wz = GetBackgroundEstimationSystematic(option, WZ, njet, flavor) * h1_wz->GetBinContent(i);
        float syst_zz = GetBackgroundEstimationSystematic(option, ZZ, njet, flavor) * h1_zz->GetBinContent(i);
        float syst_top = GetBackgroundEstimationSystematic(option, TOP, njet, flavor) * h1_top->GetBinContent(i);
        float syst_wjetsEle = GetBackgroundEstimationSystematic(option, WJETSELEDATA, njet, flavor) * h1_wjetsEle->GetBinContent(i);
        float syst_wjetsMu = GetBackgroundEstimationSystematic(option, WJETSMUDATA, njet, flavor) * h1_wjetsMu->GetBinContent(i);
        float syst_wgamma = GetBackgroundEstimationSystematic(option, WGAMMANORM, njet, flavor) * h1_wgamma->GetBinContent(i);

        // h1_dy is the sum of ZLL and ZJETS
        // get the systematic on the two parts separately, because they are different
        // (stat error for h1_dy is handled without extra complication)
        float syst_dy = GetBackgroundEstimationSystematic(option, ZLL, njet, flavor) * (h1_dy->GetBinContent(i) - h1_zjets->GetBinContent(i));
        float syst_zjets = GetBackgroundEstimationSystematic(option, ZJETS, njet, flavor) * h1_zjets->GetBinContent(i);

        float stat_qqww   = h1_qqww->GetBinError(i);
        float stat_ggww   = h1_ggww->GetBinError(i);
        float stat_wz     = h1_wz->GetBinError(i);
        float stat_zz     = h1_zz->GetBinError(i);
        float stat_top    = h1_top->GetBinError(i);
        float stat_wjetsEle = h1_wjetsEle->GetBinError(i);
        float stat_wjetsMu  = h1_wjetsMu->GetBinError(i);
        float stat_wgamma = h1_wgamma->GetBinError(i);
        float stat_dy     = h1_dy->GetBinError(i);
        
        float yield_total    = h1_wz->GetBinContent(i);
        yield_total         += h1_zz->GetBinContent(i);
        yield_total         += h1_top->GetBinContent(i);
        yield_total         += h1_wjetsEle->GetBinContent(i);
        yield_total         += h1_wjetsMu->GetBinContent(i);
        yield_total         += h1_wgamma->GetBinContent(i);
        yield_total         += h1_dy->GetBinContent(i);
        yield_total         += h1_qqww->GetBinContent(i);
        yield_total         += h1_ggww->GetBinContent(i);
        
        float syst_total2 = pow(syst_qqww, 2) + pow(syst_ggww, 2) + pow(syst_wz, 2) 
                            + pow(syst_zz, 2) + pow(syst_top, 2) + pow(syst_wjetsEle, 2) + pow(syst_wjetsMu, 2) + pow(syst_wgamma, 2) + pow(syst_dy, 2) + pow(syst_zjets, 2);
        float stat_total2 = pow(stat_qqww, 2) + pow(stat_ggww, 2) + pow(stat_wz, 2) 
                            + pow(stat_zz, 2) + pow(stat_top, 2) + pow(stat_wjetsEle, 2) + pow(syst_wjetsMu, 2) + pow(stat_wgamma, 2) + pow(stat_dy, 2);
        float err_total = sqrt(syst_total2 + stat_total2);
        
        h1_band->SetBinContent(i, yield_total);
        h1_band->SetBinError(i, err_total);
        h1_band_rel->SetBinContent(i, 1.0);
        if (yield_total > 0.0) {
            h1_band_rel->SetBinError(i, err_total/yield_total);
            h1_ratio->SetBinContent(i, h1_data->GetBinContent(i) / yield_total);
            h1_ratio->SetBinError(i, h1_data->GetBinError(i) / yield_total);
        } else {
            h1_band_rel->SetBinError(i, 0.0);
            h1_ratio->SetBinContent(i, 0.0);
            h1_ratio->SetBinError(i, 0.0);
        }
    }
    
    //
    // make the legend
    //

    TLegend *l1 = new TLegend(0.57, 0.65, 0.85, 0.92);
    l1->SetLineColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->SetFillColor(kWhite);
    if (signals.size() != 0) 
        l1->AddEntry(h1_signal, Form("mH = %u GeV", int(analysis)), "l");
    l1->AddEntry(h1_data, "Data", "lp");
    l1->AddEntry(h1_ww, "WW", "f");
    l1->AddEntry(h1_vv, "WZ/ZZ", "f");
    l1->AddEntry(h1_top, "Top", "f");
    l1->AddEntry(h1_fakes, "Fakes", "f");
    l1->AddEntry(h1_dy, "Drell-Yan", "f");
    l1->AddEntry(h1_band, "#sigma (stat. #oplus syst.)", "f");
    
    //
    // make the stack
    //
    
    THStack *st = new THStack(Form("st_%s", fullname), Form("Stack;%s;Events/bin", title));
    st->SetName(fullname);
    st->Add(h1_ww);
    st->Add(h1_dy);
    st->Add(h1_top);
    st->Add(h1_vv);
    st->Add(h1_fakes);
    
    //
    // make the canvas
    //

    TCanvas *c1 = new TCanvas(fullname);
    c1->cd();

    bool drawRatio = true;

    // Pad for plot
    TPad *pad1 = 0;
    if (drawRatio)  {
        pad1 = new TPad("p_main", "p_main", 0.0, 0.3, 1.0, 1.0);
        pad1->SetBottomMargin(0.04);
        pad1->SetRightMargin(0.07);
        pad1->Draw();
        pad1->cd();
    }
    st->Draw("HIST");
    st->GetXaxis()->SetTitle(title);
    if (signals.size() != 0)
        h1_signal->Draw("SAME HIST");
    h1_data->Draw("SAME E X0");
    h1_band->SetMarkerSize(0);
    h1_band->SetFillColor(kBlack);
    h1_band->SetLineColor(kBlack);
    h1_band->SetFillStyle(3005);
    h1_band->Draw("e2 same");
    l1->Draw();
    float textSize = 0.03;
    TLatex *tex = new TLatex(0.2,0.83,Form("#sqrt{s}=8 TeV, #int Ldt = %4.0f pb^{-1}", lumi));
    tex->SetNDC();
    tex->SetTextSize(textSize);
    tex->SetLineWidth(2);
    TLatex *tex2 = new TLatex(0.2,0.89,"CMS Preliminary 2012");
    tex2->SetNDC();
    tex2->SetTextSize(textSize);
    tex2->SetLineWidth(2);
    TLatex *tex3 = new TLatex(0.2,0.80, Form("%s (%s)", types[flav], jetbin_names[njet]));
    tex3->SetNDC();
    tex3->SetTextSize(textSize);
    tex3->SetLineWidth(2);
    tex->Draw("SAME");
    tex2->Draw("SAME");
    //tex3->Draw("SAME");
    if (drawRatio) st->GetXaxis()->SetLabelOffset(0.6);
    st->SetMaximum(TMath::Max(h1_band->GetMaximum()*2.0, h1_data->GetMaximum() + 3*sqrt(h1_data->GetMaximum())));
    if (log) {
        if (drawRatio) pad1->SetLogy(1);
        st->SetMinimum(0.1);
        st->SetMaximum(h1_data->GetMaximum() * 500);
    }

    // Pad for ratio
    if (drawRatio) {
        c1->cd();
        TPad *pad2 = new TPad("p_pull", "p_pull", 0.0, 0.0, 1.0, 0.3);
        pad2->Draw();
        pad2->cd();
        pad2->SetTopMargin(0.04);
        pad2->SetRightMargin(0.07);
        pad2->SetBottomMargin(0.4);
        h1_band_rel->Draw("E2");
        h1_ratio->Draw("E SAME X0");
        h1_band_rel->GetYaxis()->SetRangeUser(0.0, 2.0);
        h1_band_rel->SetMarkerSize(0);
        h1_band_rel->SetFillColor(kBlack);
        h1_band_rel->SetLineColor(kBlack);
        h1_band_rel->SetFillStyle(3005);
        h1_band_rel->GetXaxis()->SetTitle(title);
        h1_band_rel->GetXaxis()->SetLabelSize(0.12);
        h1_band_rel->GetXaxis()->SetTitleOffset(1.1);
        h1_band_rel->GetXaxis()->SetTitleSize(0.13);
    
        h1_band_rel->GetYaxis()->SetTitle("Data/MC");
        h1_band_rel->GetYaxis()->SetLabelSize(0.13);
        h1_band_rel->GetYaxis()->SetTitleSize(0.13);
        h1_band_rel->GetYaxis()->SetTitleOffset(0.6);
        h1_band_rel->GetYaxis()->SetNdivisions(402);
    
        float min = h1_band_rel->GetXaxis()->GetXmin();
        float max = h1_band_rel->GetXaxis()->GetXmax();
        TLine *line = new TLine(min, 1.0, max, 1.0);
        line->Draw();
    }
    
    c1->RedrawAxis();
    return c1;
}

void PrintWWYieldTableVersion1(Option option, unsigned int jetbin,
        std::vector<SmurfSample*> samples, FILE *fout)
{

    //
    // containers
    //

    double yield[kLeptonTypes][kNDataTypes];
    double stat2[kLeptonTypes][kNDataTypes];
    double syst2[kLeptonTypes][kNDataTypes];
    double nwght[kLeptonTypes][kNDataTypes];
    for (unsigned int i = 0; i < kNDataTypes; ++i) {
        for (unsigned int j = 0; j < kLeptonTypes; ++j) {
            yield[j][i] = 0.0;
            stat2[j][i] = 0.0;
            syst2[j][i] = 0.0;
            nwght[j][i] = 0.0;
        }
    }
    SetYieldsAndUncertainties(option, yield, nwght, stat2, syst2, jetbin, samples);

    bool includeSyst = false;
    if (!includeSyst) {
        for (unsigned int i = 0; i < kNDataTypes; ++i) {
            for (unsigned int j = 0; j < kLeptonTypes; ++j) {
                syst2[j][i] = 0.0;
            }
        }
    }

    fprintf(fout,"\\hline\n");
    fprintf(fout,"      &   data & all bkg. & $qq \\to \\WW$ & $gg \\to \\WW$ &  $\\ttbar/tW$    & $\\Wjets$    \\\\ \n");
    fprintf(fout,"\\hline\n");
    fprintf(fout,"\\hline\n");

    fprintf(fout," %u-jet  & ", jetbin);
    fprintf(fout,"%4.0f & ", yield[fTOTAL][DATA]);
    fprintf(fout,"%4.1f $\\pm$ %4.1f & ", yield[fTOTAL][0] + yield[fTOTAL][QQWW] + yield[fTOTAL][GGWW],
         sqrt(stat2[fTOTAL][0] + syst2[fTOTAL][0] + stat2[fTOTAL][QQWW] + syst2[fTOTAL][QQWW] + stat2[fTOTAL][GGWW] + syst2[fTOTAL][GGWW]));
    fprintf(fout,"%4.1f $\\pm$ %4.1f & ", yield[fTOTAL][QQWW], sqrt(stat2[fTOTAL][QQWW] + syst2[fTOTAL][QQWW]));
    fprintf(fout,"%4.1f $\\pm$ %4.1f & ", yield[fTOTAL][GGWW], sqrt(stat2[fTOTAL][GGWW] + syst2[fTOTAL][GGWW]));
    fprintf(fout,"%4.1f $\\pm$ %4.1f & ", yield[fTOTAL][TOP], sqrt(stat2[fTOTAL][TOP] + syst2[fTOTAL][TOP]));
    fprintf(fout,"%4.1f $\\pm$ %4.1f \\\\ \n ", yield[fTOTAL][WJETSELEDATA]+yield[fTOTAL][WJETSMUDATA], 
            sqrt(stat2[fTOTAL][WJETSMUDATA] + syst2[fTOTAL][WJETSMUDATA] + stat2[fTOTAL][WJETSELEDATA] + syst2[fTOTAL][WJETSELEDATA]));

    fprintf(fout,"\\hline\n");
    fprintf(fout,"   & $WZ/ZZ$ & \\dyll & $W+\\gamma/\\gamma^*$ & \\\\ \n");
    fprintf(fout,"\\hline\n");
    fprintf(fout,"\\hline\n");

    fprintf(fout," %u-jet  & ", jetbin);
    fprintf(fout,"%4.1f $\\pm$ %4.1f & ", yield[fTOTAL][WZ] + yield[fTOTAL][ZZ], sqrt(stat2[fTOTAL][WZ] + syst2[fTOTAL][WZ] + stat2[fTOTAL][ZZ] + syst2[fTOTAL][ZZ]));
    fprintf(fout,"%4.1f $\\pm$ %4.1f & ", yield[fTOTAL][ZLL] + yield[fTOTAL][ZJETS], sqrt(stat2[fTOTAL][ZLL] + syst2[fTOTAL][ZLL] + stat2[fTOTAL][ZJETS] + syst2[fTOTAL][ZJETS]));
    fprintf(fout,"%4.1f $\\pm$ %4.1f \\\\ \n ", yield[fTOTAL][WGAMMANORM]+yield[fTOTAL][WGAMMASTAR], sqrt(stat2[fTOTAL][WGAMMANORM] + stat2[fTOTAL][WGAMMASTAR] + syst2[fTOTAL][WGAMMASTAR] + syst2[fTOTAL][WGAMMANORM]));

}

void PrintWWYieldTableVersion2(Option option, unsigned int jetbin,
        std::vector<SmurfSample*> samples, FILE *fout)
{

    //
    // containers
    //

    double yield[kLeptonTypes][kNDataTypes];
    double stat2[kLeptonTypes][kNDataTypes];
    double syst2[kLeptonTypes][kNDataTypes];
    double nwght[kLeptonTypes][kNDataTypes];
    for (unsigned int i = 0; i < kNDataTypes; ++i) {
        for (unsigned int j = 0; j < kLeptonTypes; ++j) {
            yield[j][i] = 0.0;
            stat2[j][i] = 0.0;
            syst2[j][i] = 0.0;
            nwght[j][i] = 0.0;
        }
    }
    SetYieldsAndUncertainties(option, yield, nwght, stat2, syst2, jetbin, samples);

    bool includeSyst = false;
    if (!includeSyst) {
        for (unsigned int i = 0; i < kNDataTypes; ++i) {
            for (unsigned int j = 0; j < kLeptonTypes; ++j) {
                syst2[j][i] = 0.0;
            }
        }
    }

    fprintf(fout," %u-jet  & ", jetbin);
    fprintf(fout,"%4.0f & ", yield[fTOTAL][DATA]);
    fprintf(fout,"%4.1f $\\pm$ %4.1f & ", yield[fTOTAL][0] + yield[fTOTAL][QQWW] + yield[fTOTAL][GGWW],
         sqrt(stat2[fTOTAL][0] + syst2[fTOTAL][0] + stat2[fTOTAL][QQWW] + syst2[fTOTAL][QQWW] + stat2[fTOTAL][GGWW] + syst2[fTOTAL][GGWW]));
    fprintf(fout,"%4.1f $\\pm$ %4.1f & ", yield[fTOTAL][QQWW] + yield[fTOTAL][GGWW], 
        sqrt(stat2[fTOTAL][QQWW] + syst2[fTOTAL][QQWW] + stat2[fTOTAL][GGWW] + syst2[fTOTAL][GGWW]));
    fprintf(fout,"%4.1f $\\pm$ %4.1f & ", yield[fTOTAL][TOP], sqrt(stat2[fTOTAL][TOP] + syst2[fTOTAL][TOP]));
    fprintf(fout,"%4.1f $\\pm$ %4.1f & ", yield[fTOTAL][WJETSELEDATA]+yield[fTOTAL][WJETSMUDATA], 
            sqrt(stat2[fTOTAL][WJETSMUDATA] + syst2[fTOTAL][WJETSMUDATA] + stat2[fTOTAL][WJETSELEDATA] + syst2[fTOTAL][WJETSELEDATA]));

    fprintf(fout,"%4.1f $\\pm$ %4.1f & ", yield[fTOTAL][WZ] + yield[fTOTAL][ZZ] + yield[fTOTAL][WGAMMANORM]+yield[fTOTAL][WGAMMASTAR], 
        sqrt(stat2[fTOTAL][WZ] + syst2[fTOTAL][WZ] + stat2[fTOTAL][ZZ] + syst2[fTOTAL][ZZ] + stat2[fTOTAL][WGAMMANORM] + syst2[fTOTAL][WGAMMANORM] + stat2[fTOTAL][WGAMMASTAR] + syst2[fTOTAL][WGAMMASTAR]));
    fprintf(fout,"%4.1f $\\pm$ %4.1f & ", yield[fTOTAL][ZLL] + yield[fTOTAL][ZJETS], sqrt(stat2[fTOTAL][ZLL] + syst2[fTOTAL][ZLL] + stat2[fTOTAL][ZJETS] + syst2[fTOTAL][ZJETS]));
    fprintf(fout,"%4.1f $\\pm$ %4.1f \\\\ \n ", yield[fTOTAL][GGHWW] + yield[fTOTAL][QQHWW] + yield[fTOTAL][ZHWW] + yield[fTOTAL][WHWW], 
        sqrt(stat2[fTOTAL][GGHWW] + syst2[fTOTAL][GGHWW] + stat2[fTOTAL][QQHWW] + syst2[fTOTAL][QQHWW] + stat2[fTOTAL][ZHWW] + syst2[fTOTAL][ZHWW] + stat2[fTOTAL][WHWW] + syst2[fTOTAL][WHWW]));

}


void PrintSameSignClosureTest(Option option, unsigned int jetbin, unsigned int flavors,
        std::vector<SmurfSample*> samples, FILE *fout)
{
    
    //
    // containers
    //  
            
    double yield[kLeptonTypes][kNDataTypes];
    double stat2[kLeptonTypes][kNDataTypes];
    double syst2[kLeptonTypes][kNDataTypes];
    double nwght[kLeptonTypes][kNDataTypes];
    for (unsigned int i = 0; i < kNDataTypes; ++i) { 
        for (unsigned int j = 0; j < kLeptonTypes; ++j) {
            yield[j][i] = 0.0;
            stat2[j][i] = 0.0;
            syst2[j][i] = 0.0;
            nwght[j][i] = 0.0;
        }
    }
    SetYieldsAndUncertainties(option, yield, nwght, stat2, syst2, jetbin, samples);
        
    //  
    // print table header
    //

    fprintf(fout, "\\begin{table}[!ht]\n");
    fprintf(fout, "{\\small\n");
    fprintf(fout, "\\begin{center}\n");
    fprintf(fout, "\\begin{tabular}{|l|c|c|c|c|}\n");
    fprintf(fout, "\\hline\n");
    fprintf(fout, "Sample\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& %s\t", types[flav]);
    }
    fprintf(fout, "\\\\ \\hline\n"); 

    //
    // print table rows
    //

    fprintf(fout, "$qqWW$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][QQWW], sqrt(stat2[flav][QQWW]), sqrt(syst2[flav][QQWW]));
    }
    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$ggWW$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][GGWW], sqrt(stat2[flav][GGWW]), sqrt(syst2[flav][GGWW]));
    }
    fprintf(fout, "\\\\ \n");

    fprintf(fout, "$t\\bar{t} + tW$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][TOP], sqrt(stat2[flav][TOP]), sqrt(syst2[flav][TOP]));
    }
    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$VV$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][VV], sqrt(stat2[flav][VV]), sqrt(syst2[flav][VV]));
    }
    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$Z/\\gamma*$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
        yield[flav][ZLL] + yield[flav][ZJETS],
        sqrt(stat2[flav][ZLL] + stat2[flav][ZJETS]),
        sqrt(syst2[flav][ZLL] + syst2[flav][ZJETS]));
    }
    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$Z\\rightarrow\\tau\\tau$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
        yield[flav][ZTT],
        sqrt(stat2[flav][ZTT]),
        sqrt(syst2[flav][ZTT]));
    }
    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$W\\gamma*/W+\\gamma$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
        yield[flav][WGAMMANORM] + yield[flav][WG3L],
        sqrt(stat2[flav][WGAMMANORM] + stat2[flav][WG3L]),
        sqrt(syst2[flav][WGAMMANORM] + syst2[flav][WG3L]));
    }
    fprintf(fout, "\\\\ \\hline \n");

    fprintf(fout, "\\hline \n");
    fprintf(fout, "$W+jets$ (e)\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][WJETSELEDATA], sqrt(stat2[flav][WJETSELEDATA]), sqrt(syst2[flav][WJETSELEDATA]));
    }

    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$W+jets$ ($\\mu$)\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][WJETSMUDATA], sqrt(stat2[flav][WJETSMUDATA]), sqrt(syst2[flav][WJETSMUDATA]));
    }
    fprintf(fout, "\\\\ \\hline \n");


    //
    // print special rows at bottom of table
    //
    
    fprintf(fout, "\\hline \\hline\n");
    fprintf(fout, "Data\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%u$ \t", int(yield[flav][DATA]));
    }
    fprintf(fout, "\\\\ \n");

    fprintf(fout, "\\hline \\hline \n");
    fprintf(fout, "Data - Non W+Jets\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;

        float bg = yield[flav][QQWW] + yield[flav][GGWW] + yield[flav][TOP] + yield[flav][ZTT]
                + yield[flav][VV] + syst2[flav][ZJETS] + yield[flav][ZLL] + yield[flav][WGAMMANORM] + yield[flav][WG3L];
        
        float bgstat2 = stat2[flav][QQWW] + stat2[flav][GGWW] + stat2[flav][TOP] + stat2[flav][ZTT]
                + stat2[flav][VV] + syst2[flav][ZJETS] + stat2[flav][ZLL] + stat2[flav][WGAMMANORM] + stat2[flav][WG3L];

        float bgsyst2 = syst2[flav][QQWW] + syst2[flav][GGWW] + syst2[flav][TOP] + syst2[flav][ZTT]
                + syst2[flav][VV] + syst2[flav][ZJETS] + syst2[flav][ZLL] + syst2[flav][WGAMMANORM] + syst2[flav][WG3L];

        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][DATA] - bg, sqrt(stat2[flav][DATA]+bgstat2), sqrt(bgsyst2));
    }
    fprintf(fout, "\\\\ \n");

    fprintf(fout, "\\hline \n");
    fprintf(fout, "Total W+Jets\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][WJETSELEDATA]+yield[flav][WJETSMUDATA], sqrt(stat2[flav][WJETSELEDATA]+stat2[flav][WJETSMUDATA]),
                            sqrt(syst2[flav][WJETSELEDATA]+syst2[flav][WJETSMUDATA]));
    }
    fprintf(fout, "\\\\ \n");

    fprintf(fout, "\\hline \\hline \n");
    fprintf(fout, "Agreement\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;

        // the w+jets prediction
        float wjets = yield[flav][WJETSELEDATA]+yield[flav][WJETSMUDATA];

        // the non w+jets prediction and its uncertainties
        float bg = yield[flav][QQWW] + yield[flav][GGWW] + yield[flav][TOP] + yield[flav][ZTT]
                + yield[flav][VV] + syst2[flav][ZJETS] + yield[flav][ZLL] + yield[flav][WGAMMANORM] + yield[flav][WG3L];
        float bgstat2 = stat2[flav][QQWW] + stat2[flav][GGWW] + stat2[flav][TOP] + stat2[flav][ZTT]
                + stat2[flav][VV] + syst2[flav][ZJETS] + stat2[flav][ZLL] + stat2[flav][WGAMMANORM] + stat2[flav][WG3L];
        float bgsyst2 = syst2[flav][QQWW] + syst2[flav][GGWW] + syst2[flav][TOP] + syst2[flav][ZTT]
                + syst2[flav][VV] + syst2[flav][ZJETS] + syst2[flav][ZLL] + syst2[flav][WGAMMANORM] + syst2[flav][WG3L];

        // error on data - non w+jets
        float err_dat = sqrt(stat2[flav][DATA] + bgstat2 + bgsyst2);

        // error on the w+jets prediction
        float err_wjets = sqrt(stat2[flav][WJETSELEDATA]+stat2[flav][WJETSMUDATA]+syst2[flav][WJETSELEDATA]+syst2[flav][WJETSMUDATA]);

        // ratio of prediction to (data minus other backgrounds)
        float ratio = wjets / (yield[flav][DATA] - bg);

        // and its uncertainty
        float err_ratio = ratio * sqrt( pow(err_dat/(yield[flav][DATA] - bg), 2) + pow(err_wjets/wjets, 2));

        fprintf(fout, "& $%4.2f \\pm %4.2f $\t", ratio, err_ratio);

    }
    fprintf(fout, "\\\\ \\hline\n");

 
    //
    // print table footer
    //
    
    fprintf(fout, "\\end{tabular}\n");
    fprintf(fout, Form("\\caption{Summary of same sign closure test for %u-jet channel.}\n", jetbin));
    fprintf(fout, "\\end{center}}\n");
    fprintf(fout, "\\end{table}\n");

}



void PrintWWYieldTable(Option option, unsigned int jetbin, unsigned int flavors, 
        std::vector<SmurfSample*> samples, FILE *fout)
{

    //
    // containers
    //

    double yield[kLeptonTypes][kNDataTypes];
    double stat2[kLeptonTypes][kNDataTypes];
    double syst2[kLeptonTypes][kNDataTypes];
    double nwght[kLeptonTypes][kNDataTypes];
    for (unsigned int i = 0; i < kNDataTypes; ++i) {
        for (unsigned int j = 0; j < kLeptonTypes; ++j) {
            yield[j][i] = 0.0;
            stat2[j][i] = 0.0;
            syst2[j][i] = 0.0;
            nwght[j][i] = 0.0;
        }
    }
    SetYieldsAndUncertainties(option, yield, nwght, stat2, syst2, jetbin, samples);

    //
    // print table header
    //

    fprintf(fout, "\\begin{table}[!ht]\n");
    fprintf(fout, "{\\small\n");
    fprintf(fout, "\\begin{center}\n");
    fprintf(fout, "\\begin{tabular}{|l|c|c|c|c|}\n");
    fprintf(fout, "\\hline\n");
    fprintf(fout, "Sample\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& %s\t", types[flav]);
    }
    fprintf(fout, "\\\\ \\hline\n");

    //
    // print table rows
    //

    fprintf(fout, "$qqWW$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][QQWW], sqrt(stat2[flav][QQWW]), sqrt(syst2[flav][QQWW]));
    }
    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$qqWW$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][GGWW], sqrt(stat2[flav][GGWW]), sqrt(syst2[flav][GGWW]));
    }
    fprintf(fout, "\\\\ \n");

    //gc add higgs
    fprintf(fout, "$qqHWW$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][QQHWW], sqrt(stat2[flav][QQHWW]), sqrt(syst2[flav][QQHWW]));
    }
    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$ggHWW$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][GGHWW], sqrt(stat2[flav][GGHWW]), sqrt(syst2[flav][GGHWW]));
    }
    fprintf(fout, "\\\\ \n");
    //gc add higgs

    fprintf(fout, "$t\\bar{t} + tW$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][TOP], sqrt(stat2[flav][TOP]), sqrt(syst2[flav][TOP]));
    }

    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$W+jets$ (e)\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][WJETSELEDATA], sqrt(stat2[flav][WJETSELEDATA]), sqrt(syst2[flav][WJETSELEDATA]));
    }

    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$W+jets$ ($\\mu$)\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][WJETSMUDATA], sqrt(stat2[flav][WJETSMUDATA]), sqrt(syst2[flav][WJETSMUDATA]));
    }

    //fprintf(fout, "\\\\ \n");
    //fprintf(fout, "$W+jets$\t");
    //for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
    //    if (((1<<flav) & flavors) != (1<<flav)) continue;
    //    fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
    //        yield[flav][WJETSDATA], sqrt(stat2[flav][WJETSDATA]), sqrt(syst2[flav][WJETSDATA]));
    //}

    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$WZ$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][WZ], sqrt(stat2[flav][WZ]), sqrt(syst2[flav][WZ]));
    }
    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$ZZ$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
            yield[flav][ZZ], sqrt(stat2[flav][ZZ]), sqrt(syst2[flav][ZZ]));
    }
    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$Z/\\gamma*$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
	yield[flav][ZLL] + yield[flav][ZJETS] + yield[flav][ZTT],
        sqrt(stat2[flav][ZLL] + stat2[flav][ZJETS] + stat2[flav][ZTT]),
        sqrt(syst2[flav][ZLL] + syst2[flav][ZJETS] + syst2[flav][ZTT]));
    }
    fprintf(fout, "\\\\ \n");
    fprintf(fout, "$W\\gamma*/W+\\gamma$\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t",
        yield[flav][WGAMMANORM] + yield[flav][WGAMMASTAR],
        sqrt(stat2[flav][WGAMMANORM] + stat2[flav][WGAMMASTAR]),
        sqrt(syst2[flav][WGAMMANORM] + syst2[flav][WGAMMASTAR]));
    }
    fprintf(fout, "\\\\ \n");

    //
    // print special rows at bottom of table
    //

    fprintf(fout, "\\hline \\hline \n");
    fprintf(fout, "Total B.\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t", 
            yield[flav][0], sqrt(stat2[flav][0]), sqrt(syst2[flav][0]));
    }
    fprintf(fout, "\\\\ \\hline \\hline \n");
    fprintf(fout, "Total B.+S.\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f $\t", 
            yield[flav][0] + yield[flav][QQWW]+yield[flav][GGWW], 
            sqrt(stat2[flav][0] + stat2[flav][QQWW] + stat2[flav][GGWW]),
            sqrt(syst2[flav][0] + syst2[flav][QQWW] + syst2[flav][GGWW]));
    }
    fprintf(fout, "\\\\ \\hline \\hline\n");
    fprintf(fout, "Data\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%u$ \t", int(yield[flav][DATA]));
    }

    //
    // now for the cross section row
    //

    CrossSectionResults results[kLeptonTypes];
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        DoCrossSection(FlavorType(flav), yield[flav], nwght[flav], stat2[flav], syst2[flav], results[flav]);
    } 

    fprintf(fout, "\\\\ \\hline \\hline\n");
    fprintf(fout, "Acceptance ( \\% )\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \t$", results[flav].eff_ww_ * 100, results[flav].eff_ww_err_ * 100);        
    }
    fprintf(fout, "\\\\ \n");    
    fprintf(fout, "Cross Section ( pb )\t");
    for (unsigned int flav = 0; flav < kLeptonTypes; ++flav) {
        if (((1<<flav) & flavors) != (1<<flav)) continue;
        fprintf(fout, "& $%4.2f \\pm %4.2f \\pm %4.2f$ \t", results[flav].xs_, results[flav].stat_err_pb_, results[flav].syst_err_pt_);
    }    
    fprintf(fout, "\\\\ \\hline\n");

    //
    // print table footer
    //

    fprintf(fout, "\\end{tabular}\n");
    fprintf(fout, Form("\\caption{Summary of yields for %u-jet channel.", jetbin));
    fprintf(fout, "Uncertainties on yields and cross sections are $\\mathrm{(stat.)} \\pm \\mathrm{(syst.)}$.");
    fprintf(fout, "The systematic uncertainty on the cross section does not include the luminosity}\n");
    fprintf(fout, Form("\\label{tab:datayields_wwxsec_%uj}\n", jetbin));
    fprintf(fout, "\\end{center}}\n");
    fprintf(fout, "\\end{table}\n");

}

