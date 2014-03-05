//#include "SmurfTableWWXSec.h"
//#include "core/SmurfPlotUtilities.h"
//#include "SmurfScaleFactors.h"
//#include "core/SmurfSample.h"
//#include "core/Enums.h"

//#include "TCanvas.h"
//#include "TFile.h"
//#include "TROOT.h"
//#include "TLegend.h"
//#include "THStack.h"
//#include "TLatex.h"
//#include "TLine.h"

//#include <fstream>

static const char*          types[7]                = {"mm","me","em","ee", "sf", "of", "incl"};
static const char*          jetbin_names[3]         = { "0j", "1j", "2j"};

enum DataType {
    GGHWW=1,
    QQHWW=2,
    GGHZZ=3,
    QQHZZ=4,
    WJETS=5,
    ZZ=6,
    WZ=7,
    QQWW=8,
    GGWW=9,
    WW=10,
    TOP=11,
    ZLL=12,
    ZTT=13,
    GAMMA=14,
    DATA=15,
    ZHWW=16,
    WHWW=17,
    ZJETS=18,
    WGAMMA=19,
    VV=20,
    WJETSDATA=21, 
    ZLLLOOSEMET=22,
    TOPDATA=23,
    WWMCNLO=24,
    WWMCNLOUP=25,
    WWMCNLODOWN=26,
    WZALTER=27,
    ZZALTER=28,
    TOPALTER=29, 
    WGAMMALOOSE=30,
    WJETSMCLOOSE=31,
    GGZZ=32,
    OFDATA=33,
    WGAMMASTAR=34,
    WWTOPMC=35,
    GAMMAMC=36,
    ZVVGAMMA=37,
    ZLLGAMMA=38,
    WJETSGAMMA=39,
    TT=40,
    TW=41,
    ZLLDATA=42,
    GGHWWREF=43,
    GGHWWJHU=44,
    WG3L=45,
    WGAMMANORM=46,
    DPSWW=47,
    WJETSELEDATA=48,
    WJETSMUDATA=49,
};

enum Option {
    
    // HWW
    WW_OPT_SMURFXSECSEL = 0,

    // HWW analysis
    HWW_OPT_SMURFPRESEL = 1,
    HWW_OPT_SMURFCUTSEL = 2,
    HWW_OPT_MT2DMLL     = 3,
    
    // HWW control region 
    HWW_OPT_SSCTL       = 4,
    HWW_OPT_SSCTL2D     = 5,
    HWW_OPT_TOPTAG      = 6,
    
    // XWW  
    HWW_OPT_MT2DMLL_JCP = 7,
    XWW_OPT_MT2DMLL_JCP = 8,
    
};

float GetBackgroundEstimationSystematic_8TeV(Option option, DataType dataType, unsigned int jetbin, std::string flavor)
{

    //
    // systematics for MC derived processes
    //

    double lumi         = 0.044; 
    double ptScale      = 0.015;
    double met          = 0.020;
    double pu           = 0.023;
    double trigger      = 0.015;
    double jetVeto      = 0.047;
    double perElectron  = 0.020;
    double perMuon      = 0.015;
    double lepton = 2 * perMuon;
    if (flavor == "ee") lepton = 2 * perElectron;
    if (flavor == "em" || flavor == "me" || flavor == "of") lepton = sqrt(pow(perElectron, 2) + pow(perMuon, 2));
    double mcSyst2 = pow(lepton, 2) + pow(trigger, 2) + pow(jetVeto, 2) 
            + pow(met, 2) + pow(pu, 2) + pow(lumi, 2) + pow(ptScale, 2);
    double mcSystNoLumi2 = pow(lepton, 2) + pow(trigger, 2) + pow(jetVeto, 2)
            + pow(met, 2) + pow(pu, 2) + pow(ptScale, 2);
    // 
    // Bkg scale factors 
    // 
    Double_t TopBkgScaleFactorKappa[2]  = { 1.19259, 1.03155 }; // 8 TeV
    Double_t WWBkgScaleFactorKappa[2]   = { 1.04863, 1.09097 }; // 8 TeV
/*
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

*/
    //
    // W+Jets
    //

    

    if (dataType == WJETSDATA || dataType == WJETSELEDATA || dataType == WJETSMUDATA) {
//        double WjetsScaleFactor_[3] = {0.0, 0.0, 0.0};
//        double WjetsScaleFactorError_[3] = {0.0, 0.0, 0.0};
//        getWjetsScaleFactor(WjetsScaleFactor_, WjetsScaleFactorError_, option);
//        return WjetsScaleFactorError_[jetbin]; 
        return sqrt(0.36*0.36);
    }
    
    //
    // Top
    //

    else if (dataType == TOP) {
        //double TopScaleFactor_[3] = {0.0, 0.0, 0.0};
        //double TopScaleFactorError_[3] = {0.0, 0.0, 0.0};
        //getTopScaleFactor(TopScaleFactor_, TopScaleFactorError_, option, 0.);
        //return TopScaleFactorError_[jetbin];
        return TopBkgScaleFactorKappa[jetbin]-1;
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

    else if (dataType == WGAMMA) {
        return 0.30;
    }

    //
    // WGamma*
    //

    else if (dataType == WGAMMASTAR) {
        return sqrt(mcSyst2 + 0.40*0.40);
    }
    
    else if (dataType == ZTT) { 
        return 0.1;
    }

    //
    // qqWW
    // note - this does not include lumi systematics
    // because systematic on signal enters efficiency systematic
    // and this should not depend on lumi
    //

    else if (dataType == QQWW) {
        return sqrt( mcSystNoLumi2 + 0.015*0.015 + 0.023*0.023 + (WWBkgScaleFactorKappa[jetbin]-1)*(WWBkgScaleFactorKappa[jetbin]-1) );
    }

    //
    // ggWW
    // note - this does not include lumi systematics
    // because systematic on signal enters efficiency systematic
    // and this should not depend on lumi
    //

    else if (dataType == GGWW) {
        return sqrt( mcSystNoLumi2 + 0.008*0.008 + 0.3*0.3 + (WWBkgScaleFactorKappa[jetbin]-1)*(WWBkgScaleFactorKappa[jetbin]-1) );
    }

    return 0.0;

}

float GetBackgroundEstimationSystematic_7TeV(Option option, DataType dataType, unsigned int jetbin, std::string flavor)
{

    //
    // systematics for MC derived processes
    //

    double lumi         = 0.022; 
    double ptScale      = 0.015;
    double met          = 0.020;
    double pu           = 0.023;
    double trigger      = 0.015;
    double jetVeto      = 0.047;
    double perElectron  = 0.020;
    double perMuon      = 0.015;
    double lepton = 2 * perMuon;
    if (flavor == "ee") lepton = 2 * perElectron;
    if (flavor == "em" || flavor == "me" || flavor == "of" ) lepton = sqrt(pow(perElectron, 2) + pow(perMuon, 2));
    double mcSyst2 = pow(lepton, 2) + pow(trigger, 2) + pow(jetVeto, 2) 
            + pow(met, 2) + pow(pu, 2) + pow(lumi, 2) + pow(ptScale, 2);
    double mcSystNoLumi2 = pow(lepton, 2) + pow(trigger, 2) + pow(jetVeto, 2)
            + pow(met, 2) + pow(pu, 2) + pow(ptScale, 2);
    // 
    // Bkg scale factors 
    // 
    Double_t TopBkgScaleFactorKappa[2]  = { 1.21692, 1.06435 }; // 7 TeV
    Double_t WWBkgScaleFactorKappa[2]   = { 1.0575,  1.14452 }; // 7 TeV
/*
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

*/
    //
    // W+Jets
    //

    

    if (dataType == WJETSDATA || dataType == WJETSELEDATA || dataType == WJETSMUDATA) {
//        double WjetsScaleFactor_[3] = {0.0, 0.0, 0.0};
//        double WjetsScaleFactorError_[3] = {0.0, 0.0, 0.0};
//        getWjetsScaleFactor(WjetsScaleFactor_, WjetsScaleFactorError_, option);
//        return WjetsScaleFactorError_[jetbin]; 
        return sqrt(0.36*0.36);
    }
    
    //
    // Top
    //

    else if (dataType == TOP) {
        //double TopScaleFactor_[3] = {0.0, 0.0, 0.0};
        //double TopScaleFactorError_[3] = {0.0, 0.0, 0.0};
        //getTopScaleFactor(TopScaleFactor_, TopScaleFactorError_, option, 0.);
        //return TopScaleFactorError_[jetbin];
        return TopBkgScaleFactorKappa[jetbin]-1;
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

    else if (dataType == WGAMMA) {
        return 0.30;
    }

    //
    // WGamma*
    //

    else if (dataType == WGAMMASTAR) {
        return sqrt(mcSyst2 + 0.40*0.40);
    }
    
    else if (dataType == ZTT) { 
        return 0.1;
    }

    //
    // qqWW
    // note - this does not include lumi systematics
    // because systematic on signal enters efficiency systematic
    // and this should not depend on lumi
    //

    else if (dataType == QQWW) {
        return sqrt( mcSystNoLumi2 + 0.015*0.015 + 0.023*0.023 + (WWBkgScaleFactorKappa[jetbin]-1)*(WWBkgScaleFactorKappa[jetbin]-1) );
    }

    //
    // ggWW
    // note - this does not include lumi systematics
    // because systematic on signal enters efficiency systematic
    // and this should not depend on lumi
    //

    else if (dataType == GGWW) {
        return sqrt( mcSystNoLumi2 + 0.008*0.008 + 0.3*0.3 + (WWBkgScaleFactorKappa[jetbin]-1)*(WWBkgScaleFactorKappa[jetbin]-1) );
    }

    return 0.0;

}

// 
// Convert TGraphAsymmErrors to TH1D
// 
TH1D* makeHistogram(TGraphAsymmErrors* graph){
  
    //TH1D* hist = new TH1D(Form("h_%s",graph->GetName()),Form("%s",graph->GetName()),126,-1,1);
    TH1D* hist = new TH1D(Form("h_%s", graph->GetName()), Form("%s",graph->GetName()), 126, 0.5, 126.5);
    hist->SetDirectory(0); 
    for (int i=0; i<graph->GetN();++i){
        Int_t bin = hist->FindBin(graph->GetX()[i]);
        hist->SetBinContent(bin,graph->GetY()[i]);  
    }
    hist->SetStats(0); 
    return hist;
}
void *drawStack_78TeV(Option option, float analysis, TString file8tevname,  TString file7tevname, char* postfitfilename, 
                const unsigned int flav, const unsigned int njet, const char *name,  // name is like hww_mll
                const char *title, float lumi, float wwSF, char* region, bool splitwgamma, char* binsize) // title : x label
{

    bool doRescale = true;

    float post8_ZH      =   1.;    float post7_ZH      =   1.;
    float post8_WH      =   1.;    float post7_WH      =   1.;
    float post8_qqH     =   1.;    float post7_qqH     =   1.;
    float post8_ggH     =   1.;    float post7_ggH     =   1.;
    float post8_qqWW    =   1.;    float post7_qqWW    =   1.;
    float post8_ggWW    =   1.;    float post7_ggWW    =   1.;
    float post8_VV      =   1.;    float post7_VV      =   1.;
    float post8_Top     =   1.;    float post7_Top     =   1.;
    float post8_Zjets   =   1.;    float post7_Zjets   =   1.;
    float post8_WjetsE  =   1.;    float post7_WjetsE  =   1.;
    float post8_WjetsM  =   1.;    float post7_WjetsM  =   1.;
    float post8_Wgamma  =   1.;    float post7_Wgamma  =   1.;
    float post8_Wg3l    =   1.;    float post7_Wg3l    =   1.;
    float post8_Ztt     =   1.;    float post7_Ztt     =   1.; 

    TFile* postFile = TFile::Open(postfitfilename);
    TH1D *h1post8_ZH, *h1post8_WH, *h1post8_qqH, *h1post8_ggH, *h1post8_qqWW, *h1post8_ggWW, *h1post8_VV, *h1post8_Top, 
         *h1post8_Zjets, *h1post8_WjetsE, *h1post8_WjetsM, *h1post8_Wgamma, *h1post8_Wg3l, *h1post8_Ztt;  
    TH1D *h1post7_ZH, *h1post7_WH, *h1post7_qqH, *h1post7_ggH, *h1post7_qqWW, *h1post7_ggWW, *h1post7_VV, *h1post7_Top, 
         *h1post7_Zjets, *h1post7_WjetsE, *h1post7_WjetsM, *h1post7_Wgamma, *h1post7_Wg3l, *h1post7_Ztt;  
    h1post8_ZH     = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_ZH", njet)));      post8_ZH      =   h1post8_ZH->Integral();       
    h1post8_WH     = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_WH", njet)));      post8_WH      =   h1post8_WH->Integral();       
    h1post8_qqH    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_qqH", njet)));     post8_qqH     =   h1post8_qqH->Integral();      
    h1post8_ggH    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_ggH", njet)));     post8_ggH     =   h1post8_ggH->Integral();      
    h1post8_qqWW   = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_qqWW", njet)));    post8_qqWW    =   h1post8_qqWW->Integral();     
    h1post8_ggWW   = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_ggWW", njet)));    post8_ggWW    =   h1post8_ggWW->Integral();     
    h1post8_VV     = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_VV", njet)));      post8_VV      =   h1post8_VV->Integral();       
    h1post8_Top    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_Top", njet)));     post8_Top     =   h1post8_Top->Integral();      
    h1post8_Zjets  = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_Zjets", njet)));   post8_Zjets   =   h1post8_Zjets->Integral();    
    h1post8_WjetsE = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_WjetsE", njet)));  post8_WjetsE  =   h1post8_WjetsE->Integral();   
    h1post8_WjetsM = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_WjetsM", njet)));  post8_WjetsM  =   h1post8_WjetsM->Integral();  
    h1post8_Wgamma = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_Wgamma", njet)));  post8_Wgamma  =   h1post8_Wgamma->Integral();  
    h1post8_Wg3l   = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_Wg3l", njet)));    post8_Wg3l    =   h1post8_Wg3l->Integral();   
    h1post8_Ztt    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof8tev_Ztt", njet)));     post8_Ztt     =   h1post8_Ztt->Integral();      
    
    h1post7_ZH     = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_ZH", njet)));      post7_ZH      =   h1post7_ZH->Integral();       
    h1post7_WH     = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_WH", njet)));      post7_WH      =   h1post7_WH->Integral();       
    h1post7_qqH    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_qqH", njet)));     post7_qqH     =   h1post7_qqH->Integral();      
    h1post7_ggH    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_ggH", njet)));     post7_ggH     =   h1post7_ggH->Integral();      
    h1post7_qqWW   = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_qqWW", njet)));    post7_qqWW    =   h1post7_qqWW->Integral();     
    h1post7_ggWW   = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_ggWW", njet)));    post7_ggWW    =   h1post7_ggWW->Integral();     
    h1post7_VV     = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_VV", njet)));      post7_VV      =   h1post7_VV->Integral();       
    h1post7_Top    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_Top", njet)));     post7_Top     =   h1post7_Top->Integral();      
    h1post7_Zjets  = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_Zjets", njet)));   post7_Zjets   =   h1post7_Zjets->Integral();    
    h1post7_WjetsE = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_WjetsE", njet)));  post7_WjetsE  =   h1post7_WjetsE->Integral();   
    h1post7_WjetsM = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_WjetsM", njet)));  post7_WjetsM  =   h1post7_WjetsM->Integral();  
    h1post7_Wgamma = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_Wgamma", njet)));  post7_Wgamma  =   h1post7_Wgamma->Integral();  
    h1post7_Wg3l   = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_Wg3l", njet)));    post7_Wg3l    =   h1post7_Wg3l->Integral();   
    h1post7_Ztt    = makeHistogram( (TGraphAsymmErrors*) postFile->Get(Form("j%iof7tev_Ztt", njet)));     post7_Ztt     =   h1post7_Ztt->Integral();      


    if(1) {
        cout << "ZH     : " << post8_ZH    << " " << post7_ZH     <<endl;
        cout << "WH     : " << post8_WH    << " " << post7_WH     <<endl;
        cout << "qqH    : " << post8_qqH   << " " << post7_qqH    <<endl;
        cout << "ggH    : " << post8_ggH   << " " << post7_ggH    <<endl;
        cout << "qqWW   : " << post8_qqWW  << " " << post7_qqWW   <<endl;
        cout << "ggWW   : " << post8_ggWW  << " " << post7_ggWW   <<endl;
        cout << "VV     : " << post8_VV    << " " << post7_VV     <<endl;
        cout << "Top    : " << post8_Top   << " " << post7_Top    <<endl;
        cout << "Zjets  : " << post8_Zjets << " " << post7_Zjets  <<endl;
        cout << "WjetsE : " << post8_WjetsE<< " " << post7_WjetsE <<endl;
        cout << "WjetsM : " << post8_WjetsM<< " " << post7_WjetsM <<endl;
        cout << "Wgamma : " << post8_Wgamma<< " " << post7_Wgamma <<endl;
        cout << "Wg3l   : " << post8_Wg3l  << " " << post7_Wg3l   <<endl;
        cout << "Ztt    : " << post8_Ztt   << " " << post7_Ztt    <<endl; 
    }

    // ----------------------------------------- 
    //  pre-fit histograms
    // ----------------------------------------- 

    TFile *file8tev = TFile::Open(file8tevname);
    assert(file8tev);
    TFile *file7tev = TFile::Open(file7tevname);
    assert(file7tev);

    const char *fullname = Form("%s_%s_%s", name, jetbin_names[njet], types[flav]);
    
    TH1D *h1_ZH_8tev, *h1_WH_8tev, *h1_qqH_8tev, *h1_ggH_8tev, 
         *h1_qqWW_8tev, *h1_ggWW_8tev, *h1_VV_8tev, *h1_Top_8tev, 
         *h1_Zjets_8tev, *h1_WjetsE_8tev, *h1_WjetsM_8tev, 
         *h1_Wgamma_8tev, *h1_Wg3l_8tev, *h1_Ztt_8tev;
    TH1D *h1_ZH_7tev, *h1_WH_7tev, *h1_qqH_7tev, *h1_ggH_7tev, 
         *h1_qqWW_7tev, *h1_ggWW_7tev, *h1_VV_7tev, *h1_Top_7tev, 
         *h1_Zjets_7tev, *h1_WjetsE_7tev, *h1_WjetsM_7tev, 
         *h1_Wgamma_7tev, *h1_Wg3l_7tev, *h1_Ztt_7tev;
   
    cout << "fullname : " << fullname << endl; 
    
    //
    // get histrograms at 8 TeV
    //
    // signal 
    h1_ZH_8tev      = (TH1D*)file8tev->Get(Form("ZH_%s", fullname));    if (doRescale) h1_ZH_8tev->Scale( h1_ZH_8tev->Integral() ? post8_ZH/h1_ZH_8tev->Integral() : 1.);
    h1_WH_8tev      = (TH1D*)file8tev->Get(Form("WH_%s", fullname));    if (doRescale) h1_WH_8tev->Scale( h1_WH_8tev->Integral() ? post8_WH/h1_WH_8tev->Integral() : 1.);
    h1_qqH_8tev     = (TH1D*)file8tev->Get(Form("qqH_%s", fullname));   if (doRescale) h1_qqH_8tev->Scale( h1_qqH_8tev->Integral() ? post8_qqH/h1_qqH_8tev->Integral() : 1.);
    h1_ggH_8tev     = (TH1D*)file8tev->Get(Form("ggH_%s", fullname));   if (doRescale) h1_ggH_8tev->Scale( h1_ggH_8tev->Integral() ? post8_ggH/h1_ggH_8tev->Integral() : 1.);
    h1_signal_8tev  = (TH1D*) h1_ZH_8tev->Clone("h1_signal_8tev");
    h1_signal_8tev->Add(h1_WH_8tev);
    h1_signal_8tev->Add(h1_qqH_8tev);
    h1_signal_8tev->Add(h1_ggH_8tev);
    h1_signal_8tev->SetLineWidth(1);
    h1_signal_8tev->SetFillColor(kRed);
    h1_signal_8tev->SetLineColor(kRed);
   
    cout << "h1_ZH_8tev  : " << h1_ZH_8tev->Integral() << endl; 
    cout << "h1_WH_8tev  : " << h1_WH_8tev->Integral() << endl; 
    cout << "h1_qqH_8tev : " << h1_qqH_8tev->Integral() << endl; 
    cout << "h1_ggH_8tev : " << h1_ggH_8tev->Integral() << endl; 

    // WW
    h1_qqWW_8tev    = (TH1D*)file8tev->Get(Form("qqWW_%s", fullname));  if (doRescale) h1_qqWW_8tev->Scale( h1_qqWW_8tev->Integral() ? post8_qqWW/h1_qqWW_8tev->Integral() : 1. );
    h1_ggWW_8tev    = (TH1D*)file8tev->Get(Form("ggWW_%s", fullname));  if (doRescale) h1_ggWW_8tev->Scale( h1_ggWW_8tev->Integral() ? post8_ggWW/h1_ggWW_8tev->Integral() : 1. );
    h1_qqWW_8tev->Scale(wwSF);
    h1_ggWW_8tev->Scale(wwSF);
    TH1F *h1_WW_8tev = (TH1F*)h1_qqWW_8tev->Clone("h1_WW_8tev");
    h1_WW_8tev->Add(h1_ggWW_8tev);
    h1_WW_8tev->SetFillColor(kAzure-9);
    h1_WW_8tev->SetLineColor(kAzure-9);
   
    // VV + Wgamma(*)
    h1_VV_8tev          = (TH1D*)file8tev->Get(Form("VV_%s", fullname));    if (doRescale) h1_VV_8tev->Scale( h1_VV_8tev->Integral() ? post8_VV/h1_VV_8tev->Integral() : 1. );
    h1_Wgamma_8tev      = (TH1D*)file8tev->Get(Form("Wgamma_%s", fullname));    
    h1_Wgammanorm_8tev  = (TH1D*)file8tev->Get(Form("Wgammanorm_%s", fullname));   
    h1_Wgamma_8tev->Scale(h1_Wgammanorm_8tev->Integral(0,1000)/h1_Wgamma_8tev->Integral(0,1000)); // normalize to Wgammanorm 
                                                                                    if (doRescale) h1_Wgamma_8tev->Scale( h1_Wgamma_8tev->Integral() ? post8_Wgamma/h1_Wgamma_8tev->Integral() : 1. );
    h1_Wg3l_8tev        = (TH1D*)file8tev->Get(Form("Wg3l_%s", fullname));          if (doRescale) h1_Wg3l_8tev->Scale( h1_Wg3l_8tev->Integral() ? post8_Wg3l/h1_Wg3l_8tev->Integral() : 1. );        
    TH1F *h1_VVWgamma_8tev = (TH1F*)h1_VV_8tev->Clone("h1_VV_8tev");
    h1_VVWgamma_8tev->Add(h1_Wgamma_8tev);
    h1_VVWgamma_8tev->Add(h1_Wg3l_8tev);
    h1_VVWgamma_8tev->SetFillColor(kAzure-2);
    h1_VVWgamma_8tev->SetLineColor(kAzure-2);

    // Top 
    h1_Top_8tev     = (TH1D*)file8tev->Get(Form("Top_%s", fullname));           if (doRescale) h1_Top_8tev->Scale( h1_Top_8tev->Integral() ? post8_Top/h1_Top_8tev->Integral() : 1. );
    h1_Top_8tev->SetFillColor(kYellow);
    h1_Top_8tev->SetLineColor(kYellow);
   
    // Wjets
    h1_WjetsE_8tev  = (TH1D*)file8tev->Get(Form("WjetsE_%s", fullname));    if (doRescale) h1_WjetsE_8tev->Scale( h1_WjetsE_8tev->Integral() ? post8_WjetsE/h1_WjetsE_8tev->Integral() : 1. );
    h1_WjetsM_8tev  = (TH1D*)file8tev->Get(Form("WjetsM_%s", fullname));    if (doRescale) h1_WjetsM_8tev->Scale( h1_WjetsM_8tev->Integral() ? post8_WjetsM/h1_WjetsM_8tev->Integral() : 1. );
    TH1F *h1_Wjets_8tev = (TH1F*)h1_WjetsE_8tev->Clone("h1_Wjets_8tev");
    h1_Wjets_8tev->Add(h1_WjetsM_8tev);
    h1_Wjets_8tev->SetFillColor(kGray+1);
    h1_Wjets_8tev->SetLineColor(kGray+1);
  
    // Zjets
    h1_Ztt_8tev     = (TH1D*)file8tev->Get(Form("Ztt_%s", fullname));       if (doRescale) h1_Ztt_8tev->Scale(h1_Ztt_8tev->Integral() ?  post8_Ztt/h1_Ztt_8tev->Integral() : 1. );
    h1_Zjets_8tev    = (TH1D*)file8tev->Get(Form("Zjets_%s", fullname));    if (doRescale) h1_Zjets_8tev->Scale( h1_Zjets_8tev->Integral() ? post8_Zjets/h1_Zjets_8tev->Integral() : 1. );
    h1_Ztt_8tev->SetFillColor(kGreen+2);
    h1_Ztt_8tev->SetLineColor(kGreen+2);

    // Data
    h1_Data_8tev    = (TH1D*)file8tev->Get(Form("Data_%s", fullname));
    h1_Data_8tev->SetMarkerStyle(20);   
    h1_Data_8tev->SetLineWidth(2);
    
    //
    // get histrograms at 7 TeV
    //
    // signal 
    h1_ZH_7tev      = (TH1D*)file7tev->Get(Form("ZH_%s", fullname));    if (doRescale) h1_ZH_7tev->Scale( h1_ZH_7tev->Integral() ? post7_ZH/h1_ZH_7tev->Integral() : 1. );
    h1_WH_7tev      = (TH1D*)file7tev->Get(Form("WH_%s", fullname));    if (doRescale) h1_WH_7tev->Scale( h1_WH_7tev->Integral() ? post7_WH/h1_WH_7tev->Integral() : 1. );
    h1_qqH_7tev     = (TH1D*)file7tev->Get(Form("qqH_%s", fullname));   if (doRescale) h1_qqH_7tev->Scale( h1_qqH_7tev->Integral() ? post7_qqH/h1_qqH_7tev->Integral() : 1. );
    h1_ggH_7tev     = (TH1D*)file7tev->Get(Form("ggH_%s", fullname));   if (doRescale) h1_ggH_7tev->Scale( h1_ggH_7tev->Integral() ? post7_ggH/h1_ggH_7tev->Integral() : 1. );
    h1_signal_7tev  = (TH1D*) h1_ZH_7tev->Clone("h1_signal_7tev");
    h1_signal_7tev->Add(h1_WH_7tev);
    h1_signal_7tev->Add(h1_qqH_7tev);
    h1_signal_7tev->Add(h1_ggH_7tev);
    h1_signal_7tev->SetLineWidth(1);
    h1_signal_7tev->SetFillColor(kRed);
    h1_signal_7tev->SetLineColor(kRed);
    
    // WW
    h1_qqWW_7tev    = (TH1D*)file7tev->Get(Form("qqWW_%s", fullname));  if (doRescale) h1_qqWW_7tev->Scale( h1_qqWW_7tev->Integral() ? post7_qqWW/h1_qqWW_7tev->Integral() : 1. );
    h1_ggWW_7tev    = (TH1D*)file7tev->Get(Form("ggWW_%s", fullname));  if (doRescale) h1_ggWW_7tev->Scale( h1_ggWW_7tev->Integral() ? post7_ggWW/h1_ggWW_7tev->Integral() : 1. );
    h1_qqWW_7tev->Scale(wwSF);
    h1_ggWW_7tev->Scale(wwSF);
    TH1F *h1_WW_7tev = (TH1F*)h1_qqWW_7tev->Clone("h1_WW_7tev");
    h1_WW_7tev->Add(h1_ggWW_7tev);
    h1_WW_7tev->SetFillColor(kAzure-9);
    h1_WW_7tev->SetLineColor(kAzure-9);
   
    // VV + Wgamma(*)
    h1_VV_7tev          = (TH1D*)file7tev->Get(Form("VV_%s", fullname));        if (doRescale) h1_VV_7tev->Scale( h1_VV_7tev->Integral() ? post7_VV/h1_VV_7tev->Integral() : 1. );
    h1_Wgamma_7tev      = (TH1D*)file7tev->Get(Form("Wgamma_%s", fullname));    
    h1_Wgammanorm_7tev  = (TH1D*)file7tev->Get(Form("Wgammanorm_%s", fullname)); 
    h1_Wgamma_7tev->Scale(h1_Wgammanorm_7tev->Integral(0,1000)/h1_Wgamma_7tev->Integral(0,1000)); // normalize to Wgammanorm
                                                                                if (doRescale) h1_Wgamma_7tev->Scale( h1_Wgamma_7tev->Integral() ? post7_Wgamma/h1_Wgamma_7tev->Integral() : 1.); 
    h1_Wg3l_7tev        = (TH1D*)file7tev->Get(Form("Wg3l_%s", fullname));      if (doRescale) h1_Wg3l_7tev->Scale( h1_Wg3l_7tev->Integral() ? post7_Wg3l/h1_Wg3l_7tev->Integral() : 1. );
    TH1F *h1_VVWgamma_7tev = (TH1F*)h1_VV_7tev->Clone("h1_VV_7tev");
    h1_VVWgamma_7tev->Add(h1_Wgamma_7tev);
    h1_VVWgamma_7tev->Add(h1_Wg3l_7tev);
    h1_VVWgamma_7tev->SetFillColor(kAzure-2);
    h1_VVWgamma_7tev->SetLineColor(kAzure-2);

    // Top 
    h1_Top_7tev     = (TH1D*)file7tev->Get(Form("Top_%s", fullname));   if (doRescale) h1_Top_7tev->Scale( h1_Top_7tev->Integral() ? post7_Top/h1_Top_7tev->Integral() : 1. );
    h1_Top_7tev->SetFillColor(kYellow);
    h1_Top_7tev->SetLineColor(kYellow);
   
    // Wjets
    h1_WjetsE_7tev  = (TH1D*)file7tev->Get(Form("WjetsE_%s", fullname));    if (doRescale) h1_WjetsE_7tev->Scale( h1_WjetsE_7tev->Integral() ? post7_WjetsE/h1_WjetsE_7tev->Integral() : 1. );
    h1_WjetsM_7tev  = (TH1D*)file7tev->Get(Form("WjetsM_%s", fullname));    if (doRescale) h1_WjetsM_7tev->Scale( h1_WjetsM_7tev->Integral() ? post7_WjetsM/h1_WjetsM_7tev->Integral() : 1. );
    TH1F *h1_Wjets_7tev = (TH1F*)h1_WjetsE_7tev->Clone("h1_Wjets_7tev");
    h1_Wjets_7tev->Add(h1_WjetsM_7tev);
    h1_Wjets_7tev->SetFillColor(kGray+1);
    h1_Wjets_7tev->SetLineColor(kGray+1);
  
    // Zjets
    h1_Ztt_7tev     = (TH1D*)file7tev->Get(Form("Ztt_%s", fullname));       if (doRescale) h1_Ztt_7tev->Scale( h1_Ztt_7tev->Integral() ? post7_Ztt/h1_Ztt_7tev->Integral() : 1. );
    h1_Zjets_7tev    = (TH1D*)file7tev->Get(Form("Zjets_%s", fullname));    if (doRescale) h1_Zjets_7tev->Scale( h1_Zjets_7tev->Integral() ? post7_Zjets/h1_Zjets_7tev->Integral() : 1.);
    h1_Ztt_7tev->SetFillColor(kGreen+2);
    h1_Ztt_7tev->SetLineColor(kGreen+2);

    // Data
    h1_Data_7tev    = (TH1D*)file7tev->Get(Form("Data_%s", fullname));
    h1_Data_7tev->SetMarkerStyle(20);   
    h1_Data_7tev->SetLineWidth(2);


    // Combine 7 and 8 TeV 
    h1_ZH           = (TH1D*)h1_ZH_8tev->Clone("h1_ZH");                h1_ZH->Add(h1_ZH_7tev);
    h1_WH           = (TH1D*)h1_WH_8tev->Clone("h1_WH");                h1_WH->Add(h1_WH_7tev);
    h1_qqH          = (TH1D*)h1_qqH_8tev->Clone("h1_qqH");              h1_qqH->Add(h1_qqH_7tev);
    h1_ggH          = (TH1D*)h1_ggH_8tev->Clone("h1_ggH");              h1_ggH->Add(h1_ggH_7tev);
    h1_signal       = (TH1D*)h1_signal_8tev->Clone("h1_signal");        h1_signal->Add(h1_signal_7tev);
    h1_qqWW         = (TH1D*)h1_qqWW_8tev->Clone("h1_qqWW");            h1_qqWW->Add(h1_qqWW_7tev);
    h1_ggWW         = (TH1D*)h1_ggWW_8tev->Clone("h1_ggWW");            h1_ggWW->Add(h1_ggWW_7tev);
    h1_WW           = (TH1D*)h1_WW_8tev->Clone("h1_WW");                h1_WW->Add(h1_WW_7tev);
    h1_VV           = (TH1D*)h1_VV_8tev->Clone("h1_VV");                h1_VV->Add(h1_VV_7tev);         h1_VV->SetFillColor(kAzure-2);   h1_VV->SetLineColor(kAzure-2);
    h1_Wgamma       = (TH1D*)h1_Wgamma_8tev->Clone("h1_Wgamma");        h1_Wgamma->Add(h1_Wgamma_7tev); h1_Wgamma->SetFillColor(kOrange);   h1_Wgamma->SetLineColor(kOrange);
    h1_Wg3l         = (TH1D*)h1_Wg3l_8tev->Clone("h1_Wg3l");            h1_Wg3l->Add(h1_Wg3l_7tev);     h1_Wg3l->SetFillColor(kOrange+1);   h1_Wg3l->SetLineColor(kOrange+1);
    h1_VVWgamma     = (TH1D*)h1_VVWgamma_8tev->Clone("h1_VVWgamma");    h1_VVWgamma->Add(h1_VVWgamma_7tev);
    h1_Top          = (TH1D*)h1_Top_8tev->Clone("h1_Top");              h1_Top->Add(h1_Top_7tev);
    h1_WjetsE       = (TH1D*)h1_WjetsE_8tev->Clone("h1_WjetsE");        h1_WjetsE->Add(h1_WjetsE_7tev);
    h1_WjetsM       = (TH1D*)h1_WjetsM_8tev->Clone("h1_WjetsM");        h1_WjetsM->Add(h1_WjetsM_7tev);
    h1_Wjets        = (TH1D*)h1_Wjets_8tev->Clone("h1_Wjets");          h1_Wjets->Add(h1_Wjets_7tev);
    h1_Ztt          = (TH1D*)h1_Ztt_8tev->Clone("h1_Ztt");              h1_Ztt->Add(h1_Ztt_7tev);
    h1_Data         = (TH1D*)h1_Data_8tev->Clone("h1_Data");            h1_Data->Add(h1_Data_7tev);
    
    h1_signal_overlay = (TH1D*)h1_signal->Clone("h1_signal_overlay");
    h1_signal_overlay->SetFillColor(0);

    // note - assuming the sf uncertainty for drell-yan in the inclusive plot
    // - this is not really correct, but will make a negligible difference to
    // how the uncertainty band on the plot looks.  Really should construct
    // the inclusive plot uncertainties from the sum in quadrature of the uncertainties
    // on the exclusive channels.  But that's a bit fiddly so neglecting it for now.
    std::string flavor = "of";
  
    //
    // figure out the uncertainties
    //
    TH1F *h1_band = (TH1F*)h1_Data_8tev->Clone(Form("%s_up", fullname));
    TH1F *h1_band_rel = (TH1F*)h1_Data_8tev->Clone(Form("%s_up", fullname));
    TH1F *h1_ratio = (TH1F*)h1_Data_8tev->Clone(Form("%s_up", fullname));
    for (unsigned int i = 1; i <= h1_Data_8tev->GetNbinsX(); ++i)
    {
        float syst_qqWW_8tev    = GetBackgroundEstimationSystematic_8TeV(option, QQWW, njet, flavor) * h1_qqWW_8tev->GetBinContent(i);
        float syst_ggWW_8tev    = GetBackgroundEstimationSystematic_8TeV(option, GGWW, njet, flavor) * h1_ggWW_8tev->GetBinContent(i);
        float syst_VV_8tev      = GetBackgroundEstimationSystematic_8TeV(option, VV, njet, flavor) * h1_VV_8tev->GetBinContent(i);
        float syst_Top_8tev     = GetBackgroundEstimationSystematic_8TeV(option, TOP, njet, flavor) * h1_Top_8tev->GetBinContent(i);
        float syst_Wjets_8tev   = GetBackgroundEstimationSystematic_8TeV(option, WJETSDATA, njet, flavor) * h1_Wjets_8tev->GetBinContent(i);
        float syst_Wgamma_8tev  = GetBackgroundEstimationSystematic_8TeV(option, WGAMMA, njet, flavor) * h1_Wgamma_8tev->GetBinContent(i);
        float syst_Wg3l_8tev    = GetBackgroundEstimationSystematic_8TeV(option, WGAMMASTAR, njet, flavor) * h1_Wg3l_8tev->GetBinContent(i);
        float syst_Ztt_8tev     = GetBackgroundEstimationSystematic_8TeV(option, ZTT, njet, flavor) * h1_Ztt_8tev->GetBinContent(i);
        
        float syst_qqWW_7tev    = GetBackgroundEstimationSystematic_7TeV(option, QQWW, njet, flavor) * h1_qqWW_7tev->GetBinContent(i);
        float syst_ggWW_7tev    = GetBackgroundEstimationSystematic_7TeV(option, GGWW, njet, flavor) * h1_ggWW_7tev->GetBinContent(i);
        float syst_VV_7tev      = GetBackgroundEstimationSystematic_7TeV(option, VV, njet, flavor) * h1_VV_7tev->GetBinContent(i);
        float syst_Top_7tev     = GetBackgroundEstimationSystematic_7TeV(option, TOP, njet, flavor) * h1_Top_7tev->GetBinContent(i);
        float syst_Wjets_7tev   = GetBackgroundEstimationSystematic_7TeV(option, WJETSDATA, njet, flavor) * h1_Wjets_7tev->GetBinContent(i);
        float syst_Wgamma_7tev  = GetBackgroundEstimationSystematic_7TeV(option, WGAMMA, njet, flavor) * h1_Wgamma_7tev->GetBinContent(i);
        float syst_Wg3l_7tev    = GetBackgroundEstimationSystematic_7TeV(option, WGAMMASTAR, njet, flavor) * h1_Wg3l_7tev->GetBinContent(i);
        float syst_Ztt_7tev     = GetBackgroundEstimationSystematic_7TeV(option, ZTT, njet, flavor) * h1_Ztt_7tev->GetBinContent(i);

        if(i==3) {
            cout << " syst_qqWW_8tev   : " << syst_qqWW_8tev    << endl; 
            cout << " syst_ggWW_8tev   : " << syst_ggWW_8tev    << endl;
            cout << " syst_VV_8tev     : " << syst_VV_8tev      << endl; 
            cout << " syst_Top_8tev    : " << syst_Top_8tev     << endl; 
            cout << " syst_Wjets_8tev  : " << syst_Wjets_8tev   << endl; 
            cout << " syst_Wgamma_8tev : " << syst_Wgamma_8tev  << endl; 
            cout << " syst_Wg3l_8tev   : " << syst_Wg3l_8tev    << endl; 
            cout << " syst_Ztt_8tev    : " << syst_Ztt_8tev     << endl; 
       
//            cout << " GetBackgroundEstimationSystematic(option, qqWW, njet, flavor)  : " << GetBackgroundEstimationSystematic(option, QQWW, njet, flavor)  << endl; 
//            cout << " GetBackgroundEstimationSystematic(option, GGWW, njet, flavor)  : " << GetBackgroundEstimationSystematic(option, GGWW, njet, flavor) << endl; 
//            cout << " GetBackgroundEstimationSystematic(option, VV, njet, flavor)    : " << GetBackgroundEstimationSystematic(option, VV, njet, flavor)  << endl; 
//            cout << " GetBackgroundEstimationSystematic(option, TOP, njet, flavor)  : " << GetBackgroundEstimationSystematic(option, TOP, njet, flavor)  << endl; 
//            cout << " GetBackgroundEstimationSystematic(option, WJETSDATA, njet, flavor)  : " << GetBackgroundEstimationSystematic(option, WJETSDATA, njet, flavor) << endl; 
//            cout << " GetBackgroundEstimationSystematic(option, WGAMMA, njet, flavor)  : " << GetBackgroundEstimationSystematic(option, WGAMMA, njet, flavor) << endl; 
//            cout << " GetBackgroundEstimationSystematic(option, WGAMMASTAR, njet, flavor)  : " << GetBackgroundEstimationSystematic(option, WGAMMASTAR, njet, flavor) << endl; 
//            cout << " GetBackgroundEstimationSystematic(option, ZTT, njet, flavor)  : " << GetBackgroundEstimationSystematic(option, ZTT, njet, flavor) << endl; 
//            cout << " h1_Top_8tev->GetBinContent(i) : " << h1_Top_8tev->GetBinContent(i)  << endl; 
        }

        float stat_qqWW_8tev     = h1_qqWW_8tev->GetBinError(i);
        float stat_ggWW_8tev     = h1_ggWW_8tev->GetBinError(i);
        float stat_VV_8tev       = h1_VV_8tev->GetBinError(i);
        float stat_Top_8tev      = h1_Top_8tev->GetBinError(i);
        float stat_Wjets_8tev    = h1_Wjets_8tev->GetBinError(i);
        float stat_Wgamma_8tev   = h1_Wgamma_8tev->GetBinError(i);
        float stat_Wg3l_8tev     = h1_Wg3l_8tev->GetBinError(i);
        float stat_Ztt_8tev      = h1_Ztt_8tev->GetBinError(i);
        
        float stat_qqWW_7tev     = h1_qqWW_7tev->GetBinError(i);
        float stat_ggWW_7tev     = h1_ggWW_7tev->GetBinError(i);
        float stat_VV_7tev       = h1_VV_7tev->GetBinError(i);
        float stat_Top_7tev      = h1_Top_7tev->GetBinError(i);
        float stat_Wjets_7tev    = h1_Wjets_7tev->GetBinError(i);
        float stat_Wgamma_7tev   = h1_Wgamma_7tev->GetBinError(i);
        float stat_Wg3l_7tev     = h1_Wg3l_7tev->GetBinError(i);
        float stat_Ztt_7tev      = h1_Ztt_7tev->GetBinError(i);
        
        float yield_total    = h1_VVWgamma_8tev->GetBinContent(i);
        yield_total         += h1_Top_8tev->GetBinContent(i);
        yield_total         += h1_Wjets_8tev->GetBinContent(i);
        yield_total         += h1_Ztt_8tev->GetBinContent(i);
        yield_total         += h1_WW_8tev->GetBinContent(i);
        yield_total         += h1_VVWgamma_7tev->GetBinContent(i);
        yield_total         += h1_Top_7tev->GetBinContent(i);
        yield_total         += h1_Wjets_7tev->GetBinContent(i);
        yield_total         += h1_Ztt_7tev->GetBinContent(i);
        yield_total         += h1_WW_7tev->GetBinContent(i);
        
        float syst_total2 = pow(syst_qqWW_8tev, 2) + pow(syst_ggWW_8tev, 2) + pow(syst_VV_8tev, 2) 
                          + pow(syst_Top_8tev, 2) + pow(syst_Wjets_8tev, 2) + pow(syst_Wgamma_8tev, 2) 
                          + pow(syst_Wg3l_8tev, 2) + pow(syst_Ztt_8tev, 2);
                          + pow(syst_qqWW_7tev, 2) + pow(syst_ggWW_7tev, 2) + pow(syst_VV_7tev, 2) 
                          + pow(syst_Top_7tev, 2) + pow(syst_Wjets_7tev, 2) + pow(syst_Wgamma_7tev, 2) 
                          + pow(syst_Wg3l_7tev, 2) + pow(syst_Ztt_7tev, 2);
        float stat_total2 = pow(stat_qqWW_8tev, 2) + pow(stat_ggWW_8tev, 2) + pow(stat_VV_8tev, 2) 
                          + pow(stat_Top_8tev, 2) + pow(stat_Wjets_8tev, 2) + pow(stat_Wgamma_8tev, 2) 
                          + pow(stat_Wg3l_8tev, 2) + pow(stat_Ztt_8tev, 2);
                          + pow(stat_qqWW_7tev, 2) + pow(stat_ggWW_7tev, 2) + pow(stat_VV_7tev, 2) 
                          + pow(stat_Top_7tev, 2) + pow(stat_Wjets_7tev, 2) + pow(stat_Wgamma_7tev, 2) 
                          + pow(stat_Wg3l_7tev, 2) + pow(stat_Ztt_7tev, 2);
        float err_total = sqrt(syst_total2 + stat_total2);
        
        h1_band->SetBinContent(i, yield_total);
        h1_band->SetBinError(i, err_total);
        h1_band_rel->SetBinContent(i, 1.0);
        if (yield_total > 0.0) {
            h1_band_rel->SetBinError(i, err_total/yield_total);
            h1_ratio->SetBinContent(i, h1_Data->GetBinContent(i) / yield_total);
            h1_ratio->SetBinError(i, h1_Data->GetBinError(i) / yield_total);
        } else {
            h1_band_rel->SetBinError(i, 0.0);
            h1_ratio->SetBinContent(i, 0.0);
            h1_ratio->SetBinError(i, 0.0);
        }
    }
    
    //
    // make the legend
    //

    //TLegend *l1 = new TLegend(0.65, 0.55, 0.85, 0.85);
    TLegend* l1 = new TLegend (.16,.88-6*.035,.5,.88);
    l1->SetNColumns(2);
    l1->SetBorderSize(0);
    l1->SetFillColor(0);
    l1->SetFillStyle(0);
    l1->SetTextFont(42); 
    l1->SetTextAlign(12);
    l1->SetTextSize(0.035);
    
    l1->AddEntry(h1_Data,       "data", "lp");
    //l1->AddEntry(h1_signal_overlay, Form("m_{H} = %u GeV", int(analysis)), "l");
    //l1->AddEntry(h1_signal, Form("H%u", int(analysis)), "f");
    l1->AddEntry(h1_Wjets,      "W+jets", "f");
    if(splitwgamma) { 
        l1->AddEntry(h1_VV,         "WZ/ZZ", "f");
        l1->AddEntry(h1_Wgamma,     "W#gamma", "f");
        l1->AddEntry(h1_Wg3l,       "W#gamma*", "f"); 
    } else { 
        l1->AddEntry(h1_VVWgamma,         "VV", "f");
        l1->AddEntry(h1_Top,        "Top", "f");
        l1->AddEntry(h1_Ztt,        "Z/#gamma*", "f"); 
    }
    l1->AddEntry(h1_WW,         "WW", "f");
    l1->AddEntry(h1_band , TString(" stat. #oplus syst.") , "f");  
    
    //
    // make the stack
    //
    
    THStack *st = new THStack(Form("st_%s", fullname), Form("Stack;%s;events / %s", title, binsize));
    st->SetName(fullname);
    st->Add(h1_WW);
    st->Add(h1_Ztt);
    st->Add(h1_Top);
    if(!splitwgamma) st->Add(h1_VVWgamma); 
    if(splitwgamma) { st->Add(h1_VV);  st->Add(h1_Wgamma); st->Add(h1_Wg3l);} 
    st->Add(h1_Wjets); 
    //st->Add(h1_signal); 
    
    //
    // make the canvas
    //

    TCanvas *c1 = new TCanvas(fullname,fullname,600,600*1.2);
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
    st->GetXaxis()->SetTitleSize(0);
    st->SetTitle("");
    //h1_signal_overlay->Draw("SAME HIST");
    h1_Data->Draw("SAME E X0");
    h1_band->SetMarkerSize(0);
    h1_band->SetFillColor(kBlack);
    h1_band->SetLineColor(kWhite);
    h1_band->SetFillStyle(3005);
    h1_band->Draw("e2 same");
    l1->Draw();
    float textSize = 0.03;
    TLatex *tex = new TLatex(0.6,0.78,Form("#sqrt{s}=8 TeV, L = %.1f fb^{-1}", lumi));
    tex->SetNDC();
    tex->SetTextSize(textSize);
    tex->SetLineWidth(2);
    TLatex *tex1 = new TLatex(0.6,0.73,Form("#sqrt{s}=7 TeV, L = %.1f fb^{-1}", 4.9));
    tex1->SetNDC();
    tex1->SetTextSize(textSize);
    tex1->SetLineWidth(2);
    TLatex *tex2 = new TLatex(0.6,0.84,"CMS Preliminary");
    tex2->SetNDC();
    tex2->SetTextSize(textSize+0.02);
    tex2->SetLineWidth(2);
    TLatex *tex3 = new TLatex(0.2,0.77, Form("%s (%s)", types[flav], jetbin_names[njet]));
    tex3->SetNDC();
    tex3->SetTextSize(textSize);
    tex3->SetLineWidth(2);
    tex->Draw("SAME");
    tex1->Draw("SAME");
    tex2->Draw("SAME");
    //tex3->Draw("SAME");
    if (drawRatio) st->GetXaxis()->SetLabelOffset(1.0);
    if (drawRatio) st->GetYaxis()->SetTitleOffset(1.4);
    //st->SetMaximum(TMath::Max(h1_band->GetMaximum()*1.7, h1_Data->GetMaximum() + 3*sqrt(h1_Data->GetMaximum())));
    st->SetMaximum(h1_band->GetMaximum()*1.7);
    if (false) {
        if (drawRatio) pad1->SetLogy(1);
        st->SetMinimum(0.1);
        st->SetMaximum(h1_Data->GetMaximum() * 500);
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
        h1_band_rel->SetStats(0);
        h1_band_rel->SetTitle("");
        h1_band_rel->Draw("E2");
        h1_ratio->Draw("E SAME X0");
        h1_band_rel->GetYaxis()->SetRangeUser(0.0, 2.5);
        h1_band_rel->SetMarkerSize(0);
        h1_band_rel->SetFillColor(kBlack);
        h1_band_rel->SetLineColor(kBlack);
        h1_band_rel->SetFillStyle(3005);
        h1_band_rel->GetXaxis()->SetTitle(title);
        h1_band_rel->GetXaxis()->SetLabelSize(0.08);
        h1_band_rel->GetXaxis()->SetLabelOffset(0.03);
        h1_band_rel->GetXaxis()->SetTitleOffset(1.1);
        h1_band_rel->GetXaxis()->SetTitleSize(0.13);
        h1_band_rel->GetXaxis()->SetNdivisions(505);
    
        h1_band_rel->GetYaxis()->SetTitle("data / MC");
        h1_band_rel->GetYaxis()->SetLabelSize(0.08);
        h1_band_rel->GetYaxis()->SetTitleSize(0.08);
        h1_band_rel->GetYaxis()->SetTitleOffset(0.6);
        h1_band_rel->GetYaxis()->SetNdivisions(405);
    
        float min = h1_band_rel->GetXaxis()->GetXmin();
        float max = h1_band_rel->GetXaxis()->GetXmax();
        TLine *line = new TLine(min, 1.0, max, 1.0);
        line->Draw();
    }
    
    c1->RedrawAxis(); 
    TString var = name; var.ReplaceAll("ww_","");
    c1->Print(Form("CR_plot_2D_mH125/%s_%s_%ij.pdf", region, var.Data(), njet));
    c1->Print(Form("CR_plot_2D_mH125/%s_%s_%ij.eps", region, var.Data(), njet));
    c1->Print(Form("CR_plot_2D_mH125/%s_%s_%ij.png", region, var.Data(), njet));
    c1->Print(Form("CR_plot_2D_mH125/%s_%s_%ij.root", region, var.Data(), njet));


}


