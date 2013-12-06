#ifndef ENUMS_H
#define ENUMS_H

typedef unsigned long int Samples_t;
typedef unsigned long int Cuts_t;
typedef unsigned long int Opt_t;
typedef unsigned long int ShapeVar_t;

//
// analysis selections
//

enum Selections {

    //
    // HZZ
    //

    HZZ_PASS_2011A,

    // try and migrate to these
    HZZ_PASS_BASELINE,
    HZZ_PASS_ISZ,
    HZZ_PASS_ISNOTZ,
    HZZ_PASS_SOFTMUVETO,
    HZZ_PASS_BVETO,
    HZZ_PASS_BTAG,
    HZZ_PASS_SF,
    HZZ_PASS_OF,
    HZZ_PASS_2011ADPHI,
    HZZ_PASS_2011AMT,
    HZZ_PASS_2011AMET,
    HZZ_PASS_2011APT,

    // preselection for mva selections
    // and shape based analysis
    HZZ_PASS_MVAPRESEL,

    // mva cut based selections
    HZZ_PASS_MVASEL,
    HZZ_PASS_MESEL,

    //
    // HWW
    //

    HWW_PASS_PRESEL,
    HWW_PASS_CUTSEL,
    HWW_PASS_MVASEL,
    HWW_PASS_JCPSEL, 
    HWW_PASS_SSCTL,
    HWW_PASS_SSCTL2D,

};

//
// what analysis to do?
// these control options such as
// what scale factors need to be loaded
//

enum Option {

    //
    // HZZ
    //

    HZZ_OPT_SMURF=1,
    HZZ_OPT_SMURFPRESEL=2,
    HZZ_OPT_EPS=3,
    HZZ_OPT_SMURF42X=4,
    HZZ_OPT_SMURF42XPRESEL=5,
    HZZ_OPT_EPS42X=6,
    HZZ_OPT_EPS_LP=7,
    HZZ_OPT_EPS_POSTEPS=8,
    HZZ_OPT_2011A = 9,

    // run MVA or ME analyses
    HZZ_OPT_EPS_MVASEL=10,
    HZZ_OPT_EPS_MESEL=11,
    HZZ_OPT_EPS_MTSEL=12,
    HZZ_OPT_2011A_MTSEL = 13,
    HZZ_OPT_2011A_MESEL = 14,


    //
    // HWW
    //

    HWW_OPT_SMURFPRESEL=16,
    HWW_OPT_SMURFCUTSEL=17,
    HWW_OPT_SMURFMVASEL=18,
    HWW_OPT_SMURFMESEL=19,
    WW_OPT_SMURFXSECSEL=20,
    HWW_OPT_MT2DMLL=21,
    HWW_OPT_MT2DMLL_JCP=22,
    XWW_OPT_MT2DMLL_JCP=23,
    HWW_OPT_SSCTL=24,
    HWW_OPT_SSCTL2D=25,
    
    HWW_OPT_TOPTAG=26,

};

const Opt_t HWW_SHAPE   =  (1ll<<HWW_OPT_SMURFMVASEL) | (1ll<<HWW_OPT_SMURFMESEL) | (1ll<<HWW_OPT_MT2DMLL) | (1ll<<HWW_OPT_MT2DMLL_JCP) | (1ll<<XWW_OPT_MT2DMLL_JCP) | (1ll<<HWW_OPT_SSCTL2D); 
const Opt_t HWW_MT2DMLL =  (1ll<<HWW_OPT_MT2DMLL) | (1ll<<HWW_OPT_MT2DMLL_JCP) | (1ll<<XWW_OPT_MT2DMLL_JCP) | (1ll<<HWW_OPT_SSCTL2D); 



//
// input data samples
//

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

//
// lepton types
//

// ee, em, me, mm, sf, of, total
enum FlavorType {
    fMM,
    fME,
    fEM,
    fEE,
    fSF,
    fOF,
    fTOTAL,
};

//
// run eras
//

enum RunEra {
    RUN2011A=0,
    RUN2011B=1,
    RUN2011AB=2,
    RUN2012=3,
    RUN2012HCP=4,
    RUN2012Moriond=5,
};  

// DATA
const Samples_t data_data =    (1ll<<DATA) | (1ll<<GAMMA) | (1ll<<WJETSDATA) | (1ll<<TOPDATA) | (1ll<<OFDATA) | (1ll<<WJETSELEDATA)  | (1ll<<WJETSMUDATA) ;

// HIGGS SIGNALS
const Samples_t data_higgsww = (1ll<<GGHWW) | (1ll<<QQHWW) | (1ll<<ZHWW) | (1ll<<WHWW) | (1ll<<GGHWWREF) | (1ll<<GGHWWJHU);
const Samples_t data_gghiggs = (1ll<<GGHWW) | (1ll<<GGHWWREF) | (1ll<<GGHWWJHU);
const Samples_t data_qqhiggs = (1ll<<QQHWW);
const Samples_t data_higgs = data_higgsww;

// MC

// processes to apply the OF scaling to
// if it is applied for the analysis being done
const Samples_t data_ofscale =   (1ll<<QQWW) | (1ll<<GGWW) | (1ll<<WJETS) | (1ll<<WW) | (1ll<<TOP) | (1ll<<ZTT);

// processes not to apply the OF scaling to
const Samples_t data_nonofscale = (1ll<<ZLL) | (1ll<<ZZ) | (1ll<<GGZZ)| (1ll<<WZ);

// all MC BG processes
const Samples_t data_allmcbg = (1ll<<QQWW) | (1ll<<GGWW) | (1ll<<WJETS) | (1ll<<WW) | (1ll<<TOP) | (1ll<<ZTT) | (1ll<<ZLL) | (1ll<<ZZ) | (1ll<<GGZZ) | (1ll<<WZ) |  (1ll<<VV)  | 
  (1ll<<WGAMMA) | (1ll<<WGAMMASTAR)| (1ll<<ZJETS) | (1ll<<ZLLLOOSEMET)  |  (1ll<<WWMCNLO) | (1ll<<WWMCNLOUP) | (1ll<<WWMCNLODOWN) | (1ll<<WZALTER) | (1ll<<ZZALTER) | (1ll<<TOPALTER) | (1ll<<WGAMMALOOSE) | (1ll<<WJETSMCLOOSE) | (1ll<<ZVVGAMMA) | (1ll<<ZLLGAMMA) | (1ll<<WJETSGAMMA)| (1ll<<WG3L) | (1ll<<WGAMMANORM) | (1ll<<DPSWW);

// all MC processes
const Samples_t data_allmc = data_allmcbg | data_higgs;

// BACKGROUNDS
const Samples_t data_allbg = data_allmcbg | (1ll<<GAMMA) | (1ll<<WJETSDATA) | (1ll<<TOPDATA) | (1ll<<WJETSELEDATA) | (1ll<<WJETSMUDATA);


enum MVAShapeSyst {
  NOVAR=1,
  STATVAR=2,
  DYSHAPEVAR=3,
  QCDSCALEVAR=4,
  LEPEFFVAR=5,
  LEPRESVAR=6,
  METVAR=7,
  TOPSHAPEVAR=8,
  WWSHAPEVAR=9,
  WZSHAPEVAR=10,
  ZZSHAPEVAR=11,
  WJETSSHAPEVAR=12,
  WWTOPSHAPEVAR=13,
  JETRESVAR=14,
  WJETSELESHAPEVAR=15,
  WJETSMUSHAPEVAR=16,
  PDFSHAPEVAR=17,
};

const ShapeVar_t mva_var = (1ll<<STATVAR) | (1ll<<DYSHAPEVAR) | (1ll<<QCDSCALEVAR) | (1ll<<LEPEFFVAR) | (1ll<<LEPRESVAR) | (1ll<<METVAR) | (1ll<<TOPSHAPEVAR)| (1ll<<WWSHAPEVAR) | (1ll<<WZSHAPEVAR) | (1ll<<ZZSHAPEVAR) | (1ll<<WJETSSHAPEVAR) | (1ll<<WWTOPSHAPEVAR) | (1ll<<JETRESVAR) | (1ll<<WJETSELESHAPEVAR) | (1ll<<WJETSMUSHAPEVAR) | (1ll<<PDFSHAPEVAR); 

#endif

