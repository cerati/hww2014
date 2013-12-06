
#include "../core/Enums.h"

void doPDF() {

    //
    // the looper
    //

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gROOT->ProcessLine(".L ../../../../../Smurf/Core/SmurfTree.h+");
    gSystem->Load("/tas/dlevans/HWW2012/CMSSW_5_2_3/src/LHAPDF-5.8.92b/lib/libLHAPDF.so");
    gSystem->Load("libSmurfPDFLooper.so");

    //
    // run it
    //

    // low mass
    //doMassPoint(125.0, HWW_OPT_MT2DMLL_JCP);
    // high mass
    doMassPoint(500.0, HWW_OPT_MT2DMLL);

}


void doMassPoint(float analysis, Option option)
{

    //
    // create looper
    //

    MyScanChain *looper_cteq6ll = new MyScanChain(analysis, option, "cteq6ll.LHpdf", 0);
    //MyScanChain *looper_CT10 = new MyScanChain(analysis, option, "CT10.LHgrid", 5);

    //
    // run all pdf sets
    //

    std::vector<std::string> pdfSets;
    // CT10
    pdfSets.push_back("CT10");
    pdfSets.push_back("CT10as");
    // MSTW
    pdfSets.push_back("MSTW2008nlo68cl");
    pdfSets.push_back("MSTW2008nlo68cl_asmz+68cl");
    pdfSets.push_back("MSTW2008nlo68cl_asmz+68clhalf");
    pdfSets.push_back("MSTW2008nlo68cl_asmz-68cl");
    pdfSets.push_back("MSTW2008nlo68cl_asmz-68clhalf");
    // NNPDF
    pdfSets.push_back("NNPDF20_as_0116_100");
    pdfSets.push_back("NNPDF20_as_0117_100");
    pdfSets.push_back("NNPDF20_as_0118_100");
    pdfSets.push_back("NNPDF20_100");
    pdfSets.push_back("NNPDF20_as_0120_100");
    pdfSets.push_back("NNPDF20_as_0121_100");
    pdfSets.push_back("NNPDF20_as_0122_100");

    char *dataDir = "/smurf/yygao/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/WW/";
    SmurfSample *sample_qqww = new SmurfSample(option, QQWW, kYellow+2, "qqWW", analysis);
    SmurfSample *sample_ggww = new SmurfSample(option, GGWW, kYellow+2, "ggWW", analysis);
    SmurfSample *sample_wwmcnlo = new SmurfSample(option, WWMCNLO, kBlue+2, "WWMCNLO", analysis);
 
    for (int jetbin = 0; jetbin < 3; jetbin++) {
        sample_qqww->add(Form("%s/%ij/qqww.root", dataDir, jetbin));
        sample_ggww->add(Form("%s/%ij/ggww.root", dataDir, jetbin));
        sample_wwmcnlo->add(Form("%s/%ij/wwmcnlo.root", dataDir, jetbin));
    }
 
    // do gensets 
    looper_cteq6ll->ScanChain(sample_qqww, "cteq6ll");
    looper_cteq6ll->ScanChain(sample_ggww, "cteq6ll");
    // genset is CT10 for MC@NLO sample so no extra looping needed...

    // do other sets
    for (unsigned int i = 0; i < pdfSets.size(); ++i) {
        std::cout << "===== Doing =====> " << pdfSets[i] << std::endl;
        looper_cteq6ll->ScanChain(sample_qqww, pdfSets[i]);
        looper_cteq6ll->ScanChain(sample_ggww, pdfSets[i]);
//        looper_CT10->ScanChain(sample_wwmcnlo, pdfSets[i]);
    }

    //
    // write histograms
    // 

    const std::string outFile = Form("analysis%i_%i.root", option, int(analysis));
    saveHist(outFile.c_str());
    deleteHistos();

    //
    // tidy up
    //

    delete looper_cteq6ll;
    //delete looper_CT10;
    delete sample_qqww;
    delete sample_ggww;
    delete sample_wwmcnlo;

}

