#include "TLine.h"
#include "TFile.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include <string>
#include <math.h>
#include <iostream>

enum PDFFamily {
    PDF_CTEQ,
    PDF_CTEQ_AS,
    PDF_MSTW,
};

void FormatLegend(TLegend *lg) {
    lg->SetLineColor(kWhite);
    lg->SetShadowColor(kWhite);
    lg->SetFillColor(kWhite);
}

//
// CT10 alpha_S variation analysis
//

void AnalyseCT10as(TFile *f, std::string sampleName, TH2F *h1, TH2F *h1_up, TH2F *h1_down)
{

    TH2F *histArr[10];
    for (unsigned int i = 0; i < 10; ++i) {
        histArr[i] = (TH2F*)f->Get(Form("%s_CT10as_%i", sampleName.c_str(), i));
    }

    for (Int_t binx = 1; binx <= h1->GetNbinsX(); ++binx)
    {
        for (Int_t biny = 1; biny <= h1->GetNbinsY(); ++biny)
        {

            const unsigned int max_set = 3;
            Double_t X0 = histArr[5]->GetBinContent(binx, biny);
            //Double_t binCenter = histArr[5]->GetBinCenter(binx, biny);
            if (X0 == 0) continue;
            Double_t Xi_up = histArr[10-max_set]->GetBinContent(binx, biny);
            Double_t Xi_down = histArr[max_set]->GetBinContent(binx, biny);
            Double_t plus_max = sqrt(pow(TMath::Max(TMath::Max(Xi_up - X0, Xi_down - X0), 0.0), 2))/1.645;
            Double_t minus_max = sqrt(pow(TMath::Max(TMath::Max(X0 - Xi_up, X0 - Xi_down), 0.0), 2))/1.645;
            h1->SetBinContent(binx, biny, X0);
            h1_down->SetBinContent(binx, biny, -1*minus_max);
            h1_up->SetBinContent(binx, biny, plus_max);
        }
    }

    for (unsigned int i = 0; i < 10; ++i) {
        delete histArr[i];
    }
}

//
// MSTW/CT10 intrinsic uncertainty analysis
//

void AnalyseCTMSTW(TFile *f, std::string sampleName, std::string pdfName, const unsigned int &nsets, TH2F *h1, TH2F *h1_up, TH2F *h1_down)
{

    TH2F *histArr[nsets];
    for (unsigned int i = 0; i < nsets; ++i) {
        histArr[i] = (TH2F*)f->Get(Form("%s_%s_%i", sampleName.c_str(), pdfName.c_str(), i));
    }

    for (Int_t binx = 1; binx <= h1->GetNbinsX(); ++binx)
    {
        for (Int_t biny = 1; biny <= h1->GetNbinsY(); ++biny)
        {

            Double_t X0 = histArr[0]->GetBinContent(binx, biny);
            //Double_t binCenter = histArr[0]->GetBinCenter(binx, biny);
            Double_t plus_max = 0.0;
            Double_t minus_max = 0.0;

            if (X0 == 0) continue;

            for (unsigned int subset = 0; subset < ((nsets - 1)/2); ++subset) 
            {
                Double_t Xi_up = histArr[(subset*2) + 1]->GetBinContent(binx, biny);
                Double_t Xi_down = histArr[(subset*2) + 2]->GetBinContent(binx, biny);
                plus_max += pow(TMath::Max(TMath::Max(Xi_up - X0, Xi_down - X0), 0.0), 2);
                minus_max += pow(TMath::Max(TMath::Max(X0 - Xi_up, X0 - Xi_down), 0.0), 2);
            }
            plus_max = sqrt(plus_max);
            minus_max = sqrt(minus_max);
            h1->SetBinContent(binx, biny, X0);
            h1_down->SetBinContent(binx, biny, -1*minus_max);
            h1_up->SetBinContent(binx, biny, plus_max);
        }
    }

    for (unsigned int i = 0; i < nsets; ++i) {
        delete histArr[i];
    }

}

//
// NNPDF intrinsic uncertainty analysis
//

void AnalyseNNPDF(TFile *f, std::string sampleName, TH2F *h1, TH2F *h1_up, TH2F *h1_down)
{

    //
    // first get the average value of the set of replicas with central value alpha_S
    //

    TH2F *h1_temp = (TH2F*)f->Get(Form("%s_NNPDF20_100_0", sampleName.c_str()));
    for (Int_t binx= 1; binx <= h1->GetNbinsX(); ++binx) {
        for (Int_t biny= 1; biny <= h1->GetNbinsY(); ++biny) {
            h1->SetBinContent(binx, biny, h1_temp->GetBinContent(binx, biny));
        }
    }

    //
    // now read in the sets with larger alpha_s
    // and also smaller
    //

    std::vector<TH2F*> histArr;
    unsigned int nsets_up = 0;
    unsigned int nsets_down = 0;
    for (unsigned int i = 1; i < 101; ++i)
        histArr.push_back((TH2F*)f->Get(Form("%s_NNPDF20_100_%i", sampleName.c_str(), i)));
    for (unsigned int i = 0; i < 72; ++i) {
        histArr.push_back((TH2F*)f->Get(Form("%s_NNPDF20_as_0120_100_%i", sampleName.c_str(), i)));
        histArr.push_back((TH2F*)f->Get(Form("%s_NNPDF20_as_0118_100_%i", sampleName.c_str(), i)));
    }
    for (unsigned int i = 0; i < 27; ++i) {
        histArr.push_back((TH2F*)f->Get(Form("%s_NNPDF20_as_0121_100_%i", sampleName.c_str(), i)));
        histArr.push_back((TH2F*)f->Get(Form("%s_NNPDF20_as_0117_100_%i", sampleName.c_str(), i)));
    }
    for (unsigned int i = 0; i < 5; ++i) {
        histArr.push_back((TH2F*)f->Get(Form("%s_NNPDF20_as_0122_100_%i", sampleName.c_str(), i)));
        histArr.push_back((TH2F*)f->Get(Form("%s_NNPDF20_as_0116_100_%i", sampleName.c_str(), i)));
    }

    //
    // make the bands
    //

    for (Int_t binx = 1; binx <= h1->GetNbinsX(); ++binx)
    {
        for (Int_t biny = 1; biny <= h1->GetNbinsY(); ++biny)
        {

            float X0 = h1->GetBinContent(binx, biny);
            if (X0 == 0) continue;

            float delta_up = 0.0;
            float delta_down = 0.0;
            for (unsigned int i = 0; i < histArr.size(); ++i) {
                if (histArr[i]->GetBinContent(binx, biny) > h1->GetBinContent(binx, biny)) {
                    delta_up += pow(histArr[i]->GetBinContent(binx, biny) - h1->GetBinContent(binx, biny), 2);
                    ++nsets_up;
                }
                else {
                    delta_down += pow(h1->GetBinContent(binx, biny) - histArr[i]->GetBinContent(binx, biny), 2);
                    ++nsets_down;
                }
            }

            delta_up = sqrt( (1.0/(float(nsets_up) - 1.0)) * delta_up);
            delta_down = sqrt( (1.0/(float(nsets_down) - 1.0)) * delta_down);
            h1_up->SetBinContent(binx, biny, delta_up);
            h1_down->SetBinContent(binx, biny, -1*delta_down);
        }
    }

    //
    // clean up
    //

    delete h1_temp;
    for (unsigned int i = 0; i < histArr.size(); ++i) delete histArr[i];

    //
    // now compare to the sets with smaller and larger alpha_S
    //

}

// sampleName e.g. qqWW_DF_0j...
// genPdf e.g. cteq6ll_0 or CT10_5...
void process2D(std::string fileName, std::string sampleName, std::string genPdf)
{

    gSystem->mkdir("results2");
    gROOT->ProcessLine(".L ~/tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.3f");

    TFile *f = new TFile(Form("%s.root", fileName.c_str()), "READ");
    gROOT->cd();

    TH2F *h1_temp = (TH2F*)f->Get(Form("%s_%s", sampleName.c_str(), genPdf.c_str()))->Clone("temp");
    Int_t nbinsx    = h1_temp->GetXaxis()->GetNbins();
    Int_t nbinsy    = h1_temp->GetYaxis()->GetNbins();
    const Double_t* xbins_tmp = h1_temp->GetXaxis()->GetXbins()->GetArray();
    const Double_t* ybins_tmp = h1_temp->GetYaxis()->GetXbins()->GetArray();
    Double_t xbins[nbinsx+1];
    Double_t ybins[nbinsy+1];
    for (unsigned int i = 0; i < nbinsx+1; ++i)   xbins[i] = xbins_tmp[i];
    for (unsigned int i = 0; i < nbinsy+1; ++i)   ybins[i] = ybins_tmp[i];
    delete h1_temp;

    //
    // CTEQ
    //

    // central
    TH2F *h1_CT10   = new TH2F(Form("%s_h1_CT10", sampleName.c_str()),  "centre", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_CT10_up   = new TH2F(Form("%s_h1_CT10_up", sampleName.c_str()),  "centre up", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_CT10_down  = new TH2F(Form("%s_h1_CT10_down", sampleName.c_str()), "centre down", nbinsx, xbins, nbinsy, ybins);
    AnalyseCTMSTW(f, sampleName, "CT10", 53, h1_CT10, h1_CT10_up, h1_CT10_down);

    // alpha_s
    TH2F *h1_CT10_alpha_S   = new TH2F(Form("%s_h1_CT10_alpha_S", sampleName.c_str()),  "alpha_S", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_CT10_alpha_S_up   = new TH2F(Form("%s_h1_CT10_alpha_S_up", sampleName.c_str()),  "alpha_S up", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_CT10_alpha_S_down  = new TH2F(Form("%s_h1_CT10_alpha_S_down", sampleName.c_str()), "alpha_S down", nbinsx, xbins, nbinsy, ybins);

    AnalyseCT10as(f, sampleName, h1_CT10_alpha_S, h1_CT10_alpha_S_up, h1_CT10_alpha_S_down);
    h1_CT10_alpha_S_up->SetLineColor(kRed);
    h1_CT10_alpha_S_down->SetLineColor(kRed);
    h1_CT10_alpha_S_down->SetLineStyle(kDashed);


    TH2F *h1_CT10_env_up   = new TH2F(Form("%s_h1_CT10_env_up", sampleName.c_str()),  "CT10 up", nbinsx, xbins, nbinsy, ybins);    
    TH2F *h1_CT10_env_down   = new TH2F(Form("%s_h1_CT10_env_down", sampleName.c_str()),  "CT10 down", nbinsx, xbins, nbinsy, ybins);
    h1_CT10_env_up->SetLineColor(kRed);
    h1_CT10_env_down->SetLineColor(kRed);
    for (unsigned int binx = 0; binx < h1_CT10->GetNbinsX() + 1; ++binx) {
        for (unsigned int biny = 0; biny < h1_CT10->GetNbinsY() + 1; ++biny) {
            float up = sqrt(pow(h1_CT10_up->GetBinContent(binx, biny), 2) 
                    + pow(h1_CT10_alpha_S_up->GetBinContent(binx, biny), 2));
            float down = -1*sqrt(pow(h1_CT10_down->GetBinContent(binx, biny), 2) 
                    + pow(h1_CT10_alpha_S_down->GetBinContent(binx, biny), 2));
            h1_CT10_env_up->SetBinContent(binx, biny, up);
            h1_CT10_env_down->SetBinContent(binx, biny, down);
        }
    }

    //
    // MSTW
    //

    TH2F *h1_MSTW   = new TH2F(Form("%s_h1_MSTW", sampleName.c_str()),  "centre", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_MSTW_up   = new TH2F(Form("%s_h1_MSTW_up", sampleName.c_str()),  "centre up", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_MSTW_down  = new TH2F(Form("%s_h1_MSTW_down", sampleName.c_str()), "centre down", nbinsx, xbins, nbinsy, ybins);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl", 41, h1_MSTW, h1_MSTW_up, h1_MSTW_down);

    TH2F *h1_MSTW_asup   = new TH2F(Form("%s_h1_MSTW_asup", sampleName.c_str()),  "asup", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_MSTW_asup_up   = new TH2F(Form("%s_h1_MSTW_asup_up", sampleName.c_str()),  "asup up", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_MSTW_asup_down  = new TH2F(Form("%s_h1_MSTW_asup_down", sampleName.c_str()), "asup down", nbinsx, xbins, nbinsy, ybins);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl_asmz+68cl", 41, h1_MSTW_asup, h1_MSTW_asup_up, h1_MSTW_asup_down);
    h1_MSTW_asup_up->SetLineColor(kRed);
    h1_MSTW_asup_down->SetLineColor(kRed);

    TH2F *h1_MSTW_asdown   = new TH2F(Form("%s_h1_MSTW_asdown", sampleName.c_str()),  "asdown", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_MSTW_asdown_up   = new TH2F(Form("%s_h1_MSTW_asdown_up", sampleName.c_str()),  "asdown up", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_MSTW_asdown_down  = new TH2F(Form("%s_h1_MSTW_asdown_down", sampleName.c_str()), "asdown down", nbinsx, xbins, nbinsy, ybins);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl_asmz-68cl", 41, h1_MSTW_asdown, h1_MSTW_asdown_up, h1_MSTW_asdown_down);
    h1_MSTW_asdown_up->SetLineColor(kRed);
    h1_MSTW_asdown_up->SetLineStyle(kDashed);
    h1_MSTW_asdown_down->SetLineColor(kRed);
    h1_MSTW_asdown_down->SetLineStyle(kDashed);

    TH2F *h1_MSTW_asup05   = new TH2F(Form("%s_h1_MSTW_asup05", sampleName.c_str()),  "asup05", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_MSTW_asup05_up   = new TH2F(Form("%s_h1_MSTW_asup05_up", sampleName.c_str()),  "asup05 up", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_MSTW_asup05_down  = new TH2F(Form("%s_h1_MSTW_asup05_down", sampleName.c_str()), "asup05 down", nbinsx, xbins, nbinsy, ybins);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl_asmz+68clhalf", 41, h1_MSTW_asup05, h1_MSTW_asup05_up, h1_MSTW_asup05_down);
    h1_MSTW_asup05_up->SetLineColor(kBlue);
    h1_MSTW_asup05_down->SetLineColor(kBlue);

    TH2F *h1_MSTW_asdown05   = new TH2F(Form("%s_h1_MSTW_asdown05", sampleName.c_str()),  "asdown05", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_MSTW_asdown05_up   = new TH2F(Form("%s_h1_MSTW_asdown05_up", sampleName.c_str()),  "asdown05 up", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_MSTW_asdown05_down  = new TH2F(Form("%s_h1_MSTW_asdown05_down", sampleName.c_str()), "asdown05 down", nbinsx, xbins, nbinsy, ybins);
    AnalyseCTMSTW(f, sampleName, "MSTW2008nlo68cl_asmz-68clhalf", 41, h1_MSTW_asdown05, h1_MSTW_asdown05_up, h1_MSTW_asdown05_down);
    h1_MSTW_asdown05_up->SetLineColor(kBlue);
    h1_MSTW_asdown05_up->SetLineStyle(kDashed);
    h1_MSTW_asdown05_down->SetLineColor(kBlue);
    h1_MSTW_asdown05_down->SetLineStyle(kDashed);

    // now combine to get the PDF + alpha_S uncertainty
    // according to the envelope described in AN2011/055
    TH2F *h1_MSTW_env_up   = new TH2F(Form("%s_h1_MSTW_env_up", sampleName.c_str()),  "MSTW up", nbinsx, xbins, nbinsy, ybins);    
    TH2F *h1_MSTW_env_down   = new TH2F(Form("%s_h1_MSTW_env_down", sampleName.c_str()),  "MSTW down", nbinsx, xbins, nbinsy, ybins);
    h1_MSTW_env_up->SetLineColor(kRed);
    h1_MSTW_env_down->SetLineColor(kRed);
    for (unsigned int binx = 0; binx < h1_MSTW->GetNbinsX() + 1; ++binx) {
        for (unsigned int biny = 0; biny < h1_MSTW->GetNbinsY() + 1; ++biny) {
            float delta_A_up = h1_MSTW_up->GetBinContent(binx, biny);
            float delta_A_asup = h1_MSTW_asup_up->GetBinContent(binx, biny);
            float delta_A_asup05 = h1_MSTW_asup05_up->GetBinContent(binx, biny);
            float delta_A_down = h1_MSTW_down->GetBinContent(binx, biny);
            float delta_A_asdown = h1_MSTW_asdown_down->GetBinContent(binx, biny);
            float delta_A_asdown05 = h1_MSTW_asdown05_down->GetBinContent(binx, biny);
            h1_MSTW_env_up->SetBinContent(binx, biny, TMath::Max(TMath::Max(delta_A_up, delta_A_asup), delta_A_asup05));
            h1_MSTW_env_down->SetBinContent(binx, biny, TMath::Min(TMath::Min(delta_A_down, delta_A_asdown), delta_A_asdown05));
        }
    }

    //
    // NNPDF
    // 

    TH2F *h1_NNPDF   = new TH2F(Form("%s_h1_NNPDF", sampleName.c_str()),  "centre", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_NNPDF_up   = new TH2F(Form("%s_h1_NNPDF_up", sampleName.c_str()),  "centre up", nbinsx, xbins, nbinsy, ybins);
    TH2F *h1_NNPDF_down  = new TH2F(Form("%s_h1_NNPDF_down", sampleName.c_str()), "centre down", nbinsx, xbins, nbinsy, ybins);
    AnalyseNNPDF(f, sampleName, h1_NNPDF, h1_NNPDF_up, h1_NNPDF_down);

    //
    // ! GRAND TOTAL !
    //

    TH2F *h2_GenPDF = (TH2F*)f->Get(Form("%s_%s", sampleName.c_str(), genPdf.c_str()))->Clone("GenPDF");
    TH2F *h2_env_up   = new TH2F(Form("%s_h2_env_up", sampleName.c_str()),  "up", nbinsx, xbins, nbinsy, ybins);
    TH2F *h2_env_down   = new TH2F(Form("%s_h2_env_down", sampleName.c_str()),  "down", nbinsx, xbins, nbinsy, ybins);
    h2_env_up->SetFillColor(kGray);
    h2_env_up->SetLineColor(kWhite);
    h2_env_down->SetFillColor(10);
    h2_env_down->SetLineColor(kWhite);

    for (unsigned int binx = 0; binx < h2_GenPDF->GetNbinsX() + 1; ++binx) {
        for (unsigned int biny = 0; biny < h2_GenPDF->GetNbinsY() + 1; ++biny) {

            // get up and down variations including both
            // - intrinsic pdf variation from central value
            // - variation of central value from generated set
            Float_t up_NNPDF    = h1_NNPDF->GetBinContent(binx, biny) + h1_NNPDF_up->GetBinContent(binx, biny);
            Float_t up_CT10     = h1_CT10->GetBinContent(binx, biny) + h1_CT10_env_up->GetBinContent(binx, biny);
            Float_t up_MSTW     = h1_MSTW->GetBinContent(binx, biny) + h1_MSTW_env_up->GetBinContent(binx, biny);
            Float_t down_NNPDF  = h1_NNPDF->GetBinContent(binx, biny) + h1_NNPDF_down->GetBinContent(binx, biny);
            Float_t down_CT10   = h1_CT10->GetBinContent(binx, biny) + h1_CT10_env_down->GetBinContent(binx, biny);
            Float_t down_MSTW   = h1_MSTW->GetBinContent(binx, biny) + h1_MSTW_env_down->GetBinContent(binx, biny);

            // get the envelope of envelopes...
            h2_env_up->SetBinContent(binx, biny,    TMath::Max(TMath::Max(up_NNPDF,     up_CT10),   up_MSTW));
            h2_env_down->SetBinContent(binx, biny,  TMath::Min(TMath::Min(down_NNPDF,   down_CT10), down_MSTW));
        }
    }

    //
    // ... scale variations to central shape
    //

    // get the mid point of the up and down envelopes
    // as the central alternate pdf shape
    TH2F* h2_env_midpoint = (TH2F*)h2_env_up->Clone("h2_env_midpoint");
    h2_env_midpoint->SetTitle("midpoint");
    h2_env_midpoint->Add(h2_env_down, -1.0);
    h2_env_midpoint->Scale(0.5);
    h2_env_midpoint->Add(h2_env_down, +1.0);

    // now scale the central alternate pdf shape
    // to the generated shape in the sample
    Float_t scaleFactor = h2_GenPDF->Integral(0, h2_GenPDF->GetXaxis()->GetNbins(), 0, h2_GenPDF->GetYaxis()->GetNbins())
                            / h2_env_midpoint->Integral(0, h2_env_midpoint->GetXaxis()->GetNbins(), 0, h2_env_midpoint->GetYaxis()->GetNbins());
    h2_env_midpoint->Scale(scaleFactor);
    h2_env_up->Scale(scaleFactor);
    h2_env_down->Scale(scaleFactor);
    h1_CT10->Scale(scaleFactor);
    h1_CT10_env_up->Scale(scaleFactor);
    h1_CT10_env_down->Scale(scaleFactor);
    h1_MSTW->Scale(scaleFactor);
    h1_MSTW_env_up->Scale(scaleFactor);
    h1_MSTW_env_down->Scale(scaleFactor);
    h1_NNPDF->Scale(scaleFactor);
    h1_NNPDF_up->Scale(scaleFactor);
    h1_NNPDF_down->Scale(scaleFactor);

    // now for illustration draw collapsed version of 
    // the two shapes...  
    TH1F *h1_alternate_px           = (TH1F*)h2_env_midpoint->ProjectionX("h1_alternate_px");
    TH1F *h1_alternateUp_px         = (TH1F*)h2_env_up->ProjectionX("h1_alternateUp_px");
    TH1F *h1_alternateDown_px       = (TH1F*)h2_env_down->ProjectionX("h1_alternateDown_px");
    TH1F *h1_default_px             = (TH1F*)h2_GenPDF->ProjectionX("h1_default_px");
    TH1F *h1_CT10_px                = (TH1F*)h1_CT10->ProjectionX("h1_CT10_px");
    TH1F *h1_MSTW_px                = (TH1F*)h1_MSTW->ProjectionX("h1_MSTW_px");
    TH1F *h1_NNPDF_px                = (TH1F*)h1_NNPDF->ProjectionX("h1_NNPDF_px");
    // up
    TH1F *h1_CT10Up_px                = (TH1F*)h1_CT10_env_up->ProjectionX("h1_CT10Up_px");
    h1_CT10Up_px->Add(h1_CT10_px);
    TH1F *h1_MSTWUp_px                = (TH1F*)h1_MSTW_env_up->ProjectionX("h1_MSTWUp_px");
    h1_MSTWUp_px->Add(h1_MSTW_px);
    TH1F *h1_NNPDFUp_px                = (TH1F*)h1_NNPDF_up->ProjectionX("h1_NNPDFUp_px");
    h1_NNPDFUp_px->Add(h1_NNPDF_px);
    // down
    TH1F *h1_CT10Down_px                = (TH1F*)h1_CT10_env_down->ProjectionX("h1_CT10Down_px");
    h1_CT10Down_px->Add(h1_CT10_px);
    TH1F *h1_MSTWDown_px                = (TH1F*)h1_MSTW_env_down->ProjectionX("h1_MSTWDown_px");
    h1_MSTWDown_px->Add(h1_MSTW_px);
    TH1F *h1_NNPDFDown_px                = (TH1F*)h1_NNPDF_down->ProjectionX("h1_NNPDFDown_px");
    h1_NNPDFDown_px->Add(h1_NNPDF_px);

    // and their ratio for convenience 
    TH1F *h1_ratio = (TH1F*)h1_alternate_px->Clone("h1_ratio");
    h1_ratio->Divide(h1_default_px);
    TH1F *h1_ratioUp = (TH1F*)h1_alternateUp_px->Clone("h1_ratioUp");
    h1_ratioUp->Divide(h1_default_px);
    TH1F *h1_ratioDown = (TH1F*)h1_alternateDown_px->Clone("h1_ratioDown");
    h1_ratioDown->Divide(h1_default_px);

    TH1F *h1_ratioCT10 = (TH1F*)h1_CT10_px->Clone("h1_ratioCT10");
    h1_ratioCT10->Divide(h1_default_px);
    TH1F *h1_ratioCT10Up = (TH1F*)h1_CT10Up_px->Clone("h1_ratioCT10Up");
    h1_ratioCT10Up->Divide(h1_default_px);
    TH1F *h1_ratioCT10Down = (TH1F*)h1_CT10Down_px->Clone("h1_ratioCT10Down");
    h1_ratioCT10Down->Divide(h1_default_px);

    TH1F *h1_ratioMSTW = (TH1F*)h1_MSTW_px->Clone("h1_ratioMSTW");
    h1_ratioMSTW->Divide(h1_default_px);
    TH1F *h1_ratioMSTWUp = (TH1F*)h1_MSTWUp_px->Clone("h1_ratioMSTWUp");
    h1_ratioMSTWUp->Divide(h1_default_px);
    TH1F *h1_ratioMSTWDown = (TH1F*)h1_MSTWDown_px->Clone("h1_ratioMSTWDown");
    h1_ratioMSTWDown->Divide(h1_default_px);

    TH1F *h1_ratioNNPDF = (TH1F*)h1_NNPDF_px->Clone("h1_ratioNNPDF");
    h1_ratioNNPDF->Divide(h1_default_px);
    TH1F *h1_ratioNNPDFUp = (TH1F*)h1_NNPDFUp_px->Clone("h1_ratioNNPDFUp");
    h1_ratioNNPDFUp->Divide(h1_default_px);
    TH1F *h1_ratioNNPDFDown = (TH1F*)h1_NNPDFDown_px->Clone("h1_ratioNNPDFDown");
    h1_ratioNNPDFDown->Divide(h1_default_px);

    // and the 2-D version...
    TH2F *h2_ratio = (TH2F*)h2_env_midpoint->Clone(Form("%s_alternate", sampleName.c_str()));
    h2_ratio->Divide(h2_GenPDF);
    TH2F *h2_ratioUp = (TH2F*)h2_env_up->Clone(Form("%s_alternateUp", sampleName.c_str()));
    h2_ratioUp->Divide(h2_GenPDF);
    TH2F *h2_ratioDown = (TH2F*)h2_env_down->Clone(Form("%s_alternateDown", sampleName.c_str()));
    h2_ratioDown->Divide(h2_GenPDF);

    //
    // save it 
    //

    TFile f_out(Form("results2/PDFUncertainty_%s_%s.root",  fileName.c_str(), sampleName.c_str()), "RECREATE");
    f_out.cd();
    h2_ratio->Write();
    h2_ratioUp->Write();
    h2_ratioDown->Write();
    f_out.Close();
    

    //
    //
    //

    // sigh, this is getting tiring...
    // I need some tea
    TCanvas *c_darjeeling = new TCanvas("c_darjeeling", "c_darjeeling", 1000, 600);
    c_darjeeling->SetRightMargin(0.15);
    c_darjeeling->cd();

    TLine *line = new TLine(h1_ratio->GetXaxis()->GetXmin(), 1.0, h1_ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kBlack);

    // individual CT10 
    h1_ratioCT10Up->Draw("HIST");
    h1_ratioCT10Up->SetLineStyle(kDashed);
    h1_ratioCT10Up->SetLineColor(kRed);
    h1_ratioCT10Up->GetYaxis()->SetRangeUser(0.8, 1.2);
    h1_ratioCT10Up->GetYaxis()->SetTitle("Ratio {Alternate/Default}");
    h1_ratioCT10Down->Draw("HIST SAME");
    h1_ratioCT10Down->SetLineStyle(kDashed);
    h1_ratioCT10Down->SetLineColor(kRed);
    h1_ratioCT10->Draw("HIST SAME");
    h1_ratioCT10->SetLineWidth(2);
    h1_ratioCT10->SetLineColor(kRed);
    h1_ratioCT10->SetFillStyle(0);
    line->Draw();
    c_darjeeling->RedrawAxis();
    c_darjeeling->SaveAs(Form("results2/%s_%s_alternateShapeRatioCT10_px.png", fileName.c_str(), sampleName.c_str()));

    // individual MSTW
    h1_ratioMSTWUp->Draw("HIST");
    h1_ratioMSTWUp->SetLineStyle(kDashed);
    h1_ratioMSTWUp->SetLineColor(kBlue);
    h1_ratioMSTWUp->GetYaxis()->SetRangeUser(0.8, 1.2);
    h1_ratioMSTWUp->GetYaxis()->SetTitle("Ratio {Alternate/Default}");
    h1_ratioMSTWDown->Draw("HIST SAME");
    h1_ratioMSTWDown->SetLineStyle(kDashed);
    h1_ratioMSTWDown->SetLineColor(kBlue);
    h1_ratioMSTW->Draw("HIST SAME");
    h1_ratioMSTW->SetLineWidth(2);
    h1_ratioMSTW->SetLineColor(kBlue);
    h1_ratioMSTW->SetFillStyle(0);
    line->Draw();
    c_darjeeling->RedrawAxis();
    c_darjeeling->SaveAs(Form("results2/%s_%s_alternateShapeRatioMSTW_px.png", fileName.c_str(), sampleName.c_str()));

    // individual NNPDF
    h1_ratioNNPDFUp->Draw("HIST");
    h1_ratioNNPDFUp->SetLineStyle(kDashed);
    h1_ratioNNPDFUp->SetLineColor(kBlack);
    h1_ratioNNPDFUp->GetYaxis()->SetRangeUser(0.8, 1.2);
    h1_ratioNNPDFUp->GetYaxis()->SetTitle("Ratio {Alternate/Default}");
    h1_ratioNNPDFDown->Draw("HIST SAME");
    h1_ratioNNPDFDown->SetLineStyle(kDashed);
    h1_ratioNNPDFDown->SetLineColor(kBlack);
    h1_ratioNNPDF->Draw("HIST SAME");
    h1_ratioNNPDF->SetLineWidth(2);
    h1_ratioNNPDF->SetLineColor(kBlack);
    h1_ratioNNPDF->SetFillStyle(0);
    line->Draw();
    c_darjeeling->RedrawAxis();
    c_darjeeling->SaveAs(Form("results2/%s_%s_alternateShapeRatioNNPDF_px.png", fileName.c_str(), sampleName.c_str()));

    // alternate shape
    h2_env_midpoint->Draw("COLZ TEXT");
    c_darjeeling->SaveAs(Form("results2/%s_%s_alternateShape.png", fileName.c_str(), sampleName.c_str()));

    // and its ratio
    Int_t palette[7];
    palette[0] = kWhite;
    for (unsigned int i=0;i<7;++i){
        palette[i] = 18-i;
    }
    gStyle->SetPalette(7,palette);
    gStyle->SetPaintTextFormat("4.2f");

    h2_ratio->Draw("COLZ TEXT");
    h2_ratio->SetMaximum(1.1);
    h2_ratio->SetMinimum(0.9);
    c_darjeeling->SaveAs(Form("results2/%s_%s_alternateShapeRatio.png", fileName.c_str(), sampleName.c_str()));

    // up and down envelope
    h2_env_up->Draw("COLZ TEXT");
    c_darjeeling->SaveAs(Form("results2/%s_%s_env_up.png", fileName.c_str(), sampleName.c_str()));
    h2_env_down->Draw("COLZ TEXT");
    c_darjeeling->SaveAs(Form("results2/%s_%s_env_down.png", fileName.c_str(), sampleName.c_str()));

    // alternate projected onto x axis
    h1_alternate_px->Draw("HIST");
    h1_alternate_px->SetLineWidth(2);
    h1_alternate_px->SetLineStyle(kDashed);
    h1_alternate_px->SetLineColor(kBlack);
    h1_alternate_px->SetFillStyle(0);
    h1_default_px->Draw("SAME HIST");
    h1_default_px->SetLineWidth(2);
    c_darjeeling->SaveAs(Form("results2/%s_%s_alternateShape_px.png", fileName.c_str(), sampleName.c_str()));
    c_darjeeling->SetLogy(1);
    c_darjeeling->SaveAs(Form("results2/%s_%s_alternateShape_px_log.png", fileName.c_str(), sampleName.c_str()));

    // ratio of central and alternate...
    c_darjeeling->SetLogy(0);
    TLegend *l1 = new TLegend(0.2, 0.6, 0.6, 0.85);
    l1->SetLineColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->SetFillColor(kWhite);
    l1->AddEntry(h1_ratio, "PDF variation midpoint", "l");
    l1->AddEntry(h1_ratioUp, "PDF variation +/- 1#sigma", "f");

    h1_ratioUp->Draw("HIST");
    h1_ratioUp->GetYaxis()->SetRangeUser(0.9, 1.3);
    h1_ratioUp->GetYaxis()->SetTitle("Ratio {Alternate/Default}");
    h1_ratioDown->Draw("HIST SAME");
    h1_ratio->Draw("HIST SAME");
    h1_ratio->SetLineWidth(2);
    h1_ratio->SetLineColor(kBlack);
    h1_ratio->SetFillStyle(0);
    line->Draw();
    l1->Draw();
    c_darjeeling->RedrawAxis();
    c_darjeeling->SaveAs(Form("results2/%s_%s_alternateShapeRatio_px.png", fileName.c_str(), sampleName.c_str()));


    l1->AddEntry(h1_ratioCT10, "CT10 variation midpoint", "l");
    l1->AddEntry(h1_ratioMSTW, "MSTW variation midpoint", "l");
    l1->AddEntry(h1_ratioNNPDF, "NNPDF variation midpoint", "l");
    h1_ratioUp->Draw("HIST");
    h1_ratioDown->Draw("HIST SAME");
    h1_ratio->Draw("HIST SAME");
    h1_ratioCT10->Draw("HIST SAME");
    h1_ratioCT10->SetLineWidth(2);
    h1_ratioCT10->SetLineColor(kRed);
    h1_ratioCT10->SetLineStyle(kDashed);
    h1_ratioCT10->SetFillStyle(0);
    h1_ratioMSTW->Draw("HIST SAME");
    h1_ratioMSTW->SetLineWidth(2);
    h1_ratioMSTW->SetLineColor(kBlue);
    h1_ratioMSTW->SetLineStyle(kDashed+1);
    h1_ratioMSTW->SetFillStyle(0);
    h1_ratioNNPDF->Draw("HIST SAME");
    h1_ratioNNPDF->SetLineWidth(2);
    h1_ratioNNPDF->SetLineColor(kBlack);
    h1_ratioNNPDF->SetLineStyle(kDashed+2);
    h1_ratioNNPDF->SetFillStyle(0);
    line->Draw();
    l1->Draw();
    c_darjeeling->SaveAs(Form("results2/%s_%s_alternateShapeRatioDetailed_px.png", fileName.c_str(), sampleName.c_str()));
    
    f->Close();
}

