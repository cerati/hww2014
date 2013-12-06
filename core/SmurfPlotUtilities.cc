
#include "SmurfPlotUtilities.h"
#include "SmurfSample.h"

#include "TROOT.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TList.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TRegexp.h"
#include "TKey.h"
#include <iostream>
#include <TIterator.h>

#include <iostream>
#include <cmath>



void deleteHistos() {
    // Delete all existing histograms in memory
    TObject* obj;
    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();
    while ( (obj = (iter->Next())) ) {
        if (obj->IsA()->InheritsFrom(TH1::Class()) ||
                obj->IsA()->InheritsFrom(TH2::Class()) ) {
                delete obj;
            }
    }
}

void saveHist(const char* filename, const char* pat)
{

    TList* list = gDirectory->GetList() ;
    TIterator* iter = list->MakeIterator();
    TRegexp re(pat,kTRUE) ;

    TFile outf(filename,"RECREATE") ;
    printf("[SmurfPlotUtilities::saveHist] Saving histograms to %s\n", filename);

    int counter = 0;
    while(TObject* obj=iter->Next()) {

        // don't save TH1Keys objects
        // this is a bug fudge
        if (TString(obj->GetName()).Contains("histokeys")) continue;

        // save other stuff
        if (TString(obj->GetName()).Index(re)>=0) {
            obj->Write() ;
            ++counter;
        }

    }

    printf("[SmurfPlotUtilities::saveHist] Saved %i histograms\n", counter);
    outf.Close() ;
    delete iter ;
}

void FillHist(TH1F **hist, const unsigned int type, const float &val, const float &weight)
{
    hist[type]->Fill(val, weight);
    hist[6]->Fill(val, weight);
    if (type == 0 || type == 3) hist[4]->Fill(val, weight);
    if (type == 1 || type == 2) hist[5]->Fill(val, weight);
}

void FormatHist(TH1F **hist, SmurfSample *sample, const char *name, 
        const char *title, int n, float min, float max)
{

    gROOT->cd();
    for (unsigned int k = 0; k < kLeptonTypes; ++k) {

        char *title;
        title = Form("%s_%s_%s", sample->getName().c_str(), name, types[k]);
        hist[k] = new TH1F(title, title, n, min, max);
        hist[k]->GetXaxis()->SetTitle(title);
        hist[k]->SetFillStyle(1001);
        hist[k]->SetFillColor(sample->getColour());
        hist[k]->SetLineColor(sample->getColour());
        hist[k]->Sumw2();
        //hist[k]->SetDirectory(0);
    }

}

void FormatHist(TH1F **hist, SmurfSample *sample, const char *name,
        const char *title, int n, float bins[])
{

    gROOT->cd();
    for (unsigned int k = 0; k < kLeptonTypes; ++k) {

        char *title;
        title = Form("%s_%s_%s", sample->getName().c_str(), name, types[k]);
        hist[k] = new TH1F(title, title, n, bins);
        hist[k]->GetXaxis()->SetTitle(title);
        hist[k]->SetFillStyle(1001);
        hist[k]->SetFillColor(sample->getColour());
        hist[k]->SetLineColor(sample->getColour());
        hist[k]->Sumw2();
        //hist[k]->SetDirectory(0);
    }

}

typedef TH1F H;

H cumulate (const H &in, bool increasing)
{
    H h_out(in.GetName() + TString("tmp"), in.GetTitle(), in.GetNbinsX(),
            in.GetBinLowEdge(1), in.GetBinLowEdge(in.GetNbinsX() + 1));
    h_out.Sumw2();
    h_out.SetFillColor(in.GetFillColor());
    h_out.SetFillStyle(in.GetFillStyle());
    h_out.SetLineStyle(in.GetLineStyle());
    h_out.SetLineColor(in.GetLineColor());
    h_out.SetLineWidth(in.GetLineWidth());
    double sum = 0;
    double err2 = 0;
    if (increasing) {
        for (int j = 0; j <= in.GetNbinsX() + 1; ++j) {
            sum += in.GetBinContent(j);
            err2 += in.GetBinError(j)*in.GetBinError(j);
            h_out.SetBinContent(j, sum);
            h_out.SetBinError(j, sqrt(err2));
        }
    } else {
        for (int j = in.GetNbinsX() + 1; j >= 0; --j) {
            sum += in.GetBinContent(j);
            err2 += in.GetBinError(j)*in.GetBinError(j);
            h_out.SetBinContent(j, sum);
            h_out.SetBinError(j, sqrt(err2));
        }
    }
    return h_out;
}

TGraph eff_rej (const H &signal, H &background, bool normalize, bool increasing)
{
    H sig = *(TH1F*)signal.Clone("h_tmp_s");
    if (normalize)
        sig.Scale(1 / sig.Integral(0, sig.GetNbinsX() + 1));
    H bg = *(TH1F*)background.Clone("h_tmp_bg");
    if (normalize)
        bg.Scale(1 / bg.Integral(0, bg.GetNbinsX() + 1));
    H sig_cum = cumulate(sig, increasing);
    H bg_cum = cumulate(bg, increasing);
    TGraph ret(signal.GetNbinsX());
    for (int i = 1; i <= signal.GetNbinsX(); ++i) {
        const double x = sig_cum.GetBinCenter(i);
        const double sig = sig_cum.GetBinContent(i);
        const double bg = bg_cum.GetBinContent(bg_cum.FindBin(x));
        ret.SetPoint(i - 1, sig, bg); // gotta love offsets
        //        printf("point %d: %f sig, %f bg\n", i, sig, bg);
    }
    return ret;
}

TArrow *GetGraphArrow(TGraph *gr, int bin, Color_t colour) 
{

    Double_t x = 0.0;
    Double_t y = 0.0;
    gr->GetPoint(bin, x, y);
    TArrow *arrow = new TArrow(x, (y/2.0), x, y);
    arrow->SetLineWidth(2);
    arrow->SetLineColor(colour);
    return arrow;
}

TH1F *GetGammaJetExtrap(TFile *file, std::vector<SmurfSample*> bgSamples,
        SmurfSample* gammaSample,
        SmurfSample* dataSample,
        const char *name, const char *title)
{

    std::string histname = Form("%s_%s", gammaSample->getName().c_str(), name);
    TH1F *gamma = (TH1F*)file->Get(histname.c_str())->Clone();

    TH1F *data;
    histname = Form("%s_%s", dataSample->getName().c_str(), name);
    data = (TH1F*)file->Get(histname.c_str())->Clone();

    TH1F *background;
    TH1F *temp;
    for (unsigned int s = 0; s < bgSamples.size(); ++s) {
        histname = Form("%s_%s", bgSamples[s]->getName().c_str(), name);
        if (s == 0) background = (TH1F*)file->Get(histname.c_str())->Clone();
        else {
            temp = (TH1F*)file->Get(histname.c_str())->Clone();
            background->Add(temp);
        }
    }

    TH1F *data_bgsub = (TH1F*)data->Clone("data_bgsub");
    data_bgsub->Add(background, -1.0);

    data_bgsub->Rebin(2);
    gamma->Rebin(2);

    TH1F *ratio = (TH1F*)gamma->Clone("ratio");
    ratio->Divide(data_bgsub);
    ratio->SetLineColor(gammaSample->getColour());
    ratio->SetMarkerColor(gammaSample->getColour());
    ratio->GetXaxis()->SetRangeUser(0, 47.5);
    ratio->GetYaxis()->SetRangeUser(0.0, 1.5);
    ratio->GetXaxis()->SetTitle(title);
    ratio->GetYaxis()->SetTitle("#gamma/Data");

    //TCanvas *canvas = new TCanvas();
    //canvas->cd();
    //ratio->Draw("HIST E1");

    return ratio;
    //return canvas;

}


TGraph *GetEffRej(TFile *file, std::vector<SmurfSample*> bgSamples,
        std::vector<SmurfSample*> signalSamples,
        const char *name, const char *title, bool increasing)
{

    TH1F *temp;
    TH1F *signal;
    for (unsigned int s = 0; s < signalSamples.size(); ++s) {
        const std::string histname = Form("%s_%s", signalSamples[s]->getName().c_str(), name);
        if (s == 0) signal = (TH1F*)file->Get(histname.c_str())->Clone();
        else {
            temp = (TH1F*)file->Get(histname.c_str())->Clone();
            signal->Add(temp);
        }

    }

    std::cout << signal->FindBin(2.87979326579064354e+00) << std::endl;
    std::cout << signal->FindBin(0.28) << std::endl;


    TH1F *background;
    for (unsigned int s = 0; s < bgSamples.size(); ++s) {
        const std::string histname = Form("%s_%s", bgSamples[s]->getName().c_str(), name);
        if (s == 0) background = (TH1F*)file->Get(histname.c_str())->Clone();
        else {
            temp = (TH1F*)file->Get(histname.c_str())->Clone();
            background->Add(temp);
        }
    }

    TGraph *gr_eff_rej  = (TGraph*)eff_rej(*background, *signal, true, increasing).Clone(name);
    gr_eff_rej->GetYaxis()->SetTitle("#Signal");
    gr_eff_rej->GetXaxis()->SetTitle("#Background");
    return gr_eff_rej;

}

TH1F *GetHistogram(TFile *file, std::vector<SmurfSample*> samples, const char *name, int rebin)
{

    TH1F *temp = 0;
    TH1F *hist = 0;
    for (unsigned int s = 0; s < samples.size(); ++s) {
        const std::string histname = Form("%s_%s", samples[s]->getName().c_str(), name);
        //std::cout << histname << std::endl;
        if (s == 0) {
            hist = (TH1F*)file->Get(histname.c_str())->Clone();
        } else {
            temp = (TH1F*)file->Get(histname.c_str())->Clone();
            hist->Add(temp);
        }
    }

    float lastBinContent = hist->GetBinContent(hist->GetNbinsX());
    float lastBinError = hist->GetBinError(hist->GetNbinsX());
    float overflowBinContent = hist->GetBinContent(hist->GetNbinsX() + 1);
    float overflowBinError = hist->GetBinError(hist->GetNbinsX() + 1);
    hist->SetBinContent(hist->GetNbinsX() + 1, 0.0);
    hist->SetBinError(hist->GetNbinsX() + 1, 0.0);
    hist->SetBinContent(hist->GetNbinsX(), lastBinContent + overflowBinContent);
    hist->SetBinError(hist->GetNbinsX(), sqrt(lastBinError*lastBinError + overflowBinError*overflowBinError));
    hist->Rebin(rebin);
    
    if (temp !=0) delete temp;
    return hist;

}

TH1F *GetHistogram(TFile *file, SmurfSample* sample, const char *name, int rebin)
{
    
    TH1F *hist = 0;
    const std::string histname = Form("%s_%s", sample->getName().c_str(), name);
    //std::cout << histname << std::endl;
    hist = (TH1F*)file->Get(histname.c_str())->Clone();
    float lastBinContent = hist->GetBinContent(hist->GetNbinsX());
    float lastBinError = hist->GetBinError(hist->GetNbinsX());
    float overflowBinContent = hist->GetBinContent(hist->GetNbinsX() + 1);
    float overflowBinError = hist->GetBinError(hist->GetNbinsX() + 1);
    hist->SetBinContent(hist->GetNbinsX() + 1, 0.0);
    hist->SetBinError(hist->GetNbinsX() + 1, 0.0);
    hist->SetBinContent(hist->GetNbinsX(), lastBinContent + overflowBinContent);
    hist->SetBinError(hist->GetNbinsX(), sqrt(lastBinError*lastBinError + overflowBinError*overflowBinError));
    hist->Rebin(rebin);
    return hist;

}

TCanvas *GetStack(TFile *file, 
        std::vector<SmurfSample*> bgSamples, 
        std::vector<SmurfSample*> signalSamples,
        SmurfSample *dataSample, 
        const char *name, const char *title, float lumi, float xmin, float xmax, bool log, int rebin)
{

    TCanvas *canvas = new TCanvas();
    TLegend *legend = new TLegend(0.68, 0.55, 0.93, 0.90);
    legend->SetFillColor(10);
    legend->SetLineColor(kWhite);
    legend->SetShadowColor(kWhite);
    THStack *stack = new THStack();

    TH1F *temp = 0;
    for (unsigned int s = 0; s < bgSamples.size(); ++s) {
        const std::string histname = Form("%s_%s", bgSamples[s]->getName().c_str(), name);
        //std::cout << histname << std::endl;
        temp = (TH1F*)file->Get(histname.c_str())->Clone();
        temp->Rebin(rebin);
        if (bgSamples[s]->getDataType() == GAMMA) {
            temp->SetFillStyle(0);
            temp->SetLineColor(bgSamples[s]->getColour());
        }
        if (bgSamples[s]->getDataType() == TOPDATA) {
            temp->SetFillStyle(0);
            temp->SetLineColor(bgSamples[s]->getColour());
        }

        temp->SetFillColor(bgSamples[s]->getColour());
        temp->SetLineColor(bgSamples[s]->getColour());
        temp->SetBinContent(temp->GetNbinsX(), temp->GetBinContent(temp->GetNbinsX()) + temp->GetBinContent(temp->GetNbinsX()+1));
        temp->SetBinContent(temp->GetNbinsX()+1, 0.0);

        stack->Add(temp);
        legend->AddEntry(temp, bgSamples[s]->getName().c_str(), "f");
    }

    TH1F *signal;
    for (unsigned int s = 0; s < signalSamples.size(); ++s) {
        const std::string histname = Form("%s_%s", signalSamples[s]->getName().c_str(), name);
        std::cout << histname << std::endl;
        if (s == 0) {
            signal = (TH1F*)file->Get(histname.c_str())->Clone();
            signal->Rebin(rebin);
        }
        else {
            temp = (TH1F*)file->Get(histname.c_str())->Clone();
            temp->Rebin(rebin);
            signal->Add(temp);
        }
        signal->SetFillColor(0);
        signal->SetLineColor(kBlack);
        signal->SetBinContent(signal->GetNbinsX(), signal->GetBinContent(signal->GetNbinsX()) + signal->GetBinContent(signal->GetNbinsX()+1));
        signal->SetBinContent(signal->GetNbinsX()+1, 0.0);
    }
    
    if (signalSamples.size() > 0) legend->AddEntry(signal, "signal", "l");

    TH1F *data = 0;
    if (dataSample != 0) {
        const std::string histname = Form("%s_%s", dataSample->getName().c_str(), name);
        data = (TH1F*)file->Get(histname.c_str())->Clone();
        data->SetBinContent(data->GetNbinsX(), data->GetBinContent(data->GetNbinsX()) + data->GetBinContent(data->GetNbinsX()+1));
        data->SetBinContent(data->GetNbinsX()+1, 0.0);
        data->SetBinError(data->GetNbinsX(), sqrt(data->GetBinContent(data->GetNbinsX())) );
        data->SetMarkerStyle(20);
        data->Rebin(rebin);
        legend->AddEntry(data, dataSample->getName().c_str(), "lp");
    }

    TLatex * tex = new TLatex(0.2,0.80,Form("#sqrt{s}=7 TeV, #int Ldt = %4.0f pb^{-1}", lumi));
    tex->SetNDC();
    tex->SetTextSize(0.035);
    tex->SetLineWidth(2);
    TLatex * tex2 = new TLatex(0.2,0.86,"CMS Preliminary 2011");
    tex2->SetNDC();
    tex2->SetTextSize(0.035);
    tex2->SetLineWidth(2);

    canvas->cd();
    stack->Draw("HIST");
    if (dataSample != 0) data->Draw("SAME E1");
    if (signalSamples.size() > 0) signal->Draw("SAME HIST");
    legend->Draw("SAME");
    tex->Draw("SAME");
    tex2->Draw("SAME");

    double yMax = stack->GetMaximum();
    if (signalSamples.size() > 0) yMax = TMath::Max(stack->GetMaximum(), signal->GetMaximum());
    if (dataSample != 0) yMax = TMath::Max(stack->GetMaximum(), data->GetMaximum());
    stack->SetMaximum(yMax + TMath::Max(2*sqrt(yMax), yMax));
    if (log) {
        stack->SetMaximum(yMax*10000);
        canvas->SetLogy(1);
    }
    stack->SetMinimum(0.01);
    stack->GetXaxis()->SetRangeUser(xmin, xmax);
    stack->GetXaxis()->SetTitle(title);

    canvas->RedrawAxis();

    return canvas;

}

TCanvas *ComparePlots(TFile *f, const char *hist1, const char *hist2, const char *hist3, const char *label)
{

    TH1F *h1 = (TH1F*)f->Get(hist1)->Clone(hist1);
    TH1F *h2 = (TH1F*)f->Get(hist2)->Clone(hist2);
    TH1F *h3 = (TH1F*)f->Get(hist3)->Clone(hist3);
    h1->GetXaxis()->SetTitle(label);
    h1->SetLineWidth(2);
    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    h2->SetFillColor(kBlue);
    h2->SetFillStyle(0);
    h2->SetLineWidth(2);
    h3->SetLineColor(kBlue);
    h3->SetFillColor(kBlue);
    h3->SetLineWidth(2);
    h3->SetFillStyle(3001);

    TLegend *l1 = new TLegend(0.5, 0.5, 0.9, 0.9);
    l1->SetLineColor(kWhite);
    l1->SetFillColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(h1, "Data", "f");
    l1->AddEntry(h2, "DYee/mm 41X", "f");
    l1->AddEntry(h3, "DYee/mm 42X", "f");

    TCanvas *c1 = new TCanvas();
    c1->cd();
    h1->Draw("HIST");
    h3->Draw("SAME HIST");
    h2->Draw("SAME HIST");
    l1->Draw();
    c1->SetLogy(1);

    return c1;

}

TCanvas *ComparePlots(TFile *f, const char *hist1, const char *hist2, const char *label1, const char *label2, unsigned int rebin)
{

    // get hists
    TH1F *h1 = (TH1F*)f->Get(hist1)->Clone(hist1);
    TH1F *h2 = (TH1F*)f->Get(hist2)->Clone(hist2);

    // overflow
    h1->SetBinContent(h1->GetNbinsX(), h1->GetBinContent(h1->GetNbinsX()) + h1->GetBinContent(h1->GetNbinsX()+1));
    h1->SetBinContent(h1->GetNbinsX()+1, 0.0);
    h2->SetBinContent(h2->GetNbinsX(), h2->GetBinContent(h2->GetNbinsX()) + h2->GetBinContent(h2->GetNbinsX()+1));
    h2->SetBinContent(h2->GetNbinsX()+1, 0.0);

    // rebin
    h1->Rebin(rebin);
    h2->Rebin(rebin);

    // format
    h1->SetLineWidth(2);
    //h1->SetLineColor(kRed);
    //h2->SetLineColor(kBlue);
    //h2->SetFillColor(kBlue);
    h2->SetFillStyle(0);
    h2->SetLineWidth(2);
    //h2->SetMarkerColor(kBlue);

    // legend
    TLegend *l1 = new TLegend(0.70, 0.85, 0.98, 0.98);
    l1->SetLineColor(kWhite);
    l1->SetFillColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(h1, label1, "f");
    l1->AddEntry(h2, label2, "f");

    // draw
    TCanvas *c1 = new TCanvas();
    c1->cd();
    h1->Draw("HIST");
    h2->Draw("HIST SAME E1");
    double yMax = TMath::Max(h1->GetMaximum(), h2->GetMaximum());
    h1->GetYaxis()->SetRangeUser(0.01, yMax + 2*sqrt(yMax));

    l1->Draw();
    return c1;

}


TCanvas *ComparePlots(TFile *f1, TFile *f2, float lumi1, float lumi2, const char *hist, unsigned int rebin)
{

    // get hists
    TH1F *h1 = (TH1F*)f1->Get(hist)->Clone(hist);
    TH1F *h2 = (TH1F*)f2->Get(hist)->Clone(hist);

    // overflow
    h1->SetBinContent(h1->GetNbinsX(), h1->GetBinContent(h1->GetNbinsX()) + h1->GetBinContent(h1->GetNbinsX()+1));
    h1->SetBinContent(h1->GetNbinsX()+1, 0.0);
    h2->SetBinContent(h2->GetNbinsX(), h2->GetBinContent(h2->GetNbinsX()) + h2->GetBinContent(h2->GetNbinsX()+1));
    h2->SetBinContent(h2->GetNbinsX()+1, 0.0);

    h1->Scale(1.0/lumi1);
    h2->Scale(1.0/lumi2);

    // rebin
    h1->Rebin(rebin);
    h2->Rebin(rebin);

    // format
    h1->SetLineWidth(2);
    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    h2->SetFillColor(kBlue);
    h2->SetFillStyle(0);
    h2->SetLineWidth(2);
    h2->SetMarkerColor(kBlue);

    // legend
    TLegend *l1 = new TLegend(0.70, 0.85, 0.98, 0.98);
    l1->SetLineColor(kWhite);
    l1->SetFillColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(h1, "2011A", "f");
    l1->AddEntry(h2, "2011B", "f");

    // draw
    TCanvas *c1 = new TCanvas();
    c1->cd();
    h1->Draw("HIST");
    h2->Draw("SAME E1");
    double yMax = TMath::Max(h1->GetMaximum(), h2->GetMaximum());
    h1->GetYaxis()->SetRangeUser(0.01, yMax + 2*sqrt(yMax));

    l1->Draw();
    return c1;

}

TH1F* SmurfRebin(const TH1F *old, const unsigned int rebin)
{
  TH1F *h1_new = (TH1F*)old->Clone();
  TH1F *h1_old_tmp = (TH1F*)old->Clone("old_tmp");
  h1_old_tmp->Rebin(rebin);
  for (unsigned int i = 1; i <= h1_new->GetNbinsX(); ++i) {
      unsigned int bin = h1_old_tmp->FindBin(h1_new->GetBinCenter(i));
      h1_new->SetBinContent(i, h1_old_tmp->GetBinContent(bin) / float(rebin));
      h1_new->SetBinError  (i, h1_old_tmp->GetBinError  (bin) / float(rebin));
  }
  delete h1_old_tmp;
  return h1_new;
}

