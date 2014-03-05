#include "SmurfTable.h"

#include "SmurfScaleFactors.h"
#include "core/Enums.h"
#include "core/SmurfPlotUtilities.h"
#include <cmath>
#include <set>
#include <iostream>
#include <stdio.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "THStack.h"
#include "../Smurf/Analysis/HWWlvlv/PSUESystematics_8TeV.h"
#include "../Smurf/Analysis/HWWlvlv/HiggsQCDScaleSystematics_8TeV.h"
#include "../Smurf/Analysis/HWWlvlv/PDFgHHSystematics_8TeV.h"

//
// This function unrolls an n by m 2D histogram to an 1D (n x m bins) histogram
//
TH1F* UnrollHisto2DTo1D(TH2F* h2, const char* hname) {

    unsigned int nbinsX = h2->GetXaxis()->GetNbins();   
    unsigned int nbinsY = h2->GetYaxis()->GetNbins();   
    unsigned int nbins = nbinsX*nbinsY;     

    TString hnamestr = hname;
    hnamestr.ReplaceAll("_tmp", "");

    TH1F* histo = new TH1F(hnamestr.Data(), hnamestr.Data(), nbins, 0.5, 0.5+nbins);
    
    for(unsigned int x=1; x <= nbinsX; ++x) {
        for(unsigned int y=1; y <= nbinsY; ++y) {
            histo->SetBinContent( (x-1)*nbinsY + y, h2->GetBinContent(x, y) );
            histo->SetBinError( (x-1)*nbinsY + y, h2->GetBinError(x, y) );
        }
    }

    // storing entries is needed to calculate stat up bounding for empty bins 
    histo->SetEntries(h2->GetEntries());    

    return histo;
}

//
// print 2D templates : both central and alternate shapes  
//
void print2DShapeHistograms(std::vector<SmurfSample*> samples, Option option,
        unsigned int jetbin, float analysis, std::string cdir, unsigned int fcode, unsigned int runEra)
{

    // valid options
    if ( ! ((1<<option) & HWW_MT2DMLL ))  return;

    // hww for jet bins 0, 1, 2
    if (jetbin > 2) return;

    // -- WW PDF shape variation -------------------------
    TFile *weightPDFShapeFILE=0;   
    if(analysis < 300.) weightPDFShapeFILE = TFile::Open("/smurf/dlevans/PDFUncertainties/V00-00-01/PDFUncertainty_LowMass.root");
    else weightPDFShapeFILE = TFile::Open("/smurf/dlevans/PDFUncertainties/V00-00-01/PDFUncertainty_HighMass.root"); // TAS
    //if(analysis < 300.) weightPDFShapeFILE = TFile::Open("/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/aux/PDFUncertainty_LowMass.root"); 
    //else weightPDFShapeFILE = TFile::Open("/nfs-7/userdata/jaehyeok/smurfntuples/mitf-alljets/aux/PDFUncertainty_HighMass.root"); // UAF 
    TH2D *weightPDFShapeUp; 
    TH2D *weightPDFShapeDown; 
    // ---------------------------------------------------  

    // get the flavor
    std::string flavor;
    if (fcode == ((1<<0)|(1<<3))) {
        flavor = Form("sf");
    } else if (fcode == ((1<<1)|(1<<2))) {
        flavor = Form("of");
    } else if (fcode == ((1<<0)|(1<<1)|(1<<2)|(1<<3))) {
        flavor = Form("ll");
    } else if (fcode == (1<<0)) {
        flavor = Form("mm");
    } else if (fcode == (1<<1)) {
        flavor = Form("me");
    } else if (fcode == (1<<2)) {
        flavor = Form("em");
    } else if (fcode == (1<<3)) {
        flavor = Form("ee");
    } else {
        std::cout << "[SmurfTable::print2DShapeHistograms] Error! Must specify flavor of analysis" << std::endl;
        return;
    }

    std::string runera = "_8TeV"; 
    if ( runEra == RUN2012) 
        runera = "_8TeV";

    // create a file to save the shapes
    std::string filename = "";
    if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_MT2DMLL_JCP || option == HWW_OPT_SSCTL2D) { 
      filename = Form("%s/%i/hww%s_%ij.input%s.root", cdir.c_str(), int(analysis), flavor.c_str(), jetbin, runera.c_str());
    } 
    else if ( option == XWW_OPT_MT2DMLL_JCP ) { 
      filename = Form("%s/%i/xww%s_%ij.input%s.root", cdir.c_str(), int(analysis), flavor.c_str(), jetbin, runera.c_str());
    }
    else {
        std::cout << "[SmurfTable::print2DShapeHistograms] Error! Couldn't set output filename" << std::endl;
        return;
    }

    std::cout << "[SmurfTable::print2DShapeHistograms] " << filename << std::endl;
    TFile *f = new TFile(filename.c_str(), "RECREATE");
    f->cd();

    // save the shapes
    TH1F *temp = 0;         // unrolled 1D  ==> central 1D shape
    TH2F *temp_roll = 0;    // rolled 2D    ==> central 2D shape 
    TH2F *temp2d = 0;       // tempory for shape variation 

    TH1F *temp_dyll = 0;
    TH1F *temp_dyll_loosemet = 0;
    TH1F *temp_zjets_up = 0;
    TH1F *temp_zjets_down = 0;
    TH1F *temp_dyll_data_of = 0;
    TH1F *temp_zjets = 0;

    // - Top
    TH1F *temp_top_up = 0;
    TH1F *temp_top_down = 0;

    // - WW
    TH1F *temp_ww_up = 0;
    TH1F *temp_ww_down = 0;
    TH1F *temp_ww_mcnlo = 0;
    TH1F *temp_ww_mcnlo_up = 0;
    TH1F *temp_ww_mcnlo_down = 0;

    TH1F *temp_qqwwpdf_up = 0;
    TH1F *temp_qqwwpdf_down = 0;
    TH1F *temp_ggwwpdf_up = 0;
    TH1F *temp_ggwwpdf_down = 0;
    
    // - wjets
    TH1F *temp_wjetsE = 0;
    TH1F *temp_wjetsM = 0;
    //TH1F *temp_wjets_mc = 0;
    //TH1F *temp_wjets_mc_up = 0;
    //TH1F *temp_wjets_mc_down = 0;

    // - wgamma
    TH1F *temp_wgammanorm   = 0;

    // - qqhww
    TH1F *temp_gghww_ref = 0;
    TH1F *temp_gghww_jhu = 0;
    TH1F *temp_gghww_up = 0;
    TH1F *temp_gghww_down = 0;
    TH1F *temp_ggxww_ref = 0;

    TH2F *temp2d_gghww_ref = 0;
    TH2F *temp2d_gghww_jhu = 0;
    
    
    for (unsigned int s = 0; s < samples.size(); ++s) 
    {
        // dy shapes
        if (samples[s]->getDataType() == ZLL) {
            //const char *tmp_title     = ((TH2F*)samples[s]->get2DMVAShape(jetbin, fcode))->GetName();
            temp2d      = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone(Form("histo_%s_tmp", samples[s]->getName().c_str()));
            temp_dyll   = UnrollHisto2DTo1D(temp2d,temp2d->GetName()); 
        }
/*
        if (samples[s]->getDataType() == ZLLLOOSEMET) {
            //const char *tmp_title     = ((TH2F*)samples[s]->get2DMVAShape(jetbin, (1<<0)|(1<<3)))->GetName();
            temp2d              = (TH2F*)samples[s]->get2DMVAShape(jetbin, (1<<0)|(1<<3))->Clone(Form("histo_%s_tmp", samples[s]->getName().c_str()));
            temp_dyll_loosemet  = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
        }

        if (samples[s]->getDataType() == ZJETS) {
            //temp_zjets_up = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone(Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingUp", flavor.c_str(), jetbin));
            //temp_zjets_down = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone(Form("histo_Zjets_CMS_hww%s_%ij_MVAZBoundingDown",flavor.c_str(), jetbin));
            temp2d      = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone(Form("temp_Zjets_%ij_tmp", jetbin));
            temp_zjets  = UnrollHisto2DTo1D(temp2d,temp2d->GetName());   
        }

        if (samples[s]->getDataType() == ZLLDATA) {
            temp2d          = (TH2F*)samples[s]->get2DMVAShape(jetbin, (1<<0)|(1<<3))->Clone(Form("histo_Zjets_CMS_hwwsf_%ij_MVAZBoundingUp_tmp", jetbin));
            temp_zjets_up   = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
            temp2d          = (TH2F*)samples[s]->get2DMVAShape(jetbin, (1<<0)|(1<<3))->Clone(Form("histo_Zjets_CMS_hwwsf_%ij_MVAZBoundingDown_tmp", jetbin));
            temp_zjets_down = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
            temp2d              = (TH2F*)samples[s]->get2DMVAShape(jetbin, (1<<1)|(1<<2))->Clone(Form("histo_Zjets_data_of_%ij_tmp", jetbin)); 
            temp_dyll_data_of   = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
            temp_zjets_up->Add(temp_zjets_up, temp_dyll_data_of, -1.);
        }
*/

        if (samples[s]->getDataType() == TOPALTER) {
            temp2d      = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_Top_CMS_hww_MVATopBoundingUp_tmp");
            temp_top_up = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
            temp2d          = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_Top_CMS_hww_MVATopBoundingDown_tmp");
            temp_top_down   = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
        }

        // ww alternative shapes
        if (samples[s]->getDataType() == WWMCNLO) {
            temp2d          = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_WW_CMS_MVAWWNLO_tmp");
            temp_ww_mcnlo   =  UnrollHisto2DTo1D(temp2d,temp2d->GetName());
            temp2d          = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_qqWW_CMS_hww_MVAWWBoundingUp_tmp");
            temp_ww_up      =  UnrollHisto2DTo1D(temp2d,temp2d->GetName());
            temp2d          = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_qqWW_CMS_hww_MVAWWBoundingDown_tmp");
            temp_ww_down    =  UnrollHisto2DTo1D(temp2d,temp2d->GetName());
        }
        if (samples[s]->getDataType() == WWMCNLOUP) {
            temp2d              = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_qqWW_CMS_hww_MVAWWNLOBoundingUp_tmp");
            temp_ww_mcnlo_up    =  UnrollHisto2DTo1D(temp2d,temp2d->GetName());
        }
        if (samples[s]->getDataType() == WWMCNLODOWN) {
            temp2d              = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_qqWW_CMS_hww_MVAWWNLOBoundingDown_tmp");
            temp_ww_mcnlo_down  =  UnrollHisto2DTo1D(temp2d,temp2d->GetName());
        }
//        if (samples[s]->getDataType() == WJETSMCLOOSE) {
//            temp2d              = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_Wjets_CMS_hww_MVAWMCBoundingUp_tmp");
//            temp_wjets_mc_up    = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
//            temp2d              = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_Wjets_CMS_hww_MVAWMCBoundingDown_tmp");
//            temp_wjets_mc_down  = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
//        }
//        if (samples[s]->getDataType() == WJETS) {
//            temp2d          = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_WjetsMC_tmp"); 
//            temp_wjets_mc   = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
//        }
        if (samples[s]->getDataType() == WJETSELEDATA ) {
            temp2d          = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_WjetsE_temp_tmp");
            temp_wjetsE     = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
        }
        if (samples[s]->getDataType() == WJETSMUDATA ) {
            temp2d          = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_WjetsM_temp_tmp");
            temp_wjetsM     = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
        }
        if (samples[s]->getDataType() == WGAMMANORM ) { 
            temp2d              = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone(Form("histo_%s_tmp", samples[s]->getName().c_str()));
            temp_wgammanorm     = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
        }
        if ( option == XWW_OPT_MT2DMLL_JCP ) {
          // qqhww alternative shapes 
          if (samples[s]->getDataType() == GGHWWREF) {
            temp2d_gghww_ref    = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_ggHRef_tmp");
            temp_gghww_ref      = UnrollHisto2DTo1D(temp2d_gghww_ref,temp2d->GetName());
          }
          if (samples[s]->getDataType() == GGHWWJHU) {
            temp2d_gghww_jhu    = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_ggHJHU_tmp");
            temp_gghww_jhu      = UnrollHisto2DTo1D(temp2d_gghww_jhu,temp2d->GetName());
          }
          if (samples[s]->getDataType() == GGHWW) {
            temp2d          = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_ggH_CMS_hww_MVAggHBoundingUp_tmp");
            temp_gghww_up   = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
            temp2d          = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_ggH_CMS_hww_MVAggHBoundingDown_tmp");
            temp_gghww_down = UnrollHisto2DTo1D(temp2d,temp2d->GetName());      
            temp2d          = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_ggXRef");
            temp_ggxww_ref      = UnrollHisto2DTo1D(temp2d,temp2d->GetName());   
          }
        }
        // WW PDF shape variations
        if (samples[s]->getDataType() == QQWW) { 
            weightPDFShapeUp   = (TH2D*)(weightPDFShapeFILE->Get( Form("qqWW_DF_%ij_alternateUp", jetbin) ));
            weightPDFShapeDown = (TH2D*)(weightPDFShapeFILE->Get( Form("qqWW_DF_%ij_alternateDown", jetbin) ));

            temp2d              = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_qqWW_CMS_hww_PDFqqWWUp_tmp"); temp2d->Multiply(weightPDFShapeUp);
            temp_qqwwpdf_up     = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
            temp2d              = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_qqWW_CMS_hww_PDFqqWWDown_tmp"); temp2d->Multiply(weightPDFShapeDown);
            temp_qqwwpdf_down   = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
        }
        if (samples[s]->getDataType() == GGWW) { 
            weightPDFShapeUp   = (TH2D*)(weightPDFShapeFILE->Get( Form("ggWW_DF_%ij_alternateUp", jetbin) ));
            weightPDFShapeDown = (TH2D*)(weightPDFShapeFILE->Get( Form("ggWW_DF_%ij_alternateDown", jetbin) ));
            
            temp2d              = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_ggWW_CMS_hww_PDFggWWUp_tmp"); temp2d->Multiply(weightPDFShapeUp);
            temp_ggwwpdf_up     = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
            temp2d              = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_ggWW_CMS_hww_PDFggWWDown_tmp"); temp2d->Multiply(weightPDFShapeDown);
            temp_ggwwpdf_down   = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
        }
    }

    for (unsigned int s = 0; s < samples.size(); ++s) 
    {

        if (samples[s]->getDataType() == ZLL) continue;
//        if (samples[s]->getDataType() == ZLLLOOSEMET) continue;
//        if (samples[s]->getDataType() == ZLLDATA) continue;
        if (samples[s]->getDataType() == TOPDATA) continue;
        if (samples[s]->getDataType() == WWMCNLO) continue;
        if (samples[s]->getDataType() == WWMCNLOUP) continue;
        if (samples[s]->getDataType() == WWMCNLODOWN) continue;
        if (samples[s]->getDataType() == WZALTER) continue;
        if (samples[s]->getDataType() == ZZALTER) continue;
        if (samples[s]->getDataType() == TOPALTER) continue;
        if (samples[s]->getDataType() == GGHWWREF) continue;
        if (samples[s]->getDataType() == GGHWWJHU) continue;
        if (samples[s]->getDataType() == WGAMMANORM) continue;  
        

        // get the shape histogram
        // if it is the OF data get the keys versions
        /*
        if (samples[s]->getDataType() == OFDATA) {
            temp_roll = (TH2F*)samples[s]->getKeysMVAShape(jetbin, fcode);
            TH2F *temp_nokeys = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode);
            temp_roll->SetName("histo_WWTop");
            temp_nokeys->SetName("hist_WWTopNoKeys");
            temp_nokeys->SetDirectory(0);
            temp_nokeys->Write();
        }
       */ 
        // otherwise get the usual version
        //else {
        temp_roll = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode);   
        int tmpEntries = temp_roll->GetEntries(); 
        temp = UnrollHisto2DTo1D(temp_roll, Form("histo_%s_tmp", samples[s]->getName().c_str())); 
        temp->SetEntries(tmpEntries); 
        //}

        // Wgamma 
        // renormalize WgammaFO to Wgamma
        // Wgamma means Wgamma + Zgamma
        if (samples[s]->getDataType() == WGAMMA ) {  
            if( temp->Integral() != 0 ) temp->Scale( temp_wgammanorm->Integral()/temp->Integral() );                    
        }
/*
        // Z+Jets for HWW analysis
        // combine ZJETS and ZLL
        if (samples[s]->getDataType() == ZJETS && fcode == ((1<<0)|(1<<3))) {
            // get the z yield from the data-driven estimate
            double zSF[3] = {0.0, 0.0, 0.0};
            double zSFErr[3] = {0.0, 0.0, 0.0};
            if(option != HWW_OPT_SSCTL2D) getZScaleFactor(zSF, zSFErr, option, analysis, flavor.c_str());
            double yield_DYll;
            if  ( analysis <= 300. && option != HWW_OPT_SSCTL2D) 
                yield_DYll = zSF[jetbin]; 
            else 
                yield_DYll = temp_dyll->Integral(0,1000);
            temp_dyll_loosemet->Scale( yield_DYll /temp_dyll_loosemet->Integral(0,1000));
            temp->Add(temp, temp_dyll_loosemet);

            // variation up is the DYdata + VZ (MC)
            // add a protection against the empty bins
            double yield_Zjets = temp_zjets->Integral(0, 1000); 
            if(option == HWW_OPT_SSCTL2D) yield_Zjets = 0.;
            temp_zjets_up->Scale( yield_DYll / temp_zjets_up->Integral(0, 1000));
            temp_zjets_up->Add(temp_zjets_up, temp_zjets);
            for (int binx = 1; binx < temp_zjets_up->GetNbinsX()+1; binx++) {
                if ( temp_zjets_up->GetBinContent(binx) < 0. ) {
                    temp_zjets_up->SetBinContent(binx, 0.);
                }
            }
            temp_zjets_up->Scale( (yield_DYll+yield_Zjets) / temp_zjets_up->Integral(0, 1000));

            // write the down histogram mirroring the up-central change
            for (int binx = 1; binx < temp_zjets_down->GetNbinsX()+1; binx++) {
                double downcontent = temp->GetBinContent(binx) - ( temp_zjets_up->GetBinContent(binx) - temp->GetBinContent(binx));
                if (downcontent < 0.) 
                    downcontent = 0.0;
                temp_zjets_down->SetBinContent(binx, downcontent);
            }
            temp_zjets_down->Scale(temp->Integral(0,1000)/temp_zjets_down->Integral(0,1000));

            temp_zjets_up->SetDirectory(0);
            temp_zjets_up->Write();
            temp_zjets_down->SetDirectory(0);
            temp_zjets_down->Write();
            temp_dyll_data_of->SetDirectory(0);
            temp_dyll_data_of->Write();
        }
*/
        // 
        // Top
        // 
        if (samples[s]->getDataType() == TOP && temp_top_up != 0 ) {
            float top_yield = temp->Integral(0, temp->GetNbinsX()+1);
            float top_yield_up = temp_top_up->Integral(0, temp_top_up->GetNbinsX()+1);
            temp_top_up->Scale( top_yield / top_yield_up);

            for (int binx = 1; binx < temp_top_down->GetNbinsX() + 1; binx++) {
                double downcontent = 2*temp->GetBinContent(binx) - temp_top_up->GetBinContent(binx);
                if (downcontent < 0.) 
                    downcontent = 0.0;
                temp_top_down->SetBinContent(binx, downcontent);
            }

            temp_top_up->SetDirectory(0);
            temp_top_up->Write();
            temp_top_down->SetDirectory(0);
            temp_top_down->Write();
        }

        //
        // WW
        // 
        if ( samples[s]->getDataType() == QQWW && temp_ww_up != 0) {

            float ww_yield = temp->Integral(0,100);
            float ww_yield_mcnlo = temp_ww_up->Integral(0, 1000); 

            // deal with the central shape and MC@NLO difference
            if( ww_yield_mcnlo != 0 ) temp_ww_up->Scale( ww_yield / ww_yield_mcnlo);

            for ( int binx = 1; binx < temp_ww_up->GetNbinsX()+1; binx++) {
                double downcontent = 2*temp->GetBinContent(binx) - temp_ww_up->GetBinContent(binx);
                if ( downcontent <= 0.) downcontent = 0.0;
                temp_ww_down->SetBinContent(binx, downcontent);
            }

            // if the down histogram flucturates to 0, set the down = central shape
            if (temp_ww_down->Integral(0,1000) <= 0.0) {
                for ( int binx = 1; binx < temp_ww_up->GetNbinsX()+1; binx++) {
                    temp_ww_down->SetBinContent(binx, temp->GetBinContent(binx));
                    temp_ww_down->SetBinError(binx, temp->GetBinError(binx));
                }
            }

            // normalize to the central shape 
            temp_ww_up->Scale(temp->Integral(0,1000) / temp_ww_up->Integral(0,1000));
            temp_ww_down->Scale(temp->Integral(0,1000) / temp_ww_down->Integral(0,1000));

            temp_ww_up->SetDirectory(0);
            temp_ww_up->Write();

            temp_ww_down->SetDirectory(0);
            temp_ww_down->Write();

            // deal with the QCD scale variations
            temp_ww_mcnlo->Scale(temp->Integral(0,1000) / temp_ww_mcnlo->Integral(0,1000));
            temp_ww_mcnlo_up->Scale(temp->Integral(0,1000) / temp_ww_mcnlo_up->Integral(0,1000));
            temp_ww_mcnlo_down->Scale(temp->Integral(0,1000) / temp_ww_mcnlo_down->Integral(0,1000));

            for ( int binx = 1; binx < temp_ww_up->GetNbinsX()+1; binx++) {
                double upratio (1.0), downratio (1.0);
                if ( temp_ww_mcnlo->GetBinContent(binx) > 0 &&  temp_ww_mcnlo->GetBinError(binx) / temp_ww_mcnlo->GetBinContent(binx) < 0.5 ) {
                    upratio = temp_ww_mcnlo_up->GetBinContent(binx) / temp_ww_mcnlo->GetBinContent(binx);
                    downratio = temp_ww_mcnlo_down->GetBinContent(binx) / temp_ww_mcnlo->GetBinContent(binx);
                }
                temp_ww_mcnlo_up->SetBinContent(binx, upratio * temp->GetBinContent(binx));
                temp_ww_mcnlo_down->SetBinContent(binx, downratio * temp->GetBinContent(binx));
            }

            // normalize to the central shape 
            temp_ww_mcnlo_up->Scale(temp->Integral(0,1000) / temp_ww_mcnlo_up->Integral(0,1000));
            temp_ww_mcnlo_down->Scale(temp->Integral(0,1000) / temp_ww_mcnlo_down->Integral(0,1000));
            
            temp_ww_mcnlo_up->SetDirectory(0);
            temp_ww_mcnlo_up->Write();

            temp_ww_mcnlo_down->SetDirectory(0);
            temp_ww_mcnlo_down->Write(); 

        }


        // 
        // Wjets
        //

        if ( samples[s]->getDataType() == WJETS) continue;
        if ( samples[s]->getDataType() == WJETSMCLOOSE) continue;
        if ( samples[s]->getDataType() == WJETSELEDATA ) {

            // first do something special of the central histogram to avoid negative bin content
            float yield_wjets = temp_wjetsE->Integral(0, 1000);
            for ( int binx = 1; binx < temp->GetNbinsX()+1; binx++) {
                float bincontent = temp->GetBinContent(binx);
                if ( bincontent < 0.) {
                    temp->SetBinContent(binx, 0.);
                }
            }
            temp->Scale(temp->Integral(0, 1000) ? yield_wjets / temp->Integral(0, 1000) : 1);
/*
            // Wjets MC closure test from shape analysis note
            temp_wjets_mc->Scale( yield_wjets / temp_wjets_mc->Integral(0,1000));
            temp_wjets_mc_up->Scale( yield_wjets / temp_wjets_mc_up->Integral(0,1000));
            for ( int binx = 1; binx < temp_wjets_mc_down->GetNbinsX()+1; binx++) { 
                double downcontent = 2*temp->GetBinContent(binx) - temp_wjets_mc_up->GetBinContent(binx);
                if ( downcontent <= 0.) downcontent = 0.0;
                temp_wjets_mc_down->SetBinContent(binx, downcontent);
            }

            temp_wjets_mc_up->SetDirectory(0);
            temp_wjets_mc_up->Write();

            temp_wjets_mc_down->SetDirectory(0);
            temp_wjets_mc_down->Write();
*/
        }
        
        if ( samples[s]->getDataType() == WJETSMUDATA ) {

            // first do something special of the central histogram to avoid negative bin content
            float yield_wjets = temp_wjetsM->Integral(0, 1000);
            for ( int binx = 1; binx < temp->GetNbinsX()+1; binx++) {
                float bincontent = temp->GetBinContent(binx);
                if ( bincontent < 0.) {
                    temp->SetBinContent(binx, 0.);
                }
            }
            temp->Scale(temp->Integral(0, 1000) ? yield_wjets / temp->Integral(0, 1000) : 1);
        }

        // 
        // ggHWW
        // 
        if ( option == XWW_OPT_MT2DMLL_JCP && samples[s]->getDataType() == GGHWW ) {

          // temp_gghww_ref is the SM Higgs shape from powheg
          // temp_gghww_jhu is the sm Higgs shape from jhu
          // temp_gghww_up and temp are spin 2 out-of-box histograms
          
          // get the spin 0 reference yields and scale all the jhu samples to this yield
          float yield_gghww = temp_gghww_ref->Integral(0,1000.);
          float smto2mplusweight = yield_gghww / temp->Integral(0,1000); 
          temp->Scale( smto2mplusweight );
          temp_gghww_up->Scale( smto2mplusweight );
          temp_gghww_jhu->Scale(  yield_gghww / temp_gghww_jhu->Integral(0,1000) );
          temp_roll->Scale( yield_gghww / temp_roll->Integral(0, 1000) );
          
          // calculate the ratio between the jhu and the reference for SM Higgs
          // and apply the weight to the jhu 
          for ( int binx = 1; binx < temp_gghww_ref->GetNbinsX()+1; binx++) { 
            double refcontent = temp_gghww_ref->GetBinContent(binx);
            double jhucontent = temp_gghww_jhu->GetBinContent(binx);
            double reftojhuratio = jhucontent > 0. ? refcontent/jhucontent : 1.; 
            temp_gghww_up->SetBinContent(binx, temp->GetBinContent(binx) * reftojhuratio);      
            temp_gghww_up->SetBinError(binx, temp->GetBinError(binx) * reftojhuratio);       
          }
          
          // calculate down flucturations based on the up and central
          for ( int binx = 1; binx < temp->GetNbinsX()+1; binx++) { 
            double downcontent = 2*temp->GetBinContent(binx) - temp_gghww_up->GetBinContent(binx);
            if ( downcontent <= 0.) downcontent = 0.0;
            temp_gghww_down->SetBinContent(binx, downcontent);
          }
          
          temp_gghww_up->Scale(yield_gghww/temp_gghww_up->Integral(0,1000));
          temp_gghww_down->Scale(yield_gghww/temp_gghww_up->Integral(0,1000));

          temp_gghww_up->SetDirectory(0);
          temp_gghww_up->Write();

          temp_gghww_down->SetDirectory(0);
          temp_gghww_down->Write();
          
        }

        // PDF Shape variations for ggWW  
        if ( samples[s]->getDataType() == GGWW ) {
            temp_ggwwpdf_up->SetDirectory(0);
            temp_ggwwpdf_up->Write();
            temp_ggwwpdf_down->SetDirectory(0);
            temp_ggwwpdf_down->Write();
        }

        // PDF Shape variations for qqWW  
        if ( samples[s]->getDataType() == QQWW ) {
            temp_qqwwpdf_up->SetDirectory(0);
            temp_qqwwpdf_up->Write();
            temp_qqwwpdf_down->SetDirectory(0);
            temp_qqwwpdf_down->Write();
        }


        // add the bin-by-bin stat.
        TH1F *temp_stat_up = 0;
        TH1F *temp_stat_down = 0;

        // -- Bin-by-bin stat.
        if ( ! ((1ll<<samples[s]->getDataType()) & ((1ll<<DATA) /*| (1ll<<ZLL) | (1ll<<ZLLLOOSEMET)*/ | (1ll<<TOPDATA)
                        | (1ll<<WWMCNLO) | (1ll<<WWMCNLOUP) | (1ll<<WWMCNLODOWN) 
                        | (1ll<<WZALTER) | (1ll<<ZZALTER) | (1ll<<WWTOPMC)
                        | (1ll<<WJETSMCLOOSE) | (1ll<<WJETS) | (1ll<<TOPALTER)) ) ) {

            temp2d = (TH2F*)temp->Clone(Form("histo_%s_CMS_hww%s_%ij_MVA%sStatBounding%sUp_tmp", samples[s]->getName().c_str(),flavor.c_str(), jetbin, samples[s]->getName().c_str(),runera.c_str()));
            temp_stat_up = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
            temp2d = (TH2F*)temp->Clone(Form("histo_%s_CMS_hww%s_%ij_MVA%sStatBounding%sDown_tmp", samples[s]->getName().c_str(),flavor.c_str(), jetbin, samples[s]->getName().c_str(), runera.c_str()));
            temp_stat_down = UnrollHisto2DTo1D(temp2d,temp2d->GetName());

            for ( int binx = 1; binx < temp_stat_up->GetNbinsX()+1; binx++) {
                double upcontent = temp->GetBinContent(binx) + temp->GetBinError(binx);
                // flucturate up one event weight if the histogram is empty
                if ( upcontent == 0 && temp->GetEntries() > 0) 
                    //upcontent = temp->Integral(0,1000) / temp->GetEntries(); 
                    upcontent = 0.000001; // to be consistent with Guillelmo 
                double downcontent = temp->GetBinContent(binx) - temp->GetBinError(binx); 
                //if ( downcontent <= 0.0000001) downcontent = 0.0;
                if ( downcontent <= 0.0000001) downcontent = 0.000001; // to be consistent with Guillelmo 
                temp_stat_up->SetBinContent(binx, upcontent);
                temp_stat_down->SetBinContent(binx, downcontent);
            } 

            // if the down histogram flucturates to 0, set the down = central shape
            if (temp_stat_down->Integral(0,1000) <= 0.0) {
                for ( int binx = 1; binx < temp_stat_up->GetNbinsX()+1; binx++) {
                    temp_stat_down->SetBinContent(binx, temp->GetBinContent(binx));
                    temp_stat_down->SetBinError(binx, temp->GetBinError(binx));
                }
            }     

            temp_stat_up->SetDirectory(0);
            temp_stat_up->Write();

            temp_stat_down->SetDirectory(0);
            temp_stat_down->Write();
        }
        
        
        temp->SetDirectory(0);
        temp->Write();
        temp->SetDirectory(0);
        temp_roll->Write();

        //
        // write out any alternate shapes associated with 
        // this sample
        //

        std::set<ShapeSyst> shapeSystematics = samples[s]->getAvailableShapeSystematics();
        std::set<ShapeSyst>::const_iterator var;

        for (var = shapeSystematics.begin(); var != shapeSystematics.end(); ++var) {

            temp2d = (TH2F*)samples[s]->getShapeVariation2D(*var, true, jetbin, fcode)->Clone(Form("%s_tmp",temp2d->GetTitle()));
            TH1F *temp_up = UnrollHisto2DTo1D(temp2d,temp2d->GetTitle()); 

            std::string histname = Form("%s%s", temp_up->GetTitle(), flavor.c_str());

            // the scale variations that are 100% correlated between flavors
            // so do not include flavor in name

            if ( (1ll<<(*var)) & ( (1ll<<QCDSCALEVAR) | (1ll<<LEPEFFVAR) | (1ll<<WJETSELESHAPEVAR) | (1ll<<WJETSMUSHAPEVAR) | (1ll<<METVAR) ) | (1ll<<LEPRESVAR) ) {
                histname = Form("%s", temp_up->GetTitle());
            }
            temp_up->SetName(Form("%sUp", histname.c_str()));

            temp2d = (TH2F*)samples[s]->getShapeVariation2D(*var, false, jetbin, fcode)->Clone(Form("%s_tmp",temp2d->GetName()));
            TH1F *temp_down = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
            temp_down->SetName(Form("%sDown", histname.c_str()));

        
            if ( option == XWW_OPT_MT2DMLL_JCP && (samples[s]->getDataType() == GGHWW) ) { 
              // special treatment for the spin 2 shapes 
              float yield_gghww = temp_gghww_ref->Integral(0,1000.);
              float yield_ggxww = temp_ggxww_ref->Integral(0,1000.);
              float smto2mplusweight = yield_gghww / yield_ggxww;
              temp_up->Scale( smto2mplusweight );
              temp_down->Scale( smto2mplusweight );

            }

            // do a gymnastic on the WBouding histograms
            // (the same thing is applied to the met smearing)
            // The histo_Wjets_CMS_MVAWBoundingUp and histo_Wjets_CMS_MVAWBoundingDown are filled the same way in SmurfLooper
            // the histo_Wjets_CMS_MVAWBoundingDown should be the mirrored changed of central and histo_Wjets_CMS_MVAWBoundingUp

            if ( (*var) == WJETSELESHAPEVAR || (*var) == WJETSMUSHAPEVAR ||  (*var) == METVAR ) {
                for ( int binx = 1; binx < temp_down->GetNbinsX()+1; binx++) {
                    double downcontent = 2*temp->GetBinContent(binx) - temp_up->GetBinContent(binx); 
                    if (downcontent < 0.) downcontent = 0.0;
                    temp_down->SetBinContent(binx, downcontent);
                }
            }

            temp_up->SetDirectory(0);
            temp_down->SetDirectory(0);
            temp_up->Write(); 
            temp_down->Write();

        }

    }

    f->Close();
    delete f;

}

//
// Print Table of yields
//
void printResultsTable(std::vector<SmurfSample*> samples, Option option, bool doJetBins)
{

    printf("\n\n******************\n");
    if ( option == HWW_OPT_SMURFCUTSEL ) 
        printf("Table for Higgs Cut-Based selection\n");
    if ( option == HWW_OPT_MT2DMLL || option == HWW_OPT_MT2DMLL_JCP || option == XWW_OPT_MT2DMLL_JCP)
        printf("Table for Higgs 2D selection\n");
    if ( option == HWW_OPT_SSCTL || option == HWW_OPT_SSCTL2D)
        printf("Table for Higgs SS control region selection\n");
    if ( option == HWW_OPT_TOPTAG)
        printf("Table for Higgs Top-tagged control region selection\n");
    printf("******************\n");

    for (unsigned int jetbin = 0; jetbin < 3; ++jetbin)
    { 

        unsigned int binMinIdx = jetbin+1;
        unsigned int binMaxIdx = jetbin+1;
        if (!doJetBins) {
            binMinIdx = 0;
            binMaxIdx = 999;
            if (jetbin != 0) continue;
        } else {
            printf("\njetbin = %i\n", jetbin);
        }

        printf("\\begin{table}[!ht]\n");
        printf("\\begin{center}\n");
        printf("\\begin{tabular}{c|c|c|c|c|c}\n");
        printf("sample \t& mm \t& me \t& em \t& ee \t & TOTAL\\\\ \\hline \n");

        double s_yield[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        double s_err2[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        double data_yield[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        double data_err2[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        double b_yield[5] = {0.0, 0.0, 0.0, 0.0, 0.0}; 
        double b_err2[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

        //          
        // signal and background yields in ll
        //      

        for (unsigned int s = 0; s < samples.size(); ++s) {

            double yield[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
            double err[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

            if (samples[s]->getDataType() == WGAMMA) continue;
            for (unsigned int i = 0; i < 4; ++i) {

                // get results
                samples[s]->getResults(binMinIdx, binMaxIdx, i, yield[i], err[i]);

                if (samples[s]->getDataType() == DATA) {
                    data_yield[i] += yield[i];
                    data_err2[i] += err[i]*err[i];
                    data_yield[4] += yield[i];
                    data_err2[4] += err[i]*err[i];
                }
                else if ( 1ll<<samples[s]->getDataType() & data_higgsww ) {
                    s_yield[i] += yield[i];
                    s_err2[i] += err[i]*err[i];
                    s_yield[4] += yield[i];
                    s_err2[4] += err[i]*err[i];
                } else if ((1ll<<samples[s]->getDataType()) & data_allbg || samples[s]->getDataType() == OFDATA){
                    b_yield[i] += yield[i];
                    b_err2[i] += err[i]*err[i];
                    b_yield[4] += yield[i];
                    b_err2[4] += err[i]*err[i];
                }


            }

            printf("%s", (samples[s]->getName()).c_str());
            printf("\t& %4.2f $\\pm$ %4.2f", yield[0], err[0]);
            printf("\t& %4.2f $\\pm$ %4.2f", yield[1], err[1]);
            printf("\t& %4.2f $\\pm$ %4.2f", yield[2], err[2]);
            printf("\t& %4.2f $\\pm$ %4.2f", yield[3], err[3]);
            printf("\t& %4.2f $\\pm$ %4.2f \\\\ \\hline \n", yield[0]+yield[1]+yield[2]+yield[3], 
                    sqrt(pow(err[0], 2) + pow(err[1], 2) + pow(err[2], 2) + pow(err[3], 2)) );
        }

        printf("SIGNAL");
        printf("\t& %4.2f $\\pm$ %4.2f", s_yield[0], sqrt(s_err2[0]));
        printf("\t& %4.2f $\\pm$ %4.2f", s_yield[1], sqrt(s_err2[1]));
        printf("\t& %4.2f $\\pm$ %4.2f", s_yield[2], sqrt(s_err2[2]));
        printf("\t& %4.2f $\\pm$ %4.2f", s_yield[3], sqrt(s_err2[3]));
        printf("\t& %4.2f $\\pm$ %4.2f \\\\ \\hline \n", s_yield[0] + s_yield[1] + s_yield[2] + s_yield[3], 
                sqrt(s_err2[0] + s_err2[1] + s_err2[2] + s_err2[3]) );
        printf("BKGD");
        printf("\t& %4.2f $\\pm$ %4.2f", b_yield[0], sqrt(b_err2[0]));
        printf("\t& %4.2f $\\pm$ %4.2f", b_yield[1], sqrt(b_err2[1]));
        printf("\t& %4.2f $\\pm$ %4.2f", b_yield[2], sqrt(b_err2[2]));
        printf("\t& %4.2f $\\pm$ %4.2f", b_yield[3], sqrt(b_err2[3]));
        printf("\t& %4.2f $\\pm$ %4.2f \\\\ \\hline \n", b_yield[0] + b_yield[1] + b_yield[2] + b_yield[3], 
                sqrt(b_err2[0] + b_err2[1] + b_err2[2] + b_err2[3]) );
        printf("DATA");
        printf("\t& %4.2f $\\pm$ %4.2f", data_yield[0], sqrt(data_err2[0]));
        printf("\t& %4.2f $\\pm$ %4.2f", data_yield[1], sqrt(data_err2[1]));
        printf("\t& %4.2f $\\pm$ %4.2f", data_yield[2], sqrt(data_err2[2]));
        printf("\t& %4.2f $\\pm$ %4.2f", data_yield[3], sqrt(data_err2[3]));
        printf("\t& %4.2f $\\pm$ %4.2f \\\\ \\hline \n", data_yield[0] + data_yield[1] + data_yield[2] + data_yield[3], 
                sqrt(data_err2[0] + data_err2[1] + data_err2[2] + data_err2[3]) );

        printf("\\end{tabular}\n");
        printf("\\caption{Hello, I'm a table}\n");
        printf("\\label{tab:yield}\n");
        printf("\\end{center}\n");
        printf("\\end{table}\n");

    }

}

//
// Print cards 
//
void printCard(std::vector<SmurfSample*> samples, Option option, unsigned int jetbin, float analysis, 
        std::string cdir, unsigned int fcode, unsigned int mva_option, unsigned int runEra) 
{


    //
    // get the yields
    //

    // set up variables to get results
    // from histogram binned in njets
    unsigned int binMinIdx = jetbin+1;
    unsigned int binMaxIdx = jetbin+1;


    // find the index of each sample
    // in the vector of samples
    unsigned int i_data, i_ZH, i_WH, i_qqH, i_ggH, i_ZZ, i_WZ, i_ggWW, 
                 i_qqWW, i_Top, i_Zjets, i_DYtt, i_WjetsE, i_WjetsM, i_VV, i_DYll, 
                 i_Wgamma, i_WgammaNorm, i_Wg3l;

    // and record the yield and error
    // for that sample
    std::vector<double> yield;
    yield.reserve(samples.size());
    std::vector<double> err;
    err.reserve(samples.size());


    for (unsigned int s = 0; s < samples.size(); ++s) {

        // get the yield and error
        if (fcode == ((1<<0)|(1<<3))) {
            double y_ee, y_mm;
            double e_ee, e_mm;
            samples[s]->getResults(binMinIdx, binMaxIdx, 0, y_mm, e_mm);
            samples[s]->getResults(binMinIdx, binMaxIdx, 3, y_ee, e_ee);
            if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D) samples[s]->get2DResults(binMinIdx, binMaxIdx, 0, y_mm, e_mm);
            if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D) samples[s]->get2DResults(binMinIdx, binMaxIdx, 3, y_ee, e_ee);
            yield[s] = y_ee + y_mm;
            err[s] = sqrt(e_ee*e_ee + e_mm*e_mm);
        } else if (fcode == ((1<<1)|(1<<2))) {
            double y_em, y_me;
            double e_em, e_me;
            samples[s]->getResults(binMinIdx, binMaxIdx, 1, y_em, e_em);
            samples[s]->getResults(binMinIdx, binMaxIdx, 2, y_me, e_me);
            if ( (1ll<<option) & HWW_MT2DMLL )  samples[s]->get2DResults(binMinIdx, binMaxIdx, 1, y_em, e_em);
            if ( (1ll<<option) & HWW_MT2DMLL )  samples[s]->get2DResults(binMinIdx, binMaxIdx, 2, y_me, e_me);
            yield[s] = y_em + y_me;
            err[s] = sqrt(e_em*e_em + e_me*e_me);
        } else if (fcode == ((1<<0)|(1<<1)|(1<<2)|(1<<3))) {
            double y_ee, y_mm;
            double y_em, y_me;
            double e_em, e_me;
            double e_ee, e_mm;
            samples[s]->getResults(binMinIdx, binMaxIdx, 0, y_mm, e_mm);
            samples[s]->getResults(binMinIdx, binMaxIdx, 3, y_ee, e_ee);
            samples[s]->getResults(binMinIdx, binMaxIdx, 1, y_em, e_em);
            samples[s]->getResults(binMinIdx, binMaxIdx, 2, y_me, e_me);
            if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D)     samples[s]->get2DResults(binMinIdx, binMaxIdx, 0, y_mm, e_mm);
            if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D)     samples[s]->get2DResults(binMinIdx, binMaxIdx, 3, y_ee, e_ee);
            if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D)     samples[s]->get2DResults(binMinIdx, binMaxIdx, 1, y_em, e_em);
            if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D)     samples[s]->get2DResults(binMinIdx, binMaxIdx, 2, y_me, e_me);
            yield[s] = y_ee + y_mm + y_em + y_me;
            err[s] = sqrt(e_ee*e_ee + e_mm*e_mm + e_em*e_em + e_me*e_me); 
        } else if (fcode == (1<<0)) {
            samples[s]->getResults(binMinIdx, binMaxIdx, 0, yield[s], err[s]);
            if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D)  samples[s]->get2DResults(binMinIdx, binMaxIdx, 0, yield[s], err[s]);
        } else if (fcode == (1<<1)) {
            samples[s]->getResults(binMinIdx, binMaxIdx, 1, yield[s], err[s]);
            if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D)  samples[s]->get2DResults(binMinIdx, binMaxIdx, 1, yield[s], err[s]);
        } else if (fcode == (1<<2)) {
            samples[s]->getResults(binMinIdx, binMaxIdx, 2, yield[s], err[s]);
            if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D)  samples[s]->get2DResults(binMinIdx, binMaxIdx, 2, yield[s], err[s]);
        } else if (fcode == (1<<3)) {
            samples[s]->getResults(binMinIdx, binMaxIdx, 3, yield[s], err[s]);
            if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_SSCTL2D)     samples[s]->get2DResults(binMinIdx, binMaxIdx, 3, yield[s], err[s]);
        }

        // the yield can not be negative
        if (yield[s] <= 0)  {
            yield[s] = 0.0;
            err[s] = 0.0;
        }

        // associate this with a specific sample
        if (samples[s]->getDataType() == DATA)               i_data         = s;
        if (samples[s]->getDataType() == ZTT)                i_DYtt         = s;
        if (samples[s]->getDataType() == TOP)                i_Top          = s;
        if (samples[s]->getDataType() == WZ)                 i_WZ           = s;
        if (samples[s]->getDataType() == ZZ)                 i_ZZ           = s;
        if (samples[s]->getDataType() == VV)                 i_VV           = s;
        if (samples[s]->getDataType() == GGWW)               i_ggWW         = s;
        if (samples[s]->getDataType() == QQWW)               i_qqWW         = s;
        if (samples[s]->getDataType() == WGAMMA)             i_Wgamma       = s;
        if (samples[s]->getDataType() == WGAMMANORM)         i_WgammaNorm   = s;
        if (samples[s]->getDataType() == WG3L)               i_Wg3l         = s;
        if (samples[s]->getDataType() == GGHWW)              i_ggH          = s; 
        if (samples[s]->getDataType() == QQHWW)              i_qqH          = s; 
        if (samples[s]->getDataType() == WHWW)               i_WH           = s;
        if (samples[s]->getDataType() == ZHWW)               i_ZH           = s;
        if (samples[s]->getDataType() == WJETSELEDATA)       i_WjetsE       = s;
        if (samples[s]->getDataType() == WJETSMUDATA)        i_WjetsM       = s;
        if (samples[s]->getDataType() == ZJETS)              i_Zjets        = s; 
        if (samples[s]->getDataType() == ZLL)                i_DYll         = s;

    }

    //
    // get uncertainties
    //

    // now get the lumi scale factor uncertainty
    double lumiSF = 0.0;
    double lumiErr = 0.0;
    getLumiScaleFactor(lumiSF, lumiErr, option);

    // ueps
    Double_t ueps[3];
    for (int i = 0; i < 3; i++) {
        ueps[i] = HiggsSignalPSUESystematics(analysis, i);
    }

    // QCD Scale
    double QCDscale_ggH[3];
    double QCDscale_ggH1in[3];
    double QCDscale_ggH2in[3];
    for ( int i = 0; i < 3; i++) {
        QCDscale_ggH[i] = HiggsSignalQCDScaleKappa("QCDscale_ggH", analysis, i);
        QCDscale_ggH1in[i] = HiggsSignalQCDScaleKappa("QCDscale_ggH1in", analysis, i);
        QCDscale_ggH2in[i] = HiggsSignalQCDScaleKappa("QCDscale_ggH2in", analysis, i);
    }
    // add the higgs width uncertainty
    double theoryUncXS_HighMH = 1.0;
    //if(analysis >= 200) theoryUncXS_HighMH = 1.0+1.5*(analysis/1000.0)*(analysis/1000.0)*(analysis/1000.0);

    //
    // now print the card
    //

    // get the flavor of the analysis
    // if it is hww - this is needed for
    // the uncertainty names in the card 
    // and also the card file name

    std::string flavor;
    std::string fcardname;

    if (fcode == ((1<<0)|(1<<3))) {
        flavor = Form("sf");
    } else if (fcode == ((1<<1)|(1<<2))) {
        flavor = Form("of");
    } else if (fcode == ((1<<0)|(1<<1)|(1<<2)|(1<<3))) {
        flavor = Form("ll");
    } else if (fcode == (1<<0)) {
        flavor = Form("mm");
    } else if (fcode == (1<<1)) {
        flavor = Form("me");
    } else if (fcode == (1<<2)) {
        flavor = Form("em");
    } else if (fcode == (1<<3)) {
        flavor = Form("ee");
    } else {
        std::cout << "[printCard] Error! Invalid flavor type!" << std::endl;
    }

    std::string runera = "_8TeV"; 
    if ( runEra == RUN2012) 
        runera = "_8TeV";

    // cut based
    if (jetbin <= 2 && option == HWW_OPT_SMURFCUTSEL)
        fcardname = Form("%s/%i/hww%s_%ij_cut%s.txt", cdir.c_str(), int(analysis), flavor.c_str(), jetbin, runera.c_str());

    // 2D shape based
    if (option == HWW_OPT_MT2DMLL || option == HWW_OPT_MT2DMLL_JCP || option == HWW_OPT_SSCTL2D  )
        fcardname = Form("%s/%i/hww%s_%ij_shape%s.txt", cdir.c_str(), int(analysis), flavor.c_str(), jetbin, runera.c_str());
    if (option == XWW_OPT_MT2DMLL_JCP )
        fcardname = Form("%s/%i/xww%s_%ij_shape%s.txt", cdir.c_str(), int(analysis), flavor.c_str(), jetbin, runera.c_str());


    FILE *fcard;
    fcard = fopen (fcardname.c_str(), "w");
    printf("[SmurfTable::printCard] Writing %s\n", fcardname.c_str());

    //****************************************
    //               generic header
    //****************************************

    fprintf(fcard, "imax 1 number of channels\n");
    fprintf(fcard, "jmax * number of background\n");
    fprintf(fcard, "kmax * number of nuisance parameters\n");



    // now get the Wjets scale factor uncertainty
    double WjetsSF[3] = {0.0, 0.0, 0.0};
    double WjetsSFErr[3] = {0.0, 0.0, 0.0};
    getWjetsScaleFactor(WjetsSF, WjetsSFErr, option);

    // now get the Top scale factor uncertainty
    double topSF[3] = {0.0, 0.0, 0.0};
    double topSFErr[3] = {0.0, 0.0, 0.0};
    getTopScaleFactor(topSF, topSFErr, option, analysis);

    // now get the Z scale factor uncertainty
    double zSF[3] = {0.0, 0.0, 0.0};
    double zSFErr[3] = {0.0, 0.0, 0.0};
    if(option != HWW_OPT_SSCTL2D) getZScaleFactor(zSF, zSFErr, option, analysis, flavor.c_str());
/*    
    // calculate the relative uncertainty of DYll + VV 
    double yield_DYll_Zjets = yield[i_Zjets] + yield[i_DYll];
    double staterr_DYll_Zjets = sqrt( pow(err[i_Zjets],2) + pow(err[i_DYll],2));
    // assumes that the VV part has 10% systematic uncertainty
    double systerr_DYll_Zjets = 0.0;
    systerr_DYll_Zjets = sqrt( pow(yield[i_DYll]*zSFErr[jetbin],2) + pow(yield[i_Zjets]*0.1,2));
    if ( analysis > 0 && fcode == ((1<<0)|(1<<3)) && analysis <= 300) {
        //yield_DYll_Zjets = yield[i_Zjets] + zSF[jetbin];
        if(option != HWW_OPT_SSCTL2D) yield_DYll_Zjets = yield[i_Zjets] + zSF[jetbin];  
        staterr_DYll_Zjets = sqrt( pow(err[i_Zjets],2) );  
        systerr_DYll_Zjets = sqrt( pow(zSF[jetbin]*zSFErr[jetbin],2) + pow(yield[i_Zjets]*0.1,2)); 
    } 
*/
    // calculate the relative uncertainty of DYll
    double yield_DYll_Zjets =  yield[i_DYll];
    double staterr_DYll_Zjets = 0.;
    double systerr_DYll_Zjets = sqrt( pow(yield[i_DYll]*zSFErr[jetbin],2));
    if ( analysis > 0 && fcode == ((1<<0)|(1<<3)) && analysis <= 300) {
        if(option != HWW_OPT_SSCTL2D) yield_DYll_Zjets = zSF[jetbin];  
        systerr_DYll_Zjets = sqrt( pow(zSF[jetbin]*zSFErr[jetbin],2));    
    } 



    // for VBF analysis with mass > 300, take the value from the 300 selections
    if ( analysis > 300 && fcode == ((1<<0)|(1<<3)) && jetbin == 2 ) {
        double zSF_vbfspecial[3] = {0.0, 0.0, 0.0};
        double zSFErr_vbfspecial[3] = {0.0, 0.0, 0.0};
        getZScaleFactor(zSF_vbfspecial, zSFErr_vbfspecial, option, 300, flavor.c_str());
        yield_DYll_Zjets = yield[i_Zjets] + zSF_vbfspecial[jetbin];
        staterr_DYll_Zjets = sqrt( pow(err[i_Zjets],2) );
        systerr_DYll_Zjets = sqrt( pow(zSF_vbfspecial[jetbin]*zSFErr_vbfspecial[jetbin],2) + pow(yield[i_Zjets]*0.1,2));
    }

    // now get the WW scale factor uncertainty
    double WWSF[3] = {0.0, 0.0, 0.0};
    double WWSFErr[3] = {0.0, 0.0, 0.0}; 
    getWWScaleFactor(WWSF, WWSFErr, option, analysis);

    // get the correct yield for the gghww for JCP analysis
    double yield_gghww = yield[i_ggH]; 
    if ( option == XWW_OPT_MT2DMLL_JCP ) {
      TH1F *temp_gghww_ref = 0;
      for (unsigned int s = 0; s < samples.size(); ++s) {
        if (samples[s]->getDataType() == GGHWWREF) {
          TH2F *temp2d      = (TH2F*)samples[s]->get2DMVAShape(jetbin, fcode)->Clone("histo_ggHRef_tmp");
          temp_gghww_ref    = UnrollHisto2DTo1D(temp2d,temp2d->GetName());
        }
      }
      yield_gghww = temp_gghww_ref->Integral(0,1000);
    }
    
    fprintf(fcard, "Observation %i\n", int(yield[i_data]));

    if ( (1ll<<option) & ( (1ll<<HWW_OPT_MT2DMLL) | (1ll<<HWW_OPT_MT2DMLL_JCP) | (1ll<<HWW_OPT_SSCTL2D) ) ) {
        if  (mva_option & shape_var)
            fprintf(fcard, "shapes *   *   hww%s_%ij.input%s.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n", flavor.c_str(), jetbin, runera.c_str());
        else
        fprintf(fcard, "shapes *          * hww%s_%ij.input%s.root  histo_$PROCESS\n", flavor.c_str(), jetbin, runera.c_str());
        
        fprintf(fcard, "shapes data_obs * hww%s_%ij.input%s.root  histo_Data \n", flavor.c_str(), jetbin, runera.c_str());
    }

    if ( option == XWW_OPT_MT2DMLL_JCP ) {
        if  (mva_option & shape_var)
            fprintf(fcard, "shapes *   *   xww%s_%ij.input%s.root  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n", flavor.c_str(), jetbin, runera.c_str());
        else
        fprintf(fcard, "shapes *          * xww%s_%ij.input%s.root  histo_$PROCESS\n", flavor.c_str(), jetbin, runera.c_str());
        
        fprintf(fcard, "shapes data_obs * xww%s_%ij.input%s.root  histo_Data \n", flavor.c_str(), jetbin, runera.c_str());
    }



    if ( option == XWW_OPT_MT2DMLL_JCP  || option == HWW_OPT_MT2DMLL_JCP  ) {
        fprintf(fcard, "bin                               j%i%s%s  j%i%s%s  j%i%s%s  j%i%s%s  j%i%s%s  j%i%s%s  j%i%s%s  j%i%s%s  j%i%s%s  j%i%s%s  j%i%s%s  j%i%s%s j%i%s%s j%i%s%s\n", 
                jetbin, flavor.c_str(), runera.c_str(), jetbin, flavor.c_str(), runera.c_str(), jetbin, flavor.c_str(), runera.c_str(), 
                jetbin, flavor.c_str(), runera.c_str(), jetbin, flavor.c_str(), runera.c_str(), jetbin, flavor.c_str(), runera.c_str(),
                jetbin, flavor.c_str(), runera.c_str(), jetbin, flavor.c_str(), runera.c_str(), jetbin, flavor.c_str(), runera.c_str(),
                jetbin, flavor.c_str(), runera.c_str(), jetbin, flavor.c_str(), runera.c_str(), jetbin, flavor.c_str(), runera.c_str(), 
                jetbin, flavor.c_str(), runera.c_str(), jetbin, flavor.c_str(), runera.c_str());
    } else {
        fprintf(fcard, "bin                               j%i%s  j%i%s  j%i%s  j%i%s  j%i%s  j%i%s  j%i%s  j%i%s  j%i%s  j%i%s  j%i%s  j%i%s j%i%s j%i%s\n", 
                jetbin, flavor.c_str(), jetbin, flavor.c_str(), jetbin, flavor.c_str(), 
                jetbin, flavor.c_str(), jetbin, flavor.c_str(), jetbin, flavor.c_str(), 
                jetbin, flavor.c_str(), jetbin, flavor.c_str(), jetbin, flavor.c_str(), 
                jetbin, flavor.c_str(), jetbin, flavor.c_str(), jetbin, flavor.c_str(), 
                jetbin, flavor.c_str(), jetbin, flavor.c_str());
    }

    fprintf(fcard, "process                            ZH    WH   qqH   ggH   qqWW   ggWW  VV    Top  Zjets WjetsE Wgamma Wg3l  Ztt  WjetsM \n"); 
    fprintf(fcard, "process                            -3    -2    -1    0      1     2     3     4     5     6     7    8      9    10     \n");


    fprintf(fcard, "rate                             %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f\n",
            yield[i_ZH], yield[i_WH], yield[i_qqH], yield_gghww, yield[i_qqWW], yield[i_ggWW], yield[i_VV],
            yield[i_Top], yield_DYll_Zjets, yield[i_WjetsE], yield[i_WgammaNorm], yield[i_Wg3l], yield[i_DYtt], yield[i_WjetsM]);   

    if ( analysis <= 200 ) 
        fprintf(fcard, "lumi%s                    lnN %4.3f %4.3f %4.3f %4.3f   -     -   %4.3f   -     -     -   %4.3f %4.3f %4.3f  -\n",
                runera.c_str(), (1.0+lumiErr), (1.0+lumiErr), (1.0+lumiErr), (1.0+lumiErr), (1.0+lumiErr), (1.0+lumiErr), (1.0+lumiErr), (1.0+lumiErr));
    else 
        fprintf(fcard, "lumi%s                    lnN %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f   -     -     -   %4.3f %4.3f %4.3f  - \n",
                runera.c_str(), (1.0+lumiErr), (1.0+lumiErr), (1.0+lumiErr), (1.0+lumiErr), (1.0+lumiErr), (1.0+lumiErr), (1.0+lumiErr),(1.0+lumiErr), (1.0+lumiErr), (1.0+lumiErr));

    if  ( (mva_option & (1ll<<LEPEFFVAR)) && ( (1ll<<option) & HWW_SHAPE ) ) {
      if ( yield[i_WH] > 0 && yield[i_ZH] > 0 && yield[i_qqH] > 0 ) 
        fprintf(fcard, "CMS_hww_MVALepEffBounding  shape 1.000 1.000 1.000 1.000 1.000 1.000 1.000   -     -     -     -     -    -    - \n");
        else 
          fprintf(fcard, "CMS_hww_MVALepEffBounding  shape    -     -    -   1.000 1.000 1.000 1.000   -     -     -     -     -    -    -\n");
    } else {
        if ( analysis <= 200) {
            fprintf(fcard, "CMS_eff_m                    lnN 1.030 1.030 1.030 1.030 1.000 1.000 1.030   -     -     -   1.030 1.030 -    -\n");
            fprintf(fcard, "CMS_eff_e                    lnN 1.040 1.040 1.040 1.040 1.000 1.000 1.040   -     -     -   1.040 1.040 -    -\n");
        } else {
            fprintf(fcard, "CMS_eff_m                    lnN 1.030 1.030 1.030 1.030 1.030 1.030 1.030   -     -     -   1.030 1.030 -    -\n");
            fprintf(fcard, "CMS_eff_e                    lnN 1.040 1.040 1.040 1.040 1.040 1.040 1.040   -     -     -   1.040 1.040 -    -\n");
        }
    }

    if  ( (mva_option & (1ll<<LEPRESVAR)) && ( (1ll<<option) & HWW_SHAPE) ) {
      if ( yield[i_WH] > 0 && yield[i_ZH] > 0 && yield[i_qqH] >  0 ) 
        fprintf(fcard, "CMS_hww_MVALepResBounding  shape 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000   -     -     -     -    -    - \n");
      else 
        fprintf(fcard, "CMS_hww_MVALepResBounding  shape    -     -    -   1.000 1.000 1.000 1.000   -     -     -     -     -    -    -\n");
    } else {
        if ( analysis <= 200) {
            fprintf(fcard, "CMS_scale_m                  lnN 1.015 1.015 1.015 1.015 1.000 1.000 1.015   -     -     -   1.015 1.015 -    -\n");
            fprintf(fcard, "CMS_scale_e                  lnN 1.020 1.020 1.020 1.020 1.000 1.000 1.020   -     -     -   1.020 1.020 -    -\n");
        } else {
            fprintf(fcard, "CMS_scale_m                  lnN 1.015 1.015 1.015 1.015 1.015 1.015 1.015   -     -     -   1.015 1.015 -    -\n");
            fprintf(fcard, "CMS_scale_e                  lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020 -    -\n");
        }
    }
    
    if  ( (mva_option & (1ll<<METVAR)) && ( (1ll<<option) & HWW_SHAPE ) )  {
      if ( yield[i_WH] > 0 && yield[i_ZH] > 0 && yield[i_qqH] > 0 ) 
        fprintf(fcard, "CMS_hww_MVAMETResBounding  shape 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000   -     -     -     -    -    -\n");
      else
        fprintf(fcard, "CMS_hww_MVAMETResBounding  shape    -     -    -   1.000 1.000 1.000 1.000   -     -     -     -     -    -    -\n");
    } else {
        if ( analysis <= 200) {
            fprintf(fcard, "CMS_hww_met_resolution       lnN 1.020 1.020 1.020 1.020 1.000 1.000 1.020   -     -     -   1.020 1.020 -   -\n");
        } else {
            fprintf(fcard, "CMS_hww_met_resolution       lnN 1.020 1.020 1.020 1.020 1.000 1.000 1.020   -     -     -   1.020 1.020 -   -\n");
        }

    }

    if  ( (mva_option & (1ll<<JETRESVAR)) && ( (1ll<<option) & HWW_SHAPE ) ) {
      if ( yield[i_WH] > 0 && yield[i_ZH] > 0 && yield[i_qqH] > 0 ) 
        fprintf(fcard, "CMS_hww_MVAJESBounding     shape 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000   -     -     -     -    -    -\n");
      else 
        fprintf(fcard, "CMS_hww_MVAJESBounding     shape    -     -    -   1.000 1.000 1.000 1.000   -     -     -     -     -    -    -\n");
    } else {
        if (jetbin == 0)
            fprintf(fcard, "CMS_scale_j                  lnN 1.020 1.020 1.020 1.020 1.020 1.020 1.020   -     -     -   1.020 1.020 -    -\n");
        if (jetbin == 1)
            fprintf(fcard, "CMS_scale_j                  lnN 1.050 1.050 1.050 1.050 1.050 1.050 1.050   -     -     -   1.050 1.050 -    -\n");
        if (jetbin == 2)
            fprintf(fcard, "CMS_scale_j                  lnN 1.100 1.100 1.100 1.100 1.100 1.100 1.100   -     -     -   1.100 1.100 -    -\n");
    }

    if  ( (1ll<<option) & HWW_SHAPE ) {
        fprintf(fcard, "FakeRate_e                   lnN  -     -     -     -     -     -     -     -     -   %4.3f   -     -    -    -\n", 
                1.0+WjetsSFErr[jetbin]);
        fprintf(fcard, "FakeRate_m                   lnN   -     -     -     -     -     -     -     -     -     -     -     -    -   %4.3f\n", 
                1.0+WjetsSFErr[jetbin]);
    } else {
        fprintf(fcard, "FakeRate_cut_e               lnN   -     -     -     -     -     -     -     -     -   %4.3f   -     -    -    -\n", 
                1.0+WjetsSFErr[jetbin]);
        fprintf(fcard, "FakeRate_cut_m               lnN   -     -     -     -     -     -     -     -     -     -     -     -    -   %4.3f\n", 
                1.0+WjetsSFErr[jetbin]);
    }

    if ( (1ll<<option) & HWW_SHAPE ) {
            fprintf(fcard, "CMS_hww_MVAWEBounding        shape  -     -     -     -     -     -     -     -     -    1.0    -     -    -    -\n");
            fprintf(fcard, "CMS_hww_MVAWMBounding        shape  -     -     -     -     -     -     -     -     -     -      -     -    -   1.0\n");
    }

    fprintf(fcard, "UEPS                         lnN   -     -     -   %4.3f   -     -     -     -     -     -     -     -     -    -\n",
            ueps[jetbin]);
    //    fprintf(fcard, "theoryUncXS_HighMH           lnN %4.3f %4.3f %4.3f %4.3f   -     -     -     -     -     -     -     -    -    -\n",  
    //        theoryUncXS_HighMH, theoryUncXS_HighMH, theoryUncXS_HighMH, theoryUncXS_HighMH);
    fprintf(fcard, "interf_ggH                   lnN   -     -     -   %4.3f   -     -     -     -     -     -     -     -     -    -\n",
        theoryUncXS_HighMH);

    double gghpdfuncert = PDFgHHSystematics(int(analysis)); 

    if ( (1ll<<option) & HWW_SHAPE )  {  
        fprintf(fcard, "pdf_gg                       lnN   -     -     -   %4.3f   -     -     -     -     -     -     -     -    -    -\n", gghpdfuncert);
        fprintf(fcard, "pdf_qqbar                    lnN 1.050 1.050 1.050   -     -     -   1.040   -     -     -   1.040 1.040  -    -\n");
        fprintf(fcard, "CMS_hww_PDFggWW            shape   -     -     -     -     -   1.000   -     -     -     -     -     -    -    -\n");
        fprintf(fcard, "CMS_hww_PDFqqWW            shape   -     -     -     -   1.000   -     -     -     -     -     -     -    -    -\n");
    } else {
        fprintf(fcard, "pdf_gg                       lnN   -     -     -   %4.3f   -   1.040   -     -     -     -     -     -    -    -\n", gghpdfuncert);
        fprintf(fcard, "pdf_qqbar                    lnN 1.050 1.050 1.050   -   1.040   -   1.040   -     -     -   1.040 1.040  -    -\n");
    }

    if (QCDscale_ggH[jetbin] != 0.00)   
        fprintf(fcard, "QCDscale_ggH                 lnN   -     -     -   %4.3f   -     -     -     -     -     -     -     -    -    -\n",
                QCDscale_ggH[jetbin]);
    else                
        fprintf(fcard, "QCDscale_ggH                 lnN   -     -     -     -     -     -     -     -     -     -     -     -    -    -\n");

    if (QCDscale_ggH1in[jetbin] != 0.00) 
        fprintf(fcard, "QCDscale_ggH1in              lnN   -     -     -   %4.3f   -     -     -     -     -     -     -     -    -    -\n",
                QCDscale_ggH1in[jetbin]);
    else      
        fprintf(fcard, "QCDscale_ggH1in              lnN   -     -     -     -     -     -     -     -     -     -     -     -    -    -\n");

    if (QCDscale_ggH2in[jetbin] != 0.00)
        fprintf(fcard, "QCDscale_ggH2in              lnN   -     -     -   %4.3f   -     -     -     -     -     -     -     -     -    -\n",
                QCDscale_ggH2in[jetbin]);
    else
        fprintf(fcard, "QCDscale_ggH2in              lnN   -     -     -     -     -     -     -     -     -     -     -     -    -    -\n");

    fprintf(fcard, "QCDscale_qqH                 lnN   -     -   1.010   -     -     -     -     -     -     -     -     -    -    -\n");
    fprintf(fcard, "QCDscale_VH                  lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -    -    -\n");
    /*
    if ( analysis <= 200) {
        fprintf(fcard, "QCDscale_WW                  lnN   -     -     -     -   1.000   -     -     -     -     -     -     -    -    -\n");
        fprintf(fcard, "QCDscale_WW1in               lnN   -     -     -     -   1.000   -     -     -     -     -     -     -    -    -\n");
        fprintf(fcard, "QCDscale_WW2in               lnN   -     -     -     -   1.000   -     -     -     -     -     -     -    -    -\n");
    } else {
        if (jetbin == 0) {
            fprintf(fcard, "QCDscale_WW                  lnN   -     -     -     -   1.042   -     -     -     -     -     -     -    -    -\n");
            fprintf(fcard, "QCDscale_WW1in               lnN   -     -     -     -   0.978   -     -     -     -     -     -     -    -    -\n");
            fprintf(fcard, "QCDscale_WW2in               lnN   -     -     -     -   1.000   -     -     -     -     -     -     -    -    -\n");
        } 
        if (jetbin == 1) {
            fprintf(fcard, "QCDscale_WW                  lnN   -     -     -     -   1.000   -     -     -     -     -     -     -    -    -\n");
            fprintf(fcard, "QCDscale_WW1in               lnN   -     -     -     -   1.076   -     -     -     -     -     -     -    -    -\n");
            fprintf(fcard, "QCDscale_WW2in               lnN   -     -     -     -   0.914   -     -     -     -     -     -     -    -    -\n");
        }
        if (jetbin == 2) {
            fprintf(fcard, "QCDscale_WW                  lnN   -     -     -     -   1.000   -     -     -     -     -     -     -    -    -\n");
            fprintf(fcard, "QCDscale_WW1in               lnN   -     -     -     -   1.000   -     -     -     -     -     -     -    -    -\n");
            fprintf(fcard, "QCDscale_WW2in               lnN   -     -     -     -   1.420   -     -     -     -     -     -     -    -    -\n");
        }
    }
    */
    fprintf(fcard, "QCDscale_VV                  lnN   -     -     -     -     -     -   1.040   -     -     -     -     -    -    -\n");
    fprintf(fcard, "QCDscale_Vgamma              lnN   -     -     -     -     -     -     -     -     -     -   1.300   -    -    -\n");
    fprintf(fcard, "QCDscale_ggVV                lnN   -     -     -     -     -   1.300   -     -     -     -     -     -    -    -\n");

    if ( !((1ll<<option) & HWW_SHAPE) ) 
        fprintf(fcard, "QCDscale_WW_EXTRAP           lnN   -     -     -     -   1.060   -     -     -     -     -     -     -    -    -\n");
    fprintf(fcard, "QCDscale_ggH_ACCEPT          lnN   -     -     -   1.020   -     -     -     -     -     -     -     -    -    -\n");
    fprintf(fcard, "QCDscale_qqH_ACCEPT          lnN   -     -   1.020   -     -     -     -     -     -     -     -     -    -    -\n");
    fprintf(fcard, "QCDscale_VH_ACCEPT           lnN 1.020 1.020   -     -     -     -     -     -     -     -     -     -    -    -\n");

    // ttbar normalization error
    fprintf(fcard, "CMS_hww_%ij_ttbar%s        lnN   -     -     -     -     -     -     -   %4.3f   -     -     -     -    -    -\n",
            jetbin, runera.c_str(), (1.0+topSFErr[jetbin]) );

    // DY normalization error 
    if ( yield_DYll_Zjets > 0.)
        fprintf(fcard, "CMS_hww%s_%ij_Z%s          lnN   -     -     -     -     -     -     -     -   %4.3f   -     -     -    -    -\n",
                flavor.c_str(), jetbin, runera.c_str(), (1.0+systerr_DYll_Zjets/yield_DYll_Zjets));

    // WW normalization error   
    if ( (1ll<<option) & HWW_SHAPE ) {
        fprintf(fcard, "CMS_hww_%ij_WW%s_SHAPE     lnU   -     -     -     -   2.000   -     -     -     -     -     -     -    -    -\n",
                jetbin, runera.c_str());
    } else {    
        fprintf(fcard, "CMS_hww_%ij_WW%s           lnN   -     -     -     -   %4.3f %4.3f   -     -     -     -     -     -    -    -\n",
                jetbin, runera.c_str(), (1.0+WWSFErr[jetbin]), (1.0+WWSFErr[jetbin]) );
    }

    // Wg3l normalization error   
    fprintf(fcard, "CMS_hww_Wg3l                 lnN   -     -     -     -     -     -     -     -     -     -     -   1.400  -    -\n");
    
    // Ztt normalization error   
    fprintf(fcard, "CMS_hww_Ztt                  lnN   -     -     -     -     -     -     -     -     -     -     -    -  1.100    -\n");


    // for the QCD bounding for XWW analysis
    if ( option == XWW_OPT_MT2DMLL_JCP && ( mva_option & (1ll<<QCDSCALEVAR) ) )  {
      fprintf(fcard, "CMS_hww_MVAggHBounding      shape  -     -     -   1.000   -     -     -     -     -     -     -     -    -    -\n", 
          flavor.c_str(), jetbin);      
    }

    if ( (1ll<<option) & HWW_SHAPE )  {

        if  ( (mva_option & (1ll<<DYSHAPEVAR)) && (fcode == ((1<<0)|(1<<3)) ) ) {
            fprintf(fcard, "CMS_hww%s_%ij_MVAZBounding   shape  -     -     -     -     -     -     -     -    2.0    -     -     -    -    -\n", 
                    flavor.c_str(), jetbin);        
        }

        if  ( mva_option & (1ll<<TOPSHAPEVAR) ) {
            fprintf(fcard, "CMS_hww_MVATopBounding      shape  -     -     -     -     -     -     -    1.0    -     -     -     -    -    -\n"); 
        }
        if  ( mva_option & (1ll<<WWSHAPEVAR)) {
            fprintf(fcard, "CMS_hww_MVAWWBounding       shape   -     -     -     -    1.0    -     -     -     -    -     -     -    -    -\n");
            fprintf(fcard, "CMS_hww_MVAWWNLOBounding    shape   -     -     -     -    1.0    -     -     -     -    -     -     -    -    -\n");
        }
    }



    // Writing out the stat. uncertainty related entries
    if (yield[i_ZH] != 0.00) {
      if ( (mva_option & (1ll<<STATVAR)) && ( (1ll<<option) & HWW_SHAPE ) )
            fprintf(fcard, "CMS_hww%s_%ij_MVAZHStatBounding%s   shape 1.0 -   -     -     -     -     -     -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_ZH%s    lnN %4.3f   -     -     -     -     -     -     -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_ZH]/yield[i_ZH])));
    }

    if (yield[i_WH] != 0.00) {
                if ( mva_option & (1ll<<STATVAR) &&  ( (1ll<<option) & HWW_SHAPE ) )
            fprintf(fcard, "CMS_hww%s_%ij_MVAWHStatBounding%s   shape -  1.0  -     -     -     -     -     -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_WH%s    lnN   -   %4.3f   -     -     -     -     -     -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_WH]/yield[i_WH])));
    }

    if (yield[i_qqH] != 0.00) {
            if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) ) 
            fprintf(fcard, "CMS_hww%s_%ij_MVAqqHStatBounding%s  shape  -   -  1.0  -     -     -     -     -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_qqH%s   lnN   -     -   %4.3f   -     -     -     -     -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_qqH]/yield[i_qqH])));
    }

    if (yield[i_ggH] != 0.00) {
            if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) ) 
            fprintf(fcard, "CMS_hww%s_%ij_MVAggHStatBounding%s  shape  -   -   -  1.0    -     -     -     -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_ggH%s   lnN   -     -     -   %4.3f   -     -     -     -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_ggH]/yield[i_ggH])));
    }

    if (yield[i_qqWW] != 0.00) {
            if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) ) 
            fprintf(fcard, "CMS_hww%s_%ij_MVAqqWWStatBounding%s shape -   -   -   -   1.0     -     -     -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_WW%s    lnN   -     -     -     -   %4.3f   -     -     -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_qqWW]/yield[i_qqWW])));
    }

    if (yield[i_ggWW] != 0.00) {
            if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) ) 
            fprintf(fcard, "CMS_hww%s_%ij_MVAggWWStatBounding%s shape -   -   -   -     -    1.0    -     -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_ggWW%s  lnN   -     -     -     -     -   %4.3f   -     -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_ggWW]/yield[i_ggWW])));
    }

    if (yield[i_VV] != 0.00) {
            if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) ) 
            fprintf(fcard, "CMS_hww%s_%ij_MVAVVStatBounding%s   shape  -   -   -   -     -     -    1.0    -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_VV%s    lnN   -     -     -     -     -     -   %4.3f   -     -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_VV]/yield[i_VV])));
    }

    if (yield[i_Top] != 0.00) {
            if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) ) 
            fprintf(fcard, "CMS_hww%s_%ij_MVATopStatBounding%s  shape  -   -   -   -     -     -     -    1.0    -     -     -     -    -    -\n",
            //fprintf(fcard, "CMS_hww%s_%ij_MVATopStatBounding%s shapeStat  -   -   -   -     -     -     -    1.0    -     -     -    -    -    -\n",
                flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_ttbar%s lnN   -     -     -     -     -     -     -   %4.3f   -     -     -     -    -    -\n",
                    flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_Top]/yield[i_Top])));
    }

    if (yield_DYll_Zjets != 0.00 ) {
            //if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) &&  (fcode == ((1<<0)|(1<<3)) ))
            if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) )
            fprintf(fcard, "CMS_hww%s_%ij_MVAZjetsStatBounding%s shape -   -   -   -     -     -     -     -   1.0    -     -     -    -    -\n",
            //fprintf(fcard, "CMS_hww%s_%ij_MVAZjetsStatBounding%s shapeStat -   -   -   -     -     -     -     -   1.0    -     -     -     -    -\n",
                flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_Z%s     lnN   -     -     -     -     -     -     -     -   %4.3f   -     -     -    -    -\n",
                flavor.c_str(), jetbin, runera.c_str(), (1.0+staterr_DYll_Zjets/yield_DYll_Zjets));
    }

    if (yield[i_WjetsE] != 0.00) {
            if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) ) 
            fprintf(fcard, "CMS_hww%s_%ij_MVAWjetsEStatBounding%s  shape -   -   -   -     -     -     -     -     -   1.0    -     -    -    -\n",
            // fprintf(fcard, "CMS_hww%s_%ij_MVAWjetsStatBounding%s shapeStat -   -   -   -     -     -     -     -     -   1.0    -     -    -    -\n",
                flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_WjetsE%s lnN   -     -     -     -     -     -     -     -     -   %4.3f   -     -    -    -\n",
                flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_WjetsE]/yield[i_WjetsE])));
    }

    if (yield[i_Wgamma] != 0.00) {
            if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) ) 
            // fprintf(fcard, "CMS_hww%s_%ij_MVAWgammaStatBounding%s shapeStat -   -   -   -     -     -     -     -     -     -   1.0   -    -    -\n",
            fprintf(fcard, "CMS_hww%s_%ij_MVAWgammaStatBounding%s shape -   -   -   -     -     -     -     -     -     -   1.0   -    -    -\n",
                flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_Wgamma%s   lnN   -     -     -     -     -     -     -     -     -     -   %4.3f   -    -    -\n",
                flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_WgammaNorm]/yield[i_WgammaNorm])));
    }
    
    if (yield[i_Wg3l] != 0.00) {
            if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) ) 
            // fprintf(fcard, "CMS_hww%s_%ij_MVAWgammaStatBounding%s shapeStat -   -   -   -     -     -     -     -     -     -   1.0   -    -    -\n",
            fprintf(fcard, "CMS_hww%s_%ij_MVAWg3lStatBounding%s   shape -   -   -   -     -     -     -     -     -     -   -   1.0    -    -\n",
                flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_Wg3l%s   lnN   -     -     -     -     -     -     -     -     -     -    -   %4.3f    -    -\n",
                flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_Wg3l]/yield[i_Wg3l])));
    }

    if (yield[i_DYtt] != 0.00) {
            if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) ) 
                fprintf(fcard, "CMS_hww%s_%ij_MVAZttStatBounding%s    shape -   -   -   -     -     -     -     -     -     -     -   - 1.0    -\n",
                        flavor.c_str(), jetbin, runera.c_str());
            else
                fprintf(fcard, "CMS_hww%s_stat_%ij_Ztt%s   lnN   -     -     -     -     -     -     -     -     -     -     -    -  %4.3f    -\n",
                        flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_DYtt]/yield[i_DYtt])));
    } 
    
    if (yield[i_WjetsM] != 0.00) {
            if ( mva_option & (1ll<<STATVAR) && ( (1ll<<option) & HWW_SHAPE ) ) 
            fprintf(fcard, "CMS_hww%s_%ij_MVAWjetsMStatBounding%s  shape -   -   -   -     -     -     -     -     -   -    -     -    -    1.0\n",
            // fprintf(fcard, "CMS_hww%s_%ij_MVAWjetsMStatBounding%s shapeStat -   -   -   -     -     -     -     -     -   -    -     -    -   1.0\n",
                flavor.c_str(), jetbin, runera.c_str());
        else
            fprintf(fcard, "CMS_hww%s_stat_%ij_WjetsM%s lnN   -     -     -     -     -     -     -     -     -   -   -     -    -    %4.3f\n",
                flavor.c_str(), jetbin, runera.c_str(), (1.0+(err[i_WjetsM]/yield[i_WjetsM])));
    }

    // end of writing out the stat. uncertainty related entries



    // close file
    fclose(fcard);


}

//
// print out stacked plots 
//
TCanvas * makeHWWAnalysisStack(Option option, float analysis, std::vector<SmurfSample *> samples, DataType dyType,
        TFile *file, const unsigned int flav, const unsigned int njet, const char *dir, const char *name, const char *title, float lumi, float dyScale)
{
    // Determine to draw data or not  
    bool drawData = true;
    // Determine to draw signal or not at preselection level 
    bool drawSignal = true;

    //
    // some configuration
    //

    const char *fullname = Form("%s_%s_%s", name, jetbin_names[njet], types[flav]);

    std::string flavor = "sf";
    if (flav == fEM || flav == fME || flav == fOF) flavor = "of";

    // now get the Z scale factor uncertainty
    // needed to scale the Z yield in some cases
    double zSF[3] = {0.0, 0.0, 0.0};
    double zSFErr[3] = {0.0, 0.0, 0.0};
    getZScaleFactor(zSF, zSFErr, option, analysis, flavor);

    //
    // get histograms
    //

    // data
    TH1F *h1_data = 0;

    // signals
    TH1F *h1_gghww = 0;
    TH1F *h1_qqhww = 0;
    TH1F *h1_zhww = 0;
    TH1F *h1_whww = 0;

    // backgrounds
    TH1F *h1_qqww = 0;
    TH1F *h1_ggww = 0;
    TH1F *h1_vv = 0;
    TH1F *h1_top = 0;
    TH1F *h1_dyll = 0;
    TH1F *h1_zjets = 0;
    TH1F *h1_wjetsE = 0;
    TH1F *h1_wjetsM = 0;
    TH1F *h1_wgamma = 0;
    TH1F *h1_wgammanorm = 0;
    TH1F *h1_wg3l = 0;
    TH1F *h1_ztt = 0;

    for (unsigned int s = 0; s < samples.size(); ++s) {

        // data
        if (samples[s]->getDataType() == DATA)          h1_data         = GetHistogram(file, samples[s], fullname);

        // signals
        if (samples[s]->getDataType() == GGHWW)         h1_gghww        = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == QQHWW)         h1_qqhww        = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == ZHWW)          h1_zhww         = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == WHWW)          h1_whww         = GetHistogram(file, samples[s], fullname);

        // backgrounds
        if (samples[s]->getDataType() == QQWW)          h1_qqww         = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == GGWW)          h1_ggww         = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == VV)            h1_vv           = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == TOP)           h1_top          = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == dyType)        h1_dyll         = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == ZJETS)         h1_zjets        = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == WJETSELEDATA)  h1_wjetsE       = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == WJETSMUDATA)   h1_wjetsM       = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == WGAMMA)        h1_wgamma       = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == WGAMMANORM)    h1_wgammanorm   = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == WG3L)          h1_wg3l         = GetHistogram(file, samples[s], fullname);
        if (samples[s]->getDataType() == ZTT)           h1_ztt          = GetHistogram(file, samples[s], fullname);

    }

    //
    // group the processes (where needed)
    // and set format
    //

    // signal
    TH1F *h1_signal = 0;
    h1_signal = (TH1F*)h1_gghww->Clone("h1_signal");
    h1_signal->Add(h1_qqhww);
    h1_signal->Add(h1_zhww);
    h1_signal->Add(h1_whww);
    h1_signal->SetLineWidth(2);
    h1_signal->SetLineColor(kRed);
    h1_signal->SetFillColor(kRed);

    // wgamma
    if (h1_wgamma->Integral(0, h1_wgammanorm->GetNbinsX()+1) > 0) {
        h1_wgamma->Scale(h1_wgammanorm->Integral(0, h1_wgammanorm->GetNbinsX()+1) 
                         / h1_wgamma->Integral(0, h1_wgammanorm->GetNbinsX()+1));
    }

    // ww background
    TH1F *h1_ww = (TH1F*)h1_qqww->Clone("h1_ww");
    h1_ww->Add(h1_ggww);
    h1_ww->SetFillColor(kAzure-9);
    h1_ww->SetLineColor(kAzure-9);    
    if(analysis  == 0 && njet == 0 ) h1_ww->Scale(1.10); // this is taken from the mH115 cut based scale factor. 
    if(analysis  == 0 && njet == 1 ) h1_ww->Scale(0.93); // this is taken from the mH115 cut based scale factor. 

    // vv background
    h1_vv->SetFillColor(kAzure-2);
    h1_vv->SetLineColor(kAzure-2);

    // top background
    h1_top->SetFillColor(kYellow);
    h1_top->SetLineColor(kYellow);

    // drell-yan background
    // first scale the yield
    if (analysis > 0.0 && analysis <= 300.0) {
        double dy_yield = h1_dyll->Integral(0, h1_dyll->GetNbinsX() + 1);
        if (dy_yield > 0.0) h1_dyll->Scale(zSF[njet]/dy_yield);
    }
    h1_dyll->Scale(dyScale);

    
    // then do the addition
    TH1F *h1_dy = (TH1F*)h1_dyll->Clone("h1_dy");
    h1_dy->Add(h1_zjets);
    h1_dy->Add(h1_ztt);
    h1_dy->SetFillColor(kGreen+2);
    h1_dy->SetLineColor(kGreen+2);

    // w+jets
    TH1F *h1_WjetsE = (TH1F*)h1_wjetsE->Clone("h1_WjetsE");
    h1_WjetsE->SetFillColor(kGray);
    h1_WjetsE->SetLineColor(kGray);

    TH1F *h1_WjetsM = (TH1F*)h1_wjetsM->Clone("h1_WjetsM");
    h1_WjetsM->SetFillColor(kGray+2);
    h1_WjetsM->SetLineColor(kGray+2);

    // W + gamma
    TH1F *h1_Wgamma = (TH1F*)h1_wgamma->Clone("h1_Wgamma");
    h1_Wgamma->SetFillColor(kOrange);
    h1_Wgamma->SetLineColor(kOrange);

    TH1F *h1_Wg3l = (TH1F*)h1_wg3l->Clone("h1_Wg3l");
    h1_Wg3l->SetFillColor(kOrange+1); 
    h1_Wg3l->SetLineColor(kOrange+1);
    
    // add backgrounds to set y maximum 
    TH1F *h1_bkg = (TH1F*)h1_qqww->Clone("h1_ww");
    h1_bkg->Add(h1_vv);
    h1_bkg->Add(h1_top);
    h1_bkg->Add(h1_dy);
    h1_bkg->Add(h1_WjetsE);
    h1_bkg->Add(h1_WjetsM);
    h1_bkg->Add(h1_Wgamma);
    h1_bkg->Add(h1_Wg3l);
    
    // make the legend

    TLegend *l1 = new TLegend(0.57, 0.60, 0.85, 0.92);
    l1->SetLineColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->SetFillColor(kWhite);
    l1->AddEntry(h1_data,   "Data",         "lp");
    if (option != HWW_OPT_SMURFPRESEL && option != WW_OPT_SMURFXSECSEL && option != HWW_OPT_SSCTL && option != HWW_OPT_SSCTL2D) 
        l1->AddEntry(h1_signal, Form("mH = %u GeV", int(analysis)), "f");
    if (option == HWW_OPT_SMURFPRESEL && drawSignal) l1->AddEntry(h1_signal, Form("mH = %u GeV", 125), "f");  
    l1->AddEntry(h1_ww,     "WW",           "f");
    l1->AddEntry(h1_vv,     "WZ/ZZ",        "f");
    l1->AddEntry(h1_top,    "Top",          "f");
    l1->AddEntry(h1_WjetsM, "Wjets(#mu)",   "f");
    l1->AddEntry(h1_WjetsE, "Wjets(e)",     "f");
    l1->AddEntry(h1_Wgamma, "W#gamma",      "f");
    l1->AddEntry(h1_Wg3l,   "W#gamma*",     "f");
    l1->AddEntry(h1_dy,     "Drell-Yan",    "f");

    // make the stack

    THStack *st = new THStack(Form("st_%s", fullname), Form("Stack;%s;Events/bin", title));
    st->SetName(Form("st_%s", fullname));
    st->Add(h1_ww);
    st->Add(h1_dy);
    st->Add(h1_top);
    st->Add(h1_vv);
    st->Add(h1_Wgamma);
    st->Add(h1_Wg3l);
    st->Add(h1_WjetsM);
    st->Add(h1_WjetsE);
    if (option != HWW_OPT_SMURFPRESEL && option != WW_OPT_SMURFXSECSEL && option != HWW_OPT_SSCTL && option != HWW_OPT_SSCTL2D) st->Add(h1_signal); 
    if (option == HWW_OPT_SMURFPRESEL && drawSignal) st->Add(h1_signal);  

    // make the canvas

    TCanvas *c1 = new TCanvas(Form("c1_%s", fullname));
    TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.0, 1.0, 1.0);
    c1->cd();
    pad1->Draw();
    pad1->cd();
    st->Draw("HIST");
    st->GetXaxis()->SetTitle(title);
    //if (option != HWW_OPT_SMURFPRESEL && option != WW_OPT_SMURFXSECSEL && option != HWW_OPT_SSCTL && option != HWW_OPT_SSCTL2D) 
    //    h1_signal->Draw("SAME HIST");
    if(drawData) h1_data->Draw("SAME E X0");   
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
    tex->Draw("SAME");
    tex2->Draw("SAME");
    if(drawData) st->SetMaximum(TMath::Max(h1_data->GetMaximum() * 2.0, h1_data->GetMaximum() + 3*sqrt(h1_data->GetMaximum())));
    if(!drawData) st->SetMaximum(TMath::Max(h1_bkg->GetMaximum() * 1.7, h1_bkg->GetMaximum() +5*sqrt(h1_bkg->GetMaximum())));
    c1->RedrawAxis();

/*  
    //
    // Chi2 test 
    //  
    // chi2 = 1/Nbins * \sum_Nbins (Nbkg-Ndata)^2/(sigma_bkg^2+sigma_data^2)
    // 
    // Uncertainties include only statistical ones. 
    // It can be a better comparision if systematics are included. 
    //
    int nbins = h1_bkg->GetXaxis()->GetNbins();
    float chi2=0;
    for(unsigned int x=1; x <= nbins; ++x) { 
        float num   = (h1_bkg->GetBinContent(x) - h1_data->GetBinContent(x))*(h1_bkg->GetBinContent(x) - h1_data->GetBinContent(x));
        float deno  = (h1_bkg->GetBinError(x)*h1_bkg->GetBinError(x) + h1_data->GetBinError(x)*h1_data->GetBinError(x));
        if(h1_bkg->GetBinContent(x)*h1_data->GetBinContent(x)) chi2 = chi2 +  num/deno;
    }

    cout << "-----------------------------------------------------------------" << endl;
    cout << " Process name          : " << name                                 << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << " Chi2                  : " << chi2                                 << endl;
    cout << " Nbins                 : " << nbins                                << endl;
    cout << " Chi2/Nbins            : " << chi2/nbins                           << endl;
    cout << " From ROOT (UU)        : " << h1_bkg->Chi2Test(h1_data,"UU")       << endl;
    cout << " From ROOT (CHI2/NDF)  : " << h1_bkg->Chi2Test(h1_data,"CHI2/NDF") << endl;
    cout << "-----------------------------------------------------------------" << endl;
*/
    return c1;

}

