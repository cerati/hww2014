// 
//  ROOT script that draws a 2D mll-mT plot from 1D histogram
//  In addition to 2D plot, projections to mll and mT are drawn as well 
//
//  Main function : 
//      void Roll1DTo2D(TString infile, TString hname, bool doZoom=false) 
//
//  Arguments : 
//      TString infile : input root file that contains 1D histograms 
//      TString hname  : name of the histogram to be plotted 
//      bool doZoom    : option to zoom in the signal box for mH=125GeV
//                       signal box is defined in Zoom(TH2* h2) function
//
//  Usage example :    
//   root -q -b Roll1DTo2D.C'("ana_Moriond13_2D/112/hwwof_0j.input_8TeV.root","histo_ggH",true)' 
// 
// 
//  !! CAUTION : Rolling algorithm depends on how unrolling was done
//               This script assumes the unrolling procedure Guillelmo and Jae use 
// 
//  Author : Jae Hyeok Yoo jaehyeokyoo@gmail.com
//

// Zoom
TH2* Zoom(TH2* h2) {

	float mtbins[7]    =   {60,70,80,90,100,110,120};
	float mllbins[6]   =   {12,30,45,60,75,100};

	TH2F* temp = new TH2F(Form("%s_zoom",h2->GetName()),Form("%s_zoom",h2->GetTitle()), 6, mtbins, 5, mllbins);

	for(unsigned int ibinX=1; ibinX <= 6; ++ibinX) {
		for(unsigned int ibinY=1; ibinY <= 5; ++ibinY) {
			temp->SetBinContent(    ibinX, ibinY, h2->GetBinContent(ibinX, ibinY) );
			temp->SetBinError(  ibinX, ibinY, h2->GetBinError(ibinX, ibinY) );
		}
	}
	temp->SetStats(0);
	return temp;
}

//
// B/W color scheme for the color blind
// from Dave who got this from Dima
//
void loadColorPalette() {

	const int nColor = 6;
	Int_t palette[nColor]; 
	palette[0] = kWhite;
	for (unsigned int i=1;i<nColor;++i){
		palette[i] = 19-i;
	}
	gStyle->SetPalette(nColor,palette);

}

void Roll1DTo2D(TString infile, TString hname, bool doZoom=false) {

	// style
	gStyle->SetPaintTextFormat(".1f");
	gStyle->SetMarkerSize(2.5);

	// palette
	//loadColorPalette(); 

	TFile *File = TFile::Open(infile, "READ");
	TH1F *h1  = (TH1F*) File->Get(hname)->Clone(hname+"tmp");

	unsigned int nbins  = h1->GetXaxis()->GetNbins();

	float mtbins[15]    =   {60,70,80,90,100,110,120,140,160,180,200,220,240,260,280};
	float mllbins[10]   =   {12,30,45,60,75,100,125,150,175,200};
	TH2F* h2 = new TH2F(hname+"2d", hname+"2d", 14, mtbins, 9, mllbins);

	for(unsigned int ibin=1; ibin <= nbins; ++ibin) {

		int ibinX;
		int ibinY;
		if(ibin<=9*6) {
			ibinX = (int) (ibin-1)/9+1;
			ibinY = (int) (ibin-1)%9+1;
		} else {
			ibinX = (int) (ibin-1)/9+1;
			ibinY = (int) (ibin-1)%9+1;
		}

		h2->SetBinContent(  ibinX, ibinY, h1->GetBinContent(ibin) );
		h2->SetBinError(    ibinX, ibinY, h1->GetBinError(ibin) );
	}
	h2->SetStats(0);

	// 
	// Define histograms 
	// 
	TH2F *h2zoom = Zoom(h2);
	TH2F *h2toDraw;

	if(!doZoom) h2toDraw = (TH2F*)h2->Clone();
	else if(doZoom) h2toDraw = (TH2F*)h2zoom->Clone(); 
	
	h2toDraw->SetXTitle("M_{T} (GeV)");
	h2toDraw->SetYTitle("M_{ll} (GeV)");
	h2toDraw->SetTitleOffset(1.5,"Y");;
	h2toDraw->SetMinimum(0);

	// projection to mll and mT
	TH1F *h1mT 	= (TH1F*)h2toDraw->ProjectionX();
	h1mT->SetStats(0);
	h1mT->SetMinimum(0);
	h1mT->SetTitleOffset(1.0,"X");;
	TH1F *h1mll = (TH1F*)h2toDraw->ProjectionY();
	h1mll->SetStats(0); 
	h1mll->SetMinimum(0);
	h1mll->SetTitleOffset(1.0,"X");;

	// Draw histograms  
	TCanvas *c = new TCanvas("c", "", 1500, 500);
	c->Divide(3,1);
	c->cd(1);	
	if(!doZoom) h2toDraw->Draw("colz");
	else if(doZoom) h2toDraw->Draw("colz text");
	TGaxis *axistmp;
	if (!doZoom) axistmp    = new TGaxis(60.,200.,280.,200.,0, 2,50510,"+L");
	else if(doZoom) axistmp = new TGaxis(60.,100.,120.,100.,0, 2,50510,"+L");
    axistmp->SetNdivisions(0);
    axistmp->SetLabelFont(0);
    axistmp->Draw();
	c->cd(2);	
	h1mT->Draw("E");
	h1mT->SetXTitle("M_{T} (GeV)");
	c->cd(3);	
	h1mll->Draw("E");
	h1mll->SetXTitle("M_{ll} (GeV)"); 

	// Save histograms  
    infile.ReplaceAll(".root","");
    TObjArray* infiletoken = infile.Tokenize("/"); 
    int infiletokensize = infiletoken->GetEntries();
    outfile = hname+"_"+((TObjString*)infiletoken->At(infiletokensize-1))->GetString();
	if(!doZoom) outfile = outfile + "_unzoomed";
	else if(doZoom) outfile = outfile + "_zoomed";
    
    TString outdir = "/home/users/jaehyeok/public_html/HWW/misc/forGuillemo/xww2p/112/";
	c->SaveAs(outdir+outfile+".pdf");
	c->SaveAs(outdir+outfile+".png");
}
