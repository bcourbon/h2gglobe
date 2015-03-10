// Makes partially blinded mass distribution + fit plots for Mass-fac MVA analysis

#ifndef __CINT__
#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooGaussian.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "RooFitResult.h"
#include "RooStats/NumberCountingUtils.h"
#include "RooStats/RooStatsUtils.h"

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
using namespace RooFit;
#endif 

void doBandsFit(TGraphAsymmErrors *onesigma, TGraphAsymmErrors *twosigma, 
		RooRealVar * hmass,
		RooAbsPdf *cpdf, 
		RooCurve *nomcurve,  RooAbsData *datanorm,
		RooPlot *plot, 
		TString & catname);


//void makeBkgPlotsGeneric_2(std::string filebkg="CMS-HGG_massfacmva_legacy_80_180_haspixseed.root", std::string filesig="CMS-HGG_massfacmva_legacy_80_180_haspixseed.root", bool blind=true, bool doBands=false, bool baseline=false, bool useBinnedData=false){
int main(int argc, char *argv[])
{ 

	for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : filebkg, filesig, blind, doBands, baseline, useBinnedData" <<endl; 
                return 1;
        }	

	string output="plots"; 
	string filebkg="CMS-HGG_massfacmva_legacy_80_180_haspixseed.root";
	string filesig="CMS-HGG_massfacmva_legacy_80_180_haspixseed.root";
	bool blind=false;
	bool doBands=false;
	int ncats = 14;
	int nBins = 120;
	bool baseline=false;
	bool useBinnedData=false;
	bool drawFunction=false;

	if( argc > 1 ) output = argv[1];
	if( argc > 2 ) filebkg = argv[2];
	if( argc > 3 ) filesig = argv[3];
	if( argc > 4 ) 
        {
                std::stringstream ss ( argv[4] );
                ss >> blind;
        }
	if( argc > 5 )
        {
                std::stringstream ss ( argv[5] );
                ss >> doBands;
        }
	if( argc > 6 )
        {
                std::stringstream ss ( argv[6] );
                ss >> ncats;
        }
	if( argc > 7 )
        {
                std::stringstream ss ( argv[7] );
                ss >> nBins;
        }
	if( argc > 8 )
        {
                std::stringstream ss ( argv[8] );
                ss >> baseline;
        }
	if( argc > 9 ) 
        {
                std::stringstream ss ( argv[9] );
                ss >> useBinnedData;
        }
	if( argc > 10 ) 
        {
                std::stringstream ss ( argv[10] );
                ss >> drawFunction;
        }	

	//system(Form("rm -r %s",output.c_str()));
	system(Form("mkdir -p %s",output.c_str()));

	// Globals
	gROOT->SetStyle("Plain");
	gROOT->SetBatch(1);
	gStyle->SetOptStat(0);
	const int ncats2 = ncats;

	RooMsgService::instance().setGlobalKillBelow(RooFit::MsgLevel(RooFit::FATAL));

	std::string * labels;
	
	bool cic = false;
        if(ncats > 15 ) cic = true;
	
	std::string ciclabels[16] = { 
                "Untagged 0 (High P_{T})"
                ,"Untagged 1 (High P_{T})"
                ,"Untagged 2 (High P_{T})"
                ,"Untagged 3 (High P_{T})"
                ,"Untagged 4 (Low P_{T})"
                ,"Untagged 5 (Low P_{T})"
                ,"Untagged 6 (Low P_{T})"
                ,"Untagged 7 (Low P_{T})"
                ,"Dijet Tag 0"
                ,"Dijet Tag 1"
                ,"VH Lepton Tight"
                ,"VH Lepton Loose"
                ,"VH MET Tag"
                ,"t#bar{t}H Lepton Tag"
                ,"t#bar{t}H Multijet Tag" 
                ,"VH Dijet Tag"
        };    

	
        std::string massfactlabels[14] = { 
                "Untagged 0 btd > 0.76"
                ,"Untagged 1 btd > 0.36"
                ,"Untagged 2 btd > 0.0"
                ,"Untagged 3 btd > -0.42"
                ,"Untagged 4 bdt > -0.78"
                ,"Dijet Tag 0"
                ,"Dijet Tag 1"
                ,"Dijet Tag 2"
                ,"VH Lepton Tight"
                ,"VH Lepton Loose"
                ,"VH MET Tag"
                ,"t#bar{t}H Tag"
                ,"t#bar{t}H Multijet Tag" 
                ,"VH Dijet Tag"
        };
     
        if( cic ) { labels = ciclabels; }
        else { labels = massfactlabels; }
	
	TFile *fb = TFile::Open(filebkg.c_str());
	//fb->ls();
	TFile *fs = ( filesig.empty() ? fb : TFile::Open(filesig.c_str()) );

	cout << endl << "Coucou 1" << endl;
	
	RooWorkspace *w_bkg  = (RooWorkspace*) fb->Get("cms_hgg_workspace");
	//w_bkg->Print();
	cout << endl << "Coucou 2" << endl;

	RooRealVar *x = (RooRealVar*) w_bkg->var("CMS_hgg_mass");
	cout << endl << "Coucou 3" << endl;
	RooRealVar *intL = (RooRealVar*) w_bkg->var("IntLumi");
	cout << endl << "Coucou 4" << endl;
	double lumi = intL->getVal()/1000.;

	TLatex *latex = new TLatex();	
	latex->SetTextSize(0.025);
	latex->SetNDC();
	
	TLatex *cmslatex = new TLatex();
	cmslatex->SetTextSize(0.04);
	cmslatex->SetNDC();

	double totalGGHinDIJET = 0;
	double totalVBFinDIJET = 0;
	double totalGGHinINCL = 0;
	double totalVBFinINCL = 0;

	double totalTTHinDIJET = 0;
	double totalWZHinDIJET = 0;
	double totalWHinDIJET = 0;
	double totalZHinDIJET = 0;
	double totalTTHinINCL = 0;
	double totalWZHinINCL = 0;
	double totalWHinINCL = 0;
	double totalZHinINCL = 0;
	gROOT->Reset();
        TGaxis::SetMaxDigits(3);

	for (int cat=0;cat<ncats;cat++){
	
		cout << endl << "Coucou 5" << endl;
		TCanvas *can = new TCanvas("c","",800,800);
		//TLegend *leg = new TLegend(0.6,0.6,0.89,0.89); //FIXME
		TLegend *leg = new TLegend(.7,.75,.89,.89);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(.03);

		// Get Dataset ->
		RooAbsData *data;
		if (useBinnedData) data =  (RooDataSet*)w_bkg->data(Form("data_mass_cat%d",cat));
		else  data =  (RooDataHist*)w_bkg->data(Form("roohist_data_mass_cat%d",cat));
		data->Print();

		RooDataHist* data2 = (RooDataHist*)w_bkg->data(Form("roohist_data_mass_cat%d",cat));
	
		cout << endl << "ICI : numEntries() = " << data->numEntries() << endl;

		// Background Pdf ->
		//RooExtendPdf *bkg =  (RooExtendPdf*)w_bkg->pdf(Form("data_pol_model_8TeV_cat%d",cat));
		//RooExtendPdf *bkg =  (RooExtendPdf*)w_bkg->pdf(Form("pdf_data_pol_model_8TeV_cat%d",cat));
		RooAbsPdf *bkg =  (RooAbsPdf*)w_bkg->pdf(Form("pdf_data_pol_model_8TeV_cat%d",cat));
		//RooAbsPdf *bkg =  (RooAbsPdf*)w_bkg->pdf(Form("data_pol_model_8TeV_cat%d",cat));
		bkg->Print();
		bkg->fitTo(*data);
		RooFitResult *r = bkg->fitTo(*data,RooFit::Save(1));
    r->Print();
	
		// Get Signal pdf norms
		std::cout << "Getting Signal Components" << std::endl;
		TH1F *gghnorm = (TH1F*)fs->Get(Form("th1f_sig_ggh_mass_m125_cat%d",cat));
		TH1F *vbfnorm = (TH1F*)fs->Get(Form("th1f_sig_vbf_mass_m125_cat%d",cat));
		//TH1F *wzhnorm = (TH1F*)fs->Get(Form("th1f_sig_wzh_mass_m125_cat%d",cat)); //FIXME !!!
		TH1F *whnorm = (TH1F*)fs->Get(Form("th1f_sig_wh_mass_m125_cat%d",cat));
		TH1F *zhnorm = (TH1F*)fs->Get(Form("th1f_sig_zh_mass_m125_cat%d",cat));		
		TH1F *tthnorm = (TH1F*)fs->Get(Form("th1f_sig_tth_mass_m125_cat%d",cat));
	
		if (cat<=3){
			
			totalGGHinINCL+=gghnorm->Integral();
			totalVBFinINCL+=vbfnorm->Integral();
			totalTTHinINCL+=tthnorm->Integral();
			cout << endl << "coucou 11" << endl;
			//totalWZHinINCL+=wzhnorm->Integral();
			totalWHinINCL+=whnorm->Integral();
			totalZHinINCL+=zhnorm->Integral();
			cout << endl << "coucou 12" << endl;
		}else{
			totalGGHinDIJET+=gghnorm->Integral();
			totalVBFinDIJET+=vbfnorm->Integral();
			totalTTHinDIJET+=tthnorm->Integral();
			cout << endl << "coucou 11" << endl;
			//totalWZHinDIJET+=wzhnorm->Integral();
			totalWHinDIJET+=whnorm->Integral();
			totalZHinDIJET+=zhnorm->Integral();
		}
	
		cout << endl << "coucou 9" << endl;
	
		std::cout << "Rescaling Signal Components" << std::endl;
		gghnorm->Add(vbfnorm);
		//gghnorm->Add(wzhnorm);
		gghnorm->Add(whnorm);
		gghnorm->Add(zhnorm);
		gghnorm->Add(tthnorm);


		TH1F *allsig = (TH1F*)gghnorm->Clone();
		allsig->Rebin(2);
		allsig->SetLineColor(4);allsig->SetFillColor(38);allsig->SetFillStyle(3001) ;allsig->SetLineWidth(2);
		/// allsig->SetLineColor(1);
		/// allsig->SetFillColor(38);
		TH1F dumData("d","",80,100,180); dumData.Sumw2();dumData.SetMarkerSize(1.0);dumData.SetMarkerStyle(20);dumData.SetLineWidth(3);
		dumData.Fill(101);
		// TH1F dumSignal("s","",80,100,180); dumSignal.SetLineColor(4);dumSignal.SetFillColor(38);dumSignal.SetFillStyle(3001) ;dumSignal.SetLineWidth(2);
		TH1F dum1Sig("1s","",80,100,180); dum1Sig.SetFillColor(kYellow);dum1Sig.SetFillStyle(1001);
		TH1F dum2Sig("2s","",80,100,180); dum2Sig.SetFillColor(kGreen);dum2Sig.SetFillStyle(1001);
		TH1F dumBkg("b","",80,100,180); dumBkg.SetLineColor(kRed);dumBkg.SetLineWidth(3);
		dumBkg.Draw("P");// dumSignal.Draw("LFsame");
		dumBkg.Draw("Fsame");dum1Sig.Draw("Fsame");dum2Sig.Draw("Lsame");

		// Plot background
		RooPlot *frame = x->frame();

		std::cout << "Plotting Components" << std::endl;
		data->plotOn(frame,RooFit::Binning(nBins),RooFit::Invisible());
		/// bkg->plotOn(frame,RooFit::VisualizeError(*r,2,1),RooFit::FillColor(kGreen));
		/// bkg->plotOn(frame,RooFit::VisualizeError(*r,1,1),RooFit::FillColor(kYellow));
		//bkg->plotOn(frame,RooFit::LineColor(kRed));
		if(drawFunction) bkg->plotOn(frame,RooFit::LineColor(kRed));
		TGraphAsymmErrors *onesigma = 0, *twosigma = 0;
		if( doBands == true ) {
			onesigma = new TGraphAsymmErrors();
			twosigma = new TGraphAsymmErrors();
			TString name=Form("cat%d",cat);
      doBandsFit(onesigma, twosigma, x, bkg, dynamic_cast<RooCurve*>(frame->getObject(frame->numItems()-1)), 
				   data, frame, name );
		}
	/*	
		if( blind ) {
			//x->setRange("unblind_up",86,96);
			x->setRange("unblind_up",60,96);
			data->plotOn(frame,RooFit::Binning(nBins),RooFit::CutRange("unblind_up"));
			//x->setRange("unblind_down",150,180); 
			x->setRange("unblind_down",110,150);
			data->plotOn(frame,RooFit::Binning(nBins),RooFit::CutRange("unblind_down")); //FIXME
		} else {
			data->plotOn(frame,RooFit::Binning(nBins));
		}
	*/
		x->setRange("unblind_up",60,96);
		data->plotOn(frame,RooFit::Binning(nBins),RooFit::CutRange("unblind_up"));
		x->setRange("unblind_down",110,150);
		data->plotOn(frame,RooFit::Binning(nBins),RooFit::CutRange("unblind_down"));
		
	
		/*
		x->setRange("unblind_up",60,96);
		data->plotOn(frame,RooFit::Binning(nBins),RooFit::CutRange("unblind_up"));
		x->setRange("unblind_down",110,150);
		data->plotOn(frame,RooFit::Binning(nBins),RooFit::CutRange("unblind_down"));
		*/

		///// x->setRange("unblind_1",100,110);
		///// x->setRange("unblind_2",150,180);
		///// data->plotOn(frame,RooFit::Binning(80),RooFit::CutRange("unblind_1"));
		///// data->plotOn(frame,RooFit::Binning(80),RooFit::CutRange("unblind_2"));

		frame->SetTitle("");
		frame->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
		frame->GetXaxis()->SetNdivisions(5,5,0);
		frame->GetYaxis()->SetTitle("Events / (1 GeV)");
		frame->GetYaxis()->SetTitleOffset(1.2);
		
		leg->AddEntry(&dumData,"Data","PEL");
		//leg->AddEntry(&dumBkg,"Bkg Model","L");
		if(drawFunction == true) leg->AddEntry(&dumBkg,"Bkg Model","L");
		if(doBands == true) leg->AddEntry(&dum1Sig,"#pm 1#sigma","F");
		if(doBands == true) leg->AddEntry(&dum2Sig,"#pm 2#sigma","F");
		//leg->AddEntry(&dumSignal,"1xSM m_{H} = 125 GeV","F");
		//leg->AddEntry(allsig,"1xSM m_{H} = 125 GeV","F");
		
		frame->Draw();
		frame->SetMinimum(0.0001);
 		if( doBands == true ) {
 			twosigma->SetLineColor(kGreen);
 			twosigma->SetFillColor(kGreen);
 			twosigma->SetMarkerColor(kGreen);
 			twosigma->Draw("L3 SAME");     
 			
 			onesigma->SetLineColor(kYellow);
 			onesigma->SetFillColor(kYellow);
 			onesigma->SetMarkerColor(kYellow);
 			onesigma->Draw("L3 SAME");
 			frame->Draw("same");
 		}
		//allsig->Draw("samehistF"); //FIXME
		leg->Draw();
		//cmslatex->DrawLatex(0.15,0.8,Form("#splitline{CMS Private}{#sqrt{s} = 8TeV L = %2.1ffb^{-1}}",lumi));
		cmslatex->DrawLatex(0.36,0.80,"#scale[0.8]{#splitline{CMS Private}{#sqrt{s} = 8 TeV L = 19.7 fb^{-1}}}");
		latex->DrawLatex(0.2,0.92,labels[cat].c_str());
		can->SaveAs(Form( (baseline ? "%s/cutbasedcat%d.pdf" : "%s/massfacmvacat%d.pdf"),output.c_str(),cat));
		can->SaveAs(Form( (baseline ? "%s/cutbasedcat%d.png" : "%s/massfacmvacat%d.png"),output.c_str(),cat));
		can->SaveAs(Form( (baseline ? "%s/cutbasedcat%d.png" : "%s/massfacmvacat%d.C"),output.c_str(),cat));	
	
		if( doBands == true ) 
		{
			onesigma->Delete();
			onesigma = 0;
			twosigma->Delete();
			twosigma = 0;
		}
/*		cout<<endl<<"bkg delete"<<endl;
		bkg->Delete();
		bkg = 0;
		leg->Delete();
		leg = 0;
		delete can;
etMinimizer roofit
		cout<<endl<<"coucou end"<<endl;
	*/

		gghnorm->Delete();
		vbfnorm->Delete();
		tthnorm->Delete();
		whnorm->Delete();
		zhnorm->Delete();
		allsig->Delete();

		can->Clear();

		leg->Delete();
		leg = 0;
		delete can;
		can = 0;
	}

	// JET ID Systematics
	std::cout << "The following can be used (nick knows what they mean) as the JET Migration systematics" <<std::endl;
	std::cout << "XXX " << 1-(0.7*totalGGHinDIJET/(totalGGHinINCL))<<std::endl;
	std::cout << "YYY " << 1-(0.7*totalTTHinDIJET/(totalTTHinINCL))<<std::endl;
	std::cout << "MMM " << 1-(0.1*totalVBFinDIJET/(totalVBFinINCL))<<std::endl;
	//std::cout << "NNN " << 1-(0.1*totalWZHinDIJET/(totalWZHinINCL))<<std::endl;
	std::cout << "NNN " << 1-(0.1*totalWHinDIJET/(totalWHinINCL))<<std::endl;
	std::cout << "NNN " << 1-(0.1*totalZHinDIJET/(totalZHinINCL))<<std::endl;	

	latex->Delete();
	latex = 0;
	cmslatex->Delete();
	cmslatex = 0;


	intL->Delete();
	x->Delete();
	//w_bkg->Delete();
	fs->Delete();
	fb->Delete();


	return 0;
}


using namespace RooFit;


void doBandsFit(TGraphAsymmErrors *onesigma, TGraphAsymmErrors *twosigma, 
		RooRealVar * hmass,
		RooAbsPdf *cpdf, 
		RooCurve *nomcurve,  RooAbsData *datanorm,
		RooPlot *plot, 
		TString & catname)
{
	RooRealVar *nlim = new RooRealVar(TString::Format("nlim%s",catname.Data()),"",0.0,0.0,1e+5);


	for (int i=1; i<(plot->GetXaxis()->GetNbins()+1); ++i) {

		double lowedge = plot->GetXaxis()->GetBinLowEdge(i);
		double upedge = plot->GetXaxis()->GetBinUpEdge(i);
		double center = plot->GetXaxis()->GetBinCenter(i);
        
		double nombkg = nomcurve->interpolate(center);

		nlim->setVal(nombkg);
		hmass->setRange("errRange",lowedge,upedge);
		RooAbsPdf *epdf = 0;
		epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
        
		//RooAbsReal *nll = epdf->createNLL(*datanorm,Extended(),NumCPU(4)); //mettre 8 ou 16 CPU ou enlever commande
		RooAbsReal *nll = epdf->createNLL(*datanorm,Extended(),NumCPU(16));
		RooMinimizer minim(*nll);
		minim.setStrategy(0); // default value = 0
		minim.setPrintLevel(-1); 
		//double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
		double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);

		minim.migrad();
		//minim.minimize("Minuit");	
		minim.minos(*nlim);
		
		printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
		
		onesigma->SetPoint(i-1,center,nombkg);
		onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
        
		minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
		// eventually if cl = 0.95 this is the usual 1.92!      
		minim.migrad();
		//minim.minimize("Minuit");
		minim.minos(*nlim);
		
		twosigma->SetPoint(i-1,center,nombkg);
		twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());      
        		
		delete nll;
		delete epdf;
	}

	onesigma->Print("V");

	nlim->Delete();
	nlim = 0;
}
