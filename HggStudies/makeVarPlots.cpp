// Makes partially blinded mass distribution + fit plots for Mass-fac MVA analysis

#ifndef __CINT__
#include "TFile.h"
#include "TChain.h"
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
#include "TGaxis.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooGaussian.h"
#include "TROOT.h"
#include "TStyle.h"
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

int main(int argc, char *argv[])
{ 

	for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : fileHistos, filesig, blind, doBands, baseline, useBinnedData" <<endl; 
                return 1;
        }	

	string output="plots"; 
	string fileHistos="CMS-HGG_massfacmva_legacy_80_180_haspixseed.root";
	int ncats = 16;
	string variable = "ptoM1"; 
	int nBins = 100;
	double minBin = 0;
	double maxBin = 4;

	if( argc > 1 ) output = argv[1];
	if( argc > 2 ) fileHistos = argv[2];
	if( argc > 3 )
        {
                std::stringstream ss ( argv[3] );
                ss >> ncats;
        }
	if( argc > 4 ) variable = argv[4]; 
	if( argc > 5 ) 
        {
                std::stringstream ss ( argv[5] );
                ss >> nBins;
        }	
	if( argc > 6 )
        {
                std::stringstream ss ( argv[6] );
                ss >> minBin;
        }
	if( argc > 7 )
        {
                std::stringstream ss ( argv[7] );
                ss >> maxBin;
        }

	system(Form("rm -r %s",output.c_str()));
	system(Form("mkdir -p %s",output.c_str()));

	// Globals
	gROOT->SetStyle("Plain");
	gROOT->SetBatch(1);
	gStyle->SetOptStat(0);
	TGaxis::SetMaxDigits(3);

	RooMsgService::instance().setGlobalKillBelow(RooFit::MsgLevel(RooFit::FATAL));

	bool cic = false;
	if(ncats > 15 ) cic = true;

	std::string * labels;
		
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
		"Untagged 0"
		,"Untagged 1"
		,"Untagged 2"
		,"Untagged 3"
		,"Untagged 4"
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

	//TFile *f = TFile::Open("histo.root");
	//TTree t = (TTree)f->Get("Data"); //HERRRE !!!

	TChain * chain = new TChain("Data");
	/*string histosChain = ""; 
        fstream file;
        file.open(fileHistos.c_str(),ios::in);
	for(int i = 0; i < 100; i++)
	{
		file >> histosChain;
		cout << endl << "histosChain = " <<histosChain;
		chain->Add(histosChain.c_str());
	}
	*/
	chain->Add(fileHistos.c_str());
	
	TLatex *latex = new TLatex();	
	latex->SetTextSize(0.025);
	latex->SetNDC();
	
	TLatex *cmslatex = new TLatex();
	cmslatex->SetTextSize(0.04);
	cmslatex->SetNDC();

	//TCanvas *can = new TCanvas("can","can",800,800);

	string cut = "";

	//chain->Draw(Form("%s>>histo",variable.c_str()),cut.c_str());
	//can->Clear();

	//histo->Draw();

	//can->Print("testBlah.png");

	double binWidth = (maxBin - minBin) / nBins;

	for (int cat=0;cat<ncats;cat++){
	
		TCanvas *can = new TCanvas("c","",800,800);
		TLegend *leg = new TLegend(0.6,0.6,0.89,0.89);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);

		TH1D * histo = new TH1D("histo","histo",nBins,minBin,maxBin);

		cut = Form("dipho_mva_cat == %d",cat); //FIXME
		chain->Draw(Form("%s>>histo",variable.c_str()),cut.c_str());

		can->Clear();

		histo->Draw();
		histo->SetTitle("");
		histo->GetXaxis()->SetTitle(Form("%s",variable.c_str()));
		histo->GetXaxis()->SetNdivisions(5,5,0);
		histo->GetYaxis()->SetTitle(Form("Events / %4.1f",binWidth));
		histo->GetYaxis()->SetTitleOffset(1.2);
		histo->SetMarkerColor(kRed);
		histo->SetMarkerSize(1.5);
		histo->SetLineColor(kRed);	
	
		//frame->SetMinimum(0.0001);
		cmslatex->DrawLatex(0.50,0.8,"#splitline{CMS Preliminary}{#sqrt{s} = 8TeV L = 19.7 fb^{-1}}");
		latex->DrawLatex(0.20,0.91,labels[cat].c_str());
		can->SaveAs(Form("%s/%s_cat%d.pdf",output.c_str(),variable.c_str(),cat));
		can->SaveAs(Form("%s/%s_cat%d.png",output.c_str(),variable.c_str(),cat));	
	
		can->Clear();

		leg->Delete();
		leg = 0;
		histo->Delete();
		histo = 0;
		delete can;
		can = 0;

	}


	latex->Delete();
	latex = 0;
	cmslatex->Delete();
	cmslatex = 0;
	chain->Delete();	
	chain = 0;	

	return 0;
}

