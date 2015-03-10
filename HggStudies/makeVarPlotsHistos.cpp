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
	string variable = "all_mass"; 
	string variable2 = "";

	if( argc > 1 ) output = argv[1];
	if( argc > 2 ) fileHistos = argv[2];
	if( argc > 3 )
        {
                std::stringstream ss ( argv[3] );
                ss >> ncats;
        }
	if( argc > 4 ) variable = argv[4]; 
	if( argc > 5 ) variable2 = argv[5];

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
	TH1F ** histoTab = new TH1F*[ncats];
	TH1F ** histoTab2 = NULL;
	if(variable2 != "") histoTab2 = new TH1F*[ncats];
		
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

	//TFile *fh = TFile::Open("histograms_CMS-HGG_99.root");
	//TFile *fh = TFile::Open(fileHistos.c_str());

	//TFile *f = TFile::Open("histo.root");
	//TTree t = (TTree)f->Get("Data"); //HERRRE !!!
	
	TLatex *latex = new TLatex();	
	latex->SetTextSize(0.025);
	latex->SetNDC();
	
	TLatex *cmslatex = new TLatex();
	cmslatex->SetTextSize(0.04);
	cmslatex->SetNDC();

	string histosChain = "";
	fstream file;
        file.open(fileHistos.c_str(),ios::in);
        file >> histosChain;
        TFile *fh1 = TFile::Open(histosChain.c_str());
 
	if(variable2 != "")
	{
		for (int cat=0;cat<ncats;cat++)
		{
			histoTab[cat] = (TH1F*)fh1->Get(Form("%s_cat%d_Data",variable.c_str(),cat));
			histoTab2[cat] = (TH1F*)fh1->Get(Form("%s_cat%d_Data",variable2.c_str(),cat));
		}
	}
	else for (int cat=0;cat<ncats;cat++) histoTab[cat] = (TH1F*)fh1->Get(Form("%s_cat%d_Data",variable.c_str(),cat));

        //fh1->Delete();
        for(int i = 1; i < 100; i++)
        {
        	file >> histosChain;
                cout << endl << "histoChain = " << histosChain << endl;
                TFile *fh = TFile::Open(histosChain.c_str());
		for (int cat=0;cat<ncats;cat++)
		{
			if(variable2 != "")
			{
				TH1F *histo_temp = (TH1F*)fh->Get(Form("%s_cat%d_Data",variable.c_str(),cat));
				TH1F *histo_temp2 = (TH1F*)fh->Get(Form("%s_cat%d_Data",variable2.c_str(),cat));
				histoTab[cat]->Add(histo_temp);
				histoTab2[cat]->Add(histo_temp2);
				histo_temp->Delete();
				histo_temp2->Delete();
			}
			else
			{
                		TH1F *histo_temp = (TH1F*)fh->Get(Form("%s_cat%d_Data",variable.c_str(),cat));
                		histoTab[cat]->Add(histo_temp);
				histo_temp->Delete();
			}
		}
                fh->Delete();
        }

	for (int cat=0;cat<ncats;cat++){
	
		TCanvas *can = new TCanvas("c","",800,800);
		TLegend *leg = new TLegend(0.6,0.6,0.89,0.89);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);

		if(variable2 != "") histoTab[cat]->Divide(histoTab2[cat]);

		histoTab[cat]->Draw();
		histoTab[cat]->SetTitle("");
		histoTab[cat]->GetXaxis()->SetTitle(Form("%s (GeV)",variable.c_str()));
		histoTab[cat]->GetXaxis()->SetNdivisions(5,5,0);
		histoTab[cat]->GetYaxis()->SetTitle("Events / (2 GeV)");
		histoTab[cat]->GetYaxis()->SetTitleOffset(1.2);
		histoTab[cat]->SetMarkerColor(kRed);
		histoTab[cat]->SetMarkerSize(1.5);
		histoTab[cat]->SetLineColor(kRed);	
	
		//frame->SetMinimum(0.0001);
		cmslatex->DrawLatex(0.50,0.8,"#splitline{CMS Preliminary}{#sqrt{s} = 8TeV L = 19.6 fb^{-1}}");
		latex->DrawLatex(0.20,0.91,labels[cat].c_str());
		can->SaveAs(Form("%s/%s_%d.pdf",output.c_str(),variable.c_str(),cat));
		can->SaveAs(Form("%s/%s_%d.png",output.c_str(),variable.c_str(),cat));	
	
		can->Clear();

		leg->Delete();
		leg = 0;
		delete can;
		can = 0;

		histoTab[cat]->Delete();
		histoTab[cat] = 0;
	}

	delete [] histoTab;

	latex->Delete();
	latex = 0;
	cmslatex->Delete();
	cmslatex = 0;

	fh1->Delete();


	return 0;
}

