#include "functions.h"
#include "setTDRStyle.C"


int main(int argc, char *argv[])
{
	for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : " <<endl; 
                return 1;
        }
	

	string output = "Truthfunctions";
	string filesig = "CMS-HGG_massfit_mva_80_180.root";

	if( argc > 1 ) output = argv[1];
	if( argc > 2 ) filesig = argv[2];

	// --- style of plots --- //
        //gROOT->Reset();
        //TGaxis::SetMaxDigits(3);
        //setTDRStyle();
        //gStyle->SetOptFit(1111);  //FIXME
	//gStyle->SetOptStat(1111);

	TFile *fs = TFile::Open(filesig.c_str());

	RooWorkspace *w_bkg  = (RooWorkspace*) fs->Get("cms_hgg_workspace");

	RooRealVar *intL = (RooRealVar*) w_bkg->var("IntLumi");
	double lumi = intL->getVal()/1000.;

	cout << endl << "lumi = " << lumi;

	//intL->Delete();
	//delete w_bkg;
	//w_bkg->Delete();
	fs->Close();
	//fs->Delete("T*;*");


	return 0;

}
