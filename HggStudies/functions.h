#ifndef FUNCTIONS
#define FUNCTIONS

// --- Root libraries --- //
#include "TF1.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TBits.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TGraph.h"
#include "TText.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TROOT.h"
#include "TRint.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TPaletteAxis.h"

// --- Roofit libraries --- //
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooLandau.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "RooMsgService.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooFFTConvPdf.h"
#include "RooLognormal.h"
#include "RooGaussian.h"
#include "RooGamma.h"
#include "RooBifurGauss.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooProdPdf.h"
#include "RooDerivative.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooKeysPdf.h"
#include "RooVoigtian.h"
#include "RooBreitWigner.h"
#include "RooBinning.h"
#include "RooStringVar.h"
#include "RooBernstein.h"
#include "RooExponentialSum.h"
#include "RooPowerLawSum.h"
#include "RooProduct.h"
#include "RooExponential.h"
#include "RooMsgService.h"
#include "RooPolynomial.h"
//#include "RooCruijff.hh"
#include "RooExtendPdf.h"

// --- c++ libraries --- //
#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

// --- others --- //
//#include "CrystalBall.C"
//#include "setTDRStyle.C"
//#include "CMSStyle.C"

using namespace RooFit;
using namespace std;

	
	double Median(const TH1D * h1);
	int rowsNumberInFile(string filename);
	bool FileExists(const char* FileName);	
	//vector< pair<string,RooAbsPdf*> > readDatafile(string datafile, RooRealVar *obs_var, int cat, bool truth);
	map< pair<string, RooAbsPdf*>, pair<string*, double*> > readDatafile(string datafile, RooRealVar *obs_var, int cat, bool truth, RooRealVar* &fraction, string eVeto);
	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getPdf(string type, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto);
	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getBernstein(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getExponential(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getPowerLaw(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getLaurent(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getChebychev(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
	Double_t chiSquare(RooPlot* plot, char* pdfname, char* histname);

	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getBernsteinPlusVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto); 
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getExponentialPlusVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto); 
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getPowerLawPlusVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto); 
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getLaurentPlusVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto); 
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getChebychevPlusVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto);	
	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getPolynomialsPlusVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto);

	
	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getBernsteinCBgaus(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getExponentialCBgaus(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getPowerLawCBgaus(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getLaurentCBgaus(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getChebychevCBgaus(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
	
	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getBernsteinCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getExponentialCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getPowerLawCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getLaurentCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getChebychevCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
	
	

	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getCBgaus(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);

	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getCBpCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);
	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getCBconvCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth);

#endif
