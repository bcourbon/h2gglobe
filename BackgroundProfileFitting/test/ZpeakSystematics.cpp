#include <cmath>
#include <ctime>
#include <set>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "TKey.h"
#include "TVectorD.h"

#include "HiggsAnalysis/GBRLikelihood/interface/RooDoubleCBFast.h"   

#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>

using namespace std;
using namespace RooFit;
using namespace boost;
namespace po = program_options;


int main(int argc, char* argv[]){
 
  string fileName_DF; //double-fake (DY MC)
  string fileName_SF_allMC;  //single-fake (allMC)
  string fileName_SF_DY;  //single-fake (DY MC)
  string fileName_SF_data;  //single-fake (data)
  
  string fileNameout;

  int ncats;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilenameDF,i", po::value<string>(&fileName_DF)->default_value("CMS-HGG_DY_fit.root"),  "double-fake fit") 
    ("infilenameSFallMC,j", po::value<string>(&fileName_SF_allMC)->default_value("CMS-HGG_singleFake_allMC_fit.root"), "single-fake, allMC fit" )    
    ("infilenameSFDY,k", po::value<string>(&fileName_SF_DY)->default_value("CMS-HGG_singleFake_DY_fit.root"), "single-fake, DY fit" )    
    ("infilenameSFdata,l", po::value<string>(&fileName_SF_data)->default_value("CMS-HGG_singleFake_data_fit.root"), "single-fake, data fit" )    
    ("outfilename,o", po::value<string>(&fileNameout)->default_value("Zee_Systematics.txt"), "output text file" )    
    ("ncats,c", po::value<int>(&ncats)->default_value(4),                                       "Number of categories")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }  

  ofstream logfile_stream(fileNameout.c_str());    
 
  std::string ext = "8TeV";

  for (int cat=0; cat<ncats; cat++){ 

    TFile *inFile1 = TFile::Open(fileName_DF.c_str());
    RooWorkspace *inWS1 = (RooWorkspace*)inFile1->Get("Zpeak");
    double meanDF = ((RooRealVar*)inWS1->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_mean",ext.c_str(),cat)))->getValV();
    double errmeanDF = ((RooRealVar*)inWS1->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_mean",ext.c_str(),cat)))->getError();
    double sigmaDF = ((RooRealVar*)inWS1->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_sigma",ext.c_str(),cat)))->getValV();
    double errsigmaDF = ((RooRealVar*)inWS1->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_sigma",ext.c_str(),cat)))->getError();
    double nCB1DF = ((RooRealVar*)inWS1->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_nCB1",ext.c_str(),cat)))->getValV();
    double errnCB1DF = ((RooRealVar*)inWS1->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_nCB1",ext.c_str(),cat)))->getError();
    double nCB2DF = ((RooRealVar*)inWS1->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_nCB2",ext.c_str(),cat)))->getValV();
    double errnCB2DF = ((RooRealVar*)inWS1->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_nCB2",ext.c_str(),cat)))->getError();

    TFile *inFile2 = TFile::Open(fileName_SF_allMC.c_str());
    RooWorkspace *inWS2 = (RooWorkspace*)inFile2->Get("Zpeak");
    double meanSFallMC = ((RooRealVar*)inWS2->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_mean",ext.c_str(),cat)))->getValV();
    double errmeanSFallMC = ((RooRealVar*)inWS2->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_mean",ext.c_str(),cat)))->getError();
    double sigmaSFallMC = ((RooRealVar*)inWS2->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_sigma",ext.c_str(),cat)))->getValV();
    double errsigmaSFallMC = ((RooRealVar*)inWS2->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_sigma",ext.c_str(),cat)))->getError();

    TFile *inFile3 = TFile::Open(fileName_SF_DY.c_str());
    RooWorkspace *inWS3 = (RooWorkspace*)inFile3->Get("Zpeak");
    double meanSFDY = ((RooRealVar*)inWS3->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_mean",ext.c_str(),cat)))->getValV();
    double errmeanSFDY = ((RooRealVar*)inWS3->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_mean",ext.c_str(),cat)))->getError();
    double sigmaSFDY = ((RooRealVar*)inWS3->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_sigma",ext.c_str(),cat)))->getValV();
    double errsigmaSFDY = ((RooRealVar*)inWS3->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_sigma",ext.c_str(),cat)))->getError();

    TFile *inFile4 = TFile::Open(fileName_SF_data.c_str());
    RooWorkspace *inWS4 = (RooWorkspace*)inFile4->Get("Zpeak");
    double meanSFdata = ((RooRealVar*)inWS4->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_mean",ext.c_str(),cat)))->getValV();
    double errmeanSFdata = ((RooRealVar*)inWS4->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_mean",ext.c_str(),cat)))->getError();
    double sigmaSFdata = ((RooRealVar*)inWS4->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_sigma",ext.c_str(),cat)))->getValV();
    double errsigmaSFdata = ((RooRealVar*)inWS4->allVars().find(Form("hgg_bkg_%s_cat%d_DCB_sigma",ext.c_str(),cat)))->getError();

    double dmean_dataMCall=fabs(meanSFallMC-meanSFdata);
    double smean_dataMCall=1.5*sqrt(pow(errmeanSFallMC,2)+pow(errmeanSFdata,2));
    if (dmean_dataMCall<smean_dataMCall) dmean_dataMCall=0; //the difference is not significant

    double dmean_DYMCall=fabs(meanSFallMC-meanSFDY);
    double smean_DYMCall=1.5*sqrt(pow(errmeanSFallMC,2)+pow(errmeanSFDY,2));
    if (dmean_DYMCall<smean_DYMCall) dmean_DYMCall=0; //the difference is not significant

    double systmean=sqrt(pow(errmeanDF,2)+pow(2*dmean_dataMCall,2)+pow(2*dmean_DYMCall,2)); 

    double dsigma_dataMCall=fabs(sigmaSFallMC-sigmaSFdata);
    double ssigma_dataMCall=1.5*sqrt(pow(errsigmaSFallMC,2)+pow(errsigmaSFdata,2));
    if (dsigma_dataMCall<ssigma_dataMCall) dsigma_dataMCall=0; //the difference is not significant

    double dsigma_DYMCall=fabs(sigmaSFallMC-sigmaSFDY);
    double ssigma_DYMCall=1.5*sqrt(pow(errsigmaSFallMC,2)+pow(errsigmaSFDY,2));
    if (dsigma_DYMCall<ssigma_DYMCall) dsigma_DYMCall=0; //the difference is not significant

    double systsigma=sqrt(pow(errsigmaDF,2)+pow(2*dsigma_dataMCall,2)+pow(2*dsigma_DYMCall,2)); 

    //print results in log file

	   logfile_stream<<"cat "<<cat<<endl<<" dMean (data-MCall) "<<dmean_dataMCall<<" dMean (DY-MCall) "<<dmean_DYMCall<<endl<<" dSigma (data-MCall) "<<dsigma_dataMCall<<" dSigma (DY-MCall) "<<dsigma_DYMCall<<endl<<"Mean "<<meanDF<<" errMean (stat + syst) "<<systmean<<endl<<"Sigma "<<sigmaDF<<" errSigma (stat + syst) "<<systsigma<<endl<<"nCB1 "<<nCB1DF<<" errNCB1 (stat) "<<errnCB1DF<<endl<<"nCB2 "<<nCB2DF<<" errNCB2 (stat) "<<errnCB2DF<<endl;

       logfile_stream<<endl;

    }

}





