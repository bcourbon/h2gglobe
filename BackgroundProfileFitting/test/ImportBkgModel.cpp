#include <cmath>
#include <ctime>
#include <set>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooMsgService.h"
#include "RooDataHist.h"
#include "TLatex.h"
#include "TMacro.h"
#include "TKey.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "RooCmdArg.h"
#include "RooBernstein.h"
#include "RooVoigtian.h"


#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"

#include "../interface/PdfModelBuilderFAN.h"











using namespace std;
using namespace RooFit;
using namespace boost;
namespace po = program_options;

bool BLIND = false;

RooAbsPdf* getPdf(PdfModelBuilderFAN &pdfsModel, string type, int order, const char* ext="", float constant=0){
  
  if (type=="Bernstein") return pdfsModel.getBernstein(Form("%s_bern%d",ext,order),order); 
  else if (type=="Chebychev") return pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order); 
  else if (type=="Exponential") return pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order); 
  else if (type=="PowerLaw") return pdfsModel.getPowerLawSingle(Form("%s_pow%d",ext,order),order); 
  else if (type=="Laurent") return pdfsModel.getLaurentSeriesFAN(Form("%s_lau%d",ext,order),order,constant); 
  else {
    cerr << "ERROR -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}



pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >   getPdfSumVoigtian(PdfModelBuilderFAN &pdfsModel, string type, int orderOff, int orderVoig, const char* ext="", float constant=0, float bernDownBound=-11, float bernUpBound=11){
 
  const char* type1 = type.c_str();  
  return pdfsModel.getOfficialSumVoigtians(type, Form("%s_%s%d",ext,type1,orderOff), ext, orderOff, orderVoig, constant, bernDownBound, bernUpBound); 

  
}


pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >   getPdfSumVoigtianFix(PdfModelBuilderFAN &pdfsModel, string type, int orderOff, RooAbsPdf* pdfVoiFix, const char* ext="", float constant=0, float bernDownBound=-11, float bernUpBound=11){
 
  const char* type1 = type.c_str();  
  return pdfsModel.getOfficialSumVoigtiansFix(type, Form("%s_%s%d",ext,type1,orderOff), ext, orderOff, pdfVoiFix, constant, bernDownBound, bernUpBound); 

  
  //cerr << "ERROR -- getPdfSumVoigtian() -- type " << type << " not recognised." << endl;
  //return NULL;
}


pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >   getPdfSumVoigtianFixNew(PdfModelBuilderFAN &pdfsModel, string type, int orderOff, float voiMean, float voiMeanErrorL, float voiMeanErrorH, float voiSigma, float voiSigmaErrorL, float voiSigmaErrorH, float voiWidth, float voiWidthErrorL, float voiWidthErrorH, float ErrorRange, const char* ext="", float constant=0, float bernDownBound=-11, float bernUpBound=11){
 
  const char* type1 = type.c_str();  
  return pdfsModel.getOfficialSumVoigtiansFixNew(type, Form("%s_%s%d",ext,type1,orderOff), ext, orderOff, voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voiWidth, voiWidthErrorL, voiWidthErrorH, ErrorRange, constant, true, bernDownBound, bernUpBound); 

  
}




pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >   getPdfSumVoigtianFixNewDouleCB(PdfModelBuilderFAN &pdfsModel, string type, int orderOff, float voiMean, float voiMeanErrorL, float voiMeanErrorH, float voiSigma, float voiSigmaErrorL, float voiSigmaErrorH, float voinCB1, float voinCB1ErrorL, float voinCB1ErrorH, float voinCB2, float voinCB2ErrorL, float voinCB2ErrorH, float voialphaCB1, float voialphaCB2, float ErrorRange, const char* ext="", float constant=0, float bernDownBound=-11, float bernUpBound=11){
 
  const char* type1 = type.c_str();  
  return pdfsModel.getOfficialSumVoigtiansFixNewDoubleCB(type, Form("%s_%s%d",ext,type1,orderOff), ext, orderOff, voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voinCB1, voinCB1ErrorL, voinCB1ErrorH, voinCB2, voinCB2ErrorL, voinCB2ErrorH, voialphaCB1, voialphaCB2, ErrorRange, constant, true, bernDownBound, bernUpBound); 

  
}
















int main(int argc, char* argv[]){
 
  string fileName;
  string fileNameZee;  
  string functionName;
  string fileNameout;
  int ncats;
  int jcats;
  int bins; 
  string outfilename;
  bool is2011=false;
  bool useDoubleCB=false;  
  bool verbose=false;
  int mhLow;
  int mhHigh;



  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName),                                              "In file name")
    ("infilenameZee,I", po::value<string>(&fileNameZee),                                              "In file name Zee")   
    ("function,f", po::value<string>(&functionName),                                              "Function to use")
    ("Outfilename,o", po::value<string>(&fileNameout),                                              "Out file name")
    ("ncats,c", po::value<int>(&ncats)->default_value(5),                                       "Number of categories")
    ("jcats,j", po::value<int>(&jcats)->default_value(0),                                       "Start number of categories")
    ("mhLow,L", po::value<int>(&mhLow)->default_value(75),                                                                                                                 "Starting point for scan") 
    ("mhHigh,H", po::value<int>(&mhHigh)->default_value(120),                                                                                                               "End point for scan") 
    ("bins,B", po::value<int>(&bins)->default_value(180),                                                                                                                 "Bins for the dataset") 
    ("is2011",                                                                                  "Run 2011 config")
    ("useDoubleCB",                                                                                  "use double crystal ball function")   
    ("verbose,v",                                                                               "Run with more output")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("is2011")) is2011=true;
  if (vm.count("useDoubleCB"))  useDoubleCB=true;   
  if (vm.count("verbose")) verbose=true;

  if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
  }


  TFile *outputfile;
  //RooWorkspace *outputws = new RooWorkspace("cms_hgg_workspace");  
  RooWorkspace *outputws;  
  outputfile = new TFile(fileNameout.c_str(),"RECREATE");

  
  
  TFile *inFile = TFile::Open(fileName.c_str());
  RooWorkspace *inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
  outputws = (RooWorkspace*)inWS->Clone("cms_hgg_workspace");   
  



  vector<string> functionClasses;
  functionClasses.push_back("Chebychev");
  functionClasses.push_back("Bernstein");
  functionClasses.push_back("Exponential");
  functionClasses.push_back("PowerLaw");
  functionClasses.push_back("Laurent");
  
  map<string,string> namingMap;
  namingMap.insert(pair<string,string>("Bernstein","pol"));  
  namingMap.insert(pair<string,string>("Exponential","exp"));
  namingMap.insert(pair<string,string>("PowerLaw","pow"));
  namingMap.insert(pair<string,string>("Laurent","lau"));
  
  vector<pair<pair<string,int> ,pair<pair<int,int>, pair<float,float> > > >  fabChoice;
  int sqrts; string ext;
  if (is2011) {
    sqrts = 7;
    ext = "7TeV";
  }
  else {
    sqrts = 8;
    ext = "8TeV";
    fabChoice.push_back(pair<pair<string,int>, pair<pair<int,int>, pair<float,float> > >(make_pair("Bernstein",-3),make_pair(make_pair(5,1), make_pair(-11.0,11.0)))); //0 
    fabChoice.push_back(pair<pair<string,int>, pair<pair<int,int>, pair<float,float> > >(make_pair("Bernstein",-3),make_pair(make_pair(5,1), make_pair(-11.0,11.0)))); //1
    fabChoice.push_back(pair<pair<string,int>, pair<pair<int,int>, pair<float,float> > >(make_pair("Chebychev",-3),make_pair(make_pair(5,1), make_pair(-11.0,11.0)))); //2
    fabChoice.push_back(pair<pair<string,int>, pair<pair<int,int>, pair<float,float> > >(make_pair("Bernstein",-3),make_pair(make_pair(5,1), make_pair(-11.0,11.0)))); //3
  }

  // store results here

  PdfModelBuilderFAN pdfsModel;
  RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hgg_mass");
  mass->setRange(mhLow,mhHigh); 
  pdfsModel.setObsVar(mass);
  mass->setBins(bins); 
  ofstream outfile("Zee_Yield.log");  
 

  for (int cat=jcats; cat<ncats; cat++){ 
   
     
      
    RooDataSet *dataFull = (RooDataSet*)inWS->data(Form("data_mass_cat%d",cat));
    RooDataHist thisdataBinned(Form("roohist_data_mass_cat%d",cat),"data",*mass,*dataFull);
    RooDataSet *data = (RooDataSet*)&thisdataBinned; 
         

    string funcType = fabChoice[cat].first.first;
    float LaurentConstant = fabChoice[cat].first.second; 
    int orderOff = fabChoice[cat].second.first.first; 
    int orderBre = fabChoice[cat].second.first.second;
    float bernDownBound = fabChoice[cat].second.second.first;
    float bernUpBound = fabChoice[cat].second.second.second; 

    RooAbsPdf *pdfVoiFix; float voiMean=0.; float voiMeanErrorL=0.; float voiMeanErrorH=0.; float voiSigma=0.; float voiSigmaErrorL=0.; float voiSigmaErrorH=0.; float voiWidth=0;  float voiWidthErrorL=0.; float voiWidthErrorH=0.; float voinCB1=0.; float voinCB1ErrorL=0.; float voinCB1ErrorH=0.; float voinCB2=0.; float voinCB2ErrorL=0.; float voinCB2ErrorH=0.; float voialphaCB1=0.; float voialphaCB2=0.; float ErrorRange=1.0;
    if(orderBre != 0){
         TFile *inFileZee = TFile::Open(fileNameZee.c_str()); 
         RooWorkspace *inWS_Zee = (RooWorkspace*)inFileZee->Get("fTestVoi_Zee");
         if(!useDoubleCB)  pdfVoiFix = inWS_Zee->pdf(Form("ftest_Zee_Voi_%s_cat%d",ext.c_str(),cat));
         else pdfVoiFix = inWS_Zee->pdf(Form("ftest_Zee_DCB_%s_cat%d",ext.c_str(),cat));


         if(pdfVoiFix!=NULL){
              if(!useDoubleCB){
                   ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_Voi_%s_cat%d_mean_p0",ext.c_str(),cat)))->setConstant(true);
                   voiMean = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_Voi_%s_cat%d_mean_p0",ext.c_str(),cat)))->getValV();
                   voiMeanErrorL = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_Voi_%s_cat%d_mean_p0",ext.c_str(),cat)))->getErrorLo();
                   voiMeanErrorH = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_Voi_%s_cat%d_mean_p0",ext.c_str(),cat)))->getErrorHi();

                   ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_Voi_%s_cat%d_sigma_p0",ext.c_str(),cat)))->setConstant(true);
                   voiSigma = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_Voi_%s_cat%d_sigma_p0",ext.c_str(),cat)))->getValV();    
                   voiSigmaErrorL = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_Voi_%s_cat%d_sigma_p0",ext.c_str(),cat)))->getErrorLo();    
                   voiSigmaErrorH = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_Voi_%s_cat%d_sigma_p0",ext.c_str(),cat)))->getErrorHi();    

                   ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_Voi_%s_cat%d_width_p0",ext.c_str(),cat)))->setConstant(true);
                   voiWidth = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_Voi_%s_cat%d_width_p0",ext.c_str(),cat)))->getValV();
                   voiWidthErrorL = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_Voi_%s_cat%d_width_p0",ext.c_str(),cat)))->getErrorLo();
                   voiWidthErrorH = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_Voi_%s_cat%d_width_p0",ext.c_str(),cat)))->getErrorHi();
              }else{
                  ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_mean_p0",ext.c_str(),cat)))->setConstant(true);
                  voiMean = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_mean_p0",ext.c_str(),cat)))->getValV();
                  voiMeanErrorL = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_mean_p0",ext.c_str(),cat)))->getErrorLo();
                  voiMeanErrorH = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_mean_p0",ext.c_str(),cat)))->getErrorHi();

                  ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_sigma_p0",ext.c_str(),cat)))->setConstant(true);
                  voiSigma = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_sigma_p0",ext.c_str(),cat)))->getValV();
                  voiSigmaErrorL = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_sigma_p0",ext.c_str(),cat)))->getErrorLo();
                  voiSigmaErrorH = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_sigma_p0",ext.c_str(),cat)))->getErrorHi();

                  ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_nCB1_p0",ext.c_str(),cat)))->setConstant(true);
                  voinCB1 = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_nCB1_p0",ext.c_str(),cat)))->getValV();
                  voinCB1ErrorL = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_nCB1_p0",ext.c_str(),cat)))->getErrorLo();
                  voinCB1ErrorH = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_nCB1_p0",ext.c_str(),cat)))->getErrorHi();

                  ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_nCB2_p0",ext.c_str(),cat)))->setConstant(true);
                  voinCB2 = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_nCB2_p0",ext.c_str(),cat)))->getValV();
                  voinCB2ErrorL = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_nCB2_p0",ext.c_str(),cat)))->getErrorLo();
                  voinCB2ErrorH = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_nCB2_p0",ext.c_str(),cat)))->getErrorHi();

                  ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_alphaCB1_p0",ext.c_str(),cat)))->setConstant(true);
                  voialphaCB1 = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_alphaCB1_p0",ext.c_str(),cat)))->getValV();

                  ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_alphaCB2_p0",ext.c_str(),cat)))->setConstant(true);
                  voialphaCB2 = ((RooRealVar*)inWS_Zee->allVars().find(Form("ftest_Zee_DCB_%s_cat%d_alphaCB2_p0",ext.c_str(),cat)))->getValV();
              }
         }
    }
    else{
         pdfVoiFix = 0;
    }
 

    RooAbsPdf *bkgPdf;
    if(orderBre == 0){
        bkgPdf = getPdf(pdfsModel, funcType, orderOff, Form("pdf_data_pol_model_%dTeV_cat%d",sqrts,cat), LaurentConstant); 
        bkgPdf->SetName(Form("pdf_data_pol_model_%dTeV_cat%d",sqrts,cat));
    }
    else{ 
        if(functionName == "Voi"){
            if(!useDoubleCB){
                  bkgPdf = getPdfSumVoigtianFixNew(pdfsModel, funcType, orderOff, voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voiWidth, voiWidthErrorL, voiWidthErrorH, ErrorRange, Form("pdf_data_pol_model_%dTeV_cat%d",sqrts,cat), LaurentConstant, bernDownBound, bernUpBound).first;   
            }else{
                  bkgPdf = getPdfSumVoigtianFixNewDouleCB(pdfsModel, funcType, orderOff, voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voinCB1, voinCB1ErrorL, voinCB1ErrorH, voinCB2, voinCB2ErrorL, voinCB2ErrorH, voialphaCB1, voialphaCB2, ErrorRange, Form("pdf_data_pol_model_%dTeV_cat%d",sqrts,cat), LaurentConstant, bernDownBound, bernUpBound).first;     
            }
        }
        bkgPdf->SetName(Form("pdf_data_pol_model_%dTeV_cat%d",sqrts,cat));
    }

    RooArgSet *params = bkgPdf->getParameters(*data);
    params->Print("v");

    RooFitResult *fitRes = bkgPdf->fitTo(*data,Save(true),Range(mhLow,mhHigh));
   
    fitRes->floatParsInit().Print("v"); 
    fitRes->floatParsFinal().Print("v"); 
    


    if(voiMean != 0){
         if(!useDoubleCB){
              ((RooRealVar*)params->find(Form("pdf_data_pol_model_8TeV_cat%d_Fvoimean",cat)))->setConstant(false);
              ((RooRealVar*)params->find(Form("pdf_data_pol_model_8TeV_cat%d_Fvoisigma",cat)))->setConstant(false);
              ((RooRealVar*)params->find(Form("pdf_data_pol_model_8TeV_cat%d_Fvoiwidth",cat)))->setConstant(false);
         }else{
              ((RooRealVar*)params->find(Form("pdf_data_pol_model_8TeV_cat%d_Fdcbmean",cat)))->setConstant(false);
              ((RooRealVar*)params->find(Form("pdf_data_pol_model_8TeV_cat%d_Fdcbsigma",cat)))->setConstant(false);
              ((RooRealVar*)params->find(Form("pdf_data_pol_model_8TeV_cat%d_FdcbnCB1",cat)))->setConstant(false);
              ((RooRealVar*)params->find(Form("pdf_data_pol_model_8TeV_cat%d_FdcbnCB2",cat)))->setConstant(false);
         }

         params->Print("v");

         
         float BernFrac = ((RooRealVar*)fitRes->floatParsFinal().find(Form("pdf_data_pol_model_8TeV_cat%d_frac_sum1",cat)))->getValV();
         if(!useDoubleCB){
              outfile << Form("cat %d   ",cat) << data->sumEntries()*(1.0-BernFrac) << "   Mean " << voiMean << "   voiMeanErrorL  " << voiMeanErrorL << " voiMeanErrorH  "<< voiMeanErrorH << "  voiSigma  " << voiSigma << "  voiSigmaErrorL  " << voiSigmaErrorL << " voiSigmaErrorH  " << voiSigmaErrorH << "  voiWidth  " << voiWidth << "  voiWidthErrorL  " << voiWidthErrorL << " voiWidthErrorH  " << voiWidthErrorH << endl;
              outfile << endl;
         }else{
              outfile << Form("cat %d    ",cat) << data->sumEntries()*(1.0-BernFrac) << "    Mean " << voiMean << "   voiMeanErrorL  " << voiMeanErrorL << " voiMeanErrorH  "<< voiMeanErrorH << "  voiSigma  " << voiSigma << "  voiSigmaErrorL  " << voiSigmaErrorL << " voiSigmaErrorH  " << voiSigmaErrorH << "    nCB1  " << voinCB1 << "    nCB1ErrorL   " << voinCB1ErrorL << "   nCB1ErrorH  " << voinCB1ErrorH << "   nCB2  " << voinCB2 << "   nCB2ErrorL  " << voinCB2ErrorL << "    nCB2ErrorH   " << voinCB2ErrorH << endl;

         }
         

    }

    outputws->pdf(Form("pdf_data_pol_model_%dTeV_cat%d",sqrts,cat))->SetName(Form("pdf_data_pol_model_%dTeV_cat%d_OLD",sqrts,cat));
    outputws->import(*bkgPdf);
    outputws->pdf(Form("pdf_data_pol_model_%dTeV_cat%d",sqrts,cat))->Print();



    outputws->data(Form("data_mass_cat%d",cat))->Print("v");
   

    outputws->data(Form("roohist_data_mass_cat%d",cat))->Print("v");




 

    
  }




  
  outputfile->cd();
  outputws->Write();
  outputfile->Close();   


}






