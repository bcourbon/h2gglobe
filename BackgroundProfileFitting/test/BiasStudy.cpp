//Designed for individual categories

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMsgService.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"

#include "../interface/ProfileMultiplePdfs.h"
#include "../interface/PdfModelBuilder.h"
#include "../interface/PdfModelBuilderFAN.h" 




using namespace std;
using namespace RooFit;
using namespace boost;
namespace po = boost::program_options;

void readDatFile(string datFileName, int cat, vector<pair<int,pair<string,string> > > &toysMap, vector<pair<int,pair<string,string> > > &fabianMap, vector<pair<int,pair<string,string> > > &paulMap){
 
  ifstream datfile;
  datfile.open(datFileName.c_str());
  bool foundCat=false;
  if (datfile.fail()) return;
  while (datfile.good()){
    string line;
    getline(datfile,line);
    if (line=="\n" || line.substr(0,1)=="#" || line==" " || line.empty()) continue;
    if (line.substr(0,line.find("="))=="cat") {
      foundCat=true;
      if (atoi((line.substr(line.find("=")+1,string::npos)).c_str())!=cat) foundCat=false;
    }
    if (!foundCat) continue;
    string type = line.substr(0,line.find("="));
    string pdfs = line.substr(line.find("=")+1,string::npos);
    
    if (type=="truth" && starts_with(pdfs,"Hybrid")){
      vector<string> els;
      split(els,pdfs,boost::is_any_of(":"));
      string masses=els[1];
      string funcs=els[2];
      toysMap.push_back(pair<int,pair<string,string> >(-1,make_pair(masses,funcs)));
      continue;
    }
    if (type=="truth" && starts_with(pdfs,"KeysPdf")){
      vector<string> els;
      split(els,pdfs,boost::is_any_of(":"));
      string rho=els[1];
      string name=els[2];
      toysMap.push_back(pair<int,pair<string,string> > (-2,make_pair(rho,name)));
      continue;
    }
    if (type=="truth" && starts_with(pdfs,"File")){
      vector<string> els;
      split(els,pdfs,boost::is_any_of(":"));
      toysMap.push_back(pair<int,pair<string,string> > (-3,make_pair(els[2],els[0])));
      continue;
    }

    if (type=="truth" || type=="fabian" || type=="paul"){
      vector<string> els;
      split(els,pdfs,boost::is_any_of(":"));
      int order=atoi(els[1].c_str());
      string title=els[0];
      string name=els[2];
      if (type=="truth") {
        toysMap.push_back(pair<int,pair<string,string> >(order,make_pair(name,title)));
      }
      if (type=="fabian") {
        fabianMap.push_back(pair<int,pair<string,string> >(order,make_pair(name,title)));
      }
      if (type=="paul") {
        paulMap.push_back(pair<int,pair<string,string> >(order,make_pair(name,title)));
      }
    }
  }
}

void printOptionsMap(vector<pair<int,pair<string,string> > > opts){
  for (vector<pair<int,pair<string,string> > >::iterator it=opts.begin(); it!=opts.end(); it++){
    cout << "\t" << it->first << " " << it->second.first << " " << it->second.second << endl;
  }
}

int main(int argc, char* argv[]){

  string bkgFileName;
  string sigFileName;
  string zeeFileName; 
  string sigWSName;
  string bkgWSName;
  string outFileName;
  string datFileName;
  string outDir;
  int cat;
  int ntoys;
  int jobn;
  int seed;
  float mu_low;
  float mu_high;
  float mu_step;
  float expectSignal;
  int expectSignalMass;
  bool skipPlots=false;
  int verbosity;
  bool throwHybridToys=false;
  string runType; 
  int bins; 
  int nBinsForMass = 160; 
  int mhLow; 
  int mhHigh; 
  float LaurentConstant;  
  float bernDownBound;  
  float bernUpBound;  
  bool is2011=false;  
  bool useDoubleCB=false;  
  bool NoVoiTest=false;  
  bool FixOrSo=false;  
  vector<float> switchMass;
  vector<string> switchFunc;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("zeefilename,I", po::value<string>(&zeeFileName),                                          "Zee file name")   
    ("sigfilename,s", po::value<string>(&sigFileName),                                          "Signal file name")
    ("bkgfilename,b", po::value<string>(&bkgFileName),                                          "Background file name")
    ("sigwsname", po::value<string>(&sigWSName)->default_value("cms_hgg_workspace"),            "Signal workspace name")
    ("bkgwsname", po::value<string>(&bkgWSName)->default_value("cms_hgg_workspace"),            "Background workspace name")
    ("outfilename,o", po::value<string>(&outFileName)->default_value("BiasStudyOut.root"),      "Output file name")
    ("datfile,d", po::value<string>(&datFileName)->default_value("config.dat"),                 "Name of datfile containing pdf info")
    ("outDir,D", po::value<string>(&outDir)->default_value("./"),                               "Name of out directory for plots")
    ("cat,c", po::value<int>(&cat),                                                             "Category")
    ("ntoys,t", po::value<int>(&ntoys)->default_value(0),                                       "Number of toys to run")
    ("jobn,j", po::value<int>(&jobn)->default_value(0),                                         "Job number")
    ("seed,r", po::value<int>(&seed)->default_value(0),                                         "Set random seed")
    ("mulow,L", po::value<float>(&mu_low)->default_value(-3.),                                  "Value of mu to start scan")
    ("muhigh,H", po::value<float>(&mu_high)->default_value(3.),                                 "Value of mu to end scan")
    ("runType,R", po::value<string>(&runType)->default_value(""),                               "Official or low mass") 
    ("bins,B", po::value<int>(&bins)->default_value(45),                                                                                                                 "Bins for the plot") 
    ("nBinsForMass", po::value<int>(&nBinsForMass)->default_value(180),                                                                                               "nBinsForMass") 
    ("LaurentConstant,C", po::value<float>(&LaurentConstant)->default_value(-4.),                                                                                                                 "constant for Laurent") 
    ("bernDownBound,X", po::value<float>(&bernDownBound)->default_value(-11.),                                                                                                                 "lower bound for bern") 
    ("bernUpBound,Y", po::value<float>(&bernUpBound)->default_value(11.),                                                                                                                 "up bound for bern") 
    ("mhLow", po::value<int>(&mhLow)->default_value(75),                                                                                                                 "Starting point for scan") 
    ("mhHigh", po::value<int>(&mhHigh)->default_value(120),                                                                                                               "End point for scan") 
    ("is2011",                                                                                  "Run 2011 config")  
    ("useDoubleCB",                                                                                  "use double crystal ball function")   
    ("NoVoiTest",                                                                                  "do not add voi in test model, only pol part")   
    ("FixVoi",                                                                                  "fix voi or so")  
    ("skipBern",                                                                                  "Bern does not have logN")   
    ("mustep,S", po::value<float>(&mu_step)->default_value(0.01),                               "Value of mu step size")
    ("expectSignal", po::value<float>(&expectSignal)->default_value(0.),                        "Inject signal into toy")
    ("expectSignalMass", po::value<int>(&expectSignalMass)->default_value(125),                 "Inject signal at this mass")
    ("skipPlots",                                                                               "Skip full profile and toy plots")                        
    ("verbosity,v", po::value<int>(&verbosity)->default_value(0),                               "Verbosity level")
  ;    
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("skipPlots")) skipPlots=true;
  if (vm.count("is2011")) is2011=true;  
  if (vm.count("useDoubleCB"))  useDoubleCB=true;   
  if (vm.count("NoVoiTest"))  NoVoiTest=true;   
  if (vm.count("FixVoi")) FixOrSo=true;    


  //if (expectSignalMass!=110 && expectSignalMass!=115 && expectSignalMass!=120 && expectSignalMass!=125 && expectSignalMass!=130 && expectSignalMass!=135 && expectSignalMass!=140 && expectSignalMass!=145 && expectSignalMass!=150){
  if (expectSignalMass!=110 && expectSignalMass!=115 && expectSignalMass!=120 && expectSignalMass!=105 && expectSignalMass!=100 && expectSignalMass!=95 && expectSignalMass!=90 && expectSignalMass!=85 && expectSignalMass!=80){  
    cerr << "ERROR - expectSignalMass has to be integer in range (110,150,5)" << endl;
    exit(1);
  }

  vector<pair<int,pair<string,string> > > toysMap;
  vector<pair<int,pair<string,string> > > fabianMap;
  vector<pair<int,pair<string,string> > > paulMap;
  readDatFile(datFileName,cat,toysMap,fabianMap,paulMap);
  
  cout << "Toy vector.." << endl;
  printOptionsMap(toysMap);
  cout << "Fabian vector.." << endl;
  printOptionsMap(fabianMap);
  cout << "Paul vector.." << endl;
  printOptionsMap(paulMap);
  
  TStopwatch sw;
  sw.Start();
 
  if (verbosity<1) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
  }
 
  TFile *bkgFile = TFile::Open(bkgFileName.c_str());
  TFile *sigFile = TFile::Open(sigFileName.c_str());

  RooWorkspace *bkgWS = (RooWorkspace*)bkgFile->Get(bkgWSName.c_str());
  RooWorkspace *sigWS = (RooWorkspace*)sigFile->Get(sigWSName.c_str());

  if (!bkgWS || !sigWS){
    cerr << "ERROR - one of signal or background workspace is NULL" << endl;
    cerr << " (looked for ) signal = " << sigWSName.c_str() << ", background = " << bkgWSName.c_str() <<endl;
    exit(1);
  }

  RooRealVar *mass = (RooRealVar*)bkgWS->var("CMS_hgg_mass");
  mass->setRange(mhLow,mhHigh); 
  RooRealVar *mu = new RooRealVar("mu","mu",0.,mu_low,mu_high);


  
  ofstream outfileLogFab(Form("%s/sbFitParameters_fabian_%d.log",outDir.c_str(), jobn));
  ofstream outfileLogPau(Form("%s/sbFitParameters_paul_%d.log",outDir.c_str(), jobn));
  

  TFile *outFile = new TFile(outFileName.c_str(),"RECREATE");
  TTree *muTree = new TTree("muTree","muTree");
  int toyn;
  vector<string> truthModel;
  vector<double> muFab;
  vector<double> muPaul;
  vector<double> muChi2;
  vector<double> muAIC;
  vector<double> muFabErrLow;
  vector<double> muPaulErrLow;
  vector<double> muChi2ErrLow;
  vector<double> muAICErrLow;
  vector<double> muFabErrHigh;
  vector<double> muPaulErrHigh;
  vector<double> muChi2ErrHigh;
  vector<double> muAICErrHigh;
  muTree->Branch("jobn",&jobn);
  muTree->Branch("toyn",&toyn);
  muTree->Branch("truthModel",&truthModel);
  muTree->Branch("muFab",&muFab);
  muTree->Branch("muPaul",&muPaul);
  muTree->Branch("muChi2",&muChi2);
  muTree->Branch("muAIC",&muAIC);
  muTree->Branch("muFabErrLow",&muFabErrLow);
  muTree->Branch("muPaulErrLow",&muPaulErrLow);
  muTree->Branch("muChi2ErrLow",&muChi2ErrLow);
  muTree->Branch("muAICErrLow",&muAICErrLow);
  muTree->Branch("muFabErrHigh",&muFabErrHigh);
  muTree->Branch("muPaulErrHigh",&muPaulErrHigh);
  muTree->Branch("muChi2ErrHigh",&muChi2ErrHigh);
  muTree->Branch("muAICErrHigh",&muAICErrHigh);
  

  mass->setBins(nBinsForMass);
  std::string ext = is2011 ? "7TeV" : "8TeV";
  TFile *inFileZee; RooWorkspace *inWS_Zee; RooAbsPdf* pdfVoiFix;float voiMean=0; float voiMeanErrorL=0.; float voiMeanErrorH=0.; float voiSigma=0; float voiSigmaErrorL=0.; float voiSigmaErrorH=0.; float voiWidth=0; float voiWidthErrorL=0.; float voiWidthErrorH=0.; float voinCB1=0.; float voinCB1ErrorL=0.; float voinCB1ErrorH=0.; float voinCB2=0.; float voinCB2ErrorL=0.; float voinCB2ErrorH=0.; float voialphaCB1=0.; float voialphaCB2=0.; float ErrorRange=1.; 

  if(runType == "Voi"){
        inFileZee = TFile::Open(zeeFileName.c_str());
        inWS_Zee = (RooWorkspace*)inFileZee->Get("fTestVoi_Zee");
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
        else{
             pdfVoiFix = 0;
        }
  }
  else{
        inFileZee = 0;
        inWS_Zee = 0;
        pdfVoiFix = 0;
  }
  

 

  
  RooDataSet *data = (RooDataSet*)bkgWS->data(Form("data_mass_cat%d",cat));
  
  
  RooDataHist thisdataBinned(Form("roohist_data_mass_cat%d",cat),"data",*mass,*data);  
  RooDataSet *dataBinned = (RooDataSet*)&thisdataBinned;   
  RooDataSet *sigMC = (RooDataSet*)sigWS->data(Form("sig_ggh_mass_m%d_cat%d",expectSignalMass,cat));
  RooDataSet *sigMC_vbf = (RooDataSet*)sigWS->data(Form("sig_vbf_mass_m%d_cat%d",expectSignalMass,cat));
  RooDataSet *sigMC_wh = (RooDataSet*)sigWS->data(Form("sig_wh_mass_m%d_cat%d",expectSignalMass,cat));
  RooDataSet *sigMC_zh = (RooDataSet*)sigWS->data(Form("sig_zh_mass_m%d_cat%d",expectSignalMass,cat));
  RooDataSet *sigMC_tth = (RooDataSet*)sigWS->data(Form("sig_tth_mass_m%d_cat%d",expectSignalMass,cat));
  std::cout << "Signal Model Building " << std::endl; 
  sigMC->Print(); sigMC_vbf->Print(); sigMC_wh->Print();sigMC_zh->Print();sigMC_tth->Print();
  sigMC->append(*sigMC_vbf);
  sigMC->append(*sigMC_wh);
  sigMC->append(*sigMC_zh);
  sigMC->append(*sigMC_tth);
  
  if (!dataBinned || !sigMC){
    cerr << "ERROR -- one of data or signal is NULL" << endl;
    exit(1);
  }
  
  PdfModelBuilderFAN toysModel;  
  toysModel.setObsVar(mass);
  toysModel.setSignalModifier(mu);
  for (vector<pair<int,pair<string,string> > >::iterator it=toysMap.begin(); it!=toysMap.end(); it++){
    if (it->first==-1) { // this is a hyrbid toy
      throwHybridToys=true;
      vector<string> temp;
      split(temp,it->second.first,boost::is_any_of(","));
      split(switchFunc,it->second.second,boost::is_any_of(","));
      for (unsigned int i=0; i<temp.size(); i++){
        switchMass.push_back(atof(temp[i].c_str()));
      }
      continue; 
    }
    if (it->first==-2) { // this is a keys pdf toy
      double rho = lexical_cast<double>(it->second.first);
      toysModel.setKeysPdfAttributes(data,rho);
      toysModel.addBkgPdf("KeysPdf",0,Form("truth_%s_cat%d",it->second.second.c_str(),cat),false);
      continue;
    }
    if (it->first==-3) { // this is read pdf from file
      toysModel.addBkgPdf(it->second.second,it->first,it->second.first,false);
      continue;
    }
    
    if(runType == "Voi"){
        if(!useDoubleCB){
                 toysModel.addBkgPdfSumVoigtiansFixNew(it->second.second,it->first,Form("truth_%s_cat%d",it->second.first.c_str(),cat),voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voiWidth, voiWidthErrorL, voiWidthErrorH, ErrorRange, Form("truth_%s_cat%d",it->second.first.c_str(),cat),LaurentConstant,FixOrSo,bernDownBound,bernUpBound,false);  
        }else{
                 toysModel.addBkgPdfSumVoigtiansFixNewDoubleCB(it->second.second,it->first,Form("truth_%s_cat%d",it->second.first.c_str(),cat),voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voinCB1, voinCB1ErrorL, voinCB1ErrorH, voinCB2, voinCB2ErrorL, voinCB2ErrorH, voialphaCB1, voialphaCB2, ErrorRange, Form("truth_%s_cat%d",it->second.first.c_str(),cat),LaurentConstant,FixOrSo,bernDownBound,bernUpBound,false);  
        }
    }
    else{
        toysModel.addBkgPdfFAN(it->second.second,it->first,Form("truth_%s_cat%d",it->second.first.c_str(),cat),LaurentConstant,bernDownBound,bernUpBound,false);   
    }
    
  }
  toysModel.setSignalPdfFromMC(sigMC);
  toysModel.makeSBPdfs(true);
  map<string,RooAbsPdf*> toyBkgPdfs = toysModel.getBkgPdfs();
  map<string,RooAbsPdf*> toySBPdfs = toysModel.getSBPdfs();
  toysModel.setSeed(seed);

  PdfModelBuilderFAN fabianModel; 
  fabianModel.setObsVar(mass);
  fabianModel.setSignalModifier(mu);
  for (vector<pair<int,pair<string,string> > >::iterator it=fabianMap.begin(); it!=fabianMap.end(); it++){
    
    if(runType == "Voi" && !NoVoiTest){
         if(!useDoubleCB){
                fabianModel.addBkgPdfSumVoigtiansFixNew(it->second.second,it->first,Form("fabian_%s_cat%d",it->second.first.c_str(),cat),voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voiWidth, voiWidthErrorL, voiWidthErrorH, ErrorRange, Form("fabian_%s_cat%d",it->second.first.c_str(),cat),LaurentConstant,FixOrSo,bernDownBound,bernUpBound,false);   
         }else{
                fabianModel.addBkgPdfSumVoigtiansFixNewDoubleCB(it->second.second,it->first,Form("fabian_%s_cat%d",it->second.first.c_str(),cat),voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voinCB1, voinCB1ErrorL, voinCB1ErrorH, voinCB2, voinCB2ErrorL, voinCB2ErrorH, voialphaCB1, voialphaCB2, ErrorRange, Form("fabian_%s_cat%d",it->second.first.c_str(),cat),LaurentConstant,FixOrSo,bernDownBound,bernUpBound,false);   

         }

    }
    else{
         fabianModel.addBkgPdfFAN(it->second.second,it->first,Form("fabian_%s_cat%d",it->second.first.c_str(),cat),LaurentConstant,bernDownBound,bernUpBound,false);  
    }
    
  }
  fabianModel.setSignalPdfFromMC(sigMC);
  fabianModel.makeSBPdfs(false);
  map<string,RooAbsPdf*> fabianBkgPdfs = fabianModel.getBkgPdfs();
  map<string,RooAbsPdf*> fabianSBPdfs = fabianModel.getSBPdfs();

  PdfModelBuilderFAN paulModel;  
  paulModel.setObsVar(mass);
  paulModel.setSignalModifier(mu);
  for (vector<pair<int,pair<string,string> > >::iterator it=paulMap.begin(); it!=paulMap.end(); it++){
    
    if(runType == "Voi" && !NoVoiTest){
          if(!useDoubleCB){
                paulModel.addBkgPdfSumVoigtiansFixNew(it->second.second,it->first,Form("paul_%s_cat%d",it->second.first.c_str(),cat),voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voiWidth, voiWidthErrorL, voiWidthErrorH, ErrorRange, Form("paul_%s_cat%d",it->second.first.c_str(),cat),LaurentConstant,FixOrSo,bernDownBound,bernUpBound,false);   
          }else{  
                paulModel.addBkgPdfSumVoigtiansFixNewDoubleCB(it->second.second,it->first,Form("paul_%s_cat%d",it->second.first.c_str(),cat),voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voinCB1, voinCB1ErrorL, voinCB1ErrorH, voinCB2, voinCB2ErrorL, voinCB2ErrorH, voialphaCB1, voialphaCB2, ErrorRange, Form("paul_%s_cat%d",it->second.first.c_str(),cat),LaurentConstant,FixOrSo,bernDownBound,bernUpBound,false);  

          }

    }
    else{
          paulModel.addBkgPdfFAN(it->second.second,it->first,Form("paul_%s_cat%d",it->second.first.c_str(),cat),LaurentConstant,bernDownBound,bernUpBound,false);   
    }
    
  }
  paulModel.setSignalPdfFromMC(sigMC);
  paulModel.makeSBPdfs(false);
  map<string,RooAbsPdf*> paulBkgPdfs = paulModel.getBkgPdfs();
  map<string,RooAbsPdf*> paulSBPdfs = paulModel.getSBPdfs();

  ProfileMultiplePdfs fabianProfiler;
  for (map<string,RooAbsPdf*>::iterator pdf=fabianSBPdfs.begin(); pdf!=fabianSBPdfs.end(); pdf++){
    fabianProfiler.addPdf(pdf->second);
  }
  cout << "Fabian profiler pdfs:" << endl;
  fabianProfiler.printPdfs();

  ProfileMultiplePdfs paulProfiler;
  for (map<string,RooAbsPdf*>::iterator pdf=paulSBPdfs.begin(); pdf!=paulSBPdfs.end(); pdf++){
    paulProfiler.addPdf(pdf->second);
  }
  cout << "Paul profiler pdfs:" << endl;
  paulProfiler.printPdfs();

  if (!skipPlots) {
    system(Form("mkdir -p %s/plots/truthToData",outDir.c_str()));
    system(Form("mkdir -p %s/plots/envelopeNlls",outDir.c_str()));
    system(Form("mkdir -p %s/plots/toys",outDir.c_str()));
  }
  
  cout << "------ FITTING TRUTH TO DATA ------" << endl;
  toysModel.setSignalModifierConstant(false);
      toysModel.fitToData(dataBinned,false,false,true);
  

  // -----
  toysModel.setSignalModifierVal(0);  // Always throwing from the background Only fit (why I have no idea)
  toysModel.setSignalModifierConstant(true);
  toysModel.fitToData(dataBinned,false,true,true);
  if (!skipPlots) toysModel.plotPdfsToData(dataBinned,bins,Form("%s/plots/truthToData/datafit_mu%3.1f",outDir.c_str(),expectSignal),false); 
  toysModel.setSignalModifierVal(expectSignal);
  toysModel.setSignalModifierConstant(false);
  toysModel.saveWorkspace(outFile);
  
  for (int toy=0; toy<ntoys; toy++){
    cout << "---------------------------" << endl;
    cout << "--- RUNNING TOY " << toy << " / " << ntoys << " ----" << endl;
    cout << "---------------------------" << endl;

    // wipe stuff for tree
    truthModel.clear();
    muFab.clear();
    muPaul.clear();
    muChi2.clear();
    muAIC.clear();
    muFabErrLow.clear();
    muPaulErrLow.clear();
    muChi2ErrLow.clear();
    muAICErrLow.clear();
    muFabErrHigh.clear();
    muPaulErrHigh.clear();
    muChi2ErrHigh.clear();
    muAICErrHigh.clear();
    // throw toy
    map<string,RooAbsData*> toys; 
    if (throwHybridToys) {
      toysModel.throwHybridToy(Form("truth_job%d_toy%d",jobn,toy),dataBinned->sumEntries(),switchMass,switchFunc,false,true,true,true);
      toys = toysModel.getHybridToyData();
      if (!skipPlots) toysModel.plotToysWithPdfs(Form("%s/plots/toys/job%d_toy%d",outDir.c_str(),jobn,toy),bins,false); 
      if (!skipPlots) toysModel.plotHybridToy(Form("%s/plots/toys/job%d_toy%d",outDir.c_str(),jobn,toy),bins,switchMass,switchFunc,false); 
    }
    else {
      
      if(runType == "Voi")  toysModel.throwToyFAN(Form("truth_job%d_toy%d",jobn,toy),dataBinned->sumEntries(),false,true,true,true,nBinsForMass);
      else toysModel.throwToy(Form("truth_job%d_toy%d",jobn,toy),dataBinned->sumEntries(),false,true,true,true);
      
      toys = toysModel.getToyData();
      if (!skipPlots) toysModel.plotToysWithPdfs(Form("%s/plots/toys/job%d_toy%d",outDir.c_str(),jobn,toy),bins,false); 
    }
    for (map<string,RooAbsData*>::iterator it=toys.begin(); it!=toys.end(); it++){
      if(runType == "Voi"){
         for (map<string,RooAbsPdf*>::iterator itfab=fabianBkgPdfs.begin(); itfab!=fabianBkgPdfs.end(); itfab++){
            itfab->second->fitTo(*it->second,Save(true));
         }
      }
      
             if (!skipPlots) fabianProfiler.plotNominalFits(it->second,mass,bins,Form("%s/plots/toys/job%d_toy%d_fit_fab",outDir.c_str(),jobn,toy)); 

      if(runType == "Voi"){
         for (map<string,RooAbsPdf*>::iterator itpaul=paulBkgPdfs.begin(); itpaul!=paulBkgPdfs.end(); itpaul++){
            itpaul->second->fitTo(*it->second,Save(true));
         }
      }
      
            if (!skipPlots) paulProfiler.plotNominalFits(it->second,mass,bins,Form("%s/plots/toys/job%d_toy%d_fit_paul",outDir.c_str(),jobn,toy)); 

      //continue;
      // --------------------------------
      cout << "Fitting toy for truth model " << distance(toys.begin(),it) << "/" << toys.size() << " (" << it->first << ") " << endl;
      pair<double,map<string,TGraph*> > fabianMinNlls = fabianProfiler.profileLikelihood(it->second,mass,mu,mu_low,mu_high,mu_step);
      pair<double,map<string,TGraph*> > fabianEnvelope = fabianProfiler.computeEnvelope(fabianMinNlls,Form("fabEnvelope_job%d_%s_cat%d_toy%d",jobn,it->first.c_str(),cat,toy),0.);
      if (!skipPlots) fabianProfiler.plot(fabianEnvelope.second,Form("%s/plots/envelopeNlls/nlls_fab_%s_cat%d_toy%d",outDir.c_str(),it->first.c_str(),cat,toy));
     
      pair<double,map<string,TGraph*> > paulMinNlls = paulProfiler.profileLikelihood(it->second,mass,mu,mu_low,mu_high,mu_step);
      pair<double,map<string,TGraph*> > paulEnvelope = paulProfiler.computeEnvelope(paulMinNlls,Form("paulEnvelope_job%d_%s_cat%d_toy%d",jobn,it->first.c_str(),cat,toy),0.);
      if (!skipPlots) paulProfiler.plot(paulEnvelope.second,Form("%s/plots/envelopeNlls/nlls_paul_%s_cat%d_toy%d",outDir.c_str(),it->first.c_str(),cat,toy));
      pair<double,map<string,TGraph*> > chi2Envelope = paulProfiler.computeEnvelope(paulMinNlls,Form("chi2Envelope_job%d_%s_cat%d_toy%d",jobn,it->first.c_str(),cat,toy),1.);
      if (!skipPlots) paulProfiler.plot(chi2Envelope.second,Form("%s/plots/envelopeNlls/nlls_chi2_%s_cat%d_toy%d",outDir.c_str(),it->first.c_str(),cat,toy));
      pair<double,map<string,TGraph*> > aicEnvelope = paulProfiler.computeEnvelope(paulMinNlls,Form("aicEnvelope_job%d_%s_cat%d_toy%d",jobn,it->first.c_str(),cat,toy),2.);
      if (!skipPlots) paulProfiler.plot(aicEnvelope.second,Form("%s/plots/envelopeNlls/nlls_aic_%s_cat%d_toy%d",outDir.c_str(),it->first.c_str(),cat,toy));
     
      pair<double,pair<double,double> > muFabInfo = ProfileMultiplePdfs::getMinAndErrorAsymm(fabianEnvelope.second["envelope"],1.);
      pair<double,pair<double,double> > muPaulInfo = ProfileMultiplePdfs::getMinAndErrorAsymm(paulEnvelope.second["envelope"],1.);
      pair<double,pair<double,double> > muChi2Info = ProfileMultiplePdfs::getMinAndErrorAsymm(chi2Envelope.second["envelope"],1.);
      pair<double,pair<double,double> > muAICInfo = ProfileMultiplePdfs::getMinAndErrorAsymm(aicEnvelope.second["envelope"],1.);

      
      cout << "Pushing back results "  << endl;
      truthModel.push_back(it->first);
      muFab.push_back(muFabInfo.first);
      muPaul.push_back(muPaulInfo.first);
      muChi2.push_back(muChi2Info.first);
      muAIC.push_back(muAICInfo.first);
      muFabErrLow.push_back(muFabInfo.second.first);
      muPaulErrLow.push_back(muPaulInfo.second.first);
      muChi2ErrLow.push_back(muChi2Info.second.first);
      muAICErrLow.push_back(muAICInfo.second.first);
      muFabErrHigh.push_back(muFabInfo.second.second);
      muPaulErrHigh.push_back(muPaulInfo.second.second);
      muChi2ErrHigh.push_back(muChi2Info.second.second);
      muAICErrHigh.push_back(muAICInfo.second.second);

      cout << "Writing NLLS " << endl;
      outFile->cd();
      // Ignore poor fits
      fabianEnvelope.second["envelope"]->Write();
      paulEnvelope.second["envelope"]->Write();
      chi2Envelope.second["envelope"]->Write();
      aicEnvelope.second["envelope"]->Write();
    }
    
    toyn=toy;
    cout << "Fill tree " << endl;
    muTree->Fill();
  }

  outFile->cd();
  muTree->Write();
  cout << "Done." << endl;
  cout << "Whole process took..." << endl;
  cout << "\t "; sw.Print();
 
  outFile->Close();

  return 0;
}
