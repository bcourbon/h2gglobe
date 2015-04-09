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
#include "RooRealVar.h"
#include "RooAbsBinning.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooMsgService.h"
#include "RooDataHist.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TMacro.h"
#include "TKey.h"
#include "TVectorD.h"
#include "RooCmdArg.h"
#include "RooBernstein.h"
#include "RooVoigtian.h"


#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooDoubleCBFast.h"   

#include "../interface/PdfModelBuilder.h"
//#include "../interface/RooDoubleCBFast.h"
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>


using namespace std;
using namespace RooFit;
using namespace boost;
namespace po = program_options;

bool BLIND = true;


void plot(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name, int bins, int mhLow, int mhHigh){
 
 
  RooPlot *plot = mass->frame();
  mass->setRange("unblindReg_1",85,95);
  mass->setRange("unblindReg_2",110,120);
  if (BLIND) {
    data->plotOn(plot,Binning(bins),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(bins),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(bins),Invisible());
  }
  else data->plotOn(plot,Binning(bins));

  TCanvas *canv = new TCanvas();

  pdf->plotOn(plot);//,RooFit::NormRange("fitdata_1,fitdata_2"));
  //pdf->paramOn(plot,RooFit::Layout(0.35,0.89,0.89),RooFit::Format("NEA",AutoPrecision(1)));  
  //plot->getAttText()->SetTextSize(0.01); 
  plot->SetMaximum(plot->GetMaximum()*2);    
  if (BLIND) plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->SetYTitle(Form("Events / %.2fGeV",float(mhHigh-mhLow)/float(bins)));
  //plot->GetYaxis()->SetLabelSize(0.05);
  plot->SetTitleSize(0.045, "Y");
  plot->SetTitleOffset(1.,"Y");
  //plot->GetXaxis()->SetLabelSize(0.15);
  //plot->GetYaxis()->SetLabelSize(0.15);
  plot->SetXTitle("m_{#gamma#gamma}(GeV)");
  plot->SetTitleSize(0.15, "X");
  plot->Draw();

  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.68,0.93,"19.7fb^{-1}(8TeV)"); 
  lat->DrawLatex(0.12,0.85,"CMS Preliminary");   
  string nameTemp = name;
  int catNum = nameTemp.find("cat");
  string catNameOLD = nameTemp.replace(0, catNum, "");
  string catName = catNameOLD.replace(catNameOLD.length()-4, 4, "");
  lat->DrawLatex(0.1,0.92, Form("%s", catName.c_str()));  

     TLegend *leg = new TLegend(0.12,0.7,0.32,0.8);
     leg->SetFillColor(0);
     leg->SetLineColor(1);
     leg->AddEntry(data,"Data","lep");
     leg->Draw("same");

  canv->SaveAs(Form("%s.png",name.c_str()));
  canv->SaveAs(Form("%s.pdf",name.c_str()));
 
  delete canv;
  delete lat;
}

void plot2(RooRealVar *mass, RooAbsPdf *pdf, RooAbsPdf *pdf1, RooAbsPdf *pdf2, float frac, RooDataSet *data, string name, int bins, int mhLow, int mhHigh){
 
 
  RooPlot *plot = mass->frame();
  mass->setRange("unblindReg_1",85,95);
  mass->setRange("unblindReg_2",110,120);
  if (BLIND) {
    data->plotOn(plot,Binning(bins),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(bins),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(bins),Invisible());
  }
  else data->plotOn(plot,Binning(bins));

  TCanvas *canv = new TCanvas();

  pdf->plotOn(plot);
  pdf->plotOn(plot,Components(*pdf1),LineColor(kRed),LineStyle(2)) ;
  pdf->plotOn(plot,Components(*pdf2),LineColor(kGreen),LineStyle(2)) ;

//,RooFit::NormRange("fitdata_1,fitdata_2"));
  //pdf->paramOn(plot,RooFit::Layout(0.35,0.89,0.89),RooFit::Format("NEA",AutoPrecision(1)));  
  //plot->getAttText()->SetTextSize(0.01); 
  //plot->SetMaximum(plot->GetMaximum()*2);    
  if (BLIND) plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->SetYTitle(Form("Events / %.2fGeV",float(mhHigh-mhLow)/float(bins)));
  //plot->GetYaxis()->SetLabelSize(0.05);
  plot->SetTitleSize(0.045, "Y");
  plot->SetTitleOffset(1.,"Y");
  //plot->GetXaxis()->SetLabelSize(0.15);
  //plot->GetYaxis()->SetLabelSize(0.15);
  plot->SetXTitle("m_{#gamma#gamma}(GeV)");
  plot->SetTitleSize(0.15, "X");
  plot->Draw();

  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.68,0.93,"19.7fb^{-1}(8TeV)"); 
  lat->DrawLatex(0.12,0.85,"CMS Preliminary");   
  lat->DrawLatex(0.7,0.7,Form("Frac = %.3f",frac));  
  string nameTemp = name;
  int catNum = nameTemp.find("cat");
  string catNameOLD = nameTemp.replace(0, catNum, "");
  string catName = catNameOLD.replace(catNameOLD.length()-4, 4, "");
  lat->DrawLatex(0.1,0.92, Form("%s", catName.c_str()));  

     TLegend *leg = new TLegend(0.12,0.7,0.32,0.8);
     leg->SetFillColor(0);
     leg->SetLineColor(1);
     leg->AddEntry(data,"Data","lep");
     leg->Draw("same");

  canv->SaveAs(Form("%s.png",name.c_str()));
  canv->SaveAs(Form("%s.pdf",name.c_str()));
 
  delete canv;
  delete lat;
}



void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries, int mhLow, int mhHigh){

	int ntries=0;
	int stat=1;
	double minnll=10e8;
	while (stat!=0){
	  if (ntries>=MaxTries) break;
	  RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1),Range(mhLow,mhHigh));
	  //RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1),Range(85,110));
	  //RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1),SumW2Error(kTRUE)
          stat = fitTest->status();
	  minnll = fitTest->minNll();
	  ntries++; 
	}
	*stat_t = stat;
	*NLL = minnll;
}


string abbr(string function){

string ab;

if(function=="Bernstein") ab="ber";
if(function=="Exponential") ab="exp";
if(function=="PowerLaw") ab="pow";
if(function=="Chebychev") ab="che";
if(function=="Laurent") ab="lau";
if(function=="Polynomial") ab="pol";

return ab;

}


int main(int argc, char* argv[]){
 
  string fileName;
  string fileNameZee;  
  string fileNameout;
  string outDir;
  int ncats;
  int jcats;
  int bins; 
  string outfilename;
  string logfile;
  bool verbose=false;
  int mhLow;
  int mhHigh;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName)->default_value("CMS-HGG.root"),                                             "In file name")
    ("infilenameZee,I", po::value<string>(&fileNameZee)->default_value("CMS-HGG_DY_fit.root"),                                  "In file name Zee")   
    ("Outfilename,o", po::value<string>(&fileNameout)->default_value("CMS-HGG_finalBkg.root"),                                              "Out file name")
    ("outDir,D", po::value<string>(&outDir)->default_value("plots/FinalBackground"),                      "Out directory for plots")
    ("logfile,l", po::value<string>(&logfile)->default_value("BackgroundFit.log"),                  "log file of fit results")
    ("ncats,c", po::value<int>(&ncats)->default_value(4),                                       "Number of categories")
    ("jcats,j", po::value<int>(&jcats)->default_value(0),                                       "Start number of categories")
    ("mhLow,L", po::value<int>(&mhLow)->default_value(75),                                                                                                                 "Starting point for scan") 
    ("mhHigh,H", po::value<int>(&mhHigh)->default_value(120),                                                                                                               "End point for scan") 
    ("bins,B", po::value<int>(&bins)->default_value(45),                                                                                                                 "Bins for the dataset") 
    ("verbose,v",                                                                               "Run with more output")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }  
  if (vm.count("verbose")) verbose=true;

  if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
  }


  TFile *outputfile;
  outputfile = new TFile(fileNameout.c_str(),"RECREATE");
  RooWorkspace *outputws;  

  ofstream logfile_stream(logfile.c_str());    

  TFile *inFile = TFile::Open(fileName.c_str());
  RooWorkspace *inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
  outputws = (RooWorkspace*)inWS->Clone("cms_hgg_workspace");   

  std::string ext = "8TeV";

  vector<pair<string,int> > funChoice;   // <functionType,order>

  funChoice.push_back(pair<string,int>("Bernstein",5)); //cat0 
  funChoice.push_back(pair<string,int>("Bernstein",5)); //cat1
  funChoice.push_back(pair<string,int>("Bernstein",5)); //cat2
  funChoice.push_back(pair<string,int>("Bernstein",5)); //cat3


  PdfModelBuilder pdfsModel;
  RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hgg_mass");
  mass->setRange(mhLow,mhHigh); 
  pdfsModel.setObsVar(mass);
  //mass->setBins(bins); 
  //mass->setRange(85,110);


  for (int cat=jcats; cat<ncats; cat++){ 

    RooDataSet *data = (RooDataSet*)inWS->data(Form("data_mass_cat%d",cat));
    //RooDataHist thisdataBinned(Form("roohist_data_mass_cat%d",cat),"data",*mass,*dataFull);
    //RooDataSet *data = (RooDataSet*)&thisdataBinned; 

    RooAbsPdf *pdfZee;

    TFile *inFileZee = TFile::Open(fileNameZee.c_str()); 
    RooWorkspace *inWS_Zee = (RooWorkspace*)inFileZee->Get("Zpeak");
    
    cout<<endl<<"///////////////////////////////////"<<endl;
    inWS_Zee->Print();
    cout<<"///////////////////////////////////"<<endl<<endl;
   
    pdfZee = inWS_Zee->pdf(Form("hgg_bkg_%s_cat%d_DCB",ext.c_str(),cat));

    cout<<endl<<"///////////////////////////////////"<<endl;
    pdfZee->Print();
    cout<<"///////////////////////////////////"<<endl<<endl;

    //RooAbsPdf *pdfZeeFixed=pdfsModel.floatDoubleCB(pdfZee,data,Form("hgg_bkg_%s_cat%d_DCB",ext.c_str(),cat));
    RooAbsPdf *pdfZeeFixed=pdfsModel.fixDoubleCB(pdfZee,data,Form("hgg_bkg_%s_cat%d_DCB",ext.c_str(),cat));

    RooAbsPdf  *bkgPdf;
    RooAbsPdf  *berComponent;
    RooAbsPdf  *dcbComponent;


    string funcType = funChoice[cat].first;
    int orderOff = funChoice[cat].second; 
    
    bkgPdf = pdfsModel.getFixedDoubleCBplusContinuum(funcType,Form("hgg_bkg_%s_cat%d",ext.c_str(),cat),orderOff,pdfZeeFixed,false).first;
    berComponent = pdfsModel.getFixedDoubleCBplusContinuum(funcType,Form("hgg_bkg_%s_cat%d",ext.c_str(),cat),orderOff,pdfZeeFixed,false).second.first;
    dcbComponent = pdfsModel.getFixedDoubleCBplusContinuum(funcType,Form("hgg_bkg_%s_cat%d",ext.c_str(),cat),orderOff,pdfZeeFixed,false).second.second;

    cout<<endl<<"///////////////////////////////////"<<endl;
    bkgPdf->Print();
    cout<<"///////////////////////////////////"<<endl<<endl;

  cout<<endl<<"///////////////////////////////////"<<endl;
    berComponent->getComponents()->Print();
    cout<<"///////////////////////////////////"<<endl<<endl;


    cout<<endl<<"///////////////////////////////////"<<endl;
    bkgPdf->getParameters(data)->Print();
    cout<<"///////////////////////////////////"<<endl<<endl;


    int fitStatus = 0; 
    double thisNll=0.;
    runFit(bkgPdf,data,&thisNll,&fitStatus,3,mhLow,mhHigh);
    //RooFitResult *fitRes = bkgPdf->fitTo(*data,Save(true),Range(mhLow,mhHigh));
    //plot(mass,bkgPdf,data,Form("%s/hgg_bkg_%s_cat%d",outDir.c_str(),ext.c_str(),cat),bins,mhLow,mhHigh);

    //print results in log file


    RooArgSet *params = bkgPdf->getParameters(*data);
          float Mean = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_mean",ext.c_str(),cat)))->getValV();
          float MeanErrorL = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_mean",ext.c_str(),cat)))->getErrorLo();
          float MeanErrorH = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_mean",ext.c_str(),cat)))->getErrorHi(); 
          float Sigma = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_sigma",ext.c_str(),cat)))->getValV();
          float SigmaErrorL = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_sigma",ext.c_str(),cat)))->getErrorLo();
          float SigmaErrorH = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_sigma",ext.c_str(),cat)))->getErrorHi();
          float nCB1 = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_nCB1",ext.c_str(),cat)))->getValV();
          float nCB1ErrorL = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_nCB1",ext.c_str(),cat)))->getErrorLo();
          float nCB1ErrorH = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_nCB1",ext.c_str(),cat)))->getErrorHi();
          float nCB2 = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_nCB2",ext.c_str(),cat)))->getValV();
          float nCB2ErrorL = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_nCB2",ext.c_str(),cat)))->getErrorLo();
          float nCB2ErrorH = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_nCB2",ext.c_str(),cat)))->getErrorHi();
          float alphaCB1 = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_alphaCB1",ext.c_str(),cat)))->getValV();
          float alphaCB2 = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_DCB_alphaCB2",ext.c_str(),cat)))->getValV();

	   logfile_stream << Form("**************** cat: %d    ***************",cat) <<endl<<"    nEntries     "<<data->sumEntries() << endl << "    Mean:  " << Mean << "   MeanErrorL:  " << MeanErrorL << " MeanErrorH:  "<< MeanErrorH << "  Sigma:  " << Sigma << "  SigmaErrorL:  " << SigmaErrorL << " SigmaErrorH:  " << SigmaErrorH << "    nCB1:  " << nCB1 << "    nCB1ErrorL:   " << nCB1ErrorL << "   nCB1ErrorH:  " << nCB1ErrorH << "   nCB2:  " << nCB2 << "   nCB2ErrorL:  " << nCB2ErrorL << "    nCB2ErrorH:   " << nCB2ErrorH << "    alphaCB1:   " << alphaCB1 << "    alphaCB2:   " << alphaCB2 << endl;

          float frac = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_frac_sum1",ext.c_str(),cat)))->getValV();
          float fracErrorL = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_frac_sum1",ext.c_str(),cat)))->getErrorLo();
          float fracErrorH = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_frac_sum1",ext.c_str(),cat)))->getErrorHi();

	   logfile_stream <<"    Frac    "<<frac<<"    FracErrorL    "<<fracErrorL<<"    FracErrorH    "<<fracErrorH<<endl;


      for(int i=0;i<orderOff;i++){
          RooRealVar *coeff=(RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_%s_p%i",ext.c_str(),cat,abbr(funcType).c_str(),i));
	  logfile_stream <<Form("     %s_%i:     ",funcType.c_str(),i)<<coeff->getValV()<<Form("     %s_%i ErrorL:     ",funcType.c_str(),i)<<coeff->getErrorLo()<<Form("     %s_%i ErrorH:     ",funcType.c_str(),i)<<coeff->getErrorHi();
	}
       logfile_stream<<endl;

    plot2(mass,bkgPdf,berComponent,dcbComponent,frac,data,Form("%s/hgg_bkg_%s_cat%d",outDir.c_str(),ext.c_str(),cat),bins,mhLow,mhHigh);

    //Let the DCB parameter float within their uncertaincy range in combine ?
  
    outputws->import(*bkgPdf);

    }

    outputfile->cd();
    outputws->Write();
  
    cout<<endl<<"///////////////////////////////////"<<endl;
    //outputws->Print();
    cout<<"///////////////////////////////////"<<endl<<endl;
    
    outputfile->Close();


}





