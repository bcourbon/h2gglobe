#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

#include "TStyle.h"   
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
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooDataHist.h"
#include "RooHist.h"   
#include "TF1.h"   
#include "RooExtendPdf.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TMacro.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TArrow.h"
#include "TKey.h"

#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
//#include "HiggsAnalysis/GBRLikelihood/interface/RooDoubleCB.h"   

#include "../interface/RooDoubleCB.h"
#include "../interface/PdfModelBuilder.h"
#include "../interface/PdfModelBuilderFAN.h" 
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>

using namespace std;
using namespace RooFit;
using namespace boost;

namespace po = program_options;

TRandom3 *RandomGen = new TRandom3();

  bool BLIND = false;

void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries, int mhLow, int mhHigh){

	int ntries=0;
	int stat=1;
	double minnll=10e8;
	while (stat!=0){
	  if (ntries>=MaxTries) break;
	  RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1),Range(mhLow,mhHigh));
	  //RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1),SumW2Error(kTRUE)
          stat = fitTest->status();
	  minnll = fitTest->minNll();
	  ntries++; 
	}
	*stat_t = stat;
	*NLL = minnll;
}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name, int bins=0){  

  double prob;
  int ntoys = 500;
  //int ntoys = 1000;

  name+="_gofTest.png";  
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);

  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);

  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(bins),Name("data"));  

  pdf->plotOn(plot_chi2,Name("pdf"));
  int np = pdf->getParameters(*data)->getSize() - 2;   

  double chi2 = plot_chi2->chiSquare("pdf","data",np);
  std::cout << "Calculating GOF for pdf " << pdf->GetName() << ", using " <<np << " fitted parameters" <<std::endl;

 
  if ((double)data->sumEntries()/bins < 5 ){

    std::cout << "Running toys for GOF test " << std::endl;
    RooArgSet *params = pdf->getParameters(*data);
    RooArgSet preParams;   
    params->snapshot(preParams);
    int ndata = data->sumEntries();

    int npass =0;
    std::vector<double> toy_chi2;
    for (int itoy = 0 ; itoy < ntoys ; itoy++){
      std::cout << Form("\t.. %.1f %% complete\r",100*float(itoy)/ntoys) << std::flush;
      params->assignValueOnly(preParams);
      int nToyEvents = RandomGen->Poisson(ndata);
      RooDataHist *binnedtoy = pdf->generateBinned(RooArgSet(*mass),nToyEvents,0,1);
      pdf->fitTo(*binnedtoy,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Strategy(0));

      RooPlot *plot_t = mass->frame();
      if(bins==0)  binnedtoy->plotOn(plot_t);
      else  binnedtoy->plotOn(plot_t,Binning(bins));  
      pdf->plotOn(plot_t);//,RooFit::NormRange("fitdata_1,fitdata_2"));

      double chi2_t = plot_t->chiSquare(np);
      if( chi2_t>=chi2) npass++;
      toy_chi2.push_back(chi2_t*(bins-np));
      delete plot_t;
    }
    std::cout << "complete" << std::endl;
    prob = (double)npass / ntoys;
 
    TCanvas *can = new TCanvas();
    double medianChi2 = toy_chi2[(int)(((float)ntoys)/2)];
    double rms = TMath::Sqrt(medianChi2);

    TH1F toyhist(Form("gofTest_%s.pdf",pdf->GetName()),";Chi2;",50,medianChi2-5*rms,medianChi2+5*rms);
    for (std::vector<double>::iterator itx = toy_chi2.begin();itx!=toy_chi2.end();itx++){
      toyhist.Fill((*itx));
    }
    toyhist.Draw();

    TArrow lData(chi2*(bins-np),toyhist.GetMaximum(),chi2*(bins-np),0);
    lData.SetLineWidth(2);
    lData.Draw();
    can->SaveAs(name.c_str());

    params->assignValueOnly(preParams);
  } else {
    if(bins == 0)   prob = TMath::Prob(chi2*(bins-np),bins-np);
    else   prob = TMath::Prob(chi2*(bins-np),bins-np); 
    
  }
  if(bins == 0)   std::cout << "Chi2 in Observed =  " << chi2*(bins-np) << std::endl;
  else    std::cout << "Chi2 in Observed =  " << chi2*(bins-np) << std::endl; 
  std::cout << "p-value  =  " << prob << std::endl;
  delete pdf;
  return prob;

}

void GoodnessOfFit (RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, int nFitParam, double &chi2, int &ndf, double &pval, int bins){

RooPlot *plot = mass->frame();//to convert pdf to curve
data->plotOn(plot,Name("data"),Binning(bins));
pdf->plotOn(plot,Name("pdf"));

  // Calculate the chi^2/NDOF of this curve with respect to the histogram
  // 'hist' accounting nFitParam floating parameters in case the curve
  // was the result of a fit
  // Find curve object
  RooCurve* curve = (RooCurve*) plot->findObject("pdf", RooCurve::Class());
  // Find histogram object
  RooHist* hist = (RooHist*) plot->findObject("data", RooHist::Class()) ;
  Int_t i,np = hist->GetN();
  Double_t x,y,xl,xh ;
  // Find starting and ending bin of histogram based on range of RooCurve
  Double_t xstart,xstop ;
#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  curve->GetPoint(0,xstart,y) ;
  curve->GetPoint(curve->GetN()-1,xstop,y) ;
#else
  const_cast<RooCurve*>(curve)->GetPoint(0,xstart,y) ;
  const_cast<RooCurve*>(curve)->GetPoint(curve->GetN() - 1,xstop,y) ;
#endif
  Int_t nbin(0) ;
  Int_t non0bin(0) ;
  chi2=0;
  for (i=0 ; i<np ; i++) {   
    // Retrieve histogram contents
    hist->GetPoint(i,x,y) ;
    xl = x - hist->GetEXlow()[i] ;
    xh = x + hist->GetEXhigh()[i] ;
    // Check if the whole bin is in range of curve
    if (xl < xstart || xstop < xh) continue ;
    nbin++ ;
    // Integrate function over this bin.
    // Start a hack to work around a bug in RooCurve::interpolate
    // that sometimes gives a wrong result.
    Double_t avg = curve->average(xl, xh);
    Double_t avg2 = 0.5 * (curve->average(xl, x) + curve->average(x, xh));
    if (avg + avg2 > 0 &&
	(avg2 - avg) / (avg2 + avg) > 0.1) {
      avg = curve->interpolate(x);
    }
    // End of hack around the bug in RooCurve::interpolate
    // JV: Adjust observed and expected number of events for bin width to represent
    // number of events.
    Double_t norm = (xh - xl) / plot->getFitRangeBinW();
    y *= norm;
    avg *= norm;
    double dy=hist->GetErrorY(i);
    // JV: Use the expected number of events for the y uncertainty,
    // See (33.34) of http://pdg.lbl.gov/2011/reviews/rpp2011-rev-statistics.pdf
    // Add pull^2 to chisq
    if (y != 0) {      
      Double_t resid = y - avg;
      chi2 += pow(resid/dy,2) ;
      non0bin++;
  //cout<<"bin "<<i<<" x "<<x<<" y "<<y<<" pdf "<<avg<<" dy "<<dy<<" (y-pdf)Â²/dyÂ² "<<pow(resid/dy,2)<<" chi2 "<<chi2<<endl;
    }
  }
  ndf=non0bin -nFitParam;
  pval=TMath::Prob(chi2,ndf);
   
}



void plot(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name, int status, int bins, int mhLow, int mhHigh, bool runGOF=false, int nFitParams=1, float alphacb1=1, float alphacb2=1, bool isData=false){
 
 
  double chi2=0;
  double pval=0;
  int ndf=0;

  if(runGOF)   GoodnessOfFit(mass,pdf,data,nFitParams,chi2,ndf,pval,bins);    
 
  RooPlot *plot = mass->frame();
  mass->setRange("unblindReg_1",85,95);
  mass->setRange("unblindReg_2",mhHigh,mhHigh);
  if (BLIND) {
    data->plotOn(plot,Binning(bins),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(bins),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(bins),Invisible());
  }
  else data->plotOn(plot,Binning(bins));
  //else data->plotOn(plot,Binning(bins),DataError(RooAbsData::SumW2));  

 // data->plotOn(plot,Binning(bins));
  TCanvas *canv = new TCanvas();
  TPad *pad =new TPad("haut","haut",0,0.25,1,1);
  pad->SetNumber(1);
  pad->SetBottomMargin(0);
  pad->Draw();
  TPad *pad2 =new TPad("bas","bas",0,0,1,0.23);
  pad2->SetNumber(2);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();

  pdf->plotOn(plot);//,RooFit::NormRange("fitdata_1,fitdata_2"));
  pdf->paramOn(plot,RooFit::Layout(0.35,0.89,0.89),RooFit::Format("NEA",AutoPrecision(1)));  
  canv->cd(1);
  plot->SetMaximum(plot->GetMaximum()*2);    
  if (BLIND) plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->SetXTitle("");
  plot->SetYTitle(Form("Events / %.2fGeV",float(mhHigh-mhLow)/float(bins)));
  //plot->GetYaxis()->SetLabelSize(0.05);
  plot->SetTitleSize(0.045, "Y");
  plot->SetTitleOffset(1.,"Y");
  plot->SetLabelColor(0,"X");
  plot->Draw();

  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  if(runGOF)   lat->DrawLatex(0.2,0.93,Form("chi2 / ndf = %.1f / %i , pval = %.3f",chi2,ndf,pval));  
  lat->DrawLatex(0.68,0.93,"19.7fb^{-1}(8TeV)"); 
  lat->DrawLatex(0.12,0.85,"CMS Preliminary");   
  string nameTemp = name;
  int catNum = nameTemp.find("cat");
  string catNameOLD = nameTemp.replace(0, catNum, "");
  string catName = catNameOLD.replace(catNameOLD.length()-4, 4, "");
  lat->DrawLatex(0.1,0.92, Form("%s", catName.c_str()));  
  cout << Form("Fit Status = %d",status) << endl; 

     TLegend *leg = new TLegend(0.12,0.7,0.32,0.8);
     leg->SetFillColor(0);
     leg->SetLineColor(1);
     if(!isData)   leg->AddEntry(data,"Simulation","lep");
     else    leg->AddEntry(data,"Data","lep");
     leg->Draw("same");

  canv->cd(2);
  RooHist* hPull = plot->pullHist();
  TLine *line = new TLine(mhLow, 0., mhHigh, 0.);  
  line->SetLineWidth(1.5);
  RooPlot *plotPull = mass->frame();
  plotPull->addPlotable(hPull,"P");
  plotPull->addObject(line);
  plotPull->GetXaxis()->SetLabelSize(0.15);
  plotPull->GetYaxis()->SetLabelSize(0.15);
  plotPull->SetXTitle("m_{#gamma#gamma}(GeV)");
  plotPull->SetTitleSize(0.15, "XY");
  plotPull->SetTitleOffset(0.3,"Y");
  plotPull->SetYTitle("Pull");
  //plotPull->SetMarkerSize(0.05);
  plotPull->Draw();


  //pull plots
  TCanvas *canv_Pull = new TCanvas();
  canv_Pull->cd(1);

  const int histoBins = bins;
  int yMaxIndex[histoBins];
  TMath::Sort(histoBins, hPull->GetY(), yMaxIndex, true);
  int MaxIndexNum = 0;
  if(TMath::Abs(hPull->GetY()[yMaxIndex[0]]) < TMath::Abs(hPull->GetY()[yMaxIndex[bins-1]]))  MaxIndexNum = bins-1;
  TH1F *h_pull = new TH1F("h_pull", "h_pull", int(2*TMath::Abs(hPull->GetY()[yMaxIndex[MaxIndexNum]])/0.5), -TMath::Abs(hPull->GetY()[yMaxIndex[MaxIndexNum]]), TMath::Abs(hPull->GetY()[yMaxIndex[MaxIndexNum]]));
  for(int pd=0; pd<=bins-1; pd++){
       double xp_d=0; double yp_d=0;
       hPull->GetPoint(pd, xp_d, yp_d);
       h_pull->Fill(yp_d);
  }
  h_pull->SetXTitle("pull");
  //h_pull->Rebin(2);
  //h_pull->Fit("gaus");
  TF1 *f1_gau = new TF1("gaus","gaus",-TMath::Abs(hPull->GetY()[yMaxIndex[MaxIndexNum]]), TMath::Abs(hPull->GetY()[yMaxIndex[MaxIndexNum]]));
  h_pull->Fit(f1_gau);
  h_pull->SetMaximum(h_pull->GetMaximum()*1.5);     
  f1_gau->SetLineColor(4);
  h_pull->Draw();
  f1_gau->Draw("same");

  TLatex *lat_pull = new TLatex();
  lat_pull->SetNDC();
  lat_pull->SetTextFont(42);
  lat_pull->DrawLatex(0.32,0.82,""); 
  

    canv->SaveAs(Form("%s.png",name.c_str()));
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv_Pull->SaveAs(Form("%s_pull.png", name.c_str())); 
  canv_Pull->SaveAs(Form("%s_pull.pdf", name.c_str())); 
 
  delete canv;
  delete lat;
}

void transferMacros(TFile *inFile, TFile *outFile){
  
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    if (string(key->ReadObj()->ClassName())=="TMacro") {
      TMacro *macro = (TMacro*)inFile->Get(key->GetName());
      outFile->cd();
      macro->Write();
    }
  }
}

int main(int argc, char* argv[]){

gStyle->SetOptFit(111);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

  string fileName;
  int ncats; 
  string logfile;//log
  string outDir;//plots
  string outfilename;//workspace
  int bins; 
  int mhLow; 
  int mhHigh; 
  bool runGOF=true;  
  bool isData=false;  
  bool addExp=false;  
  bool verbose=true;


po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName)->default_value("CMS-HGG_DY.root"),                                              "In file name")
    ("ncats,c", po::value<int>(&ncats)->default_value(4),                                       "Number of categories")
    ("logfile,l", po::value<string>(&logfile)->default_value("Zee_Yield_Zfit.log"),                  "log file of fit results")
    ("outDir,D", po::value<string>(&outDir)->default_value("plots/Zpeak"),                      "Out directory for plots")
    ("outfile,o", po::value<string>(&outfilename)->default_value("CMS-HGG_DY_fit.root"),         		"output file with workspace")
    ("mhLow,L", po::value<int>(&mhLow)->default_value(75),                                                                                                                 "Starting point for scan") 
    ("mhHigh,H", po::value<int>(&mhHigh)->default_value(120),                                                                                                               "End point for scan") 
    ("bins,B", po::value<int>(&bins)->default_value(20),                                                                                                                 "Bins for the plot") 
    ("runGOF",                                                                                  "Run GOF")   
    ("isData",                                                                                  "Run zee data")   
    ("addExp",                                                                                  "add Exponential")   
    ("unblind",  									        "Dont blind plots")
    ("verbose,v",                                                                               "Run with more output")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("runGOF")) runGOF=true;  
  if (vm.count("isData")) isData=true;   
  if (vm.count("addExp")) addExp=true;   
  if (vm.count("unblind")) BLIND=false;

  if (vm.count("verbose")) verbose=true;

 if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
    gErrorIgnoreLevel=kWarning;
  }

  int startingCategory=0;

  std::string ext = "8TeV";

  TFile *outputfile;
  outputfile = new TFile(outfilename.c_str(),"RECREATE");
  RooWorkspace *outputws;
  outputws = new RooWorkspace(); 
  outputws->SetName("Zpeak");
  
  system(Form("mkdir -p %s",outDir.c_str()));
  TFile *inFile = TFile::Open(fileName.c_str());
  RooWorkspace *inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");


	transferMacros(inFile,outputfile);
	RooRealVar *intL  = (RooRealVar*)inWS->var("IntLumi");
	RooRealVar *sqrts = (RooRealVar*)inWS->var("Sqrts");
	outputws->import(*intL);
	outputws->import(*sqrts);

  vector<pair<float,float> > alphaCB12;   
  alphaCB12.push_back(pair<float,float>(1, 1.)); //cat 0
  alphaCB12.push_back(pair<float,float>(0.7, 1.3));   //cat 1
  alphaCB12.push_back(pair<float,float>(0.8, 0.8));   //cat 2
  alphaCB12.push_back(pair<float,float>(0.9, 1.8)); //cat 3

  // store results here
  PdfModelBuilder pdfsModel;
  RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hgg_mass");
  mass->setRange(mhLow,mhHigh); 
  pdfsModel.setObsVar(mass);
 
  ofstream logfile_stream(logfile.c_str());  

for (int cat=startingCategory; cat<ncats; cat++){
    cout << "processing cat " << cat << endl; 
    
    string source = "bkg";
    if(isData) source = "data";
    RooDataSet *data = (RooDataSet*)inWS->data(Form("%s_mass_cat%d",source.c_str(),cat));
    
    //mass->setBins(bins);
    //RooDataHist thisdataBinned(Form("roohist_%s_mass_cat%d",source.c_str(),cat),"data",*mass,*dataFull);
    //RooDataSet *data = (RooDataSet*)&thisdataBinned;


    if(data->sumEntries() <= 1e-6)  continue;

    RooAbsPdf  *bkgPdf;

    float alphacb1Comm = 0;
    float alphacb2Comm = 0;
    alphacb1Comm = alphaCB12[cat].first;
    alphacb2Comm = alphaCB12[cat].second;


    if(!addExp)
	bkgPdf = pdfsModel.getDoubleCB(Form("hgg_bkg_%s_cat%d",ext.c_str(),cat), alphacb1Comm, alphacb2Comm);//fix alpha1 and alpha2
	
    else 
	bkgPdf = pdfsModel.getDoubleCBplusContinuum("Exponential",Form("hgg_bkg_%s_cat%d",ext.c_str(),cat),1,true).first;//float all parameters
	
   
    int fitStatus = 0; 
    double thisNll=0.;
    runFit(bkgPdf,data,&thisNll,&fitStatus,3,mhLow,mhHigh);

    if(!addExp) plot(mass,bkgPdf,data,Form("%s/Zee_DCB_cat%d",outDir.c_str(),cat), fitStatus,bins,mhLow,mhHigh,runGOF,4,alphacb1Comm,alphacb2Comm,isData);
else plot(mass,bkgPdf,data,Form("%s/Zee_DCB_exp_cat%d",outDir.c_str(),cat), fitStatus,bins,mhLow,mhHigh,runGOF,6,alphacb1Comm,alphacb2Comm,isData);


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

          logfile_stream << Form("cat: %d    ",cat) <<"    nEntries     "<<data->sumEntries() << "    Mean:  " << Mean << "   MeanErrorL:  " << MeanErrorL << " MeanErrorH:  "<< MeanErrorH << "  Sigma:  " << Sigma << "  SigmaErrorL:  " << SigmaErrorL << " SigmaErrorH:  " << SigmaErrorH << "    nCB1:  " << nCB1 << "    nCB1ErrorL:   " << nCB1ErrorL << "   nCB1ErrorH:  " << nCB1ErrorH << "   nCB2:  " << nCB2 << "   nCB2ErrorL:  " << nCB2ErrorL << "    nCB2ErrorH:   " << nCB2ErrorH << "    alphaCB1:   " << alphaCB1 << "    alphaCB2:   " << alphaCB2 << endl;

if(addExp){
              float expParam = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_exp_p1",ext.c_str(),cat)))->getValV(); 
              float expParamErrorL = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_exp_p1",ext.c_str(),cat)))->getErrorLo();
              float expParamErrorH = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_exp_p1",ext.c_str(),cat)))->getErrorHi();
              float expCBFrac = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_frac_sum1",ext.c_str(),cat)))->getValV();
              float expCBFracErrorL = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_frac_sum1",ext.c_str(),cat)))->getErrorLo();
              float expCBFracErrorH = ((RooRealVar*)params->find(Form("hgg_bkg_%s_cat%d_frac_sum1",ext.c_str(),cat)))->getErrorHi();
 
              logfile_stream << Form("cat: %d    ",cat) <<"   expParam:  " << expParam << "   expParamErrorL:  " << expParamErrorL << "   expParamErrorH: " << expParamErrorH << "   expCBFrac:  " << expCBFrac << "   expCBFracErrorL:  " << expCBFracErrorL << "   expCBFracErrorH:  " << expCBFracErrorH << endl;
}


    logfile_stream << endl;
    params->Print("v");

    outputws->import(*bkgPdf);
    outputws->import(*data);

  }


    outputfile->cd();
    outputws->Write();

    cout<<endl<<"///////////////////////////////////"<<endl;
    outputws->Print();
    cout<<"///////////////////////////////////"<<endl<<endl;
    
    outputfile->Close();


}

