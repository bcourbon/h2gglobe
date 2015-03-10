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

#include "../interface/PdfModelBuilder.h"
#include "../interface/PdfModelBuilderFAN.h" 
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>

using namespace std;
using namespace RooFit;
using namespace boost;

namespace po = program_options;

bool BLIND = true;
bool runFtestCheckWithToys=false;
int nBinsForMass = 320;

TRandom3 *RandomGen = new TRandom3();

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext=""){
  
  if (type=="Bernstein") return pdfsModel.getBernstein(Form("%s_bern%d",ext,order),order); 
  else if (type=="Chebychev") return pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order); 
  else if (type=="Exponential") return pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order); 
  else if (type=="PowerLaw") return pdfsModel.getPowerLawSingle(Form("%s_pow%d",ext,order),order); 
  else if (type=="Laurent") return pdfsModel.getLaurentSeries(Form("%s_lau%d",ext,order),order); 
  else {
    cerr << "ERROR -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}





pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getPdfLowVoi(PdfModelBuilderFAN &pdfsModel, string type, int orderOff, int orderVoig, const char* ext="", float constant=0, float bernDownBound=-11., float bernUpBound=11.){

  const char* type1 = type.c_str();
  return pdfsModel.getOfficialSumVoigtiansForZeeFit(type, ext, Form("%s_plus_%s_%d",ext,type1,orderOff), orderOff, orderVoig, constant, bernDownBound, bernUpBound);

}





pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getPdfLowDoubleCB(PdfModelBuilderFAN &pdfsModel, string type, int orderOff, int orderVoig, const char* ext="", float alphacb1=1., float alphacb2=1., float constant=0, float bernDownBound=-11., float bernUpBound=11.){

  const char* type1 = type.c_str();
  return pdfsModel.getOfficialSumVoigtiansForZeeFitDoubleCB(type, ext, Form("%s_plus_%s_%d",ext,type1,orderOff), orderOff, orderVoig, constant, bernDownBound, bernUpBound, alphacb1, alphacb2);

}


pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getPdfLowDoubleCB_Float(PdfModelBuilderFAN &pdfsModel, string type, int orderOff, int orderVoig, const char* ext="", float alphacb1=1., float alphacb2=1., float constant=0, float bernDownBound=-11., float bernUpBound=11.){

  const char* type1 = type.c_str();
  return pdfsModel.getOfficialSumVoigtiansForZeeFitDoubleCB_Float(type, ext, Form("%s_plus_%s_%d",ext,type1,orderOff), orderOff, orderVoig, constant, bernDownBound, bernUpBound, alphacb1, alphacb2);

}


void runFitFAN(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries, int mhLow, int mhHigh){

	int ntries=0;
  	RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
	//std::cout << " BEFORE ITERATIONS-------------------------------" << std::endl;
	//params_test->Print("v");
	int stat=1;
	double minnll=10e8;
	while (stat!=0){
	  if (ntries>=MaxTries) break;
	  //std::cout << "----------------------------- BEFORE FIT-------------------------------" << std::endl;
	  //params_test->Print("v");
	  //std::cout << "-----------------------------------------------------------------------" << std::endl;
	  RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1)
	  //RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1),SumW2Error(kTRUE)
 	  	,Range(mhLow,mhHigh));   
          stat = fitTest->status();
	  minnll = fitTest->minNll();
	  if (stat!=0) params_test->assignValueOnly(fitTest->randomizePars());
	  ntries++; 
	}
	*stat_t = stat;
	*NLL = minnll;
}



void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries){

	int ntries=0;
  	RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
	int stat=1;
	double minnll=10e8;
	while (stat!=0){
	  if (ntries>=MaxTries) break;
	  RooFitResult *fitTest = pdf->fitTo(*data,RooFit::Save(1)
		,RooFit::Minimizer("Minuit2","minimize"));
          stat = fitTest->status();
	  minnll = fitTest->minNll();
	  if (stat!=0) params_test->assignValueOnly(fitTest->randomizePars());
	  ntries++; 
	}
	*stat_t = stat;
	*NLL = minnll;
}
double getProbabilityFtest(double chi2, int ndof,RooAbsPdf *pdfNull, RooAbsPdf *pdfTest, RooRealVar *mass, RooDataSet *data, std::string name){
 
  double prob_asym = TMath::Prob(chi2,ndof);
  if (!runFtestCheckWithToys) return prob_asym;

  int ndata = data->sumEntries();
  
  RooFitResult *fitNullData = pdfNull->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1));
  RooFitResult *fitTestData = pdfTest->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1));

  RooArgSet *params_null = pdfNull->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_null;
  params_null->snapshot(preParams_null);
  RooArgSet *params_test = pdfTest->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_test;
  params_test->snapshot(preParams_test);
 
  int ntoys =5000;
  TCanvas *can = new TCanvas();
  can->SetLogy();
  TH1F toyhist(Form("toys_fTest_%s.pdf",pdfNull->GetName()),";Chi2;",60,-2,10);
  TH1I toyhistStatN(Form("Status_%s.pdf",pdfNull->GetName()),";FitStatus;",8,-4,4);
  TH1I toyhistStatT(Form("Status_%s.pdf",pdfTest->GetName()),";FitStatus;",8,-4,4);

  TGraph *gChi2 = new TGraph();
  gChi2->SetLineColor(kGreen+2);
  double w = toyhist.GetBinWidth(1);

  int ipoint=0;

  for (int b=0;b<toyhist.GetNbinsX();b++){
	double x = toyhist.GetBinCenter(b+1);
	if (x>0){
	  gChi2->SetPoint(ipoint,x,(ROOT::Math::chisquared_pdf(x,ndof)));
	  ipoint++;
	}
  }
  int npass =0; int nsuccesst =0;
  mass->setBins(nBinsForMass);
  for (int itoy = 0 ; itoy < ntoys ; itoy++){

        params_null->assignValueOnly(preParams_null);
        params_test->assignValueOnly(preParams_test);
  	RooDataHist *binnedtoy = pdfNull->generateBinned(RooArgSet(*mass),ndata,0,1);

	int stat_n=1;
        int stat_t=1;
	int ntries = 0;
	double nllNull,nllTest;
	// Iterate on the fit 
	int MaxTries = 2;
	while (stat_n!=0){
	  if (ntries>=MaxTries) break;
	  RooFitResult *fitNull = pdfNull->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
		//,RooFit::Optimize(0));

	  nllNull = fitNull->minNll();
          stat_n = fitNull->status();
	  if (stat_n!=0) params_null->assignValueOnly(fitNullData->randomizePars());
	  ntries++; 
	}
	
	ntries = 0;
	while (stat_t!=0){
	  if (ntries>=MaxTries) break;
	  RooFitResult *fitTest = pdfTest->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
	  nllTest = fitTest->minNll();
          stat_t = fitTest->status();
	  if (stat_t!=0) params_test->assignValueOnly(fitTestData->randomizePars()); 
	  ntries++; 
	}
       
	toyhistStatN.Fill(stat_n);
	toyhistStatT.Fill(stat_t);

        if (stat_t !=0 || stat_n !=0) continue;
	nsuccesst++;
	double chi2_t = 2*(nllNull-nllTest);
	if (chi2_t >= chi2) npass++;
        toyhist.Fill(chi2_t);
  }

  double prob=0;
  if (nsuccesst!=0)  prob = (double)npass / nsuccesst;
  toyhist.Scale(1./(w*toyhist.Integral()));
  toyhist.Draw();
  TArrow lData(chi2,toyhist.GetMaximum(),chi2,0);
  lData.SetLineWidth(2);
  lData.Draw();
  gChi2->Draw("L");
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.91,Form("Prob (asymptotic) = %.4f (%.4f)",prob,prob_asym));
  can->SaveAs(name.c_str());

  TCanvas *stas =new TCanvas();
  toyhistStatN.SetLineColor(2);
  toyhistStatT.SetLineColor(1); 
  TLegend *leg = new TLegend(0.2,0.6,0.4,0.89); leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(&toyhistStatN,"Null Hyp","L");
  leg->AddEntry(&toyhistStatT,"Test Hyp","L");
  toyhistStatN.Draw();
  toyhistStatT.Draw("same");
  leg->Draw();
  //stas->SaveAs(Form("%s_fitstatus.pdf",name.c_str()));
  stas->SaveAs(Form("%s_fitstatus.png",name.c_str()));  
  //reassign params
  params_null->assignValueOnly(preParams_null);
  params_test->assignValueOnly(preParams_test);

  delete can; delete stas;
  delete gChi2;
  delete leg;
  delete lat;

  // Still return the asymptotic prob (usually its close to the toys one)
  return prob_asym;
  //return prob; 

}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name, int bins=0){  

  double prob;
  int ntoys = 500;
  //int ntoys = 1000;

  name+="_gofTest.png";  
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);

  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);

  RooPlot *plot_chi2 = mass->frame();
  if(bins == 0)   data->plotOn(plot_chi2,Binning(nBinsForMass),Name("data"));
  else   data->plotOn(plot_chi2,Binning(bins),Name("data"));  

  pdf->plotOn(plot_chi2,Name("pdf"));
  int np = pdf->getParameters(*data)->getSize() - 2;   

  double chi2 = plot_chi2->chiSquare("pdf","data",np);
  std::cout << "Calculating GOF for pdf " << pdf->GetName() << ", using " <<np << " fitted parameters" <<std::endl;

 
  if ((double)data->sumEntries()/nBinsForMass < 5 ){

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
      toy_chi2.push_back(chi2_t*(nBinsForMass-np));
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

    TArrow lData(chi2*(nBinsForMass-np),toyhist.GetMaximum(),chi2*(nBinsForMass-np),0);
    lData.SetLineWidth(2);
    lData.Draw();
    can->SaveAs(name.c_str());

    params->assignValueOnly(preParams);
  } else {
    if(bins == 0)   prob = TMath::Prob(chi2*(nBinsForMass-np),nBinsForMass-np);
    else   prob = TMath::Prob(chi2*(bins-np),bins-np); 
    
  }
  if(bins == 0)   std::cout << "Chi2 in Observed =  " << chi2*(nBinsForMass-np) << std::endl;
  else    std::cout << "Chi2 in Observed =  " << chi2*(bins-np) << std::endl; 
  std::cout << "p-value  =  " << prob << std::endl;
  delete pdf;
  return prob;

}


void plot(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name, int status, double *prob){
  
  // Chi2 taken from full range fit
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(nBinsForMass));
  pdf->plotOn(plot_chi2);

  int np = pdf->getParameters(*data)->getSize()+1; //Because this pdf has no extend
  double chi2 = plot_chi2->chiSquare(np);
 
  *prob = getGoodnessOfFit(mass,pdf,data,name);
  RooPlot *plot = mass->frame();
  mass->setRange("unblindReg_1",100,110);
  mass->setRange("unblindReg_2",150,180);
  if (BLIND) {
    data->plotOn(plot,Binning(80),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(80),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(80),Invisible());
  }
  else data->plotOn(plot,Binning(80));

  TCanvas *canv = new TCanvas();
  pdf->plotOn(plot);//,RooFit::NormRange("fitdata_1,fitdata_2"));
  pdf->paramOn(plot,RooFit::Layout(0.34,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->Draw();

  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.92,Form("#chi^{2} = %.3f, Prob = %.2f, Fit Status = %d ",chi2*(nBinsForMass-np),*prob,status));
  canv->SaveAs(name.c_str());
  
  delete canv;
  delete lat;
}
void plot(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooDataSet *data, string name, int cat, int bestFitPdf=-1){
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.5,0.55,0.92,0.92);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",100,110);
  mass->setRange("unblindReg_2",150,180);
  if (BLIND) {
    data->plotOn(plot,Binning(80),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(80),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(80),Invisible());
  }
  else data->plotOn(plot,Binning(80)); 

  int currentIndex = catIndex->getIndex();
  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,Form("Data - cat%d",cat),"LEP");
  int style=1;
  for (int icat=0;icat<catIndex->numTypes();icat++){
    int col;
    if (icat<=6) col=color[icat];
    else {col=kBlack; style++;}
    catIndex->setIndex(icat);
    pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"));	
    pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style));//,RooFit::NormRange("fitdata_1,fitdata_2"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==icat) ext=" (Best Fit Pdf) ";
    leg->AddEntry(pdfLeg,Form("%s%s",pdfs->getCurrentPdf()->GetName(),ext.c_str()),"L");
  }
  plot->SetTitle(Form("Category %d",cat));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  canv->SaveAs(Form("%s.png",name.c_str()));
  catIndex->setIndex(currentIndex);
  delete canv;
}

void plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooDataSet *data, string name, int cat, int bestFitPdf=-1){
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.6,0.65,0.89,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",100,110);
  mass->setRange("unblindReg_2",150,180);
  if (BLIND) {
    data->plotOn(plot,Binning(80),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(80),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(80),Invisible());
  }
  else data->plotOn(plot,Binning(80));

  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,Form("Data - cat%d",cat),"LEP");
  int i=0;
  int style=1;
  for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
    int col;
    if (i<=6) col=color[i];
    else {col=kBlack; style++;}
    it->second->plotOn(plot,LineColor(col),LineStyle(style));//,RooFit::NormRange("fitdata_1,fitdata_2"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==i) ext=" (Best Fit Pdf) ";
    leg->AddEntry(pdfLeg,Form("%s%s",it->first.c_str(),ext.c_str()),"L");
    i++;
  }
  plot->SetTitle(Form("Category %d",cat));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  canv->SaveAs(Form("%s.png",name.c_str()));
  delete canv;
}



void plotFAN(RooRealVar *mass, RooAbsPdf *pdf, RooAbsPdf *pdf_gof, RooDataSet *data, string name, int status, double *prob, int bins=0, int mhLow=0, int mhHigh=0, bool runGOF=false, bool useDoubleCB=false, float alphacb1=1, float alphacb2=1, bool isData=false, bool runWrong=false, bool runTDR=false){
 
 
  RooPlot *plot_chi2 = mass->frame();
  if(!runWrong)  data->plotOn(plot_chi2,Binning(nBinsForMass));
  else  data->plotOn(plot_chi2,Binning(bins));   
  pdf->plotOn(plot_chi2);

  int np = pdf->getParameters(*data)->getSize()+1-2; 
  double chi2 = plot_chi2->chiSquare(np);

  if(!runWrong){
      if(runGOF)   *prob = getGoodnessOfFit(mass,pdf_gof,data,name);    
  }else{
      if(runGOF)   *prob = getGoodnessOfFit(mass,pdf_gof,data,name,bins);    
  }
  

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
  if(!runTDR)   pdf->paramOn(plot,RooFit::Layout(0.28,0.96,0.95),RooFit::Format("NEA",AutoPrecision(1)));  
  canv->cd(1);
  if(!runTDR)   plot->SetMaximum(plot->GetMaximum()*2);    
  else  plot->SetMaximum(plot->GetMaximum()*1.2);    
  if (BLIND) plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->SetXTitle("");
  plot->SetYTitle(Form("Events / %.2fGeV",float(mhHigh-mhLow)/float(bins)));
  plot->SetTitleSize(0.045, "Y");
  plot->SetTitleOffset(1.5,"Y");
  plot->SetLabelColor(0,"X");
  plot->Draw();

  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  if(runGOF && useDoubleCB && !runTDR)   lat->DrawLatex(0.55,0.42,Form("alpha1 = %.1f, alpha2 = %.1f",alphacb1,alphacb2));  
  if(runTDR)   lat->DrawLatex(0.68,0.93,"19.7fb^{-1}(8TeV)"); 
  if(runTDR)   lat->DrawLatex(0.18,0.85,"CMS Preliminary");   
  string nameTemp = name;
  int catNum = nameTemp.find("cat");
  string catNameOLD = nameTemp.replace(0, catNum, "");
  string catName = catNameOLD.replace(catNameOLD.length()-4, 4, "");
  if(!runTDR)  lat->DrawLatex(0.1,0.92, Form("%s", catName.c_str()));  
  if(!runWrong)  cout << "chi2  " << chi2*(nBinsForMass-np) << endl;    
  else  cout << "chi2  " << chi2*(bins-np) << endl;      
  cout << Form("Fit Status = %d",status) << endl; 

  if(runTDR){
     TLegend *leg = new TLegend(0.55,0.55,0.85,0.75);
     leg->SetFillColor(0);
     leg->SetLineColor(1);
     if(!isData)   leg->AddEntry(data,"Simulation","lep");
     else    leg->AddEntry(data,"Data","lep");
     leg->Draw("same");
  }


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
  

  canv->SaveAs(name.c_str());
  if(runTDR)   canv->SaveAs(Form("%s.pdf", name.substr(0, name.length()-4).c_str()));
  canv_Pull->SaveAs(Form("%s_pull.png", name.substr(0, name.length()-4).c_str())); 
  if(runTDR)   canv_Pull->SaveAs(Form("%s_pull.pdf", name.substr(0, name.length()-4).c_str())); 
 
  delete canv;
  delete lat;
}

void plotFAN(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooDataSet *data, string name, int cat, int bestFitPdf=-1, int bins=0, int mhLow=0, int mhHigh=0){
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.5,0.55,0.92,0.92);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",85,95);
  mass->setRange("unblindReg_2",mhHigh,mhHigh);
  if (BLIND) {
    data->plotOn(plot,Binning(bins),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(bins),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(bins),Invisible());
  }
  else data->plotOn(plot,Binning(bins)); 

  int currentIndex = catIndex->getIndex();
  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,Form("Data - cat%d",cat),"LEP");
  int style=1;
  for (int icat=0;icat<catIndex->numTypes();icat++){
    int col;
    if (icat<=6) col=color[icat];
    else {col=kBlack; style++;}
    catIndex->setIndex(icat);
    pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"));	
    pdfs->getCurrentPdf()->plotOn(plot,LineColor(col),LineStyle(style));//,RooFit::NormRange("fitdata_1,fitdata_2"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==icat) ext=" (Best Fit Pdf) ";
    leg->AddEntry(pdfLeg,Form("%s%s",pdfs->getCurrentPdf()->GetName(),ext.c_str()),"L");
  }
  plot->SetTitle(Form("Category %d",cat));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  canv->SaveAs(Form("%s.png",name.c_str()));
  catIndex->setIndex(currentIndex);
  delete canv;
}

void plotFAN(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooDataSet *data, string name, int cat, int bestFitPdf=-1, int bins=0, int mhHigh=0){
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.6,0.65,0.89,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",85,95);
  mass->setRange("unblindReg_2",mhHigh,mhHigh);
  if (BLIND) {
    data->plotOn(plot,Binning(bins),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(bins),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(bins),Invisible());
  }
  else data->plotOn(plot,Binning(bins));

  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,Form("Data - cat%d",cat),"LEP");
  int i=0;
  int style=1;
  for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
    int col;
    if (i<=6) col=color[i];
    else {col=kBlack; style++;}
    it->second->plotOn(plot,LineColor(col),LineStyle(style));//,RooFit::NormRange("fitdata_1,fitdata_2"));
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==i) ext=" (Best Fit Pdf) ";
    leg->AddEntry(pdfLeg,Form("%s%s",it->first.c_str(),ext.c_str()),"L");
    i++;
  }
  plot->SetTitle(Form("Category %d",cat));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  canv->SaveAs(Form("%s.png",name.c_str()));
  delete canv;
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
int getBestFitFunction(RooMultiPdf *bkg, RooDataSet *data, RooCategory *cat, bool silent=false){


	double global_minNll = 1E10;
	int best_index = 0;
	int number_of_indeces = cat->numTypes();
		
	RooArgSet snap,clean;
	RooArgSet *params = bkg->getParameters((const RooArgSet*)0);
	params->remove(*cat);
	params->snapshot(snap);
	params->snapshot(clean);
	if (!silent) {
		std::cout << "CLEAN SET OF PARAMETERS" << std::endl;
		//params->Print("V");
		std::cout << "-----------------------" << std::endl;
	}
	
	
	for (int id=0;id<number_of_indeces;id++){		
		params->assignValueOnly(clean);
		cat->setIndex(id);


		if (!silent) {
		}
		
		//minim.minimize("Minuit2","minimize");
		double minNll=0; //(nllm->getVal())+bkg->getCorrection();
		int fitStatus=1;		
		runFit(bkg->getCurrentPdf(),data,&minNll,&fitStatus,/*max iterations*/3);
		// Add the penalty

		minNll=minNll+bkg->getCorrection();

		if (!silent) {
			/*
			std::cout << "After Minimization ------------------  " <<std::endl;
			std::cout << bkg->getCurrentPdf()->GetName() << " " << minNll <<std::endl;
			bkg->Print("v");
			bkg->getCurrentPdf()->getParameters(*data)->Print("V");
			std::cout << " ------------------------------------  " << std::endl;
	
			params->Print("V");
			*/
			std::cout << "AFTER FITTING" << std::endl;
			std::cout << " Function was " << bkg->getCurrentPdf()->GetName() <<std::endl;
			std::cout << " Correction Applied is " << bkg->getCorrection() <<std::endl;
			std::cout << " NLL + c = " <<  minNll << std::endl;
			std::cout << "-----------------------" << std::endl;
		}
			
		if (minNll < global_minNll){
        		global_minNll = minNll;
			snap.assignValueOnly(*params);
        		best_index=id;
		}
	}
    	cat->setIndex(best_index);
	params->assignValueOnly(snap);
	
	if (!silent) {
		std::cout << "Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<std::endl;
		//bkg->getCurrentPdf()->getParameters(*data)->Print("v");
	}
	return best_index;
}




int main(int argc, char* argv[]){

  //bool runTDR=true;
  bool runTDR=false;

  gStyle->SetOptFit(111);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

 
  string fileName;
  int ncats;
  int singleCategory;
  int forLoop;   
  float alphacb1;   
  float alphacb2;   
  string datfile;
  string outDir;
  string outfilename;
  string runType; 
  int bins; 
  int mhLow; 
  int mhHigh; 
  bool is2011=false;
  bool runGOF=false;  
  bool isData=false;  
  bool isFloat=false;  
  bool useDoubleCB=false;  
  bool addExp=false;  
  bool verbose=false;
  bool saveZeePdf=false;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName),                                              "In file name")
    ("ncats,c", po::value<int>(&ncats)->default_value(5),                                       "Number of categories")
    ("singleCat", po::value<int>(&singleCategory)->default_value(-1),                           "Run A single Category")
    ("datfile,d", po::value<string>(&datfile)->default_value("Zee_Yield_Zfit.log"),                  "log file of fit results")
    ("outDir,D", po::value<string>(&outDir)->default_value("plots/fTest"),                      "Out directory for plots")
    ("saveZeePdf", po::value<string>(&outfilename),         					"Save a ZeePdf model with the appropriate pdfs")
    ("runType,R", po::value<string>(&runType)->default_value(""),                               "Official or low mass") 
    ("mhLow,L", po::value<int>(&mhLow)->default_value(100),                                                                                                                 "Starting point for scan") 
    ("mhHigh,H", po::value<int>(&mhHigh)->default_value(180),                                                                                                               "End point for scan") 
    ("bins,B", po::value<int>(&bins)->default_value(80),                                                                                                                 "Bins for the plot") 
    ("nBinsForMass,b", po::value<int>(&nBinsForMass)->default_value(320),                                                                                               "nBinsForMass") 
    ("forLoop", po::value<int>(&forLoop)->default_value(-1),                           "loop index when choose alphacb1 and alphacb2")
    ("alphacb1", po::value<float>(&alphacb1)->default_value(1.),                           "alphacb1 for doule CB")
    ("alphacb2", po::value<float>(&alphacb2)->default_value(1.),                           "alphacb2 for double CB")
    ("runFtestCheckWithToys", 									"When running the F-test, use toys to calculate pvals (and make plots) ")
    ("is2011",                                                                                  "Run 2011 config")
    ("runGOF",                                                                                  "Run GOF")   
    ("isData",                                                                                  "Run zee data")   
    ("isFloat",                                                                                  "Float alpha1 and alpha2 or so")   
    ("useDoubleCB",                                                                                  "use double crystal ball function")   
    ("addExp",                                                                                  "add Exponential")   
    ("unblind",  									        "Dont blind plots")
    ("verbose,v",                                                                               "Run with more output")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("is2011")) is2011=true;
  if (vm.count("runGOF")) runGOF=true;  
  if (vm.count("isData")) isData=true;   
  if (vm.count("isFloat")) isFloat=true;   
  if (vm.count("useDoubleCB"))  useDoubleCB=true;   
  if (vm.count("addExp")) addExp=true;   
  if (vm.count("unblind")) BLIND=false;
  saveZeePdf = vm.count("saveZeePdf");

  if (vm.count("verbose")) verbose=true;
  if (vm.count("runFtestCheckWithToys")) runFtestCheckWithToys=true;

  if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
    gErrorIgnoreLevel=kWarning;
  }
  int startingCategory=0;
  if (singleCategory >-1){
	ncats=singleCategory+1;	
	startingCategory=singleCategory;
  }

  std::cout << "SaveZeePdf? " << saveZeePdf << std::endl;
  TFile *outputfile;
  RooWorkspace *outputws;

  outputfile = new TFile(outfilename.c_str(),"RECREATE");
  outputws = new RooWorkspace(); outputws->SetName("fTestVoi_Zee");

  system(Form("mkdir -p %s",outDir.c_str()));
  TFile *inFile = TFile::Open(fileName.c_str());
  RooWorkspace *inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");

  if (saveZeePdf){
	transferMacros(inFile,outputfile);
	RooRealVar *intL  = (RooRealVar*)inWS->var("IntLumi");
	RooRealVar *sqrts = (RooRealVar*)inWS->var("Sqrts");
	outputws->import(*intL);
	outputws->import(*sqrts);
  }
  
  vector<string> functionClasses;
  functionClasses.push_back("Bernstein");
  functionClasses.push_back("Exponential");
  functionClasses.push_back("PowerLaw");
  functionClasses.push_back("Laurent");
  map<string,string> namingMap;
  namingMap.insert(pair<string,string>("Bernstein","pol"));
  namingMap.insert(pair<string,string>("Exponential","exp"));
  namingMap.insert(pair<string,string>("PowerLaw","pow"));
  namingMap.insert(pair<string,string>("Laurent","lau"));
  vector<pair<string,int> > fabChoice;
  vector<pair<float,float> > alphaCB12;   

  if (is2011) {
    fabChoice.push_back(pair<string,int>("Bernstein",4));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
  }
  else {
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",4));
    fabChoice.push_back(pair<string,int>("Bernstein",4));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",5));
    fabChoice.push_back(pair<string,int>("Bernstein",4));
    fabChoice.push_back(pair<string,int>("Bernstein",4));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
    fabChoice.push_back(pair<string,int>("Bernstein",3));
  }

  
  if(useDoubleCB && alphacb1 == -1 && alphacb1 == -1 ){
           alphaCB12.push_back(pair<float,float>(1, 1.)); //0
           alphaCB12.push_back(pair<float,float>(0.7, 1.3));   //1
           alphaCB12.push_back(pair<float,float>(0.8, 0.8));   //2
           alphaCB12.push_back(pair<float,float>(0.9, 1.8)); //3
           alphaCB12.push_back(pair<float,float>(1, 1));   //4
           //alphaCB12.push_back(pair<float,float>(0.7, 1.4));   //4
           alphaCB12.push_back(pair<float,float>(1, 1));   //5
           alphaCB12.push_back(pair<float,float>(1, 1));   //6
           alphaCB12.push_back(pair<float,float>(1, 1));   //7
           alphaCB12.push_back(pair<float,float>(1, 1));   //8
           alphaCB12.push_back(pair<float,float>(1, 1));   //9
           alphaCB12.push_back(pair<float,float>(1, 1));   //10
           alphaCB12.push_back(pair<float,float>(1, 1));   //11
           alphaCB12.push_back(pair<float,float>(1, 1));   //12
           alphaCB12.push_back(pair<float,float>(1, 1));   //13
           //alphaCB12.push_back(pair<float,float>(1, 1.4));

  }
  


  // store results here
  PdfModelBuilder pdfsModel;
  PdfModelBuilderFAN pdfsModelFAN; 
  RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hgg_mass");
  mass->setRange(mhLow,mhHigh); 
  pdfsModel.setObsVar(mass);
  pdfsModelFAN.setObsVar(mass); 
 
 
  std::string ext = is2011 ? "7TeV" : "8TeV";
  ofstream outfile(datfile.c_str());  

  for (int cat=startingCategory; cat<ncats; cat++){
    cout << "FAN processing cat " << cat << "..." << "  forLoop   " << forLoop << endl; 
    
    
    string dataorMC = "bkg";
    if(isData) dataorMC = "data";
    RooDataSet *dataFull = (RooDataSet*)inWS->data(Form("%s_mass_cat%d",dataorMC.c_str(),cat));
    

    mass->setBins(nBinsForMass);
    RooDataHist thisdataBinned(Form("roohist_%s_mass_cat%d",dataorMC.c_str(),cat),"data",*mass,*dataFull);
    RooDataSet *data = (RooDataSet*)&thisdataBinned;

    if(data->sumEntries() <= 1e-6)  continue;
  
 
    //fitting start
    RooAbsPdf  *bkgPdf;
    RooAbsPdf  *bkgPdf_gof;
    float alphacb1Comm = 0; float alphacb2Comm = 0;
    if(!useDoubleCB){ 
           bkgPdf = pdfsModelFAN.getSumOfVoigtians(Form("ftest_Zee_Voi_%s_cat%d",ext.c_str(),cat), 1, false);  
           bkgPdf_gof = pdfsModelFAN.getSumOfVoigtians(Form("ftest_Zee_Voi_gof_%s_cat%d",ext.c_str(),cat), 1, false);  
           if(addExp)   bkgPdf = getPdfLowVoi(pdfsModelFAN, "Exponential", 1, 1, Form("ftest_Zee_Voi_%s_cat%d",ext.c_str(),cat)).first;  
           if(addExp)   bkgPdf_gof = getPdfLowVoi(pdfsModelFAN, "Exponential", 1, 1, Form("ftest_Zee_Voi_gof_%s_cat%d",ext.c_str(),cat)).first;  
    }else{
           if(alphacb1 == -1 && alphacb2 == -1){
                 alphacb1Comm = alphaCB12[cat].first;
                 alphacb2Comm = alphaCB12[cat].second;
           }else{
                 alphacb1Comm = alphacb1;
                 alphacb2Comm = alphacb2;
           } 

           if(!isFloat){
                bkgPdf = pdfsModelFAN.getConvDoubleCBnew(Form("ftest_Zee_DCB_%s_cat%d",ext.c_str(),cat), alphacb1Comm, alphacb2Comm);  
                bkgPdf_gof = pdfsModelFAN.getConvDoubleCBnew(Form("ftest_Zee_DCB_gof_%s_cat%d",ext.c_str(),cat), alphacb1Comm, alphacb2Comm);  
                if(addExp)   bkgPdf = getPdfLowDoubleCB(pdfsModelFAN, "Exponential", 1, 1, Form("ftest_Zee_DCB_%s_cat%d",ext.c_str(),cat), alphacb1Comm, alphacb2Comm).first;   
                if(addExp)   bkgPdf_gof = getPdfLowDoubleCB(pdfsModelFAN, "Exponential", 1, 1, Form("ftest_Zee_DCB_gof_%s_cat%d",ext.c_str(),cat), alphacb1Comm, alphacb2Comm).first;
            }else{
                bkgPdf = pdfsModelFAN.getConvDoubleCBnew_Float(Form("ftest_Zee_DCB_%s_cat%d",ext.c_str(),cat), alphacb1Comm, alphacb2Comm);  
                bkgPdf_gof = pdfsModelFAN.getConvDoubleCBnew_Float(Form("ftest_Zee_DCB_gof_%s_cat%d",ext.c_str(),cat), alphacb1Comm, alphacb2Comm);  
                if(addExp)   bkgPdf = getPdfLowDoubleCB_Float(pdfsModelFAN, "Exponential", 1, 1, Form("ftest_Zee_DCB_%s_cat%d",ext.c_str(),cat), alphacb1Comm, alphacb2Comm).first;  
                if(addExp)   bkgPdf_gof = getPdfLowDoubleCB_Float(pdfsModelFAN, "Exponential", 1, 1, Form("ftest_Zee_DCB_gof_%s_cat%d",ext.c_str(),cat), alphacb1Comm, alphacb2Comm).first;
           } 
   }

 


 
    int fitStatus = 0; double thisNll=0.;double gofProb=0;
    runFitFAN(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/3,mhLow,mhHigh);

    
    string outFileName;   
    if(!useDoubleCB){
        if(forLoop == -1)  outFileName = Form("%s/ftest_Zee_Voi_cat%d.png",outDir.c_str(),cat);
        else  outFileName = Form("%s/ftest_Zee_Voi_cat%d_loop%d.png",outDir.c_str(),cat,forLoop);
    }else{
        if(forLoop == -1)  outFileName = Form("%s/ftest_Zee_DCB_cat%d.png",outDir.c_str(),cat);
        else  outFileName = Form("%s/ftest_Zee_DCB_cat%d_loop%d.png",outDir.c_str(),cat,forLoop);
    }
    
    plotFAN(mass,bkgPdf,bkgPdf_gof,data,outFileName,fitStatus,&gofProb,bins,mhLow,mhHigh,runGOF,useDoubleCB,alphacb1Comm,alphacb2Comm,isData,true,runTDR);

   
    RooArgSet *paramsVoi = bkgPdf->getParameters(*data);
    if(!useDoubleCB){
          float voiMean = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_mean_p0",ext.c_str(),cat)))->getValV();
          float voiMeanErrorL = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_mean_p0",ext.c_str(),cat)))->getErrorLo();
          float voiMeanErrorH = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_mean_p0",ext.c_str(),cat)))->getErrorHi(); 
          float voiSigma = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_sigma_p0",ext.c_str(),cat)))->getValV();
          float voiSigmaErrorL = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_sigma_p0",ext.c_str(),cat)))->getErrorLo();
          float voiSigmaErrorH = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_sigma_p0",ext.c_str(),cat)))->getErrorHi();
          float voiWidth = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_width_p0",ext.c_str(),cat)))->getValV();
          float voiWidthErrorL = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_width_p0",ext.c_str(),cat)))->getErrorLo();
          float voiWidthErrorH = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_width_p0",ext.c_str(),cat)))->getErrorHi();


          if(addExp){
              float expParam = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_Off_p1",ext.c_str(),cat)))->getValV(); 
              float expParamErrorL = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_Off_p1",ext.c_str(),cat)))->getErrorLo();
              float expParamErrorH = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_Off_p1",ext.c_str(),cat)))->getErrorHi();
              float expVoiFrac = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_plus_Exponential_1_frac_sum1",ext.c_str(),cat)))->getValV();
              float expVoiFracErrorL = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_plus_Exponential_1_frac_sum1",ext.c_str(),cat)))->getErrorLo();
              float expVoiFracErrorH = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_Voi_%s_cat%d_plus_Exponential_1_frac_sum1",ext.c_str(),cat)))->getErrorHi();
 
              outfile << Form("cat: %d    ",cat) << data->sumEntries()*(1.0-expVoiFrac) << "    Mean: " << voiMean << "   voiMeanErrorL:  " << voiMeanErrorL << " voiMeanErrorH:  "<< voiMeanErrorH << "  voiSigma:  " << voiSigma << "  voiSigmaErrorL:  " << voiSigmaErrorL << " voiSigmaErrorH:  " << voiSigmaErrorH << "  voiWidth:  " << voiWidth << "  voiWidthErrorL:  " << voiWidthErrorL << " voiWidthErrorH:  " << voiWidthErrorH << "   expParam:  " << expParam << "   expParamErrorL:  " << expParamErrorL << "   expParamErrorH:  " << expParamErrorH << "   expVoiFrac:  " << expVoiFrac << "   expVoiFracErrorL:  " << expVoiFracErrorL << "   expVoiFracErrorH:  " << expVoiFracErrorH << endl;
          }
          else{
              outfile << Form("cat: %d    ",cat) << data->sumEntries() << "    Mean: " << voiMean << "   voiMeanErrorL:  " << voiMeanErrorL << " voiMeanErrorH:  "<< voiMeanErrorH << "  voiSigma:  " << voiSigma << "  voiSigmaErrorL:  " << voiSigmaErrorL << " voiSigmaErrorH:  " << voiSigmaErrorH << "  voiWidth:  " << voiWidth << "  voiWidthErrorL:  " << voiWidthErrorL << " voiWidthErrorH:  " << voiWidthErrorH << endl;

          }
    }else{
          float voiMean = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_mean_p0",ext.c_str(),cat)))->getValV();
          float voiMeanErrorL = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_mean_p0",ext.c_str(),cat)))->getErrorLo();
          float voiMeanErrorH = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_mean_p0",ext.c_str(),cat)))->getErrorHi(); 
          float voiSigma = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_sigma_p0",ext.c_str(),cat)))->getValV();
          float voiSigmaErrorL = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_sigma_p0",ext.c_str(),cat)))->getErrorLo();
          float voiSigmaErrorH = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_sigma_p0",ext.c_str(),cat)))->getErrorHi();
          float nCB1 = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_nCB1_p0",ext.c_str(),cat)))->getValV();
          float nCB1ErrorL = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_nCB1_p0",ext.c_str(),cat)))->getErrorLo();
          float nCB1ErrorH = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_nCB1_p0",ext.c_str(),cat)))->getErrorHi();
          float nCB2 = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_nCB2_p0",ext.c_str(),cat)))->getValV();
          float nCB2ErrorL = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_nCB2_p0",ext.c_str(),cat)))->getErrorLo();
          float nCB2ErrorH = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_nCB2_p0",ext.c_str(),cat)))->getErrorHi();
          float alphaCB1 = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_alphaCB1_p0",ext.c_str(),cat)))->getValV();
          float alphaCB2 = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_alphaCB2_p0",ext.c_str(),cat)))->getValV();




          if(addExp){
              float expParam = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_Off_p1",ext.c_str(),cat)))->getValV(); 
              float expParamErrorL = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_Off_p1",ext.c_str(),cat)))->getErrorLo();
              float expParamErrorH = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_Off_p1",ext.c_str(),cat)))->getErrorHi();
              float expCBFrac = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_plus_Exponential_1_frac_sum1",ext.c_str(),cat)))->getValV();
              float expCBFracErrorL = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_plus_Exponential_1_frac_sum1",ext.c_str(),cat)))->getErrorLo();
              float expCBFracErrorH = ((RooRealVar*)paramsVoi->find(Form("ftest_Zee_DCB_%s_cat%d_plus_Exponential_1_frac_sum1",ext.c_str(),cat)))->getErrorHi();
 
              outfile << Form("cat: %d    ",cat) << data->sumEntries()*(1.0-expCBFrac) << "    Mean: " << voiMean << "   voiMeanErrorL:  " << voiMeanErrorL << " voiMeanErrorH:  "<< voiMeanErrorH << "  voiSigma:  " << voiSigma << "  voiSigmaErrorL:  " << voiSigmaErrorL << " voiSigmaErrorH:  " << voiSigmaErrorH << "    nCB1:  " << nCB1 << "    nCB1ErrorL:   " << nCB1ErrorL << "   nCB1ErrorH:  " << nCB1ErrorH << "   nCB2:  " << nCB2 << "   nCB2ErrorL:  " << nCB2ErrorL << "    nCB2ErrorH:   " << nCB2ErrorH << "    alphaCB1:   " << alphaCB1 << "    alphaCB2:   " << alphaCB2 << "   expParam:  " << expParam << "   expParamErrorL:  " << expParamErrorL << "   expParamErrorH: " << expParamErrorH << "   expCBFrac:  " << expCBFrac << "   expCBFracErrorL:  " << expCBFracErrorL << "   expCBFracErrorH:  " << expCBFracErrorH << endl;
          }
          else{
              outfile << Form("cat: %d    ",cat) << data->sumEntries() << "    Mean:  " << voiMean << "   voiMeanErrorL:  " << voiMeanErrorL << " voiMeanErrorH:  "<< voiMeanErrorH << "  voiSigma:  " << voiSigma << "  voiSigmaErrorL:  " << voiSigmaErrorL << " voiSigmaErrorH:  " << voiSigmaErrorH << "    nCB1:  " << nCB1 << "    nCB1ErrorL:   " << nCB1ErrorL << "   nCB1ErrorH:  " << nCB1ErrorH << "   nCB2:  " << nCB2 << "   nCB2ErrorL:  " << nCB2ErrorL << "    nCB2ErrorH:   " << nCB2ErrorH << "    alphaCB1:   " << alphaCB1 << "    alphaCB2:   " << alphaCB2 << endl;

          }

    }



    outfile << endl;
    paramsVoi->Print("v");

 

    //outputws->Print();

 
    outputfile->cd();
    outputws->Write();


  }

    outputfile->Close();



}
