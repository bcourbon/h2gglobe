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
  return pdfsModel.getOfficialSumVoigtians(type, Form("%s_%s%d",ext,type1,orderOff), Form("%s_%s%d",ext,type1,orderOff), orderOff, orderVoig, constant, bernDownBound, bernUpBound);

}



pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getPdfLowVoiFix(PdfModelBuilderFAN &pdfsModel, string type, int orderOff, RooAbsPdf* pdfVoiFix, const char* ext="", float constant=0, float bernDownBound=-11., float bernUpBound=11.){

  const char* type1 = type.c_str();
  return pdfsModel.getOfficialSumVoigtiansFix(type, Form("%s_%s%d",ext,type1,orderOff), Form("%s_%s%d",ext,type1,orderOff), orderOff, pdfVoiFix, constant, bernDownBound, bernUpBound);

}



pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getPdfLowVoiFixNew(PdfModelBuilderFAN &pdfsModel, string type, int orderOff, float voiMean, float voiMeanErrorL, float voiMeanErrorH, float voiSigma, float voiSigmaErrorL, float voiSigmaErrorH, float voiWidth, float voiWidthErrorL, float voiWidthErrorH, float ErrorRange, const char* ext="", float constant=0, bool fixVoi=false, float bernDownBound=-11., float bernUpBound=11.){

  const char* type1 = type.c_str();
  return pdfsModel.getOfficialSumVoigtiansFixNew(type, Form("%s_%s%d",ext,type1,orderOff), Form("%s_%s%d",ext,type1,orderOff), orderOff, voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voiWidth, voiWidthErrorL, voiWidthErrorH, ErrorRange, constant, fixVoi, bernDownBound, bernUpBound);

}


pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getPdfLowVoiFixNewDouleCB(PdfModelBuilderFAN &pdfsModel, string type, int orderOff, float voiMean, float voiMeanErrorL, float voiMeanErrorH, float voiSigma, float voiSigmaErrorL, float voiSigmaErrorH, float voinCB1, float voinCB1ErrorL, float voinCB1ErrorH, float voinCB2, float voinCB2ErrorL, float voinCB2ErrorH, float voialphaCB1, float voialphaCB2, float ErrorRange, const char* ext="", float constant=0, bool fixVoi=false, float bernDownBound=-11., float bernUpBound=11.){

  const char* type1 = type.c_str();
  return pdfsModel.getOfficialSumVoigtiansFixNewDoubleCB(type, Form("%s_%s%d",ext,type1,orderOff), Form("%s_%s%d",ext,type1,orderOff), orderOff, voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voinCB1, voinCB1ErrorL, voinCB1ErrorH, voinCB2, voinCB2ErrorL, voinCB2ErrorH, voialphaCB1, voialphaCB2, ErrorRange, constant, fixVoi, bernDownBound, bernUpBound);

}





void runFitFAN(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries, int mhLow, int mhHigh, bool fitMethod){

	int ntries=0;
  	RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
	std::cout << " BEFORE ITERATIONS-------------------------------" << std::endl;
	params_test->Print("v");
	int stat=1;
	double minnll=10e8;
	while (stat!=0){
	  if (ntries>=MaxTries) break;
          RooFitResult *fitTest; 
          fitTest = pdf->fitTo(*data,RooFit::Save(1), Range(mhLow,mhHigh));  
 
          cout << "****************************************************before" << endl; 
          fitTest->floatParsInit().Print("v");
          cout << "****************************************************after" << endl;
 	  fitTest->floatParsFinal().Print("v");  
          stat = fitTest->status();
	  minnll = fitTest->minNll();
	  if (stat!=0) params_test->assignValueOnly(fitTest->randomizePars());
          cout << "========================= fit ========= " << stat << "   =====   " << ntries << endl;
	  ntries++; 
	}
	*stat_t = stat;
	*NLL = minnll;
}



void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries){

	int ntries=0;
  	RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
	//std::cout << " BEFORE ITERATIONS-------------------------------" << std::endl;
	//params_test->Print("v");
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
  
  // fit the pdfs to the data and keep this fit Result (for randomizing)
  RooFitResult *fitNullData = pdfNull->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1));
  RooFitResult *fitTestData = pdfTest->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::PrintLevel(-1));

  // Ok we want to check the distribution in toys then 
  // Step 1, cache the parameters of each pdf so as not to upset anything 
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



double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name){

  double prob;
  int ntoys = 500;

  // Routine to calculate the goodness of fit. 
  //name+="_gofTest.pdf";
  name+="_gofTest.png";  
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);
  //norm.removeRange();

  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);

  // get The Chi2 value from the data
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(nBinsForMass),Name("data"));

  pdf->plotOn(plot_chi2,Name("pdf"));
  int np = pdf->getParameters(*data)->getSize();

  double chi2 = plot_chi2->chiSquare("pdf","data",np);
  std::cout << "Calculating GOF for pdf " << pdf->GetName() << ", using " <<np << " fitted parameters" <<std::endl;

  // The first thing is to check if the number of entries in any bin is < 5 
  // if so, we don't rely on asymptotic approximations
 
  if ((double)data->sumEntries()/nBinsForMass < 5 ){

    std::cout << "Running toys for GOF test " << std::endl;
    // store pre-fit params 
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
      binnedtoy->plotOn(plot_t);
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

    // back to best fit 	
    params->assignValueOnly(preParams);
  } else {
    prob = TMath::Prob(chi2*(nBinsForMass-np),nBinsForMass-np);
  }
  std::cout << "Chi2 in Observed =  " << chi2*(nBinsForMass-np) << std::endl;
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

 // data->plotOn(plot,Binning(80));
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
  //canv->SaveAs(Form("%s.pdf",name.c_str()));  
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
  //canv->SaveAs(Form("%s.pdf",name.c_str()));  
  canv->SaveAs(Form("%s.png",name.c_str()));
  delete canv;
}



void plotFAN(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name, int status, double *prob, int bins=0, int mhHigh=0){
  
  // Chi2 taken from full range fit
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(nBinsForMass));
  pdf->plotOn(plot_chi2);

  int np = pdf->getParameters(*data)->getSize()+1; //Because this pdf has no extend
  double chi2 = plot_chi2->chiSquare(np);
 
  *prob = getGoodnessOfFit(mass,pdf,data,name);
  RooPlot *plot = mass->frame();
  mass->setRange("unblindReg_1",85,95);
  mass->setRange("unblindReg_2",mhHigh,mhHigh);
  if (BLIND) {
    data->plotOn(plot,Binning(bins),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(bins),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(bins),Invisible());
  }
  else data->plotOn(plot,Binning(bins));

 // data->plotOn(plot,Binning(bins));
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
  //canv->SaveAs(Form("%s.pdf",name.c_str()));  
  canv->SaveAs(Form("%s.png",name.c_str()));
  catIndex->setIndex(currentIndex);
  delete canv;
}

void plotFAN(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooDataSet *data, string name, int cat, int bestFitPdf=-1, int bins=0, int mhLow=0, int mhHigh=0, float LaurentConstant=0, bool runTDR=false){
  
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.6,0.65,0.89,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",86,96);     
  mass->setRange("unblindReg_2",110,mhHigh);  
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
    //if(it->first.find("Laurent")!=string::npos)  leg->AddEntry(pdfLeg,Form("%s%s(%.f)",it->first.c_str(),ext.c_str(),LaurentConstant),"L");  no need
    //else  leg->AddEntry(pdfLeg,Form("%s%s",it->first.c_str(),ext.c_str()),"L");  no need
    i++;
  }
  plot->SetTitle(Form("Category %d",cat));
  plot->SetXTitle("m_{#gamma#gamma}(GeV)"); 
  plot->SetYTitle(Form("Events / %.2fGeV",float(mhHigh-mhLow)/float(bins))); 
  if(runTDR){
     plot->SetTitleSize(0.045, "XY");
     plot->SetTitleOffset(1.8,"Y");
  }
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  //canv->SaveAs(Form("%s.pdf",name.c_str()));  
  canv->SaveAs(Form("%s.png",name.c_str()));
  if(runTDR)   canv->SaveAs(Form("%s.pdf",name.c_str())); 
  delete canv;
}



void transferMacros(TFile *inFile, TFile *outFile){
  
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    if (string(key->ReadObj()->ClassName())=="TMacro") {
      //cout << key->ReadObj()->ClassName() << " : " << key->GetName() << endl;
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
	
	//bkg->setDirtyInhibit(1);
	//RooAbsReal *nllm = bkg->createNLL(*data);
	//RooMinimizer minim(*nllm);
	//minim.setStrategy(1);
	
	for (int id=0;id<number_of_indeces;id++){		
		params->assignValueOnly(clean);
		cat->setIndex(id);

		//RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);

		if (!silent) {
			/*
			std::cout << "BEFORE  MAKING FIT" << std::endl;
			params->Print("V");
			std::cout << "-----------------------" << std::endl;		
			*/
		}
		
		//minim.minimize("Minuit2","minimize");
		double minNll=0; //(nllm->getVal())+bkg->getCorrection();
		int fitStatus=1;		
		runFit(bkg->getCurrentPdf(),data,&minNll,&fitStatus,/*max iterations*/3);
		// Add the penalty

		minNll=minNll+bkg->getCorrection();

		if (!silent) {
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

  bool runTDR=true;
  //bool runTDR=false;

  if(runTDR){
     gStyle->SetOptFit(111);
     gStyle->SetOptStat(0);
     gStyle->SetOptTitle(0);
  }
  
 
  string fileName;
  string fileNameZee;  
  int ncats;
  int singleCategory;
  string datfile;
  string outDir;
  string outfilename;
  string runType; 
  int bins; 
  int mhLow; 
  int mhHigh; 
  bool FixOrSo=true;  
  float LaurentConstant;  
  float bernDownBound;  
  float bernUpBound;  
  bool fitMethod=false;   
  bool is2011=false;
  bool addExp=false;  
  bool useDoubleCB=false;  
  bool verbose=false;
  bool saveMultiPdf=false;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName),                                              "In file name")
    ("infilenameZee,I", po::value<string>(&fileNameZee),                                              "In file name Zee")   
    ("ncats,c", po::value<int>(&ncats)->default_value(5),                                       "Number of categories")
    ("singleCat", po::value<int>(&singleCategory)->default_value(-1),                           "Run A single Category")
    ("datfile,d", po::value<string>(&datfile)->default_value("dat/fTest.dat"),                  "Right results to datfile for BiasStudy")
    ("outDir,D", po::value<string>(&outDir)->default_value("plots/fTest"),                      "Out directory for plots")
    ("saveMultiPdf", po::value<string>(&outfilename),         					"Save a MultiPdf model with the appropriate pdfs")
    ("runType,R", po::value<string>(&runType)->default_value(""),                               "Official or low mass") 
    ("mhLow,L", po::value<int>(&mhLow)->default_value(75),                                                                                                                 "Starting point for scan") 
    ("mhHigh,H", po::value<int>(&mhHigh)->default_value(120),                                                                                                               "End point for scan") 
    ("bins,B", po::value<int>(&bins)->default_value(45),                                                                                                                 "Bins for the plot") 
    ("nBinsForMass,b", po::value<int>(&nBinsForMass)->default_value(180),                                                                                               "nBinsForMass") 
    ("LaurentConstant,C", po::value<float>(&LaurentConstant)->default_value(-4.),                                                                                                                 "constant for Laurent") 
    ("bernDownBound,X", po::value<float>(&bernDownBound)->default_value(-11.),                                                                                                                 "lower bound for bern") 
    ("bernUpBound,Y", po::value<float>(&bernUpBound)->default_value(11.),                                                                                                                 "up bound for bern") 
    ("runFtestCheckWithToys", 									"When running the F-test, use toys to calculate pvals (and make plots) ")
    ("is2011",                                                                                  "Run 2011 config")
    ("addExp",                                                                                  "add Exponential")   
    ("useDoubleCB",                                                                                  "use double crystal ball function")   
    ("FixVoi",                                                                                  "Fix voi parameter")   
    ("skipBern",                                                                                  "Bern does not have logN")   
    ("unblind",  									        "Dont blind plots")
    ("verbose,v",                                                                               "Run with more output")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("is2011")) is2011=true;
  if (vm.count("addExp")) addExp=true;   
  if (vm.count("useDoubleCB"))  useDoubleCB=true;   
  if (vm.count("FixVoi")) FixOrSo=true;   	
  if (vm.count("unblind")) BLIND=false;
  saveMultiPdf = vm.count("saveMultiPdf");

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

  std::cout << "SaveMultiPdf? " << saveMultiPdf << std::endl;
  TFile *outputfile;
  RooWorkspace *outputws;

  if (saveMultiPdf){
	outputfile = new TFile(outfilename.c_str(),"RECREATE");
	outputws = new RooWorkspace(); outputws->SetName("multipdf");
  }

  system(Form("mkdir -p %s",outDir.c_str()));
  TFile *inFile = TFile::Open(fileName.c_str());
  RooWorkspace *inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");

  if (saveMultiPdf){
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

  // store results here
  
  FILE *resFile ;
  if  (singleCategory >-1) resFile = fopen(Form("%s/fTestResults_cat%d.txt",outDir.c_str(),singleCategory),"w");
  else resFile = fopen(Form("%s/fTestResults.txt",outDir.c_str()),"w");
  vector<map<string,int> > choices_vec;
  vector<map<string,std::vector<int> > > choices_envelope_vec;
  vector<map<string,RooAbsPdf*> > pdfs_vec;

  PdfModelBuilder pdfsModel;
  PdfModelBuilderFAN pdfsModelFAN; 
  RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hgg_mass");
  mass->setRange(mhLow,mhHigh); 
  pdfsModel.setObsVar(mass);
  pdfsModelFAN.setObsVar(mass); 
  double upperEnvThreshold = 0.1; // upper threshold on delta(chi2) to include function in envelope (looser than truth function)
  

  std::string ext = is2011 ? "7TeV" : "8TeV";
  for (int cat=startingCategory; cat<ncats; cat++){
    
    map<string,int> choices;
    map<string,std::vector<int> > choices_envelope;
    map<string,RooAbsPdf*> pdfs;
    map<string,RooAbsPdf*> allPdfs;
    
    RooDataSet *dataFull = (RooDataSet*)inWS->data(Form("data_mass_cat%d",cat));
    

    mass->setBins(nBinsForMass);
    RooDataHist thisdataBinned(Form("roohist_data_mass_cat%d",cat),"data",*mass,*dataFull);
    RooDataSet *data = (RooDataSet*)&thisdataBinned;

    RooArgList storedPdfs("store");

    fprintf(resFile,"\\centerline{\\begin{tabular}[c]{|c|c|c|c|}\n");  
    fprintf(resFile,"\\hline\n");  
    fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Category %d}} \\\\\n",cat);
    fprintf(resFile,"\\hline\n");
    fprintf(resFile,"Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\\\n");    
    fprintf(resFile,"\\hline\n");    

    double MinimimNLLSoFar=1e10;
    int simplebestFitPdfIndex = 0;

    
    TFile *inFileZee; RooWorkspace *inWS_Zee; RooAbsPdf* pdfVoiFix;float voiMean=0.; float voiMeanErrorL=0.; float voiMeanErrorH=0.; float voiSigma=0.; float voiSigmaErrorL=0.; float voiSigmaErrorH=0.; float voiWidth=0.;  float voiWidthErrorL=0.; float voiWidthErrorH=0.; float voinCB1=0.; float voinCB1ErrorL=0.; float voinCB1ErrorH=0.; float voinCB2=0.; float voinCB2ErrorL=0.; float voinCB2ErrorH=0.; float voialphaCB1=0.; float voialphaCB2=0.; float ErrorRange=1.0;

    if(runType == "Voi"){
        inFileZee = TFile::Open(fileNameZee.c_str());
        inWS_Zee = (RooWorkspace*)inFileZee->Get("fTestVoi_Zee");
        if(!useDoubleCB) pdfVoiFix = inWS_Zee->pdf(Form("ftest_Zee_Voi_%s_cat%d",ext.c_str(),cat));  
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
    

        cout << "I am here  ...===================   " << voiMean << endl;        

    // Standard F-Test to find the truth functions
    int resBern_order=0; int resExp_order=0; int resPow_order=0; int resLau_order=0; 
    for (vector<string>::iterator funcType=functionClasses.begin(); 
	 funcType!=functionClasses.end(); funcType++){
      
      double thisNll=0.; double prevNll=0.; double chi2=0.; double prob=0.; 
      int order=1; int prev_order=0; int cache_order=0;

      RooAbsPdf *prev_pdf=NULL;
      RooAbsPdf *cache_pdf=NULL;
      std::vector<int> pdforders;

      while (prob<0.05){

        
            string TypeFAN = *funcType;
            if((TypeFAN == "Exponential" || TypeFAN == "PowerLaw") && order%2 == 0 ) order++;
        

        
        RooAbsPdf *bkgPdf = 0;
        if(runType == "Voi"){
            if(!useDoubleCB){
                      bkgPdf = getPdfLowVoiFixNew(pdfsModelFAN, *funcType, order, voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voiWidth, voiWidthErrorL, voiWidthErrorH, ErrorRange, Form("ftest_pdf_%d_%s",cat,ext.c_str()), LaurentConstant, FixOrSo, bernDownBound, bernUpBound).first; //fix voi new
             }else{
                      bkgPdf = getPdfLowVoiFixNewDouleCB(pdfsModelFAN, *funcType, order, voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voinCB1, voinCB1ErrorL, voinCB1ErrorH, voinCB2, voinCB2ErrorL, voinCB2ErrorH, voialphaCB1, voialphaCB2, ErrorRange, Form("ftest_pdf_%d_%s",cat,ext.c_str()), LaurentConstant, FixOrSo, bernDownBound, bernUpBound).first; 
                      if(addExp)  bkgPdf = getPdfLowVoiFixNewDouleCB(pdfsModelFAN, *funcType, order, voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voinCB1, voinCB1ErrorL, voinCB1ErrorH, voinCB2, voinCB2ErrorL, voinCB2ErrorH, voialphaCB1, voialphaCB2, ErrorRange, Form("ftest_pdf_%d_%s",cat,ext.c_str()), LaurentConstant, FixOrSo,bernDownBound, bernUpBound).first; 
            }
        }
        else{
            bkgPdf = getPdf(pdfsModel,*funcType,order,Form("ftest_pdf_%d_%s",cat,ext.c_str()));
        }
        
        
        
        if (!bkgPdf){
          // assume this order is not allowed
          order++;
        }
        else {
	
	  int fitStatus = 0;
          
          
          if(runType == "Voi"){
               runFitFAN(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/1,mhLow,mhHigh, fitMethod);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));  
          }
          else{
             runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/3);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));  
          }


 	  if (fitStatus!=0) std::cout << "Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
          chi2 = 2.*(prevNll-thisNll);
          if (chi2<0. && order>1) chi2=0.;
	  if (prev_pdf!=NULL){
	    prob = getProbabilityFtest(chi2,order-prev_order,prev_pdf,bkgPdf,mass,data
		   ,Form("%s/Ftest_from_%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,cat));
	    std::cout << " F-test Prob(chi2>chi2(data) == " << prob << std::endl;
	  } else {
	    prob = 0;
	  }
	  double gofProb=0;
          
          if (!saveMultiPdf){
             if(runType == "Voi"){
                plotFAN(mass,bkgPdf,data,Form("%s/%s%d_cat%d.png",outDir.c_str(),funcType->c_str(),order,cat),fitStatus,&gofProb,bins,mhHigh); 
             }
             else{
                plotFAN(mass,bkgPdf,data,Form("%s/%s%d_cat%d.png",outDir.c_str(),funcType->c_str(),order,cat),fitStatus,&gofProb,bins,mhHigh); 
             }
          }
          

          cout << "\t " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
          prevNll=thisNll;
          cache_order=prev_order;
          cache_pdf=prev_pdf;
          prev_order=order;
          prev_pdf=bkgPdf;
          order++;
        }
      }

      if(!useDoubleCB)   fprintf(resFile,"%d%s+Voi & %d + 1 & %5.2f & %5.2f \\\\\n",cache_order,funcType->c_str(),cache_order+1,chi2,prob);   
      else  fprintf(resFile,"%d%s+DCB & %d + 1 & %5.2f & %5.2f \\\\\n",cache_order,funcType->c_str(),cache_order+1,chi2,prob);   
            
      string TypeRes = *funcType;
      if(TypeRes == "Bernstein")  resBern_order = cache_order; 
      if(TypeRes == "Exponential")   resExp_order = cache_order; 
      if(TypeRes == "PowerLaw")  resPow_order = cache_order; 
      if(TypeRes == "Laurent")   resLau_order = cache_order; 
      
      choices.insert(pair<string,int>(*funcType,cache_order));
      pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));

      int truthOrder = cache_order;

      if (saveMultiPdf){
        chi2=0.;
        thisNll=0.;
        prevNll=0.;
        prob=0.;
        order=1;
        prev_order=0;
        cache_order=0;
	std::cout << "Determining Envelope Functions for Family " << *funcType << ", cat " << cat << std::endl;
	std::cout << "Upper end Threshold for highest order function " << upperEnvThreshold <<std::endl;

        while (prob<upperEnvThreshold){
          
              string TypeFAN = *funcType;
              if((TypeFAN == "Exponential" || TypeFAN == "PowerLaw") && order%2 == 0 ) order++;

          RooAbsPdf *bkgPdf = 0;
          if(runType == "Voi"){
                    bkgPdf = getPdfLowVoiFixNew(pdfsModelFAN, *funcType, order, voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voiWidth, voiWidthErrorL, voiWidthErrorH, ErrorRange, Form("ftest_Epdf_%d_%s",cat,ext.c_str()), LaurentConstant, FixOrSo, bernDownBound, bernUpBound).first;   //fix voi new
          } 
          else{
               bkgPdf = getPdf(pdfsModel,*funcType,order,Form("ftest_Epdf_%d_%s",cat,ext.c_str()));
          }
                    

          if (!bkgPdf ){
          // assume this order is not allowed
          order++;
          }

          else {
	   int fitStatus=0;
 
           
	   if(runType == "Voi"){
               runFitFAN(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/1,mhLow,mhHigh, fitMethod);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));  
           }
           else{
               runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/3);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));  
           }
           

           //thisNll = fitRes->minNll();
 	   if (fitStatus!=0) std::cout << "Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
	   double myNll = 2.*thisNll;
           chi2 = 2.*(prevNll-thisNll);
           if (chi2<0. && order>1) chi2=0.;
           
           if (prev_pdf!=NULL){
              prob = getProbabilityFtest(chi2,order-prev_order,prev_pdf,bkgPdf,mass,data 
                                 ,Form("%s/Ftest_from_%s%d_cat%d.png",outDir.c_str(),funcType->c_str(),order,cat)); 
              std::cout << " F-test Envelope Prob(chi2>chi2(data) == " << prob << std::endl;
           } else {
              prob = 0;
           }
             

           cout << "\t " << *funcType << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
           prevNll=thisNll;
           cache_order=prev_order;
           cache_pdf=prev_pdf;

	   double gofProb =0; 
           
           if(runType == "Voi"){
                plotFAN(mass,bkgPdf,data,Form("%s/%s%d_cat%d.png",outDir.c_str(),funcType->c_str(),order,cat),fitStatus,&gofProb,bins,mhHigh); 
           }
           else{
                plot(mass,bkgPdf,data,Form("%s/%s%d_cat%d.png",outDir.c_str(),funcType->c_str(),order,cat),fitStatus,&gofProb); 
           }
           
           

	   if ((prob < upperEnvThreshold) ) { // Looser requirements for the envelope
		
	     if (gofProb > 0.01 || order == truthOrder ) {  // Good looking fit or one of our regular truth functions

		std::cout << "Adding to Envelope " << bkgPdf->GetName() << " "<< gofProb 
			  << " 2xNLL + c is " << myNll + bkgPdf->getVariables()->getSize() <<  std::endl;
      		allPdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),order),bkgPdf));
		storedPdfs.add(*bkgPdf);
		pdforders.push_back(order);
 
		// Keep track but we shall redo this later
		if ((myNll + bkgPdf->getVariables()->getSize()) < MinimimNLLSoFar) {
		 simplebestFitPdfIndex = storedPdfs.getSize()-1;
		 MinimimNLLSoFar = myNll + bkgPdf->getVariables()->getSize();
		}
	     }
	   }

           prev_order=order;
           prev_pdf=bkgPdf;
           order++;
        }
      }
      
      //fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
      fprintf(resFile,"%d%s+Voi & %d + 3 & %5.2f & %5.2f \\\\\n",cache_order,funcType->c_str(),cache_order+1,chi2,prob);   
      choices_envelope.insert(pair<string,std::vector<int> >(*funcType,pdforders));
      }
    }

    fprintf(resFile,"\\hline\n");
    if(!useDoubleCB)   fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Fit Model: FinalFitFunc+Voi }} \\\\\n");   
    else  fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Fit Model: FinalFitFunc+DCB }} \\\\\n");   
    if(!useDoubleCB)   fprintf(resFile,"%dBernstein+Voi & %dExponential+Voi & %dPowerLaw+Voi & %dLaurent+Voi \\\\\n", resBern_order,resExp_order,resPow_order,resLau_order);   
    else   fprintf(resFile,"%dBernstein+DCB & %dExponential+DCB & %dPowerLaw+DCB & %dLaurent+DCB \\\\\n", resBern_order,resExp_order,resPow_order,resLau_order);   
    choices_vec.push_back(choices);
    choices_envelope_vec.push_back(choices_envelope);
    pdfs_vec.push_back(pdfs);

    
    if(runType == "Voi"){
        plotFAN(mass,pdfs,data,Form("%s/truths_cat%d",outDir.c_str(),cat),cat,-1,bins,mhLow,mhHigh,LaurentConstant,runTDR);
    }
    else{
        plotFAN(mass,pdfs,data,Form("%s/truths_cat%d",outDir.c_str(),cat),cat,-1,bins,mhLow,mhHigh,LaurentConstant,runTDR);  
    }
    
 
    if (saveMultiPdf){


	  RooCategory catIndex(Form("pdfindex_%d_%s",cat,ext.c_str()),"c");
	  RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_hgg_cat%d_%s_bkgshape",cat,ext.c_str()),"all pdfs",catIndex,storedPdfs);
	  RooRealVar nBackground(Form("CMS_hgg_cat%d_%s_bkgshape_norm",cat,ext.c_str()),"nbkg",data->sumEntries(),0,10E8);
	  int bestFitPdfIndex;
              bestFitPdfIndex = getBestFitFunction(pdf,data,&catIndex,!verbose);
          
          
	  catIndex.setIndex(bestFitPdfIndex);
	  std::cout << "// ------------------------------------------------------------------------- //" <<std::endl; 
	  std::cout << "Created MultiPdf " << pdf->GetName() << ", in Category " << cat << " with a total of " << catIndex.numTypes() << " pdfs"<< std::endl;
	  storedPdfs.Print();
	  std::cout << "Best Fit Pdf = " << bestFitPdfIndex << ", " << storedPdfs.at(bestFitPdfIndex)->GetName() << std::endl;
	  std::cout << "// ------------------------------------------------------------------------- //" <<std::endl;
	  std::cout << " Simple check of index "<< simplebestFitPdfIndex <<std::endl;

	  mass->setBins(nBinsForMass);
	  RooDataHist dataBinned(Form("roohist_data_mass_cat%d",cat),"data",*mass,*dataFull);

	  outputws->import(*pdf);
	  outputws->import(nBackground);
	  outputws->import(catIndex);
	  outputws->import(dataBinned);
	  outputws->import(*data);

          plot(mass,pdf,&catIndex,data,Form("%s/multipdf_cat%d",outDir.c_str(),cat),cat,bestFitPdfIndex);

    }
    
  }
  if (saveMultiPdf){
	outputfile->cd();
	outputws->Write();
	outputfile->Close();	
  }

  FILE *dfile = fopen(datfile.c_str(),"w");
  cout << "Recommended options" << endl;
  
  for (int cat=startingCategory; cat<ncats; cat++){
    cout << "Cat " << cat << endl;
    fprintf(dfile,"cat=%d\n",cat); 
    for (map<string,int>::iterator it=choices_vec[cat-startingCategory].begin(); it!=choices_vec[cat-startingCategory].end(); it++){
      cout << "\t" << it->first << " - " << it->second << endl;
      fprintf(dfile,"truth=%s:%d:%s%d\n",it->first.c_str(),it->second,namingMap[it->first].c_str(),it->second);
    }
    fprintf(dfile,"fabian=%s:%d:%s%d\n",fabChoice[cat-startingCategory].first.c_str()
	,fabChoice[cat-startingCategory].second,namingMap[fabChoice[cat-startingCategory].first].c_str(),fabChoice[cat-startingCategory].second);
    for (map<string,std::vector<int> >::iterator it=choices_envelope_vec[cat-startingCategory].begin(); it!=choices_envelope_vec[cat-startingCategory].end(); it++){
	std::vector<int> ords = it->second;
        for (std::vector<int>::iterator ordit=ords.begin(); ordit!=ords.end(); ordit++){
          fprintf(dfile,"paul=%s:%d:%s%d\n",it->first.c_str(),*ordit,namingMap[it->first].c_str(),*ordit);
        }
    }
    fprintf(dfile,"\n");
  }
  inFile->Close();

  return 0;
}
