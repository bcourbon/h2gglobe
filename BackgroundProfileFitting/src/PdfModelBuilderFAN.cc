#include "TCanvas.h"

#include "RooPlot.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "RooPowerLaw.h"
#include "RooPowerLawSum.h"
#include "RooKeysPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooGaussian.h" 
#include "RooVoigtian.h" 
#include "RooBreitWigner.h" 
#include "RooCBShape.h" 
#include "RooFFTConvPdf.h" 


#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "../interface/PdfModelBuilderFAN.h"

#include "HiggsAnalysis/CombinedLimit/interface/HGGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooDoubleCBFast.h"   



using namespace std;
using namespace RooFit;
using namespace boost;

PdfModelBuilderFAN::PdfModelBuilderFAN():
  obs_var_set(false),
  signal_modifier_set(false),
  signal_set(false),
  bkgHasFit(false),
  sbHasFit(false),
  keysPdfAttributesSet(false),
  verbosity(0)
{
  
  recognisedPdfTypes.push_back("Bernstein");
  recognisedPdfTypes.push_back("Exponential");
  recognisedPdfTypes.push_back("PowerLaw");
  recognisedPdfTypes.push_back("Laurent");
  recognisedPdfTypes.push_back("Chebychev");  
  recognisedPdfTypes.push_back("Polynomial");  
  recognisedPdfTypes.push_back("KeysPdf");
  recognisedPdfTypes.push_back("File");

  wsCache = new RooWorkspace("PdfModelBuilderFANCache");

};

PdfModelBuilderFAN::~PdfModelBuilderFAN(){};

void PdfModelBuilderFAN::setObsVar(RooRealVar *var){
  obs_var=var;
  obs_var_set=true;
}

void PdfModelBuilderFAN::setSignalModifier(RooRealVar *var){
  signalModifier=var;
  signal_modifier_set=true;
}

void PdfModelBuilderFAN::setSignalModifierVal(float val){
  signalModifier->setVal(val);
}

void PdfModelBuilderFAN::setSignalModifierConstant(bool val){
  signalModifier->setConstant(val);
}

RooAbsPdf* PdfModelBuilderFAN::getChebychev(string prefix, int order){
  
  RooArgList *coeffList = new RooArgList();
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
    //RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0*@0",RooArgList(*param));
    params.insert(pair<string,RooRealVar*>(name,param));
    //prods.insert(pair<string,RooFormulaVar*>(name,form));
    coeffList->add(*params[name]);
  }
  RooChebychev *cheb = new RooChebychev(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);   
  //RooPolynomial *cheb = new RooPolynomial(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  return cheb;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(bern->GetName(),bern));

}



RooAbsPdf* PdfModelBuilderFAN::getChebychevFAN(string prefix, int order, float bernDownBound, float bernUpBound){
  
  RooArgList *coeffList = new RooArgList();
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01, bernDownBound, bernUpBound); 
    //RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0*@0",RooArgList(*param));
    params.insert(pair<string,RooRealVar*>(name,param));
    //prods.insert(pair<string,RooFormulaVar*>(name,form));
    coeffList->add(*params[name]);
  }
  RooChebychev *cheb = new RooChebychev(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);   
  //RooPolynomial *cheb = new RooPolynomial(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  return cheb;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(bern->GetName(),bern));

}


RooAbsPdf* PdfModelBuilderFAN::getPolynomial(string prefix, int order){
  
  RooArgList *coeffList = new RooArgList();
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
    //RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0*@0",RooArgList(*param));
    params.insert(pair<string,RooRealVar*>(name,param));
    //prods.insert(pair<string,RooFormulaVar*>(name,form));
    coeffList->add(*params[name]);
  }
  //RooChebychev *cheb = new RooChebychev(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);   
  RooPolynomial *cheb = new RooPolynomial(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  return cheb;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(bern->GetName(),bern));

}



RooAbsPdf* PdfModelBuilderFAN::getBernstein(string prefix, int order){
  
  RooArgList *coeffList = new RooArgList();
  //coeffList->add(RooConst(1.0)); // no need for cnstant in this interface
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.1*(i+1),-11.,11.); 
    RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0*@0",RooArgList(*param));
    params.insert(pair<string,RooRealVar*>(name,param));
    prods.insert(pair<string,RooFormulaVar*>(name,form));
    coeffList->add(*prods[name]);
  }
  //RooBernstein *bern = new RooBernstein(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  if (order==1) {
	RooBernsteinFast<1> *bern = new RooBernsteinFast<1>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==2) {
	RooBernsteinFast<2> *bern = new RooBernsteinFast<2>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==3) {
	RooBernsteinFast<3> *bern = new RooBernsteinFast<3>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==4) {
	RooBernsteinFast<4> *bern = new RooBernsteinFast<4>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==5) {
	RooBernsteinFast<5> *bern = new RooBernsteinFast<5>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==6) {
	RooBernsteinFast<6> *bern = new RooBernsteinFast<6>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==7) {
	RooBernsteinFast<7> *bern = new RooBernsteinFast<7>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else {
	return NULL;
  }
  //return bern;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(bern->GetName(),bern));

}



RooAbsPdf* PdfModelBuilderFAN::getBernsteinFAN(string prefix, int order, float bernDownBound, float bernUpBound){
  
  RooArgList *coeffList = new RooArgList();
  //coeffList->add(RooConst(1.0)); // no need for cnstant in this interface
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.1*(i+1),bernDownBound,bernUpBound); 
    RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0*@0",RooArgList(*param));
    params.insert(pair<string,RooRealVar*>(name,param));
    prods.insert(pair<string,RooFormulaVar*>(name,form));
    coeffList->add(*prods[name]);
  }
  //RooBernstein *bern = new RooBernstein(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  if (order==1) {
	RooBernsteinFast<1> *bern = new RooBernsteinFast<1>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==2) {
	RooBernsteinFast<2> *bern = new RooBernsteinFast<2>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==3) {
	RooBernsteinFast<3> *bern = new RooBernsteinFast<3>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==4) {
	RooBernsteinFast<4> *bern = new RooBernsteinFast<4>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==5) {
	RooBernsteinFast<5> *bern = new RooBernsteinFast<5>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==6) {
	RooBernsteinFast<6> *bern = new RooBernsteinFast<6>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else if (order==7) {
	RooBernsteinFast<7> *bern = new RooBernsteinFast<7>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  	return bern;
  } else {
	return NULL;
  }
  //return bern;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(bern->GetName(),bern));

}



RooAbsPdf* PdfModelBuilderFAN::getPowerLawGeneric(string prefix, int order){
  
  if (order%2==0){
    cerr << "ERROR -- addPowerLaw -- only odd number of params allowed" << endl;
    return NULL;
  }
  else {
    int nfracs=(order-1)/2;
    int npows=order-nfracs;
    assert(nfracs==npows-1);
    string formula="";
    RooArgList *dependents = new RooArgList();
    dependents->add(*obs_var);
    // first do recursive fraction
    if (order>1) {
      formula += "(1.-";
      for (int i=1; i<=nfracs; i++){
        if (i<nfracs) formula += Form("@%d-",i);
        else formula += Form("@%d)*",i);
        string name =  Form("%s_f%d",prefix.c_str(),i);
        params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.1,0.,1.)));
        dependents->add(*params[name]);
      }
    }
    for (int i=1; i<=npows; i++){
      string pname =  Form("%s_p%d",prefix.c_str(),i);
      string fname =  Form("%s_f%d",prefix.c_str(),i-1);
      params.insert(pair<string,RooRealVar*>(pname, new RooRealVar(pname.c_str(),pname.c_str(),TMath::Max(-10.,-2.*(i+1)),-10.,0.)));
      if (i==1) {
        formula += Form("TMath::Power(@0,@%d)",nfracs+i);
        dependents->add(*params[pname]);
      }
      else {
        formula += Form(" + @%d*TMath::Power(@0,@%d)",i-1,nfracs+i);
        dependents->add(*params[pname]);
      }
    }
    cout << "FORM -- " << formula << endl;
    dependents->Print("v");
    RooGenericPdf *pow = new RooGenericPdf(prefix.c_str(),prefix.c_str(),formula.c_str(),*dependents);
    pow->Print("v");
    return pow;
    //bkgPdfs.insert(pair<string,RooAbsPdf*>(pow->GetName(),pow));

  }
}

RooAbsPdf* PdfModelBuilderFAN::getPowerLaw(string prefix, int order){
  
  RooArgList coefList;
  for (int i=0; i<order; i++){
    double start=-2.;
    double low=-10.;
    double high=0.;
    if (order>0){
      start=-0.001/double(i);
      low=-0.01;
      high=0.01;
    }
    RooRealVar *var = new RooRealVar(Form("%s_p%d",prefix.c_str(),i),Form("%s_p%d",prefix.c_str(),i),start,low,high);
    coefList.add(*var);
  }
  RooPowerLawSum *pow = new RooPowerLawSum(prefix.c_str(),prefix.c_str(),*obs_var,coefList);
  return pow;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(pow->GetName(),pow));

}

RooAbsPdf* PdfModelBuilderFAN::getExponential(string prefix, int order){
  
  RooArgList coefList;
  for (int i=0; i<order; i++){
    double start=-1.;
    double low=-2.;
    double high=0.;
    if (order>0){
      start=-0.001/double(i);
      low=-0.01;
      high=0.01;
    }
    RooRealVar *var = new RooRealVar(Form("%s_p%d",prefix.c_str(),i),Form("%s_p%d",prefix.c_str(),i),start,low,high);
    coefList.add(*var);
  }
  RooPowerLawSum *exp = new RooPowerLawSum(prefix.c_str(),prefix.c_str(),*obs_var,coefList);
  return exp;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(exp->GetName(),exp));

}

RooAbsPdf* PdfModelBuilderFAN::getPowerLawSingle(string prefix, int order){
  
  if (order%2==0){
    cerr << "ERROR -- addPowerLaw -- only odd number of params allowed" << endl;
    return NULL;
  }
  else {
    int nfracs=(order-1)/2;
    int npows=order-nfracs;
    assert(nfracs==npows-1);
    RooArgList *fracs = new RooArgList();
    RooArgList *pows = new RooArgList();
    for (int i=1; i<=nfracs; i++){
      string name =  Form("%s_f%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.9-float(i-1)*1./nfracs,0.,1.)));
      //params[name]->removeRange();
      fracs->add(*params[name]);
    }
    for (int i=1; i<=npows; i++){
      string name =  Form("%s_p%d",prefix.c_str(),i);
      string ename =  Form("%s_e%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),TMath::Max(-9.,-1.*(i+1)),-9.,1.)));
      //params[name]->removeRange();
      utilities.insert(pair<string,RooAbsPdf*>(ename, new RooPower(ename.c_str(),ename.c_str(),*obs_var,*params[name])));
      pows->add(*utilities[ename]);
    }
    //cout << "RooArgLists..." << endl;
    //fracs->Print("v");
    //pows->Print("v");
    //cout << "Function..." << endl;
    RooAbsPdf *pow = new RooAddPdf(prefix.c_str(),prefix.c_str(),*pows,*fracs,true); 
    //pow->Print("v");
    return pow;
    //bkgPdfs.insert(pair<string,RooAbsPdf*>(pow->GetName(),pow));
  }
}

RooAbsPdf* PdfModelBuilderFAN::getLaurentSeries(string prefix, int order){
 
  int nlower=int(ceil(order/2.));
  int nhigher=order-nlower;
  // first do 0th order
  RooArgList *pows = new RooArgList();
  RooArgList *plist = new RooArgList();
  string pname =  Form("%s_pow0",prefix.c_str());
  utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.))));
  pows->add(*utilities[pname]);

  // even terms
  for (int i=1; i<=nlower; i++){
    string name = Form("%s_l%d",prefix.c_str(),i);
    params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.25/order,0.000001,0.999999)));
    plist->add(*params[name]);
    string pname =  Form("%s_powl%d",prefix.c_str(),i);
    utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.-i))));
    pows->add(*utilities[pname]);
  }
  // odd terms
  for (int i=1; i<=nhigher; i++){
    string name = Form("%s_h%d",prefix.c_str(),i);
    params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.25/order,0.000001,0.999999)));
    plist->add(*params[name]);
    string pname =  Form("%s_powh%d",prefix.c_str(),i);
    utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.+i))));
    pows->add(*utilities[pname]);
  }
  RooAddPdf *pdf = new RooAddPdf(prefix.c_str(),prefix.c_str(),*pows,*plist,true);
  return pdf;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
}



RooAbsPdf* PdfModelBuilderFAN::getLaurentSeriesFAN(string prefix, int order, float constant){
 
  int nlower=int(ceil(order/2.));
  int nhigher=order-nlower;
  // first do 0th order
  RooArgList *pows = new RooArgList();
  RooArgList *plist = new RooArgList();
  string pname =  Form("%s_pow0",prefix.c_str());
  utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(constant))));  
  pows->add(*utilities[pname]);

  // even terms
  for (int i=1; i<=nlower; i++){
    string name = Form("%s_l%d",prefix.c_str(),i);
    params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.25/order,0.000001,0.999999)));
    plist->add(*params[name]);
    string pname =  Form("%s_powl%d",prefix.c_str(),i);
    utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(constant-i))));  
    pows->add(*utilities[pname]);
  }
  // odd terms
  for (int i=1; i<=nhigher; i++){
    string name = Form("%s_h%d",prefix.c_str(),i);
    params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.25/order,0.000001,0.999999)));
    plist->add(*params[name]);
    string pname =  Form("%s_powh%d",prefix.c_str(),i);
    utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(constant+i))));  
    pows->add(*utilities[pname]);
  }
  RooAddPdf *pdf = new RooAddPdf(prefix.c_str(),prefix.c_str(),*pows,*plist,true);
  return pdf;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
}


RooAbsPdf* PdfModelBuilderFAN::getSumOfVoigtians(string prefix, int nVoigtians, bool recursive){
  RooArgList *voigtians = new RooArgList();
  RooArgList *coeffs = new RooArgList();

  for (int g=0; g<nVoigtians; g++){

    RooRealVar *mean = new RooRealVar(Form("%s_mean_p%d",prefix.c_str(),g),Form("%s_mean_p%d",prefix.c_str(),g),91.,80.,100.);
    RooRealVar *sigma = new RooRealVar(Form("%s_sigma_p%d",prefix.c_str(),g),Form("%s_sigma_p%d",prefix.c_str(),g),2.,0.1,20.);
    RooRealVar *width = new RooRealVar(Form("%s_width_p%d",prefix.c_str(),g),Form("%s_width_p%d",prefix.c_str(),g),2.,0.1,20.);
    RooVoigtian *voig = new RooVoigtian(Form("%s_voig_g%d",prefix.c_str(),g),Form("%s_voig_g%d",prefix.c_str(),g),*obs_var,*mean,*sigma, *width);

    params.insert(pair<string,RooRealVar*>(string(sigma->GetName()),sigma));
    params.insert(pair<string,RooRealVar*>(string(mean->GetName()),mean));
    params.insert(pair<string,RooRealVar*>(string(width->GetName()),width));
    voigtians->add(*voig);

    if (g<nVoigtians-1) {
      RooRealVar *frac = new RooRealVar(Form("%s_frac_g%d",prefix.c_str(),g),Form("%s_frac_g%d",prefix.c_str(),g),0.01,0.000001,0.999999);
      params.insert(pair<string,RooRealVar*>(string(frac->GetName()),frac));
      coeffs->add(*frac);
    }

    forceFracUnity_=true;
    if (g==nVoigtians-1 && forceFracUnity_) {
      string formula="1.";
      for (int i=0; i<nVoigtians-1; i++) formula += Form("-@%d",i);
      RooFormulaVar *recFrac = new RooFormulaVar(Form("%s_frac_g%d",prefix.c_str(),g),Form("%s_frac_g%d",prefix.c_str(),g),formula.c_str(),*coeffs);
      prods.insert(pair<string,RooFormulaVar*>(string(recFrac->GetName()),recFrac));
      coeffs->add(*recFrac);
    }
  }
  assert(voigtians->getSize()==nVoigtians && coeffs->getSize()==nVoigtians-(1*!forceFracUnity_));
  RooAbsPdf *temp = new RooAddPdf(prefix.c_str(),prefix.c_str(),*voigtians,*coeffs,recursive);
  return temp;
}





RooAbsPdf* PdfModelBuilderFAN::getConvDoubleCB(string prefix, RooRealVar *alphaCB1, RooRealVar *alphaCB2){

    RooRealVar *nCB1 = new RooRealVar(Form("%s_nCB1_p%d",prefix.c_str(),0),Form("%s_nCB1_p%d",prefix.c_str(),0), 2.,0.01,5000.);
    RooRealVar *nCB2 = new RooRealVar(Form("%s_nCB2_p%d",prefix.c_str(),0),Form("%s_nCB2_p%d",prefix.c_str(),0), 2.,0.01,5000);
    RooRealVar *meanCB = new RooRealVar(Form("%s_mean_p%d",prefix.c_str(),0),Form("%s_mean_p%d",prefix.c_str(),0), 91.,80.,100.);
    RooRealVar *sigmaCB = new RooRealVar(Form("%s_sigma_p%d",prefix.c_str(),0),Form("%s_sigma_p%d",prefix.c_str(),0), 2., 0.1, 20.);

    RooAbsPdf *temp = new RooDoubleCBFast(prefix.c_str(),prefix.c_str(), *obs_var,*meanCB,*sigmaCB, *alphaCB1, *nCB1, *alphaCB2, *nCB2);


    return temp;
}



RooAbsPdf* PdfModelBuilderFAN::getConvDoubleCBnew(string prefix, float alphacb1, float alphacb2){

    RooRealVar *nCB1 = new RooRealVar(Form("%s_nCB1_p%d",prefix.c_str(),0),Form("%s_nCB1_p%d",prefix.c_str(),0), 2.,0.01,5000.);
    RooRealVar *nCB2 = new RooRealVar(Form("%s_nCB2_p%d",prefix.c_str(),0),Form("%s_nCB2_p%d",prefix.c_str(),0), 2.,0.01,5000);
    RooRealVar *meanCB = new RooRealVar(Form("%s_mean_p%d",prefix.c_str(),0),Form("%s_mean_p%d",prefix.c_str(),0), 91.,80.,100.);
    RooRealVar *sigmaCB = new RooRealVar(Form("%s_sigma_p%d",prefix.c_str(),0),Form("%s_sigma_p%d",prefix.c_str(),0), 2., 0.1, 20.);
    RooRealVar *alphaCB1 = new RooRealVar(Form("%s_alphaCB1_p%d",prefix.c_str(),0),Form("%s_alphaCB1_p%d",prefix.c_str(),0), alphacb1);   
    RooRealVar *alphaCB2 = new RooRealVar(Form("%s_alphaCB2_p%d",prefix.c_str(),0),Form("%s_alphaCB2_p%d",prefix.c_str(),0), alphacb2);   
    alphaCB1->setConstant(true);
    alphaCB2->setConstant(true);

    RooAbsPdf *temp = new RooDoubleCBFast(prefix.c_str(),prefix.c_str(), *obs_var,*meanCB,*sigmaCB, *alphaCB1, *nCB1, *alphaCB2, *nCB2);


    return temp;
}


RooAbsPdf* PdfModelBuilderFAN::getConvDoubleCBnew_Float(string prefix, float alphacb1, float alphacb2){

    RooRealVar *nCB1 = new RooRealVar(Form("%s_nCB1_p%d",prefix.c_str(),0),Form("%s_nCB1_p%d",prefix.c_str(),0), 2.,0.01,5000.);
    RooRealVar *nCB2 = new RooRealVar(Form("%s_nCB2_p%d",prefix.c_str(),0),Form("%s_nCB2_p%d",prefix.c_str(),0), 2.,0.01,5000);
    RooRealVar *meanCB = new RooRealVar(Form("%s_mean_p%d",prefix.c_str(),0),Form("%s_mean_p%d",prefix.c_str(),0), 91.,80.,100.);
    RooRealVar *sigmaCB = new RooRealVar(Form("%s_sigma_p%d",prefix.c_str(),0),Form("%s_sigma_p%d",prefix.c_str(),0), 2., 0.1, 20.);
    RooRealVar *alphaCB1 = new RooRealVar(Form("%s_alphaCB1_p%d",prefix.c_str(),0),Form("%s_alphaCB1_p%d",prefix.c_str(),0), alphacb1,0,2);    
    RooRealVar *alphaCB2 = new RooRealVar(Form("%s_alphaCB2_p%d",prefix.c_str(),0),Form("%s_alphaCB2_p%d",prefix.c_str(),0), alphacb2,0,2);   
    RooAbsPdf *temp = new RooDoubleCBFast(prefix.c_str(),prefix.c_str(), *obs_var,*meanCB,*sigmaCB, *alphaCB1, *nCB1, *alphaCB2, *nCB2);


    return temp;
}



pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >   PdfModelBuilderFAN::getOfficialSumVoigtians(string type, string name, string prefix, int nOfficial, int nVoigtians, float constant, float bernDownBound, float bernUpBound, bool recursive){
  RooArgList *OfficialSumVoigtians = new RooArgList();
  RooArgList *coeffs = new RooArgList();

  RooAbsPdf *pdfOfficial;
  RooAbsPdf *pdfSumVoigtians;

  string nameOff = name.append("_Off");
  if (type=="Bernstein")  pdfOfficial =  getBernsteinFAN(nameOff, nOfficial, bernDownBound, bernUpBound);  
  if (type=="Exponential")  pdfOfficial = getExponentialSingle(nameOff, nOfficial);
  if (type=="PowerLaw")  pdfOfficial = getPowerLawSingle(nameOff, nOfficial);
  if (type=="Laurent")   pdfOfficial = getLaurentSeriesFAN(nameOff, nOfficial, constant);  
  if (type=="Chebychev") pdfOfficial = getChebychevFAN(nameOff, nOfficial, bernDownBound, bernUpBound);
  if (type=="Polynomial") pdfOfficial = getPolynomial(nameOff, nOfficial);  

  string nameVoi = name.append("_Voi");
  pdfSumVoigtians = getSumOfVoigtians(nameVoi, nVoigtians, recursive);

  OfficialSumVoigtians->add(*pdfOfficial);
  OfficialSumVoigtians->add(*pdfSumVoigtians);

  RooRealVar *frac = new RooRealVar(Form("%s_frac_sum1",prefix.c_str()),Form("%s_frac_sum1",prefix.c_str()),0.01,0.000001,0.999999);
  coeffs->add(*frac);

  forceFracUnity_=true;
  if (forceFracUnity_) {
    string formula="1.";
    formula += Form("-@%d", 0);
    RooFormulaVar *recFrac = new RooFormulaVar(Form("%s_frac_sum2",prefix.c_str()),Form("%s_frac_sum2",prefix.c_str()),formula.c_str(),*coeffs);
    coeffs->add(*recFrac);
  }

  RooAbsPdf *temp = new RooAddPdf(prefix.c_str(),prefix.c_str(),*OfficialSumVoigtians,*coeffs,recursive);
  return pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >(temp, make_pair(pdfOfficial,pdfSumVoigtians));
}



pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >   PdfModelBuilderFAN::getOfficialSumVoigtiansForZeeFit(string type, string name, string prefix, int nOfficial, int nVoigtians, float constant, float bernDownBound, float bernUpBound, bool recursive){
  RooArgList *OfficialSumVoigtians = new RooArgList();
  RooArgList *coeffs = new RooArgList();

  RooAbsPdf *pdfOfficial;
  RooAbsPdf *pdfSumVoigtians;

  string nameOLD = name;
  string nameOff = name.append("_Off");
  if (type=="Bernstein")  pdfOfficial =  getBernsteinFAN(nameOff, nOfficial, bernDownBound, bernUpBound);  
  if (type=="Exponential")  pdfOfficial = getExponentialSingle(nameOff, nOfficial);
  if (type=="PowerLaw")  pdfOfficial = getPowerLawSingle(nameOff, nOfficial);
  if (type=="Laurent")   pdfOfficial = getLaurentSeriesFAN(nameOff, nOfficial, constant);  
  if (type=="Chebychev") pdfOfficial = getChebychevFAN(nameOff, nOfficial, bernDownBound, bernUpBound);
  if (type=="Polynomial") pdfOfficial = getPolynomial(nameOff, nOfficial);  

  string nameVoi = nameOLD;  
  pdfSumVoigtians = getSumOfVoigtians(nameVoi, nVoigtians, recursive);

  OfficialSumVoigtians->add(*pdfOfficial);
  OfficialSumVoigtians->add(*pdfSumVoigtians);

  RooRealVar *frac = new RooRealVar(Form("%s_frac_sum1",prefix.c_str()),Form("%s_frac_sum1",prefix.c_str()),0.01,0.000001,0.999999);
  coeffs->add(*frac);

  forceFracUnity_=true;
  if (forceFracUnity_) {
    string formula="1.";
    formula += Form("-@%d", 0);
    RooFormulaVar *recFrac = new RooFormulaVar(Form("%s_frac_sum2",prefix.c_str()),Form("%s_frac_sum2",prefix.c_str()),formula.c_str(),*coeffs);
    coeffs->add(*recFrac);
  }

  RooAbsPdf *temp = new RooAddPdf(prefix.c_str(),prefix.c_str(),*OfficialSumVoigtians,*coeffs,recursive);
  return pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >(temp, make_pair(pdfOfficial,pdfSumVoigtians));
}




pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >   PdfModelBuilderFAN::getOfficialSumVoigtiansForZeeFitDoubleCB(string type, string name, string prefix, int nOfficial, int nVoigtians, float constant, float bernDownBound, float bernUpBound, float alphacb1, float alphacb2, bool recursive){
  RooArgList *OfficialSumVoigtians = new RooArgList();
  RooArgList *coeffs = new RooArgList();

  RooAbsPdf *pdfOfficial;
  RooAbsPdf *pdfSumVoigtians;

  string nameOLD = name;
  string nameOff = name.append("_Off");
  if (type=="Bernstein")  pdfOfficial =  getBernsteinFAN(nameOff, nOfficial, bernDownBound, bernUpBound);  
  if (type=="Exponential")  pdfOfficial = getExponentialSingle(nameOff, nOfficial);
  if (type=="PowerLaw")  pdfOfficial = getPowerLawSingle(nameOff, nOfficial);
  if (type=="Laurent")   pdfOfficial = getLaurentSeriesFAN(nameOff, nOfficial, constant);  
  if (type=="Chebychev") pdfOfficial = getChebychevFAN(nameOff, nOfficial, bernDownBound, bernUpBound);  
  if (type=="Polynomial") pdfOfficial = getPolynomial(nameOff, nOfficial);  

  string nameVoi = nameOLD;  
  pdfSumVoigtians = getConvDoubleCBnew(nameVoi, alphacb1, alphacb2);

  OfficialSumVoigtians->add(*pdfOfficial);
  OfficialSumVoigtians->add(*pdfSumVoigtians);

  RooRealVar *frac = new RooRealVar(Form("%s_frac_sum1",prefix.c_str()),Form("%s_frac_sum1",prefix.c_str()),0.01,0.000001,0.999999);
  coeffs->add(*frac);

  forceFracUnity_=true;
  if (forceFracUnity_) {
    string formula="1.";
    formula += Form("-@%d", 0);
    RooFormulaVar *recFrac = new RooFormulaVar(Form("%s_frac_sum2",prefix.c_str()),Form("%s_frac_sum2",prefix.c_str()),formula.c_str(),*coeffs);
    coeffs->add(*recFrac);
  }

  RooAbsPdf *temp = new RooAddPdf(prefix.c_str(),prefix.c_str(),*OfficialSumVoigtians,*coeffs,recursive);
  return pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >(temp, make_pair(pdfOfficial,pdfSumVoigtians));
}


pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >   PdfModelBuilderFAN::getOfficialSumVoigtiansForZeeFitDoubleCB_Float(string type, string name, string prefix, int nOfficial, int nVoigtians, float constant, float bernDownBound, float bernUpBound, float alphacb1, float alphacb2, bool recursive){
  RooArgList *OfficialSumVoigtians = new RooArgList();
  RooArgList *coeffs = new RooArgList();

  RooAbsPdf *pdfOfficial;
  RooAbsPdf *pdfSumVoigtians;

  string nameOLD = name;
  string nameOff = name.append("_Off");
  if (type=="Bernstein")  pdfOfficial =  getBernsteinFAN(nameOff, nOfficial, bernDownBound, bernUpBound);  
  if (type=="Exponential")  pdfOfficial = getExponentialSingle(nameOff, nOfficial);
  if (type=="PowerLaw")  pdfOfficial = getPowerLawSingle(nameOff, nOfficial);
  if (type=="Laurent")   pdfOfficial = getLaurentSeriesFAN(nameOff, nOfficial, constant);  
  if (type=="Chebychev") pdfOfficial = getChebychevFAN(nameOff, nOfficial, bernDownBound, bernUpBound);  
  if (type=="Polynomial") pdfOfficial = getPolynomial(nameOff, nOfficial);  

  string nameVoi = nameOLD;  
  pdfSumVoigtians = getConvDoubleCBnew_Float(nameVoi, alphacb1, alphacb2);

  OfficialSumVoigtians->add(*pdfOfficial);
  OfficialSumVoigtians->add(*pdfSumVoigtians);

  RooRealVar *frac = new RooRealVar(Form("%s_frac_sum1",prefix.c_str()),Form("%s_frac_sum1",prefix.c_str()),0.01,0.000001,0.999999);
  coeffs->add(*frac);

  forceFracUnity_=true;
  if (forceFracUnity_) {
    string formula="1.";
    formula += Form("-@%d", 0);
    RooFormulaVar *recFrac = new RooFormulaVar(Form("%s_frac_sum2",prefix.c_str()),Form("%s_frac_sum2",prefix.c_str()),formula.c_str(),*coeffs);
    coeffs->add(*recFrac);
  }

  RooAbsPdf *temp = new RooAddPdf(prefix.c_str(),prefix.c_str(),*OfficialSumVoigtians,*coeffs,recursive);
  return pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >(temp, make_pair(pdfOfficial,pdfSumVoigtians));
}



pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >   PdfModelBuilderFAN::getOfficialSumVoigtiansFix(string type, string name, string prefix, int nOfficial, RooAbsPdf* pdfVoiFix, float constant, float bernDownBound, float bernUpBound, bool recursive){
  RooArgList *OfficialSumVoigtians = new RooArgList();
  RooArgList *coeffs = new RooArgList();

  RooAbsPdf *pdfOfficial;
  RooAbsPdf *pdfSumVoigtians;

  string nameOff = name.append("_Off");
  if (type=="Bernstein")  pdfOfficial =  getBernsteinFAN(nameOff, nOfficial, bernDownBound, bernUpBound);  
  if (type=="Exponential")  pdfOfficial = getExponentialSingle(nameOff, nOfficial);
  if (type=="PowerLaw")  pdfOfficial = getPowerLawSingle(nameOff, nOfficial);
  if (type=="Laurent")   pdfOfficial = getLaurentSeriesFAN(nameOff, nOfficial, constant);  
  if (type=="Chebychev") pdfOfficial = getChebychevFAN(nameOff, nOfficial, bernDownBound, bernUpBound);
  if (type=="Polynomial") pdfOfficial = getPolynomial(nameOff, nOfficial);  

  pdfSumVoigtians = pdfVoiFix;
 
  if( pdfSumVoigtians == 0)  return pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >(pdfOfficial, make_pair(pdfOfficial,pdfOfficial));

  OfficialSumVoigtians->add(*pdfOfficial);
  OfficialSumVoigtians->add(*pdfSumVoigtians);

  RooRealVar *frac = new RooRealVar(Form("%s_frac_sum1",prefix.c_str()),Form("%s_frac_sum1",prefix.c_str()),0.01,0.000001,0.999999);
  coeffs->add(*frac);

  forceFracUnity_=true;
  if (forceFracUnity_) {
    string formula="1.";
    formula += Form("-@%d", 0);
    RooFormulaVar *recFrac = new RooFormulaVar(Form("%s_frac_sum2",prefix.c_str()),Form("%s_frac_sum2",prefix.c_str()),formula.c_str(),*coeffs);
    coeffs->add(*recFrac);
  }

  RooAbsPdf *temp = new RooAddPdf(prefix.c_str(),prefix.c_str(),*OfficialSumVoigtians,*coeffs,recursive);
  return pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >(temp, make_pair(pdfOfficial,pdfSumVoigtians));
}


pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >   PdfModelBuilderFAN::getOfficialSumVoigtiansFixNew(string type, string name, string prefix, int nOfficial, float voiMean, float voiMeanErrorL, float voiMeanErrorH, float voiSigma, float voiSigmaErrorL, float voiSigmaErrorH, float voiWidth, float voiWidthErrorL, float voiWidthErrorH, float ErrorRange, float constant, bool constOrSo, float bernDownBound, float bernUpBound, bool recursive){
  RooArgList *OfficialSumVoigtians = new RooArgList();
  RooArgList *coeffs = new RooArgList();

  RooAbsPdf *pdfOfficial;
  RooAbsPdf *pdfSumVoigtians;

  string nameOff = name.append("_Off");
  if (type=="Bernstein")  pdfOfficial =  getBernsteinFAN(nameOff, nOfficial, bernDownBound, bernUpBound);  
  if (type=="Exponential")  pdfOfficial = getExponentialSingle(nameOff, nOfficial);
  if (type=="PowerLaw")  pdfOfficial = getPowerLawSingle(nameOff, nOfficial);
  if (type=="Laurent")   pdfOfficial = getLaurentSeriesFAN(nameOff, nOfficial, constant);  
  if (type=="Chebychev") pdfOfficial = getChebychevFAN(nameOff, nOfficial, bernDownBound, bernUpBound);  
  if (type=="Polynomial") pdfOfficial = getPolynomial(nameOff, nOfficial);  

  //define voi
  if(voiMean != 0){
        RooRealVar *mean = new RooRealVar(Form("%s_Fvoimean",prefix.c_str()),Form("%s_Fvoimean",prefix.c_str()), voiMean, voiMean+ErrorRange*voiMeanErrorL, voiMean+ErrorRange*voiMeanErrorH);
        mean->setConstant(constOrSo);
        RooRealVar *sigma = new RooRealVar(Form("%s_Fvoisigma",prefix.c_str()),Form("%s_Fvoisigma",prefix.c_str()), voiSigma, voiSigma+ErrorRange*voiSigmaErrorL, voiSigma+ErrorRange*voiSigmaErrorH);
        sigma->setConstant(constOrSo);
        RooRealVar *width = new RooRealVar(Form("%s_Fvoiwidth",prefix.c_str()),Form("%s_Fvoiwidth",prefix.c_str()), voiWidth, voiWidth+ErrorRange*voiWidthErrorL, voiWidth+ErrorRange*voiWidthErrorH);
        width->setConstant(constOrSo);
        RooVoigtian *voig = new RooVoigtian(Form("%s_Fvoig",prefix.c_str()),Form("%s_Fvoig",prefix.c_str()),*obs_var,*mean,*sigma, *width); 
        pdfSumVoigtians = voig;
  }

  if( pdfSumVoigtians == 0)  return pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >(pdfOfficial, make_pair(pdfOfficial,pdfOfficial));

  OfficialSumVoigtians->add(*pdfOfficial);
  OfficialSumVoigtians->add(*pdfSumVoigtians);

  RooRealVar *frac = new RooRealVar(Form("%s_frac_sum1",prefix.c_str()),Form("%s_frac_sum1",prefix.c_str()),0.01,0.000001,0.999999);
  coeffs->add(*frac);

  forceFracUnity_=true;
  if (forceFracUnity_) {
    string formula="1.";
    formula += Form("-@%d", 0);
    RooFormulaVar *recFrac = new RooFormulaVar(Form("%s_frac_sum2",prefix.c_str()),Form("%s_frac_sum2",prefix.c_str()),formula.c_str(),*coeffs);
    coeffs->add(*recFrac);
  }

  RooAbsPdf *temp = new RooAddPdf(prefix.c_str(),prefix.c_str(),*OfficialSumVoigtians,*coeffs,recursive);
  return pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >(temp, make_pair(pdfOfficial,pdfSumVoigtians));
}


pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >   PdfModelBuilderFAN::getOfficialSumVoigtiansFixNewDoubleCB(string type, string name, string prefix, int nOfficial, float voiMean, float voiMeanErrorL, float voiMeanErrorH, float voiSigma, float voiSigmaErrorL, float voiSigmaErrorH, float voinCB1, float voinCB1ErrorL, float voinCB1ErrorH, float voinCB2, float voinCB2ErrorL, float voinCB2ErrorH, float voialphaCB1, float voialphaCB2, float ErrorRange, float constant, bool constOrSo, float bernDownBound, float bernUpBound, bool recursive){
  RooArgList *OfficialSumVoigtians = new RooArgList();
  RooArgList *coeffs = new RooArgList();

  RooAbsPdf *pdfOfficial;
  RooAbsPdf *pdfSumVoigtians;

  string nameOff = name.append("_Off");
  if (type=="Bernstein")  pdfOfficial =  getBernsteinFAN(nameOff, nOfficial, bernDownBound, bernUpBound);  
  if (type=="Exponential")  pdfOfficial = getExponentialSingle(nameOff, nOfficial);
  if (type=="PowerLaw")  pdfOfficial = getPowerLawSingle(nameOff, nOfficial);
  if (type=="Laurent")   pdfOfficial = getLaurentSeriesFAN(nameOff, nOfficial, constant);  
  if (type=="Chebychev") pdfOfficial = getChebychevFAN(nameOff, nOfficial, bernDownBound, bernUpBound);
  if (type=="Polynomial") pdfOfficial = getPolynomial(nameOff, nOfficial);  


  //define voi
  if(voiMean != 0){
        RooRealVar *mean = new RooRealVar(Form("%s_Fdcbmean",prefix.c_str()),Form("%s_Fdcbmean",prefix.c_str()), voiMean, voiMean+ErrorRange*voiMeanErrorL, voiMean+ErrorRange*voiMeanErrorH);
        mean->setConstant(constOrSo);
        RooRealVar *sigma = new RooRealVar(Form("%s_Fdcbsigma",prefix.c_str()),Form("%s_Fdcbsigma",prefix.c_str()), voiSigma, voiSigma+ErrorRange*voiSigmaErrorL, voiSigma+ErrorRange*voiSigmaErrorH);
        sigma->setConstant(constOrSo);
        RooRealVar *nCB1 = new RooRealVar(Form("%s_FdcbnCB1",prefix.c_str()),Form("%s_FdcbnCB1",prefix.c_str()), voinCB1, voinCB1+ErrorRange*voinCB1ErrorL,  voinCB1+ErrorRange*voinCB1ErrorH);
        nCB1->setConstant(constOrSo);
        RooRealVar *nCB2 = new RooRealVar(Form("%s_FdcbnCB2",prefix.c_str()),Form("%s_FdcbnCB2",prefix.c_str()), voinCB2, voinCB2+ErrorRange*voinCB2ErrorL,  voinCB2+ErrorRange*voinCB2ErrorH);
        nCB2->setConstant(constOrSo);
        RooRealVar *alphaCB1 = new RooRealVar(Form("%s_FdcbalphaCB1",prefix.c_str()),Form("%s_FdcbalphaCB1",prefix.c_str()),voialphaCB1);
        alphaCB1->setConstant(true);
        RooRealVar *alphaCB2 = new RooRealVar(Form("%s_FdcbalphaCB2",prefix.c_str()),Form("%s_FdcbalphaCB2",prefix.c_str()),voialphaCB2);
        alphaCB2->setConstant(true);
        RooAbsPdf *dcb = new RooDoubleCBFast(Form("%s_Fdcb",prefix.c_str()),Form("%s_Fdcb",prefix.c_str()), *obs_var,*mean,*sigma, *alphaCB1, *nCB1, *alphaCB2, *nCB2);
        pdfSumVoigtians = dcb;
  }

  if( pdfSumVoigtians == 0)  return pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >(pdfOfficial, make_pair(pdfOfficial,pdfOfficial));

  OfficialSumVoigtians->add(*pdfOfficial);
  OfficialSumVoigtians->add(*pdfSumVoigtians);

  RooRealVar *frac = new RooRealVar(Form("%s_frac_sum1",prefix.c_str()),Form("%s_frac_sum1",prefix.c_str()),0.01,0.000001,0.999999);
  coeffs->add(*frac);

  forceFracUnity_=true;
  if (forceFracUnity_) {
    string formula="1.";
    formula += Form("-@%d", 0);
    RooFormulaVar *recFrac = new RooFormulaVar(Form("%s_frac_sum2",prefix.c_str()),Form("%s_frac_sum2",prefix.c_str()),formula.c_str(),*coeffs);
    coeffs->add(*recFrac);
  }

  RooAbsPdf *temp = new RooAddPdf(prefix.c_str(),prefix.c_str(),*OfficialSumVoigtians,*coeffs,recursive);
  return pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >(temp, make_pair(pdfOfficial,pdfSumVoigtians));
}




RooAbsPdf* PdfModelBuilderFAN::getKeysPdf(string prefix){
  if (!keysPdfAttributesSet){
    cerr << "ERROR -- keysPdf attributes not set" << endl;
    exit(1);
  }
  return new RooKeysPdf(prefix.c_str(),prefix.c_str(),*obs_var,*keysPdfData,RooKeysPdf::MirrorBoth,keysPdfRho);
}

RooAbsPdf* PdfModelBuilderFAN::getPdfFromFile(string &prefix){
  vector<string> details;
  split(details,prefix,boost::is_any_of(","));

  string fname = details[2];
  string wsname = details[1];
  string pdfname = details[0];

  TFile *tempFile = TFile::Open(fname.c_str());
  if (!tempFile){
    cerr << "PdfModelBuilderFAN::getPdfFromFile -- file not found " << fname << endl;
    assert(0);
  }
  RooWorkspace *tempWS = (RooWorkspace*)tempFile->Get(wsname.c_str());
  if (!tempWS){
    cerr << "PdfModelBuilderFAN::getPdfFromFile -- workspace not found " << wsname << endl;
    assert(0);
  }
  RooAbsPdf *tempPdf = (RooAbsPdf*)tempWS->pdf(pdfname.c_str());
  if (!tempPdf){
    cerr << "PdfModelBuilderFAN::getPdfFromFile -- pdf not found " << pdfname << endl;
    assert(0);
  }
  prefix = pdfname;
  RooAbsPdf *pdf = (RooAbsPdf*)tempPdf->Clone(prefix.c_str());
  tempFile->Close();
  delete tempFile;
  return pdf;
}



RooAbsPdf* PdfModelBuilderFAN::getExponentialSingle(string prefix, int order){
  
  if (order%2==0){
    cerr << "ERROR -- addExponential -- only odd number of params allowed" << endl;
    return NULL;
  }
  else {
    int nfracs=(order-1)/2;
    int nexps=order-nfracs;
    assert(nfracs==nexps-1);
    RooArgList *fracs = new RooArgList();
    RooArgList *exps = new RooArgList();
    for (int i=1; i<=nfracs; i++){
      string name =  Form("%s_f%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.9-float(i-1)*1./nfracs,0.0001,0.9999)));
      fracs->add(*params[name]);
    }
    for (int i=1; i<=nexps; i++){
      string name =  Form("%s_p%d",prefix.c_str(),i);
      string ename =  Form("%s_e%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),TMath::Max(-1.,-0.04*(i+1)),-1.,0.)));
      utilities.insert(pair<string,RooAbsPdf*>(ename, new RooExponential(ename.c_str(),ename.c_str(),*obs_var,*params[name])));
      exps->add(*utilities[ename]);
    }
    //fracs->Print("v");
    //exps->Print("v");
    RooAbsPdf *exp = new RooAddPdf(prefix.c_str(),prefix.c_str(),*exps,*fracs,true);
    //exp->Print("v");
    cout << "--------------------------" << endl;
    return exp;
    //bkgPdfs.insert(pair<string,RooAbsPdf*>(exp->GetName(),exp));

  }
}



void PdfModelBuilderFAN::addBkgPdf(string type, int nParams, string name, bool cache){
 
  if (!obs_var_set){
    cerr << "ERROR -- obs Var has not been set!" << endl;
    exit(1);
  }
  bool found=false;
  for (vector<string>::iterator it=recognisedPdfTypes.begin(); it!=recognisedPdfTypes.end(); it++){
    if (*it==type) found=true;
  }
  if (!found){
    cerr << "Pdf of type " << type << " is not recognised!" << endl;
    exit(1);
  }
  RooAbsPdf *pdf;
  if (type=="Bernstein") pdf = getBernstein(name,nParams);
  if (type=="Exponential") pdf = getExponentialSingle(name,nParams);
  if (type=="PowerLaw") pdf = getPowerLawSingle(name,nParams);
  if (type=="Laurent") pdf = getLaurentSeries(name,nParams);
  if (type=="KeysPdf") pdf = getKeysPdf(name);
  if (type=="File") pdf = getPdfFromFile(name);

  if (cache) {
    wsCache->import(*pdf);
    RooAbsPdf *cachePdf = wsCache->pdf(pdf->GetName());
    bkgPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
  }
  else {
    bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
  }

}



void PdfModelBuilderFAN::addBkgPdfFAN(string type, int nParams, string name, float constant, float bernDownBound, float bernUpBound, bool cache){
 
  if (!obs_var_set){
    cerr << "ERROR -- obs Var has not been set!" << endl;
    exit(1);
  }
  bool found=false;
  for (vector<string>::iterator it=recognisedPdfTypes.begin(); it!=recognisedPdfTypes.end(); it++){
    if (*it==type) found=true;
  }
  if (!found){
    cerr << "Pdf of type " << type << " is not recognised!" << endl;
    exit(1);
  }
  RooAbsPdf *pdf;
  if (type=="Bernstein") pdf = getBernsteinFAN(name,nParams, bernDownBound, bernUpBound);   
  if (type=="Exponential") pdf = getExponentialSingle(name,nParams);
  if (type=="PowerLaw") pdf = getPowerLawSingle(name,nParams);
  if (type=="Laurent") pdf = getLaurentSeriesFAN(name,nParams, constant);  
  if (type=="Chebychev") pdf = getChebychevFAN(name,nParams,bernDownBound, bernUpBound);  
  if (type=="Polynomial") pdf = getPolynomial(name,nParams);  
  if (type=="KeysPdf") pdf = getKeysPdf(name);
  if (type=="File") pdf = getPdfFromFile(name);

  if (cache) {
    wsCache->import(*pdf);
    RooAbsPdf *cachePdf = wsCache->pdf(pdf->GetName());
    bkgPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
  }
  else {
    bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
  }

}


void PdfModelBuilderFAN::addBkgPdfSumVoigtians(string type, int nParams, string name, int nBre, string ext, float constant, float bernDownBound, float bernUpBound, bool cache, bool recursive){

  if (!obs_var_set){
    cerr << "ERROR -- obs Var has not been set!" << endl;
    exit(1);
  }
  bool found=false;
  for (vector<string>::iterator it=recognisedPdfTypes.begin(); it!=recognisedPdfTypes.end(); it++){
    if (*it==type) found=true;
  }
  if (!found){
    cerr << "Pdf of type " << type << " is not recognised!" << endl;
    exit(1);
  }
  RooAbsPdf *pdf;
  if(nBre == 0){
       //if (type=="Bernstein") pdf = getBernstein(name,nParams);
       if (type=="Bernstein") pdf = getBernsteinFAN(name,nParams, bernDownBound, bernUpBound);   
       if (type=="Exponential") pdf = getExponentialSingle(name,nParams);
       if (type=="PowerLaw") pdf = getPowerLawSingle(name,nParams);
       //if (type=="Laurent") pdf = getLaurentSeries(name,nParams);
       if (type=="Laurent") pdf = getLaurentSeriesFAN(name,nParams, constant);  
       //if (type=="Chebychev") pdf = getChebychev(name,nParams);  
       if (type=="Chebychev") pdf = getChebychevFAN(name,nParams,bernDownBound, bernUpBound);  
       if (type=="Polynomial") pdf = getPolynomial(name,nParams);  
       if (type=="KeysPdf") pdf = getKeysPdf(name);
       if (type=="File") pdf = getPdfFromFile(name);
  }  else{
       pdf = getOfficialSumVoigtians(type, name, ext, nParams, nBre, constant, bernDownBound, bernUpBound, recursive).first;
  }

  if (cache) {
    wsCache->import(*pdf);
    RooAbsPdf *cachePdf = wsCache->pdf(pdf->GetName());
    bkgPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
  }
  else {
    bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
  }
}



void PdfModelBuilderFAN::addBkgPdfSumVoigtiansFix(string type, int nParams, string name, RooAbsPdf* pdfVoiFix, string ext, float constant, float bernDownBound, float bernUpBound, bool cache, bool recursive){

  if (!obs_var_set){
    cerr << "ERROR -- obs Var has not been set!" << endl;
    exit(1);
  }
  bool found=false;
  for (vector<string>::iterator it=recognisedPdfTypes.begin(); it!=recognisedPdfTypes.end(); it++){
    if (*it==type) found=true;
  }
  if (!found){
    cerr << "Pdf of type " << type << " is not recognised!" << endl;
    exit(1);
  }
  RooAbsPdf *pdf;
  if(pdfVoiFix == 0){
       if (type=="Bernstein") pdf = getBernsteinFAN(name,nParams, bernDownBound, bernUpBound);   
       if (type=="Exponential") pdf = getExponentialSingle(name,nParams);
       if (type=="PowerLaw") pdf = getPowerLawSingle(name,nParams);
       if (type=="Laurent") pdf = getLaurentSeriesFAN(name,nParams, constant);  
       if (type=="Chebychev") pdf = getChebychevFAN(name,nParams,bernDownBound, bernUpBound);  
       if (type=="Polynomial") pdf = getPolynomial(name,nParams); 
       if (type=="KeysPdf") pdf = getKeysPdf(name);
       if (type=="File") pdf = getPdfFromFile(name);
  }  else{
       pdf = getOfficialSumVoigtiansFix(type, name, ext, nParams, pdfVoiFix, constant, bernDownBound, bernUpBound, recursive).first;
  }

  if (cache) {
    wsCache->import(*pdf);
    RooAbsPdf *cachePdf = wsCache->pdf(pdf->GetName());
    bkgPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
  }
  else {
    bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
  }
}


void PdfModelBuilderFAN::addBkgPdfSumVoigtiansFixNew(string type, int nParams, string name, float voiMean, float voiMeanErrorL, float voiMeanErrorH, float voiSigma, float voiSigmaErrorL, float voiSigmaErrorH, float voiWidth, float voiWidthErrorL, float voiWidthErrorH, float ErrorRange, string ext, float constant, bool constOrSo, float bernDownBound, float bernUpBound, bool cache, bool recursive){

  if (!obs_var_set){
    cerr << "ERROR -- obs Var has not been set!" << endl;
    exit(1);
  }
  bool found=false;
  for (vector<string>::iterator it=recognisedPdfTypes.begin(); it!=recognisedPdfTypes.end(); it++){
    if (*it==type) found=true;
  }
  if (!found){
    cerr << "Pdf of type " << type << " is not recognised!" << endl;
    exit(1);
  }
  RooAbsPdf *pdf;
  if(voiMean == 0){
       if (type=="Bernstein") pdf = getBernsteinFAN(name,nParams,bernDownBound, bernUpBound);  
       if (type=="Exponential") pdf = getExponentialSingle(name,nParams);
       if (type=="PowerLaw") pdf = getPowerLawSingle(name,nParams);
       if (type=="Laurent") pdf = getLaurentSeriesFAN(name,nParams, constant);  
       if (type=="Chebychev") pdf = getChebychevFAN(name,nParams,bernDownBound, bernUpBound); 
       if (type=="Polynomial") pdf = getPolynomial(name,nParams); 
       if (type=="KeysPdf") pdf = getKeysPdf(name);
       if (type=="File") pdf = getPdfFromFile(name);
  }  else{
       pdf = getOfficialSumVoigtiansFixNew(type, name, ext, nParams, voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voiWidth, voiWidthErrorL, voiWidthErrorH, ErrorRange, constant, constOrSo, bernDownBound, bernUpBound, recursive).first;
  }

  if (cache) {
    wsCache->import(*pdf);
    RooAbsPdf *cachePdf = wsCache->pdf(pdf->GetName());
    bkgPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
  }
  else {
    bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
  }
}


void PdfModelBuilderFAN::addBkgPdfSumVoigtiansFixNewDoubleCB(string type, int nParams, string name, float voiMean, float voiMeanErrorL, float voiMeanErrorH, float voiSigma, float voiSigmaErrorL, float voiSigmaErrorH, float voinCB1, float voinCB1ErrorL, float voinCB1ErrorH, float voinCB2, float voinCB2ErrorL, float voinCB2ErrorH, float voialphaCB1, float voialphaCB2, float ErrorRange, string ext, float constant, bool constOrSo, float bernDownBound, float bernUpBound, bool cache, bool recursive){

  if (!obs_var_set){
    cerr << "ERROR -- obs Var has not been set!" << endl;
    exit(1);
  }
  bool found=false;
  for (vector<string>::iterator it=recognisedPdfTypes.begin(); it!=recognisedPdfTypes.end(); it++){
    if (*it==type) found=true;
  }
  if (!found){
    cerr << "Pdf of type " << type << " is not recognised!" << endl;
    exit(1);
  }
  RooAbsPdf *pdf;
  if(voiMean == 0){
       if (type=="Bernstein") pdf = getBernsteinFAN(name,nParams,bernDownBound, bernUpBound);  
       if (type=="Exponential") pdf = getExponentialSingle(name,nParams);
       if (type=="PowerLaw") pdf = getPowerLawSingle(name,nParams);
       if (type=="Laurent") pdf = getLaurentSeriesFAN(name,nParams, constant);  
       if (type=="Chebychev") pdf = getChebychevFAN(name,nParams,bernDownBound, bernUpBound); 
       if (type=="Polynomial") pdf = getPolynomial(name,nParams); 
       if (type=="KeysPdf") pdf = getKeysPdf(name);
       if (type=="File") pdf = getPdfFromFile(name);
  }  else{
       pdf = getOfficialSumVoigtiansFixNewDoubleCB(type, name, ext, nParams, voiMean, voiMeanErrorL, voiMeanErrorH, voiSigma, voiSigmaErrorL, voiSigmaErrorH, voinCB1, voinCB1ErrorL, voinCB1ErrorH, voinCB2, voinCB2ErrorL, voinCB2ErrorH, voialphaCB1, voialphaCB2, ErrorRange, constant, constOrSo, bernDownBound, bernUpBound, recursive).first;
  }

  if (cache) {
    wsCache->import(*pdf);
    RooAbsPdf *cachePdf = wsCache->pdf(pdf->GetName());
    bkgPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
  }
  else {
    bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
  }
}




void PdfModelBuilderFAN::setKeysPdfAttributes(RooDataSet *data, double rho){
  keysPdfData = data;
  keysPdfRho = rho;
  keysPdfAttributesSet=true;
}

void PdfModelBuilderFAN::setSignalPdf(RooAbsPdf *pdf, RooRealVar *norm){
  sigPdf=pdf;
  sigNorm=norm;
  signal_set=true;
}

void PdfModelBuilderFAN::setSignalPdfFromMC(RooDataSet *data){
  
  RooDataHist *sigMCBinned = new RooDataHist(Form("roohist_%s",data->GetName()),Form("roohist_%s",data->GetName()),RooArgSet(*obs_var),*data);
  sigPdf = new RooHistPdf(Form("pdf_%s",data->GetName()),Form("pdf_%s",data->GetName()),RooArgSet(*obs_var),*sigMCBinned);
  sigNorm = new RooConstVar(Form("sig_events_%s",data->GetName()),Form("sig_events_%s",data->GetName()),data->sumEntries());
  signal_set=true;
}

void PdfModelBuilderFAN::makeSBPdfs(bool cache){
  
  if (!signal_set){
    cerr << "ERROR - no signal model set!" << endl;
    exit(1);
  }
  if (!signal_modifier_set){
    cerr << "ERROR - no signal modifier set!" << endl;
    exit(1);
  }
 
  if (sigNorm) {
    sigYield = new RooProduct("sig_yield","sig_yield",RooArgSet(*signalModifier,*sigNorm));
  }
  else {
    sigYield = signalModifier;
  }
  bkgYield = new RooRealVar("bkg_yield","bkg_yield",1000.,0.,1.e6);

  for (map<string,RooAbsPdf*>::iterator bkg=bkgPdfs.begin(); bkg!=bkgPdfs.end(); bkg++){
    RooAbsPdf *sbMod = new RooAddPdf(Form("sb_%s",bkg->first.c_str()),Form("sb_%s",bkg->first.c_str()),RooArgList(*(bkg->second),*sigPdf),RooArgList(*bkgYield,*sigYield));
    if (cache) {
      wsCache->import(*sbMod,RecycleConflictNodes());
      RooAbsPdf *cachePdf = (RooAbsPdf*)wsCache->pdf(sbMod->GetName());
      signalModifier = (RooRealVar*)wsCache->var(signalModifier->GetName());
      sbPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
    }
    else {
      sbPdfs.insert(pair<string,RooAbsPdf*>(sbMod->GetName(),sbMod));
    }
  }
}

map<string,RooAbsPdf*> PdfModelBuilderFAN::getBkgPdfs(){
  return bkgPdfs;
}

map<string,RooAbsPdf*> PdfModelBuilderFAN::getSBPdfs(){
  return sbPdfs;
}

RooAbsPdf* PdfModelBuilderFAN::getSigPdf(){
  return sigPdf;
}

void PdfModelBuilderFAN::plotPdfsToData(RooAbsData *data, int binning, string name, bool bkgOnly,string specificPdfName){
  
  TCanvas *canv = new TCanvas();
  bool specPdf=false;
  if (specificPdfName!="") specPdf=true;

  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) pdfSet = bkgPdfs;
  else pdfSet = sbPdfs;
  
  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    if (specPdf && it->first!=specificPdfName && specificPdfName!="NONE") continue;
    RooPlot *plot = obs_var->frame();
    data->plotOn(plot,Binning(binning));
    if (specificPdfName!="NONE") {
	 it->second->plotOn(plot);
	 it->second->paramOn(plot,RooFit::Layout(0.34,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
    }	
    plot->Draw();
    //canv->Print(Form("%s_%s.pdf",name.c_str(),it->first.c_str()));  
    canv->Print(Form("%s_%s.png",name.c_str(),it->first.c_str()));
  }
  delete canv;
}

void PdfModelBuilderFAN::fitToData(RooAbsData *data, bool bkgOnly, bool cache, bool print){
  
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) pdfSet = bkgPdfs;
  else pdfSet = sbPdfs;

  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    RooFitResult *fit = (RooFitResult*)it->second->fitTo(*data,Save(true));
    if (print){
      cout << "Fit Res Before: " << endl;
      fit->floatParsInit().Print("v");
      cout << "Fit Res After: " << endl;
      fit->floatParsFinal().Print("v");
    }
    if (cache) {
      RooArgSet *fitargs = (RooArgSet*)it->second->getParameters(*obs_var);
      // remove the signal strength since this will be set AFTER fitting the background 
      fitargs->remove(*signalModifier); 
      wsCache->defineSet(Form("%s_params",it->first.c_str()),*fitargs);
      wsCache->defineSet(Form("%s_observs",it->first.c_str()),*obs_var);
      wsCache->saveSnapshot(it->first.c_str(),*fitargs,true);
      if (print) {
        cout << "Cached values: " << endl;
        fitargs->Print("v");
      }
    }
  }
  if (bkgOnly) bkgHasFit=true;
  else sbHasFit=true;
}



void PdfModelBuilderFAN::fitToDataFAN(RooAbsData *data, bool bkgOnly, bool cache, bool print){
 
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) pdfSet = bkgPdfs;
  else pdfSet = sbPdfs;

  
  map<string,RooArgSet> PdfInitial;
  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
       RooArgSet *paraInitial = (RooArgSet*)it->second->getParameters((const RooArgSet*)(0));
       RooArgSet preParaInitial;
       paraInitial->snapshot(preParaInitial);
       PdfInitial.insert(pair<string,RooArgSet>(it->first, preParaInitial));     
  }
  

  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    RooArgSet *paraInitial = (RooArgSet*)it->second->getParameters((const RooArgSet*)(0)); 
    cout << "I am here +++++++++++" << endl;
    paraInitial->Print("v");
    for (map<string,RooArgSet>::iterator im=PdfInitial.begin(); im!=PdfInitial.end(); im++){
       cout << im->first << endl;
       im->second.Print("v");
       if(it->first == im->first)    paraInitial->assignValueOnly(im->second);    
       cout << "I am here ***************" << endl;
       paraInitial->Print("v");
       cout << "I am here ***************" << endl;
    }
   

    RooFitResult *fit = (RooFitResult*)it->second->fitTo(*data,Save(true));
    if (print){
      cout << "Fit Res Before: " << endl;
      fit->floatParsInit().Print("v");
      cout << "Fit Res After: " << endl;
      fit->floatParsFinal().Print("v");
    }
    if (cache) {
      RooArgSet *fitargs = (RooArgSet*)it->second->getParameters(*obs_var);
      // remove the signal strength since this will be set AFTER fitting the background 
      fitargs->remove(*signalModifier); 
      wsCache->defineSet(Form("%s_params",it->first.c_str()),*fitargs);
      wsCache->defineSet(Form("%s_observs",it->first.c_str()),*obs_var);
      wsCache->saveSnapshot(it->first.c_str(),*fitargs,true);
      if (print) {
        cout << "Cached values: " << endl;
        fitargs->Print("v");
      }
    }
  }
  if (bkgOnly) bkgHasFit=true;
  else sbHasFit=true;
}



void PdfModelBuilderFAN::setSeed(int seed){
  RooRandom::randomGenerator()->SetSeed(seed);
}

RooDataSet* PdfModelBuilderFAN::makeHybridDataset(vector<float> switchOverMasses, vector<RooDataSet*> dataForHybrid){
  
  assert(switchOverMasses.size()==dataForHybrid.size()-1);

  vector<string> cut_strings;
  cut_strings.push_back("cutstring0");
  obs_var->setRange("cutstring0",obs_var->getMin(),switchOverMasses[0]);
  for (unsigned int i=1; i<switchOverMasses.size(); i++){
    cut_strings.push_back(Form("cutstring%d",i));
    obs_var->setRange(Form("cutstring%d",i),switchOverMasses[i-1],switchOverMasses[i]);
  }
  cut_strings.push_back(Form("cutstring%d",int(switchOverMasses.size())));
  obs_var->setRange(Form("cutstring%d",int(switchOverMasses.size())),switchOverMasses[switchOverMasses.size()-1],obs_var->getMax());
  
  obs_var->Print("v");
  assert(cut_strings.size()==dataForHybrid.size());
  
  RooDataSet *data;
  for (unsigned int i=0; i<dataForHybrid.size(); i++){
    RooDataSet *cutData = (RooDataSet*)dataForHybrid[i]->reduce(Name("hybridToy"),Title("hybridToy"),CutRange(cut_strings[i].c_str()));
    //RooDataSet *cutData = new RooDataSet("hybridToy","hybridToy",RooArgSet(*obs_var),Import(*dataForHybrid[i]),CutRange(cut_strings[i].c_str()));
    if (i==0) data=cutData;
    else data->append(*cutData);
  }
  return data;
}

void PdfModelBuilderFAN::throwHybridToy(string postfix, int nEvents, vector<float> switchOverMasses, vector<string> functions, bool bkgOnly, bool binned, bool poisson, bool cache){
  
  assert(switchOverMasses.size()==functions.size()-1);
  toyHybridData.clear();

  // have to throw unbinned for the hybrid
  throwToy(postfix,nEvents,bkgOnly,false,poisson,cache);

  vector<RooDataSet*> dataForHybrid;
  string hybridName = "hybrid";
  for (vector<string>::iterator func=functions.begin(); func!=functions.end(); func++){
    hybridName += "_"+*func;
    for (map<string,RooDataSet*>::iterator it=toyDataSet.begin(); it!=toyDataSet.end(); it++){
      if (it->first.find(*func)!=string::npos){
        dataForHybrid.push_back(it->second);
      }
    }
  }
  if (dataForHybrid.size()!=functions.size()){
    cerr << "One of the requested hybrid functions has not been found" << endl;
    exit(1);
  }

  RooDataSet *hybridData = makeHybridDataset(switchOverMasses,dataForHybrid);
  toyHybridData.clear();
  if (binned) {
    RooDataHist *hybridDataHist = hybridData->binnedClone();
    hybridDataHist->SetName(Form("%s_%s",hybridName.c_str(),postfix.c_str()));
    toyHybridData.insert(pair<string,RooAbsData*>(hybridDataHist->GetName(),hybridDataHist));
  }
  else {
    hybridData->SetName(Form("%s_%s",hybridName.c_str(),postfix.c_str()));
    toyHybridData.insert(pair<string,RooAbsData*>(hybridData->GetName(),hybridData));
  }
}

void PdfModelBuilderFAN::throwToy(string postfix, int nEvents, bool bkgOnly, bool binned, bool poisson, bool cache){

  toyData.clear();
  toyDataSet.clear();
  toyDataHist.clear();
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
    if (!bkgHasFit) cerr << "WARNING -- bkg has not been fit to data. Are you sure this is wise?" << endl; 
  }
  else {
    pdfSet = sbPdfs;
    if (!sbHasFit) cerr << "WARNING -- sb has not been fit to data. Are you sure this is wise?" << endl;
  }
  
  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    if (cache) {
      wsCache->loadSnapshot(it->first.c_str());
      cout << "Loaded snapshot, params at.." << endl;
      it->second->getVariables()->Print("v");
    }
    RooAbsData *toy;
    if (binned){
      RooDataHist *toyHist;
      if (poisson) toyHist = it->second->generateBinned(RooArgSet(*obs_var),nEvents,Extended(),Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      else toyHist = it->second->generateBinned(RooArgSet(*obs_var),nEvents,Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      toyDataHist.insert(pair<string,RooDataHist*>(toyHist->GetName(),toyHist));
      toy=toyHist;
    }
    else {
      RooDataSet *toySet;
      if (poisson) toySet = it->second->generate(RooArgSet(*obs_var),nEvents,Extended(),Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      else toySet = it->second->generate(RooArgSet(*obs_var),nEvents,Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      toyDataSet.insert(pair<string,RooDataSet*>(toySet->GetName(),toySet));
      toy=toySet;
    }
    toyData.insert(pair<string,RooAbsData*>(toy->GetName(),toy));
  }
  
}




void PdfModelBuilderFAN::throwToyFAN(string postfix, int nEvents, bool bkgOnly, bool binned, bool poisson, bool cache, int bins){

  toyData.clear();
  toyDataSet.clear();
  toyDataHist.clear();
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
    if (!bkgHasFit) cerr << "WARNING -- bkg has not been fit to data. Are you sure this is wise?" << endl;
  }
  else {
    pdfSet = sbPdfs;
    if (!sbHasFit) cerr << "WARNING -- sb has not been fit to data. Are you sure this is wise?" << endl;
  }

  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    if (cache) {
      wsCache->loadSnapshot(it->first.c_str());
      cout << "Loaded snapshot, params at.." << endl;
      it->second->getVariables()->Print("v");
    }
    RooAbsData *toy;
    if (binned){
      RooDataHist *toyHist;
      if (poisson) toyHist = it->second->generateBinned(RooArgSet(*obs_var),nEvents,Extended(),Name(Form("%s_%s",it->first.c_str(),postfix.c_str())),Binning(bins));
      else toyHist = it->second->generateBinned(RooArgSet(*obs_var),nEvents,Name(Form("%s_%s",it->first.c_str(),postfix.c_str())),Binning(bins));
      toyDataHist.insert(pair<string,RooDataHist*>(toyHist->GetName(),toyHist));
      toy=toyHist;
    }
    else {
      RooDataSet *toySet;
      if (poisson) toySet = it->second->generate(RooArgSet(*obs_var),nEvents,Extended(),Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      else toySet = it->second->generate(RooArgSet(*obs_var),nEvents,Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      toyDataSet.insert(pair<string,RooDataSet*>(toySet->GetName(),toySet));
      toy=toySet;
    }
    toyData.insert(pair<string,RooAbsData*>(toy->GetName(),toy));
  }

}



map<string,RooAbsData*> PdfModelBuilderFAN::getToyData(){
  return toyData;
}

map<string,RooAbsData*> PdfModelBuilderFAN::getHybridToyData(){
  return toyHybridData;
}

void PdfModelBuilderFAN::plotHybridToy(string prefix, int binning, vector<float> switchOverMasses, vector<string> functions, bool bkgOnly){

  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
  }
  else {
    pdfSet = sbPdfs;
  }
  
  int tempColors[4] = {kBlue,kRed,kGreen+2,kMagenta};

  vector<string> cut_strings;
  cut_strings.push_back("cutstring0");
  obs_var->setRange("cutstring0",obs_var->getMin(),switchOverMasses[0]);
  for (unsigned int i=1; i<switchOverMasses.size(); i++){
    cut_strings.push_back(Form("cutstring%d",i));
    obs_var->setRange(Form("cutstring%d",i),switchOverMasses[i-1],switchOverMasses[i]);
  }
  cut_strings.push_back(Form("cutstring%d",int(switchOverMasses.size())));
  obs_var->setRange(Form("cutstring%d",int(switchOverMasses.size())),switchOverMasses[switchOverMasses.size()-1],obs_var->getMax());
  
  RooPlot *plot = obs_var->frame();
  TCanvas *canv = new TCanvas();
  int i=0;
  for (vector<string>::iterator func=functions.begin(); func!=functions.end(); func++){
    for (map<string,RooAbsPdf*>::iterator pdfIt = pdfSet.begin(); pdfIt != pdfSet.end(); pdfIt++){
      // check if in list of hybrid functions
      if (pdfIt->first.find(*func)!=string::npos) {
        for (map<string,RooAbsData*>::iterator toyIt = toyData.begin(); toyIt != toyData.end(); toyIt++){
          //cout << "pdf: " << pdfIt->first << " - toy: " << toyIt->first << endl; 
          if (toyIt->first.find(pdfIt->first)!=string::npos){
            RooAbsData *data = toyIt->second->reduce(CutRange(cut_strings[i].c_str()));
            data->plotOn(plot,Binning(binning),MarkerColor(tempColors[i]),LineColor(tempColors[i]),CutRange(cut_strings[i].c_str()));
            pdfIt->second->plotOn(plot,LineColor(tempColors[i]),Range(cut_strings[i].c_str()));
            i++;
          }
        }
      }
    }
  }
  for (map<string,RooAbsData*>::iterator hybrid=toyHybridData.begin(); hybrid!=toyHybridData.end(); hybrid++){
    hybrid->second->plotOn(plot,Binning(binning),MarkerSize(0.8),MarkerStyle(kFullSquare));
    plot->SetMinimum(0.0001);
    plot->Draw();
    canv->Print(Form("%s_%s.pdf",prefix.c_str(),hybrid->first.c_str()));
  }
  delete canv;
}

void PdfModelBuilderFAN::plotToysWithPdfs(string prefix, int binning, bool bkgOnly){
  
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
  }
  else {
    pdfSet = sbPdfs;
  }
  TCanvas *canv = new TCanvas();
  for (map<string,RooAbsPdf*>::iterator pdfIt = pdfSet.begin(); pdfIt != pdfSet.end(); pdfIt++){
    for (map<string,RooAbsData*>::iterator toyIt = toyData.begin(); toyIt != toyData.end(); toyIt++){
      //cout << "pdf: " << pdfIt->first << " - toy: " << toyIt->first << endl; 
      if (toyIt->first.find(pdfIt->first)!=string::npos){
        RooPlot *plot = obs_var->frame();
        toyIt->second->plotOn(plot,Binning(binning));
        pdfIt->second->plotOn(plot,LineColor(kRed));
        pdfIt->second->paramOn(plot,LineColor(kRed),RooFit::Layout(0.34,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
        plot->Draw();
        //canv->Print(Form("%s_%s.pdf",prefix.c_str(),pdfIt->first.c_str()));  
        canv->Print(Form("%s_%s.png",prefix.c_str(),pdfIt->first.c_str()));
      }
    }
  }
  delete canv;

}

void PdfModelBuilderFAN::saveWorkspace(string filename){
  
  TFile *outFile = new TFile(filename.c_str(),"RECREATE");
  outFile->cd();
  wsCache->Write();
  outFile->Close();
  delete outFile;
}

void PdfModelBuilderFAN::saveWorkspace(TFile *file){
  file->cd();
  wsCache->Write();
}
