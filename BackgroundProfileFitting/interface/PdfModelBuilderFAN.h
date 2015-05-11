#ifndef PdfModelBuilderFAN_h 
#define PdfModelBuilderFAN_h

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooProduct.h"
#include "RooWorkspace.h"
#include "RooLognormal.h"  //FAN
#include "TFile.h"
#include "Math/SMatrix.h"  //FAN
#include "Math/SVector.h"  //FAN 



using namespace std;
using namespace RooFit;

class PdfModelBuilderFAN {
  
  public:
    PdfModelBuilderFAN();
    ~PdfModelBuilderFAN();

    void setObsVar(RooRealVar *var);
    void setSignalModifier(RooRealVar *var);
    void setSignalModifierVal(float val);
    void setSignalModifierConstant(bool val);
    void setKeysPdfAttributes(RooDataSet *data, double rho=2);

    void addBkgPdf(string type, int nParams, string name, bool cache=true);
    void addBkgPdfFAN(string type, int nParams, string name, float constant, float bernDownBound, float bernUpBound, bool cache=true);   //FAN

    void setSignalPdf(RooAbsPdf *pdf, RooRealVar *norm=NULL);
    void setSignalPdfFromMC(RooDataSet *data);
    void makeSBPdfs(bool cache=true);

    map<string,RooAbsPdf*> getBkgPdfs();
    map<string,RooAbsPdf*> getSBPdfs();
    RooAbsPdf *getSigPdf();

    void fitToData(RooAbsData *data, bool bkgOnly=true, bool cache=true, bool print=false);
    void fitToDataFAN(RooAbsData *data, bool bkgOnly=true, bool cache=true, bool print=false);  //FAN
    void plotPdfsToData(RooAbsData *data, int binning, string name, bool bkgOnly=true, string specificPdfName="");
    void plotToysWithPdfs(string prefix, int binning, bool bkgOnly=true);
    void plotHybridToy(string prefix, int binning, vector<float> switchOverMasses, vector<string> functions, bool bkgOnly=true);
    
    void setSeed(int seed);
    void throwToy(string name, int nEvents, bool bkgOnly=true, bool binned=true, bool poisson=true, bool cache=true);
    void throwToyFAN(string name, int nEvents, bool bkgOnly=true, bool binned=true, bool poisson=true, bool cache=true, int bins=160); //FAN
    RooDataSet *makeHybridDataset(vector<float> switchOverMasses, vector<RooDataSet*> dataForHybrid);
    void throwHybridToy(string name, int nEvents, vector<float> switchOverMasses, vector<string> functions, bool bkgOnly=true, bool binned=true, bool poisson=true, bool cache=true);
    map<string,RooAbsData*> getToyData();
    map<string,RooAbsData*> getHybridToyData();

    void saveWorkspace(TFile* file);
    void saveWorkspace(string filename);

    RooAbsPdf* getBernstein(string prefix, int order);
    RooAbsPdf* getBernsteinFAN(string prefix, int order, float bernDownBound, float bernUpBound);   //FAN
    RooAbsPdf* getChebychev(string prefix, int order);
    RooAbsPdf* getChebychevFAN(string prefix, int order, float bernDownBound, float bernUpBound);  //FAN
    RooAbsPdf* getPolynomial(string prefix, int order);  //FAN
    RooAbsPdf* getPowerLaw(string prefix, int order);
    RooAbsPdf* getPowerLawSingle(string prefix, int order);
    RooAbsPdf* getPowerLawGeneric(string prefix, int order);
    RooAbsPdf* getExponential(string prefix, int order);
    RooAbsPdf* getExponentialSingle(string prefix, int order);
    RooAbsPdf* getLaurentSeries(string prefix, int order);
    RooAbsPdf* getLaurentSeriesFAN(string prefix, int order, float constant);  //FAN
    RooAbsPdf* getKeysPdf(string prefix);
    RooAbsPdf* getPdfFromFile(string &prefix);
    RooAbsPdf* getConvDoubleCB(string prefix, RooRealVar *alphaCB1, RooRealVar *alphaCB2);    //FAN
    RooAbsPdf* getConvDoubleCBnew(string prefix,float alphacb1, float alphacb2);    //FAN
    RooAbsPdf* getConvDoubleCBnew_Float(string prefix,float alphacb1, float alphacb2);    //FAN
    RooAbsPdf* getSumOfVoigtians(string prefix, int order, bool recursive=false); //FAN
    pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getOfficialSumVoigtians(string type, string name, string prefix, int order1, int order2, float constant, float bernDownBound, float bernUpBound, bool recursive=false); //FAN
    pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getOfficialSumVoigtiansForZeeFit(string type, string name, string prefix, int order1, int order2, float constant, float bernDownBound, float bernUpBound, bool recursive=false); //FAN
    pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getOfficialSumVoigtiansForZeeFitDoubleCB(string type, string name, string prefix, int order1, int order2, float constant, float bernDownBound, float bernUpBound, float alphacb1, float alphacb2, bool recursive=false); //FAN
    pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getOfficialSumVoigtiansForZeeFitDoubleCB_Float(string type, string name, string prefix, int order1, int order2, float constant, float bernDownBound, float bernUpBound, float alphacb1, float alphacb2, bool recursive=false); //FAN
    pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getOfficialSumVoigtiansFix(string type, string name, string prefix, int order1, RooAbsPdf* pdfVoiFix, float constant, float bernDownBound, float bernUpBound, bool recursive=false); //FAN
    pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getOfficialSumVoigtiansFixNew(string type, string name, string prefix, int order1, float voiMean, float voiMeanErrorL, float voiMeanErrorH, float voiSigma, float voiSigmaErrorL, float voiSigmaErrorH, float voiWidth, float voiWidthErrorL, float voiWidthErrorH, float ErrorRange, float constant, bool constOrSo, float bernDownBound, float bernUpBound, bool recursive=false); //FAN
    pair<RooAbsPdf*,pair<RooAbsPdf*,RooAbsPdf*> >  getOfficialSumVoigtiansFixNewDoubleCB(string type, string name, string prefix, int order1, float voiMean, float voiMeanErrorL, float voiMeanErrorH, float voiSigma, float voiSigmaErrorL, float voiSigmaErrorH, float voinCB1, float voinCB1ErrorL, float voinCB1ErrorH, float voinCB2, float voinCB2ErrorL, float voinCB2ErrorH, float voialphaCB1, float voialphaCB2, float ErrorRange, float constant, bool constOrSo, float bernDownBound, float bernUpBound, bool recursive=true); //FAN
    void addBkgPdfSumVoigtians(string type, int nParams, string name, int nBre, string ext, float constant, float bernDownBound, float bernUpBound, bool cache=true, bool recursive=false); //FAN
    void addBkgPdfSumVoigtiansFix(string type, int nParams, string name, RooAbsPdf* pdfVoiFix, string ext, float constant, float bernDownBound, float bernUpBound, bool cache=true, bool recursive=false); //FAN
    void addBkgPdfSumVoigtiansFixNew(string type, int nParams, string name, float voiMean, float voiMeanErrorL, float voiMeanErrorH, float voiSigma, float voiSigmaErrorL, float voiSigmaErrorH, float voiWidth, float voiWidthErrorL, float voiWidthErrorH, float ErrorRange, string ext, float constant, bool constOrSo, float bernDownBound, float bernUpBound, bool cache=true, bool recursive=false); //FAN
    void addBkgPdfSumVoigtiansFixNewDoubleCB(string type, int nParams, string name, float voiMean, float voiMeanErrorL, float voiMeanErrorH, float voiSigma, float voiSigmaErrorL, float voiSigmaErrorH, float voinCB1, float voinCB1ErrorL, float voinCB1ErrorH, float voinCB2, float voinCB2ErrorL, float voinCB2ErrorH, float voialphaCB1, float voialphaCB2, float ErrorRange, string ext, float constant, bool constOrSo, float bernDownBound, float bernUpBound, bool cache=true, bool recursive=false); //FAN
    



  private:
   
    vector<string> recognisedPdfTypes;

    map<string,RooAbsPdf*> bkgPdfs;
    map<string,RooAbsPdf*> sbPdfs;
    RooAbsPdf* sigPdf;
    RooAbsReal* sigNorm;
    RooRealVar *bkgYield;
    RooAbsReal *sigYield;
    RooDataSet *keysPdfData;
    double keysPdfRho;

    map<string,RooAbsData*> toyData;
    map<string,RooAbsData*> toyHybridData;
    map<string,RooDataSet*> toyDataSet;
    map<string,RooDataHist*> toyDataHist;

    map<string,RooRealVar*> params;
    map<string,RooFormulaVar*> prods;
    map<string,RooAbsPdf*> utilities;

    RooRealVar *obs_var;
    bool obs_var_set;
    RooRealVar *signalModifier;
    bool signal_modifier_set;
    bool signal_set;
    bool bkgHasFit;
    bool sbHasFit;
    bool keysPdfAttributesSet;
    vector<string> cut_strings;

    RooWorkspace *wsCache;

    int verbosity;
    bool forceFracUnity_; //FAN


};
#endif
