//Continue to correct chi2 -> if less 5 entries/bin -> we don't use this bin to compute chi2/pval.
//for low stat categories, don't use the chi2 criteria !!

#include "functions.h"
#include "setTDRStyle.C"


int main(int argc, char *argv[])
{
	for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 ) 
        {
                cerr << "arguments should be passed : " <<endl; 
                return 1;
        }
	

	string output = "Truthfunctions";
	string filebkg = "CMS-HGG_massfit_mva_80_180.root"; 
	string filesig = "CMS-HGG_massfit_mva_80_180.root";
	bool blind = false;
	int nbin = 80;
	string eVeto = "defaultE";

	if( argc > 1 ) output = argv[1];
	if( argc > 2 ) filebkg = argv[2];
	if( argc > 3 ) filesig = argv[3];
	if( argc > 4 ) 
        {
                std::stringstream ss ( argv[4] );
                ss >> blind;
        }   
	if( argc > 5 ) 
	{   
                std::stringstream ss ( argv[5] );
                ss >> nbin;
        }	
	if( argc > 6 ) eVeto = argv[6];

	int ncats = 14;

	// --- style of plots --- //
        gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
        //gStyle->SetOptFit(1111);  //FIXME
	//gStyle->SetOptStat(1111);
	TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);

	TFile *fb = TFile::Open(filebkg.c_str());
	TFile *fs = TFile::Open(filesig.c_str());

	RooWorkspace *w_bkg  = (RooWorkspace*) fb->Get("cms_hgg_workspace");

	RooRealVar *x = (RooRealVar*) w_bkg->var("CMS_hgg_mass");
	RooRealVar *intL = (RooRealVar*) w_bkg->var("IntLumi");
	double lumi = intL->getVal()/1000.;

	double thisNll=0.0;
      	double prevNll=0.0;
      	double chi2=0.0;
      	double prob=0.0;
	double pvalueGF=0.0;
	double fractionValue=0.0;
	double fractionValueError=0.0;
	
	
	//int nbin = 100;
	int nFitParam = 0;


	//system(Form("rm -rf %s",output.c_str()));	

	system(Form("mkdir -p %s",output.c_str()));

	vector <string> functions;


	functions.push_back("BernsteinPlusVoigtian");	
	functions.push_back("ChebychevPlusVoigtian");
	//functions.push_back("ExponentialPlusVoigtian");
	//functions.push_back("LaurentPlusVoigtian");
	//functions.push_back("PowerLawPlusVoigtian");

/*
	functions.push_back("ExponentialPlusCB");
        functions.push_back("BernsteinPlusCB");
        functions.push_back("ChebychevPlusCB");
        functions.push_back("PowerLawPlusCB");
        functions.push_back("LaurentPlusCB");	
	

	functions.push_back("Bernstein");
        functions.push_back("Chebychev");		
	functions.push_back("Exponential");
        functions.push_back("Laurent");	
	functions.push_back("PowerLaw");
*/

	//functions.push_back("PolynomialsPlusVoigtian");

	//functions.push_back("ExponentialPlusCBgaus");
        //functions.push_back("BernsteinPlusCBgaus");
        //functions.push_back("ChebychevPlusCBgaus");
        //functions.push_back("PowerLawPlusCBgaus");
        //functions.push_back("LaurentPlusCBgaus");

	for(vector<string>::iterator iter = functions.begin(); iter != functions.end(); iter++)
	{

		int order = 4;
        	int prev_order=0;

		cout << endl << "// ----- Function : " << *iter <<" ----- //"<< endl;		

		for (int cat=0; cat<ncats; cat++)
		{
			cout << "--- Category " << cat << " ---" << endl;
			
				cout << "*** Order " << order << " ***" << endl;
				
				RooAbsData *data;
				data = (RooDataHist*)w_bkg->data(Form("roohist_data_mass_cat%d",cat));
				data->Print();
				
				bool truth = true;	
				RooRealVar * fraction = NULL;
				pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > pdf = getPdf((*iter).c_str(),order,x,fraction,cat,true,eVeto); 
                                RooAbsPdf *bkgPdf = pdf.first.first;
                                //RooAbsPdf *bkgPdf = getPdf((*iter).c_str(),order,x,fraction,cat,true).first.first; 
                                string* paramNames = pdf.second.first;

				//RooAbsPdf *bkgPdf = getPdf((*iter).c_str(),order,x,fraction,cat,true).first.first; 
				//RooRealVar * fraction = getPdf((*iter).c_str(),order,x).second;
	
	
				RooFitResult *fitRes = NULL;
				//if(cat == 3 || cat == 2)fitRes = bkgPdf->fitTo(*data,Save(true));
				//else fitRes = bkgPdf->fitTo(*data,Save(true));
				RooRealVar temp_var("pdf_temp_var_norm","temp_var",100.,0.0,1e6);
				RooExtendPdf temp_extend("temp_extend","temp_extend",*bkgPdf,temp_var);
				//fitRes = bkgPdf->fitTo(*data,RooFit::Save(true),RooFit::Extended(true),RooFit::Strategy(1));		
				fitRes = temp_extend.fitTo(*data,RooFit::Save(true),RooFit::Extended(true),RooFit::Strategy(1));	

				
				//RooRealVar* frac_temp = (RooRealVar*) bkgPdf->getVariables()->FindObject("fraction"); //FIXME
				//fractionValue = frac_temp->getVal();
				fractionValue = fraction->getVal();
				fractionValueError = fraction->getError();

				cout << endl << "fractionValue =  "<<fractionValue<< endl;
				RooPlot *frame = x->frame();
				data->plotOn(frame,RooFit::Binning(nbin),Name("myhist"));
				bkgPdf->plotOn(frame,RooFit::LineColor(kRed),Name("mycurve"));
				frame->Draw();

	
				double chi2GF = chiSquare(frame,(char *)"mycurve",(char *)"myhist");

				nFitParam = order;

				int np = bkgPdf->getParameters(*data)->getSize();
	

				TLatex latexLabel;
				latexLabel.SetNDC();
				latexLabel.SetTextSize(0.030);
				latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2012, #sqrt{s} = 8 TeV");
				latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %2.1f fb^{-1}",lumi));
				//latexLabel.DrawLatex(0.60,0.88,Form("#chi^{2} = %.2f, Prob = %.2f ",chi2*(105-np),prob));
				latexLabel.DrawLatex(0.60,0.83,Form("fraction = %.2f #pm %.2f",fractionValue,fractionValueError));	
	
			
				RooArgSet *params = bkgPdf->getParameters(*x);			
	
				int loopEnd = pdf.second.second[0];
                                cout << endl << "loopEnd = " <<loopEnd << endl;

                                fstream summaryFile_params;
                                summaryFile_params.open(Form("%s/%s_cat%d_order%d.txt",output.c_str(),(*iter).c_str(),cat,order),ios::out);    
                                //double temp_leg = 0.83;
                                for(int n = 1; n < loopEnd; n++) //FIXME
                                {
                                        //temp_leg -= 0.05;
                                        ///cout << endl << "temp_leg = " <<temp_leg << endl;
                                        RooRealVar * var_temp = (RooRealVar*) params->find(paramNames[n].c_str());
                                        //cout << endl << "var_temp->getVal() = "<<var_temp->getVal() <<", var_temp->getError() = " <<var_temp->getError() << endl;        
     
                                        //latexLabel.DrawLatex(0.60,temp_leg,Form("par_{%d} = %.2f #pm %.2f",n,var_temp->getVal(),var_temp->getError()));
                                        summaryFile_params << Form("par_{%d} = %.2f #pm %.2f",n,var_temp->getVal(),var_temp->getError()) << endl;
					if(n == loopEnd - 1) summaryFile_params << Form("fraction = %.2f #pm %.2f",fractionValue,fractionValueError);
                                }
                                summaryFile_params.close();





				c1->Print(Form("%s/%s_cat%d_order%d.png",output.c_str(),(*iter).c_str(),cat,order));
	
	
				cout << endl << "coucou 1 " << endl;
				bkgPdf->Delete();
				bkgPdf = 0;
				cout << endl << "coucou 2 " << endl;
				//data->Delete();
				//data = 0;
				cout << endl << "coucou 3 " << endl;
				c1->Clear();
		
				
		}
	}

	cout << endl << "coucou 4 " << endl;
	x->Delete();
        x = 0;
	cout << endl << "coucou 5 " << endl;
        intL->Delete();
        intL = 0;
        w_bkg->Delete();
        w_bkg = 0;
	cout << endl << "coucou 6 " << endl;
        fs->Delete();
        fs = 0;
	cout << endl << "coucou 7 " << endl;
        fb->Delete();   
        fb = 0;
	cout << endl << "coucou 8 " << endl;
	delete c1;
	c1 = 0;
	return 0;
}
