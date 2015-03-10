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
	string eVeto = "defaultE";

	if( argc > 1 ) output = argv[1];
	if( argc > 2 ) filebkg = argv[2];
	if( argc > 3 ) filesig = argv[3];
	if( argc > 4 ) blind = argv[4];
	if( argc > 5 ) eVeto = argv[5];

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
	
	
	int nbin = 240; //FIXME
	int nFitParam = 0;


	system(Form("rm -rf %s",output.c_str()));	

	system(Form("mkdir -p %s",output.c_str()));
	system(Form("mkdir -p %s/Trial/",output.c_str()));
	system(Form("mkdir -p %s/Chosen/",output.c_str()));

	vector <string> functions;

	functions.push_back("Voigtian");
	//functions.push_back("CB");
	//functions.push_back("CBgaus");	
	//functions.push_back("CBpCB");
	//functions.push_back("CBconvCB");

	int order = 0;


	for(vector<string>::iterator iter = functions.begin(); iter != functions.end(); iter++)
	{
		cout << endl << "// ----- Function : " << *iter <<" ----- //"<< endl;		

		for (int cat=0; cat<ncats; cat++)
		{
			cout << "--- Category " << cat << " ---" << endl;
			//while ( pvalueGF < 0.01)
			//{
			
				RooAbsData *data;
				data = (RooDataHist*)w_bkg->data(Form("roohist_bkg_mass_cat%d",cat));

				data->Print();
	
				bool truth = true;	
				RooRealVar * fraction = NULL;
			
				pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > pdf = getPdf((*iter).c_str(),order,x,fraction,cat,true,eVeto); 
				RooAbsPdf *bkgPdf = pdf.first.first;
				//RooAbsPdf *bkgPdf = getPdf((*iter).c_str(),order,x,fraction,cat,true).first.first; 
				string* paramNames = pdf.second.first;	
				
				RooFitResult *fitRes = NULL;
				if(cat == 4 || cat == 3 || cat == 2) fitRes = bkgPdf->fitTo(*data,Save(true),Range(86,93)); //FIXME
				else fitRes = bkgPdf->fitTo(*data,Save(true),Range(85,95));
	          		thisNll = fitRes->minNll();
			
				cout << endl << "thisNll = " << thisNll << endl;
				
				RooPlot *frame = x->frame();
				data->plotOn(frame,RooFit::Binning(nbin),Name("myhist"),Range(85,95)); //FIXME
				if(cat == 4 || cat == 3 || cat == 2) bkgPdf->plotOn(frame,RooFit::LineColor(kRed),Name("mycurve"),Range(86,93)); //FIXME
				else bkgPdf->plotOn(frame,RooFit::LineColor(kRed),Name("mycurve"),Range(85,95)); //FIXME	
				frame->Draw();

				double chi2GF = chiSquare(frame,(char *)"mycurve",(char *)"myhist");

				nFitParam = pdf.second.second[0] - 1;
				//pvalueGF = TMath::Prob(chi2GF, nbin - nFitParam);
				pvalueGF = TMath::Prob(chi2GF, 40 - nFitParam); //FIXME !

				int np = bkgPdf->getParameters(*data)->getSize();

				RooArgSet *params = bkgPdf->getParameters(*x);  
     
                                int loopEnd = pdf.second.second[0];
                                cout << endl << "loopEnd = " <<loopEnd << endl;
     

				TLatex latexLabel;
				latexLabel.SetNDC();
				latexLabel.SetTextSize(0.030);
				latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2012, #sqrt{s} = 8 TeV");
				latexLabel.DrawLatex(0.17, 0.88, Form("Data, #int L = %2.1f fb^{-1}",lumi));
				latexLabel.DrawLatex(0.60,0.88,Form("NLL_{N} = %.2f",thisNll));
				latexLabel.DrawLatex(0.60,0.83,Form("#chi^{2}_{GF} = %.2f, Prob_{GF} = %.2f ",chi2GF / (nbin - nFitParam) ,pvalueGF));
				//latexLabel.DrawLatex(0.60,0.83,Form("#chi^{2}_{GF} = %.2f, Prob_{GF} = %.2f ",));			
	
				fstream summaryFile_params;
				summaryFile_params.open(Form("%s/Trial/%s_cat%d.txt",output.c_str(),(*iter).c_str(),cat),ios::out);	
				double temp_leg = 0.83;
				for(int n = 1; n < loopEnd; n++) //FIXME
                                {
                                        temp_leg -= 0.05;
                                        cout << endl << "temp_leg = " <<temp_leg << endl;
					RooRealVar * var_temp = (RooRealVar*) params->find(paramNames[n].c_str());
					cout << endl << "var_temp->getVal() = "<<var_temp->getVal() <<", var_temp->getError() = " <<var_temp->getError() << endl;				
	
                                        latexLabel.DrawLatex(0.60,temp_leg,Form("par_{%d} = %.2f #pm %.2f",n,var_temp->getVal(),var_temp->getError()));
        				summaryFile_params << Form("par_{%d} = %.2f #pm %.2f",n,var_temp->getVal(),var_temp->getError()) << endl;                       
	
				}
				summaryFile_params.close();

				frame->GetXaxis()->SetLimits(80,100);
				//frame->GetXaxis()->SetLimits(60,180);

				c1->Print(Form("%s/Trial/%s_cat%d.png",output.c_str(),(*iter).c_str(),cat));
				c1->Print(Form("%s/Trial/%s_cat%d.C",output.c_str(),(*iter).c_str(),cat));	
	
				if(!(pvalueGF < 0.01) ) 
				{
					c1->Print(Form("%s/Chosen/Chosen_%s_cat%d.png",output.c_str(),(*iter).c_str(),cat)); 
					fstream summaryFile;
					summaryFile.open(Form("%s/Chosen/Chosen_cat%d.txt",output.c_str(),cat),ios::out | ios::app);
					summaryFile << Form("%.2f & %.2f & %.2f",thisNll,chi2GF / (nbin - nFitParam),pvalueGF);
					summaryFile << endl;
					summaryFile.close();  

				}
				

				cout << endl << "coucou 1 " << endl;
				bkgPdf->Delete();
				bkgPdf = 0;
				cout << endl << "coucou 2 " << endl;
				//data->Delete();
				//data = 0;
				cout << endl << "coucou 3 " << endl;
				c1->Clear();
			//}
			
			prob = 0.0;
			pvalueGF = 0.0;
		}
	}

	cout << endl << "coucou 4 " << endl;
	x->Delete();
        x = 0;
	cout << endl << "coucou 5 " << endl;
        intL->Delete();
        intL = 0;
        //w_bkg->Delete();
        //w_bkg = 0;
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
