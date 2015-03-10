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

	cout << endl << "coucou" << endl;	

	string output = "plotsBias";
	string filebkg = "CMS-HGG_massfit_80_180_mva_haspixseed.root"; 
	string filesig = "CMS-HGG_massfit_80_180_mva_haspixseed.root";
	string datafile = "config.dat";
	int cat = 0;
	double expectSignal = 0.0;
	int expectSignalMass = 125;
	double muLow = -3.0;
	double muHigh = 3.0;
	double muStep = 0.01;
	bool blind = false;
	int ntoys = 10;
	int jobID = 0;

	cout << endl << "Coucou A" << endl;

	if( argc > 1 ) output = argv[1];
	if( argc > 2 ) filebkg = argv[2];
	if( argc > 3 ) filesig = argv[3];
	if( argc > 4 ) datafile = argv[4];
	if( argc > 5 ) 
        {
                std::stringstream ss ( argv[5] );
                ss >> cat;
        }
	if( argc > 6 )
        {
                std::stringstream ss ( argv[6] );
                ss >> expectSignal;
        }
	if( argc > 7 )
        {
                std::stringstream ss ( argv[7] );
                ss >> expectSignalMass;
        }
	if( argc > 8 )
        {
                std::stringstream ss ( argv[8] );
                ss >> muLow;
        }
	if( argc > 9 )
        {
                std::stringstream ss ( argv[9] );
                ss >> muHigh;
        }
	if( argc > 10 )
        {
                std::stringstream ss ( argv[10] );
                ss >> muStep;
        }
	if( argc > 11 )
        {
                std::stringstream ss ( argv[11] );
                ss >> blind;
        }
	if( argc > 12 )
        {
                std::stringstream ss ( argv[12] );
                ss >> ntoys;
        }
	if( argc > 13 )
        {
                std::stringstream ss ( argv[13] );
                ss >> jobID;
        }

	
	srand(time(0));

	cout << endl << "Coucou B" << endl;

	if(expectSignalMass!=90 && expectSignalMass!=95 && expectSignalMass!=100 && expectSignalMass!=105 && expectSignalMass!=110 && expectSignalMass!=115 && expectSignalMass!=120 && expectSignalMass!=125 && expectSignalMass!=130 && expectSignalMass!=135 && expectSignalMass!=140 && expectSignalMass!=145 && expectSignalMass!=150)
	{
    		cerr << "ERROR - expectSignalMass has to be integer in range (90,150,5)" << endl;
    		exit(1);
  	}

	system(Form("mkdir -p %s",output.c_str()));

	int ncats = 9;
	int binning = 105;

	// --- style of plots --- //
        gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
        //gStyle->SetOptFit(1111);  //FIXME
	//gStyle->SetOptStat(1111);
	//TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);

	cout << endl << "Coucou C" << endl;


	cout << endl << "output = " <<output << ", filebkg = " <<filebkg <<", filesig = " <<filesig <<", datafile = " <<datafile << endl;


	TFile *fb = TFile::Open(filebkg.c_str());
	TFile *fs = TFile::Open(filesig.c_str());

	cout << endl << "Coucou D" << endl;

	RooWorkspace *w_bkg  = (RooWorkspace*) fb->Get("cms_hgg_workspace");
	RooWorkspace *w_sig  = (RooWorkspace*) fs->Get("cms_hgg_workspace");

	cout << endl << "Coucou E" << endl;

	RooRealVar *mass = (RooRealVar*) w_bkg->var("CMS_hgg_mass");
	RooRealVar *mu = new RooRealVar("mu","mu",0.0,muLow,muHigh); 
	RooRealVar *intL = (RooRealVar*) w_bkg->var("IntLumi");
	double lumi = intL->getVal()/1000.;

	cout << endl << "Coucou 0" << endl;

	RooDataSet *data = (RooDataSet*)w_bkg->data(Form("data_mass_cat%d",cat));
	RooDataHist *dataBinned = new RooDataHist(Form("roohist_data_mass_cat%d",cat),Form("roohist_data_mass_cat%d",cat),RooArgSet(*mass),*data);

	RooDataSet *sigMC = (RooDataSet*)w_sig->data(Form("sig_ggh_mass_m%d_cat%d",expectSignalMass,cat));
  	RooDataSet *sigMC_vbf = (RooDataSet*)w_sig->data(Form("sig_wzh_mass_m%d_cat%d",expectSignalMass,cat));
  	RooDataSet *sigMC_wzh = (RooDataSet*)w_sig->data(Form("sig_vbf_mass_m%d_cat%d",expectSignalMass,cat));
  	RooDataSet *sigMC_tth = (RooDataSet*)w_sig->data(Form("sig_tth_mass_m%d_cat%d",expectSignalMass,cat));
  	sigMC->append(*sigMC_vbf);
  	sigMC->append(*sigMC_wzh);
  	sigMC->append(*sigMC_tth);

	if (!dataBinned || !sigMC)
	{
    		cerr << "ERROR -- one of data or signal is NULL" << endl;
		cerr << "You probably have to change the name of the workspace..." << endl;
    		exit(1);
  	}

	cout << endl << "Coucou 1" << endl;

	// ----- Signal Pdf ----- //
	RooDataHist *sigMCBinned = new RooDataHist(Form("roohist_%s",sigMC->GetName()),Form("roohist_%s",sigMC->GetName()),RooArgSet(*mass),*sigMC);
  	RooAbsPdf * sigPdf = new RooHistPdf(Form("pdf_%s",sigMC->GetName()),Form("pdf_%s",sigMC->GetName()),RooArgSet(*mass),*sigMCBinned);
  	RooAbsReal * sigNorm = new RooConstVar(Form("sig_events_%s",sigMC->GetName()),Form("sig_events_%s",sigMC->GetName()),sigMC->sumEntries());
		
	// ----- Signal + Background Pdf : Truth ----- //

	RooAbsReal * sigYield = new RooProduct("sig_yield","sig_yield",RooArgSet(*mu,*sigNorm));
	RooRealVar * bkgYield = new RooRealVar("bkg_yield","bkg_yield",1000.,0.,1.e6);

	cout << endl << "Coucou 2" << endl;
	
	RooRealVar * fraction = NULL;

	map<string,RooAbsPdf*> bkgPdfs = readDatafile(datafile, mass, cat, true,fraction); 
	//map<string,RooAbsPdf*> testPdfs = readDatafile(datafile, mass, cat, false); //FIXME
	
	map<string,RooAbsPdf*> sbPdfs;

	cout << endl << "Truth functions : " << endl;
	
	for(map<string,RooAbsPdf*>::iterator bkg=bkgPdfs.begin(); bkg!=bkgPdfs.end(); bkg++)
	{
		cout << endl <<  bkg->first << endl;
    		RooAbsPdf *sbMod = new RooAddPdf(Form("sb_%s",bkg->first.c_str()),Form("sb_%s",bkg->first.c_str()),RooArgList(*(bkg->second),*sigPdf),RooArgList(*bkgYield,*sigYield));
	
		sbPdfs.insert(pair<string,RooAbsPdf*>(sbMod->GetName(),sbMod));
	}

	
	// ----- Signal + Background Pdf : Test ----- //

        map<string,RooAbsPdf*> testPdfs = readDatafile(datafile, mass, cat, false,fraction); //FIXME
     
        map<string,RooAbsPdf*> test_sbPdfs;

	cout << endl << "Test functions : " << endl;

        for(map<string,RooAbsPdf*>::iterator test=testPdfs.begin(); test!=testPdfs.end(); test++)
        {
                cout << endl <<  test->first << endl;
                RooAbsPdf *sbMod = new RooAddPdf(Form("sb_%s",test->first.c_str()),Form("sb_%s",test->first.c_str()),RooArgList(*(test->second),*sigPdf),RooArgList(*bkgYield,*sigYield));
     
                test_sbPdfs.insert(pair<string,RooAbsPdf*>(sbMod->GetName(),sbMod));
        }
		

	// ----- Fitting Truth to data ----- //

	mu->setVal(expectSignal);
	mu->setConstant(true);


	TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);

	for (map<string,RooAbsPdf*>::iterator it=sbPdfs.begin(); it!=sbPdfs.end(); it++)
	{
		RooPlot * frame = mass->frame();
    		RooFitResult *fit = (RooFitResult*)it->second->fitTo(*data,Save(true));
		cout << endl << "Fit Res Before: " << endl;
      		fit->floatParsInit().Print("v");
      		cout << endl << "Fit Res After: " << endl;
      		fit->floatParsFinal().Print("v");

		data->plotOn(frame,Binning(binning));
		it->second->plotOn(frame);
		frame->Draw();
		system(Form("mkdir -p %s/cat%d/mass_%d/",output.c_str(),cat,expectSignalMass));
		system(Form("mkdir -p %s/cat%d/mass_%d/truthToData/",output.c_str(),cat,expectSignalMass));
		system(Form("mkdir -p %s/cat%d/mass_%d/toyHisto/",output.c_str(),cat,expectSignalMass));
		system(Form("mkdir -p %s/cat%d/mass_%d/toyHistoPlusTest/",output.c_str(),cat,expectSignalMass));
		c1->Print(Form("%s/cat%d/mass_%d/truthToData/%s_truthToData_jobID%d.pdf",output.c_str(),cat,expectSignalMass,it->first.c_str(),jobID));
		//frame->Clear();
		c1->Clear();
		frame->Delete();

	}

	c1->Clear();

	// ----- Generation of toys ----- //

	cout << endl << "------------------------------";
	cout << endl << "----- Generation of toys -----";
	cout << endl << "------------------------------" << endl; 

	mu->setVal(expectSignal);
	mu->setConstant(true);

	int iseed = 0;

	for (int toy=0; toy<ntoys; toy++)
	{

		mu->setVal(expectSignal);
        	mu->setConstant(true);
		
		map<string,RooAbsData*> toys;

		for (map<string,RooAbsPdf*>::iterator it=sbPdfs.begin(); it!=sbPdfs.end(); it++)
		{
			cout << endl <<  it->first << endl;

			iseed=rand();
			RooRandom::randomGenerator()->SetSeed(iseed) ; 
			//RooDataHist *toyHist = it->second->generateBinned(RooArgSet(*mass),dataBinned->sumEntries(),Extended(),Name(Form("%s_%d",it->first.c_str(),toy)));
			RooDataHist *toyHist = it->second->generateBinned(RooArgSet(*mass),dataBinned->sumEntries(),Name(Form("%s_%d",it->first.c_str(),toy)));
			toys.insert(pair<string,RooDataHist*>(toyHist->GetName(),toyHist));	

			cout << endl << toyHist->GetName() << endl;
	
			RooPlot * frame = mass->frame();
			toyHist->plotOn(frame,Binning(binning));
			it->second->plotOn(frame,LineColor(kBlue),LineStyle(kDashed));
			frame->Draw();
			c1->Print(Form("%s/cat%d/mass_%d/toyHisto/%s_toyHisto%d_jobID%d.pdf",output.c_str(),cat,expectSignalMass,it->first.c_str(),toy,jobID));
			c1->Clear();	


		}

		mu->setConstant(false);
		
		for (map<string,RooAbsData*>::iterator it=toys.begin(); it!=toys.end(); it++)
		{

			double globalMinNLL=1.e6;
			double bestFitValError = 1.0;
			double bestFitVal = 0.0;
			double bias = 0.0;
			for (map<string,RooAbsPdf*>::iterator m=test_sbPdfs.begin(); m!=test_sbPdfs.end(); m++) 
			{
				RooAbsPdf *pdf = m->second;
				RooFitResult *nom = pdf->fitTo(*(it->second),PrintLevel(-1),PrintEvalErrors(-1),Warnings(false),Save(true));

				RooArgSet *params = pdf->getParameters(*mass);
    				/*if(2*nom->minNll()<globalMinNLL)
				{
      					globalMinNLL = 2*nom->minNll();
      					bestFitVal = mu->getVal();
    				}
				*/
				
				globalMinNLL = 2*nom->minNll();
				bestFitVal = mu->getVal();
				bestFitValError = mu->getError();

				cout << endl << "bestFitVal = " << bestFitVal << endl;
				cout << endl << "bestFitValError = " << bestFitValError << endl;
				cout << endl << "globalMinNLL = " << globalMinNLL << endl;
				
				bias = (expectSignal - bestFitVal) / bestFitValError;

				cout << endl << "bias" << endl;

				string correspondingTruth = "";

				for (map<string,RooAbsPdf*>::iterator it2=sbPdfs.begin(); it2!=sbPdfs.end(); it2++)
				{
					if(it->first.find(it2->first)) 
					{
						correspondingTruth = it2->first;
						break;
					}
				}


				ofstream fileOut(Form("%s/cat%d/mass_%d/Truth_%s_Test_%s.txt",output.c_str(),cat,expectSignalMass,correspondingTruth.c_str(),m->first.c_str()),ios::app);
				//fileOut.open(Form("/afs/cern.ch/work/l/lsgandur/CMSSW_5_3_9_patch1_h2gglobe_60_180/src/h2gglobe/BiasStudy_15_11_13/%s/Cat%d_Truth_%s_Test_%s_mass_%d.txt",output.c_str(),cat,it->first.c_str(),m->first.c_str(),expectSignalMass),ios::out | ios::app);

				//fileOut.open(Form("%s/cat%d/mass_%d/Truth_%s_Test_%s.txt",output.c_str(),cat,expectSignalMass,it->first.c_str(),m->first.c_str()),ios::out | ios::app); //FIXME

				cout << endl << Form("%s/cat%d/mass_%d/Truth_%s_Test_%s.txt",output.c_str(),cat,expectSignalMass,it->first.c_str(),m->first.c_str()) <<endl;


				//fileOut << bestFitVal << endl;
				fileOut << bias << endl;			


				fileOut.close();
	
				RooPlot * frame = mass->frame();
				it->second->plotOn(frame);
				m->second->plotOn(frame);
				frame->Draw();
				c1->Print(Form("%s/cat%d/mass_%d/toyHistoPlusTest//%s_toyHisto%d_%s_jobID%d.pdf",output.c_str(),cat,expectSignalMass,it->first.c_str(),toy,m->first.c_str(),jobID));
				c1->Clear();

			}

		}


	}

	return 0;
}







