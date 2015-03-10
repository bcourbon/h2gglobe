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
	string filebkg = "CMS-HGG_massfacmva_legacy_notHasPixelSeed_80_120_PTl036_PTs03_DynPtCut45.root"; 
	string filesig = "CMS-HGG_massfacmva_legacy_notHasPixelSeed_80_120_PTl036_PTs03_DynPtCut45.root";
	string datafile = "config_massfacmva_legacy_notHasPixelSeed_80_120_PTl036_PTs03_DynPtCut45_noChi2GF_v3.dat";
	int cat = 11;
	double expectSignal = 0.0;
	int expectSignalMass = 110;
	double muLow = -1000.0;
	double muHigh = 1000.0;
	double muStep = 0.01;
	bool blind = false;
	int ntoys = 1000;
	int jobID = 0;
	int binning = 40;
	string eVeto = "defaultE";
	
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
	if( argc > 14 )
        {
                std::stringstream ss ( argv[14] );
                ss >> binning;
        }
	if( argc > 15 ) eVeto = argv[15];

	
	srand(time(0));

	cout << endl << "Coucou B" << endl;

	if(expectSignalMass!=90 && expectSignalMass!=95 && expectSignalMass!=100 && expectSignalMass!=105 && expectSignalMass!=110 && expectSignalMass!=115 && expectSignalMass!=120 && expectSignalMass!=125 && expectSignalMass!=130 && expectSignalMass!=135 && expectSignalMass!=140 && expectSignalMass!=145 && expectSignalMass!=150)
	{
    		cerr << "ERROR - expectSignalMass has to be integer in range (90,150,5)" << endl;
    		exit(1);
  	}

	system(Form("mkdir -p %s",output.c_str()));

	int ncats = 9;
	//int binning = 105;

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
	cout << endl << "Coucou 1" << endl;
	RooDataHist *dataBinned = new RooDataHist(Form("roohist_data_mass_cat%d",cat),Form("roohist_data_mass_cat%d",cat),RooArgSet(*mass),*data);
	cout << endl << "Coucou 2" << endl;

	RooDataSet *sigMC = (RooDataSet*)w_sig->data(Form("sig_ggh_mass_m%d_cat%d",expectSignalMass,cat));
  	RooDataSet *sigMC_wh = (RooDataSet*)w_sig->data(Form("sig_wh_mass_m%d_cat%d",expectSignalMass,cat));
	RooDataSet *sigMC_zh = (RooDataSet*)w_sig->data(Form("sig_zh_mass_m%d_cat%d",expectSignalMass,cat));	
	cout << endl << "Coucou 3" << endl;
  	RooDataSet *sigMC_vbf = (RooDataSet*)w_sig->data(Form("sig_vbf_mass_m%d_cat%d",expectSignalMass,cat));
  	cout << endl << "Coucou 4" << endl;
	RooDataSet *sigMC_tth = (RooDataSet*)w_sig->data(Form("sig_tth_mass_m%d_cat%d",expectSignalMass,cat));
  	cout << endl << "Coucou 4b" << endl;
	sigMC->append(*sigMC_vbf);
	cout << endl << "Coucou 5" << endl;
  	sigMC->append(*sigMC_wh);
	 cout << endl << "Coucou 6" << endl;
	sigMC->append(*sigMC_zh);
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

	bool truth = true;

	RooAbsReal * sigYield = new RooProduct("sig_yield","sig_yield",RooArgSet(*mu,*sigNorm));
	RooRealVar * bkgYield = new RooRealVar("bkg_yield","bkg_yield",1000.,0.,1.e6);

	cout << endl << "Coucou 2" << endl;
	
	RooRealVar * fraction = NULL;

	//map<string,RooAbsPdf*> bkgPdfs = readDatafile(datafile, mass, cat, truth, fraction); 
	//map<string,RooAbsPdf*> testPdfs = readDatafile(datafile, mass, cat, false); //FIXME
	map< pair<string, RooAbsPdf*>, pair<string*, double*> > bkgPdfs = readDatafile(datafile, mass, cat, truth, fraction, eVeto);	

	//map<string,RooAbsPdf*> sbPdfs;
	map< pair<string, RooAbsPdf*>, pair<string*, double*> > sbPdfs;

	cout << endl << "Truth functions : " << endl;
	
	//for(map<string,RooAbsPdf*>::iterator bkg=bkgPdfs.begin(); bkg!=bkgPdfs.end(); bkg++)
	for(map< pair<string, RooAbsPdf*>, pair<string*, double*> >::iterator bkg=bkgPdfs.begin(); bkg!=bkgPdfs.end(); bkg++)
	{
		cout << endl <<  bkg->first.first << endl;
    		RooAbsPdf *sbMod = new RooAddPdf(Form("sb_%s",bkg->first.first.c_str()),Form("sb_%s",bkg->first.first.c_str()),RooArgList(*(bkg->first.second),*sigPdf),RooArgList(*bkgYield,*sigYield));
	
		pair<string, RooAbsPdf*> pair1 = make_pair(sbMod->GetName(),sbMod);
		pair<string*, double*> pair2 = bkg->second; 

		sbPdfs.insert(make_pair(pair1,pair2));
	}

	
	// ----- Signal + Background Pdf : Test ----- //

	truth = false;

        //map<string,RooAbsPdf*> testPdfs = readDatafile(datafile, mass, cat, truth, fraction); //FIXME
	map< pair<string, RooAbsPdf*>, pair<string*, double*> > testPdfs = readDatafile(datafile, mass, cat, truth, fraction, eVeto);    

 
        //map<string,RooAbsPdf*> test_sbPdfs;
	map< pair<string, RooAbsPdf*>, pair<string*, double*> > test_sbPdfs;

	cout << endl << "Test functions : " << endl;

        for(map< pair<string, RooAbsPdf*>, pair<string*, double*> >::iterator test=testPdfs.begin(); test!=testPdfs.end(); test++)
        {
                cout << endl <<  test->first.first << endl;
                RooAbsPdf *sbMod = new RooAddPdf(Form("sb_%s",test->first.first.c_str()),Form("sb_%s",test->first.first.c_str()),RooArgList(*(test->first.second),*sigPdf),RooArgList(*bkgYield,*sigYield));
    
		pair<string, RooAbsPdf*> pair1 = make_pair(sbMod->GetName(),sbMod);
                pair<string*, double*> pair2 = test->second;
 
                test_sbPdfs.insert(make_pair(pair1,pair2));
        }
		

	// ----- Fitting Truth to data ----- //

	mu->setConstant(false);
	mu->setVal(expectSignal);
	mu->setConstant(true);


	TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);

	for (map< pair<string, RooAbsPdf*>, pair<string*, double*> >::iterator it=sbPdfs.begin(); it!=sbPdfs.end(); it++)
	{
		RooPlot * frame = mass->frame();
    		RooFitResult *fit = (RooFitResult*)it->first.second->fitTo(*data,Save(true));
		cout << endl << "Fit Res Before: " << endl;
      		fit->floatParsInit().Print("v");
      		cout << endl << "Fit Res After: " << endl;
      		fit->floatParsFinal().Print("v");

		data->plotOn(frame,Binning(binning));
		it->first.second->plotOn(frame);
		frame->Draw();
		system(Form("mkdir -p %s/cat%d/mass_%d/",output.c_str(),cat,expectSignalMass));
		system(Form("mkdir -p %s/cat%d/mass_%d/truthToData/",output.c_str(),cat,expectSignalMass));
		system(Form("mkdir -p %s/cat%d/mass_%d/toyHisto/",output.c_str(),cat,expectSignalMass));
		system(Form("mkdir -p %s/cat%d/mass_%d/toyHistoPlusTest/",output.c_str(),cat,expectSignalMass));
		c1->Print(Form("%s/cat%d/mass_%d/truthToData/%s_truthToData_jobID%d.png",output.c_str(),cat,expectSignalMass,it->first.first.c_str(),jobID));
		
		//frame->Clear();
		c1->Clear();
		cout << endl << "coucou fit";
		fit->Delete();
		cout << endl << "coucou fit2";
		frame->Delete();

	}

	c1->Clear();

	// ----- Generation of toys ----- //

	cout << endl << "------------------------------";
	cout << endl << "----- Generation of toys -----";
	cout << endl << "------------------------------" << endl; 

	mu->setConstant(false);
	mu->setVal(expectSignal);
	mu->setConstant(true);

	int iseed = 0;
	bool fitOk = true;

	for (int toy=0; toy<ntoys; toy++)
	{
		cout << endl << "// ----- Toy " << toy << " ----- //" << endl;

		mu->setConstant(false);
		mu->setVal(expectSignal);
        	mu->setConstant(true);
		
		map<string,RooAbsData*> toys;

		for (map< pair<string, RooAbsPdf*>, pair<string*, double*> >::iterator it=sbPdfs.begin(); it!=sbPdfs.end(); it++)
		{
			cout << endl <<  it->first.first<< endl;

			iseed=rand();
			RooRandom::randomGenerator()->SetSeed(iseed) ; 
			//RooDataHist *toyHist = it->second->generateBinned(RooArgSet(*mass),dataBinned->sumEntries(),Extended(),Name(Form("%s_%d",it.c_str(),toy)));
			RooDataHist *toyHist = 0;
			//if(dataBinned->sumEntries() > 200) toyHist = it->first.second->generateBinned(RooArgSet(*mass),dataBinned->sumEntries(),Name(Form("%s_%d",it->first.first.c_str(),toy)));
			//else toyHist = it->first.second->generateBinned(RooArgSet(*mass),2000,Name(Form("%s_%d",it->first.first.c_str(),toy)));
			//toyHist = it->first.second->generateBinned(RooArgSet(*mass),5000,Name(Form("%s_%d",it->first.first.c_str(),toy)));	
			toyHist = it->first.second->generateBinned(RooArgSet(*mass),dataBinned->sumEntries(),Name(Form("%s_%d",it->first.first.c_str(),toy)));		
	
			toys.insert(pair<string,RooDataHist*>(toyHist->GetName(),toyHist));	

			cout << endl << toyHist->GetName() << endl;
	
			RooPlot * frame = mass->frame();
			toyHist->plotOn(frame,Binning(binning));
			it->first.second->plotOn(frame,LineColor(kBlue),LineStyle(kDashed));
			frame->Draw();
			//c1->Print(Form("%s/cat%d/mass_%d/toyHisto/%s_toyHisto%d_jobID%d.png",output.c_str(),cat,expectSignalMass,it->first.first.c_str(),toy,jobID));
			c1->Clear();	
			frame->Delete();

		}

		mu->setConstant(false);
		
		for (map<string,RooAbsData*>::iterator it=toys.begin(); it!=toys.end(); it++)
		{
			cout << endl << " -----  " << it->first << " ----- " << endl;

			double globalMinNLL=1.e6;
			double bestFitValError = 1.0;
			double bestFitVal = 0.0;
			double bias = 0.0;
			for (map< pair<string, RooAbsPdf*>, pair<string*, double*> >::iterator m=test_sbPdfs.begin(); m!=test_sbPdfs.end(); m++) 
			{
				cout << endl << " -----  " << m->first.first << " ----- " << endl;
				
				RooAbsPdf *pdf = m->first.second;
				

				RooArgSet *params = pdf->getParameters(*mass);	
			
				int loopEnd = m->second.second[0];
				cout << endl << "loopEnd = " <<loopEnd << endl;
				 
				//Reset all the parameters to their initial values	
				for(int n = 1; n <= loopEnd; n++) //FIXME
				{
					cout << endl << "m->second.first[n] = " <<m->second.first[n];
					cout << endl << "m->second.second[n] = " << m->second.second[n];
					RooRealVar * var_temp = (RooRealVar*) params->find(m->second.first[n].c_str());
					var_temp->setVal(m->second.second[n]);
				}


				RooFitResult *nom = pdf->fitTo(*(it->second),PrintLevel(-1),PrintEvalErrors(-1),Warnings(false),Save(true));
				nom->floatParsInit().Print("v");
				nom->floatParsFinal().Print("v");
				//RooArgSet *params = pdf->getParameters(*mass);
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

				cout << endl << "bias = " << bias << endl;

				string correspondingTruth = "";

				for (map< pair<string, RooAbsPdf*>, pair<string*, double*> >::iterator it2=sbPdfs.begin(); it2!=sbPdfs.end(); it2++)
				{
					
					if(it->first.find(it2->first.first) == 0)
					{
						correspondingTruth = it2->first.first;
						cout << endl << "correspondingTruth = "<< correspondingTruth; 
						break;
					}
				}

				cout << endl;

				string chainOut = Form("%s/cat%d/mass_%d/Truth_%s_Test_%s.txt",output.c_str(),cat,expectSignalMass,correspondingTruth.c_str(),m->first.first.c_str()); 			

				ofstream fileOut(chainOut.c_str(),ios::app);
				//ofstream fileOut(Form("%s/cat%d/mass_%d/Truth_%s_Test_%s.txt",output.c_str(),cat,expectSignalMass,correspondingTruth.c_str(),m.c_str()),ios::app);



				cout << endl << "ICIIIIIIII !!!! ";
				cout << endl << chainOut << endl;

				//if(bias != 0.0) fileOut << bias << endl; //FIXME			
				if(mu->getVal() > muLow && mu->getVal() < muHigh && mu->getError() > 0.001 && mu->getError() < 99.0 && nom->edm() < 0.01) //FIXME
				{
					fitOk = true;
				}
				else 
				{
					fitOk = false;
				}
				
				if(fitOk == true) fileOut << bias << endl;

				fileOut.close();

				//nom->Delete();
	
				RooPlot * frame = mass->frame();
				it->second->plotOn(frame);
				m->first.second->plotOn(frame);
				frame->Draw();
				//c1->Print(Form("%s/cat%d/mass_%d/toyHistoPlusTest/%s_toyHisto%d_%s_jobID%d.png",output.c_str(),cat,expectSignalMass,it->first.c_str(),toy,m->first.first.c_str(),jobID));
				c1->Clear();
				cout << endl << "before nom";
				nom->Delete();
				cout << endl << "after nom";
				frame->Delete();
	
				/*
				RooRealVar * var_temp = (RooRealVar*) params->find("mean");
				var_temp->setVal(91.187);

				cout << endl << "LAAAAA : var_temp->getVal() = " << var_temp->getVal() << endl;
				*/
			}

		}


	}


	
	system(Form("touch %s/cat_%d_expectSignalMass_%d_job_%d.done",output.c_str(),cat,expectSignalMass,jobID));
	system(Form("rm %s/cat_%d_expectSignalMass_%d_job_%d.fail",output.c_str(),cat,expectSignalMass,jobID));


	delete c1;

	cout << endl << "coucou 1";

	mass->Delete();
	cout << endl << "coucou 1b";
	mu->Delete();
	cout << endl << "coucou 2";
	intL->Delete();
	cout << endl << "coucou 2b";
        sigPdf->Delete();
	cout << endl << "coucou 3";
        sigNorm->Delete();
	cout << endl << "coucou 4";
        sigYield->Delete();
	cout << endl << "coucou 5";
        bkgYield->Delete();
	cout << endl << "coucou 6";

	for(map< pair<string, RooAbsPdf*>, pair<string*, double*> >::iterator m=sbPdfs.begin(); m!=sbPdfs.end(); m++)
        {
                m->first.second->Delete();
		delete [] m->second.first;
		delete [] m->second.second;
        }

	cout << endl << "coucou 7";
	for(map< pair<string, RooAbsPdf*>, pair<string*, double*> >::iterator m=test_sbPdfs.begin(); m!=test_sbPdfs.end(); m++)
	{
		m->first.second->Delete();
		delete [] m->second.first;
		delete [] m->second.second;
	}

	cout << endl << "coucou 7b";

	map<string,RooAbsData*> toys;

        for(map<string,RooAbsData*>::iterator m=toys.begin(); m!=toys.end(); m++)
        {
                m->second->Delete();
        }	

	cout << endl << "coucou 8";
	dataBinned->Delete();
	cout << endl << "coucou 9";
	sigMCBinned->Delete();
	cout << endl << "coucou 10";

	data->Delete();
	cout << endl << "coucou 11";

	sigMC_wh->Delete();
	cout << endl << "coucou 12";
	sigMC_zh->Delete();
	cout << endl << "coucou 13";
	sigMC_vbf->Delete();
	cout << endl << "coucou 14";
	sigMC_tth->Delete();
	cout << endl << "coucou 15";
	//sigMC->Delete();
	cout << endl << "coucou 16";
	

	w_bkg->Delete();
	cout << endl << "coucou 17";
	w_sig->Delete();
	cout << endl << "coucou 18";
	fb->Delete();
	cout << endl << "coucou 19";
	fs->Delete();
	cout << endl << "coucou 20";

	return 0;
}







