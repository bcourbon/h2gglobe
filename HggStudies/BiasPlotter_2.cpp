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
                cerr << "arguments should be passed : output, datafile, startCat, massLow, massHigh, muLow, muHight" <<endl; 
                return 1;
        }
	

	string output = "plotsBias";
	string datafile = "config.dat";
	int startCat = 11;
	double massLow = 95;
	double massHigh = 110;
	double muLow = -5.0;
	double muHigh = 5.0;

	if( argc > 1 ) output = argv[1];
	if( argc > 2 ) datafile = argv[2];
	if( argc > 3 ) 
	{
		std::stringstream ss ( argv[3] );
		ss >> startCat;
	}
	if( argc > 4 )
        {
                std::stringstream ss ( argv[4] );
                ss >> massLow;
        }
	if( argc > 5 ) 
        {
                std::stringstream ss ( argv[5] );
                ss >> massHigh;
        }	
	if( argc > 6 )
        {
                std::stringstream ss ( argv[6] );
                ss >> muLow;
        }
	if( argc > 7 )
        {
                std::stringstream ss ( argv[7] );
                ss >> muHigh;
        }

	system(Form("mkdir -p %s/CombinedPlots",output.c_str()));

	//int startCat = 0;
	int ncats = startCat + 1;
	int binning = 100;

	int nMassPoints = (massHigh - massLow) / 5 + 1;

	// --- style of plots --- //
        gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
        //gStyle->SetOptFit(1111);  //FIXME
	//gStyle->SetOptStat(1111);
	//TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);

	double mu_temp = 0;

	vector <string> truthFunctions;
	vector <string> testFunctions;

	int nbPoints = (massHigh - massLow) / nMassPoints + 1;	

	double biasVector[nbPoints];
	double massVector[nbPoints];
	for(int i = 0; i < nbPoints; i++)
	{
		 biasVector[i] = 0;
		massVector[i] = massLow + i * 5;
	}

	string chain_temp = "";
	int order_temp = 0;

	int nbTruth = 0;
	int nbTest = 0;

	TLatex latexLabel;

	for(int i = startCat; i < ncats; i++)
	{
		cout << "// ----- cat " << i << " ----- //" << endl;

		

		fstream datFile;
        	datFile.open(datafile.c_str(),ios::in);
		
		while(chain_temp != Form("cat%d",i)) datFile >> chain_temp;
		datFile >> chain_temp;
		while(true)
		{ 
			datFile >> chain_temp;
			if(chain_temp == "TestFunctions") break;
			datFile >> order_temp;
			truthFunctions.push_back(Form("%s_%d",chain_temp.c_str(),order_temp));
		}
		while(! datFile.eof())
                {	
			datFile >> chain_temp;
                        if(chain_temp == Form("cat%d",i+1)) break;
                        datFile >> order_temp;
                        testFunctions.push_back(Form("%s_%d",chain_temp.c_str(),order_temp));
		}
		
        	datFile.close();

		nbTruth = truthFunctions.size();
		nbTest = testFunctions.size();

		cout << endl << "coucou A" << endl;
		vector < vector<double> > maxBiasVector;
	
		cout << endl << "coucou B" << endl;
		for(vector<string>::iterator truth = truthFunctions.begin(); truth != truthFunctions.end(); truth++)
                {
			int it_test = 0;
			vector<double> maxBiasVectorRow;
                	for(vector<string>::iterator test = testFunctions.begin(); test != testFunctions.end(); test++)
                        {
				
				for(int j = massLow; j <= massHigh; j+=5)
				{
					cout << endl << "coucou C" << endl;
					TCanvas * c1 = new TCanvas("c1", "c1",0,0,600,600);
					TH1D * muTh1 = new TH1D("muTh1","muTh1",binning,muLow,muHigh);
			
					string temp_chain1 = "";
					string temp_chain2 = "";


					if(truth->find("Bernstein") !=std::string::npos) temp_chain1 = "bern";		
					if(truth->find("Chebychev") !=std::string::npos) temp_chain1 = "che";
					if(truth->find("Exponential") !=std::string::npos) temp_chain1 = "exp";
					if(truth->find("Laurent") !=std::string::npos) temp_chain1 = "lau";
					if(truth->find("PowerLaw") !=std::string::npos) temp_chain1 = "pow";

					if(test->find("Bernstein") !=std::string::npos) temp_chain2 = "bern";    
                                        if(test->find("Chebychev") !=std::string::npos) temp_chain2 = "che";
                                        if(test->find("Exponential") !=std::string::npos) temp_chain2 = "exp";
                                        if(test->find("Laurent") !=std::string::npos) temp_chain2 = "lau";
                                        if(test->find("PowerLaw") !=std::string::npos) temp_chain2 = "pow";	
/*
					for(int order = 1; order <= 11; order++)
					{
						cout << endl <<Form("Bernstein_%d",order);
						cout << endl << "truth->c_str() = "<<truth->c_str();
						if(truth->c_str() == Form("Bernstein_%d",order)) 
						{
							cout << endl <<"Yooo";
							temp_chain1 = Form("bern_%d",order);
						}
						if(truth->c_str() == Form("Chebychev_%d",order)) temp_chain1 = Form("che_%d",order);	
						if(truth->c_str() == Form("Exponential_%d",order)) temp_chain1 = Form("exp_%d",order);
						if(truth->c_str() == Form("Laurent_%d",order)) temp_chain1 = Form("lau_%d",order);
						if(truth->c_str() == Form("PowerLaw_%d",order)) temp_chain1 = Form("pow_%d",order);				

						if(test->c_str() == Form("Bernstein_%d",order)) temp_chain2 = Form("bern_%d",order);
                                                if(test->c_str() == Form("Chebychev_%d",order)) temp_chain2 = Form("che_%d",order);    
                                                if(test->c_str() == Form("Exponential_%d",order)) temp_chain2 = Form("exp_%d",order);
                                                if(test->c_str() == Form("Laurent_%d",order)) temp_chain2 = Form("lau_%d",order);
                                                if(test->c_str() == Form("PowerLaw_%d",order)) temp_chain2 = Form("pow_%d",order); 	


					}
*/
					cout << endl << "truth->c_str() = "<<truth->c_str() << ", test->c_str() = " << test->c_str();
					cout << endl << "temp_chain1 = " << temp_chain1 <<", temp_chain2 = "<<temp_chain2; 

					string chainIn = Form("%s/cat%d/mass_%d/Truth_sb_%s_Test_sb_%s.txt",output.c_str(),i,j,temp_chain1.c_str(),temp_chain2.c_str());
					cout << endl << "chainIn = " << chainIn << endl;
					//system(Form("source suppressNullCar.sh %s",chainIn.c_str()));

					fstream file;
	
					file.open(chainIn.c_str(),ios::in);	

					if(file) cout << endl << "file "<<chainIn << " open" << endl;
					else cout << endl << "file doesn't exist" << endl;

					while(! file.eof())
					{
						file >> mu_temp;
						//cout << endl << mu_temp;
						muTh1->Fill(mu_temp);	
					}
					file.close();
				
		
					muTh1->GetXaxis()->SetTitle("(#mu_{Fit} - #mu_{Gen}) / #sigma_{Fit}");
					muTh1->GetYaxis()->SetTitle("Events");
					
					RooRealVar x ("x","x",0.0,muLow,muHigh);    
                                        RooDataHist rooTh1("rooTh1", "rooTh1", x, muTh1);


					double xMaxHisto = muTh1->GetXaxis()->GetBinCenter(muTh1->GetMaximumBin());
                                        double sigH = muTh1->GetRMS();


                                        RooRealVar meanGaus("meanGaus","meanGaus",xMaxHisto,-10.0,10.0);
                                        RooRealVar sigmaGaus("sigmaGaus","sigmaGaus",sigH,0.0,10.0);
                                        //RooGaussian pdf_fit("pdf_fit","pdf_fit",x,meanGaus,sigmaGaus);

                                        RooRealVar meanV("meanV","meanV",xMaxHisto,-10.0,10.0);
                                        RooRealVar sigmaV("sigmav","sigmaV",0.5,0.0,sigH);
                                        RooRealVar widthV("widthV","widthV",0.5,0.0,sigH);
                                        //RooVoigtian * pdf = new RooVoigtian("pdf","Voigtian",variable,mean,sigma,width);
                                        RooVoigtian pdf_fit("pdf_fit","pdf_fit",x,meanV,sigmaV,widthV);


					double median = Median(muTh1);
                                        RooPlot * frame = x.frame();
                                        RooFitResult * fit = 0;
					if(fabs(median) < 5.0) fit = (RooFitResult*) pdf_fit.fitTo(rooTh1,Save(true),Range(median-1.6,median+1.6));

                                        rooTh1.plotOn(frame);
                                        pdf_fit.plotOn(frame);
                                        frame->Draw();

					/*	
					TF1 * pdf = new TF1("pdf","pdf_fit",muLow,muHigh);
					muTh1->Fit("pdf");
					muTh1->Draw("E1");
					pdf->Draw("SAMES");
					
					double mean = pdf->GetParameter(1);	
					double sigma = pdf->GetParameter(2);
					*/
					//double mean = meanGaus.getVal();
                                        //double sigma = sigmaGaus.getVal();
                                        double mean = meanV.getVal();                         
                                        double sigma = sigmaV.getVal();

					//double bias = (expectSignal - mean) / sigma;
					// vv FIXME vv
					double bias = mean;
					if(fabs(median) < 5.0) bias = mean;	
					else bias = median;
					//double bias = muTh1->GetMean();
					//double median = Median(muTh1);	

					//biasVector.push_back(bias);
					int temp_bV = (j-massLow) / 5;
					biasVector[temp_bV] = bias;
					

					latexLabel.SetTextSize(0.030);
                			latexLabel.SetNDC();
					latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2012, #sqrt{s} = 8 TeV");
					//latexLabel.DrawLatex(0.17, 0.88,Form("mean = %3.1f #pm %3.1f",mean,pdf->GetParError(1)));
					//latexLabel.DrawLatex(0.17, 0.83,Form("#sigma = %3.1f #pm %3.1f",sigma,pdf->GetParError(2)));
					latexLabel.DrawLatex(0.17, 0.88,Form("mean = %3.2f #pm %3.2f",mean,meanGaus.getError()));
                                        latexLabel.DrawLatex(0.17, 0.83,Form("#sigma = %3.2f #pm %3.2f",sigma,sigmaGaus.getError()));
					latexLabel.DrawLatex(0.17, 0.78,Form("median = %3.2f",median));
					//latexLabel.DrawLatex(0.17, 0.78,Form("#bias = %3.1f",bias));	

					c1->Print(Form("%s/cat%d/mass_%d/Bias_Truth_sb_%s_Test_sb_%s.pdf",output.c_str(),i,j,truth->c_str(),test->c_str()));
					c1->Print(Form("%s/cat%d/mass_%d/Bias_Truth_sb_%s_Test_sb_%s.png",output.c_str(),i,j,truth->c_str(),test->c_str()));	
					c1->Clear();

					muTh1->Delete();
					muTh1 = 0;

					delete c1;
					c1 = 0;
				}

				double maxBias = 0.0;
				double maxBiasAbs = 0.0;

				for(int k = 0; k < nbPoints; k++)
				{
					if(abs(biasVector[k]) > maxBiasAbs)
					{
						maxBiasAbs = abs(biasVector[k]);
						maxBias = biasVector[k];
					}		
				}

				cout << endl << "coucou 1" << endl;
				cout << "maxBias = "<<maxBias << endl;
				cout << "it_test = "<<it_test << endl;
				//maxBiasVector[it_test].push_back(maxBias); //FIXME	
				maxBiasVectorRow.push_back(maxBias);
				
				cout << endl << "coucou 2" << endl;
				TCanvas * c2 = new TCanvas("c2", "c2",0,0,600,600);

				TGraph * graph = new TGraph(nbPoints,massVector,biasVector);
				
				TLine lineD(massLow - 5,-0.14,massHigh + 5,-0.14);
				TLine lineU(massLow - 5,0.14,massHigh + 5,0.14);
				lineD.SetLineStyle(kDashed);
				lineD.SetLineColor(kRed);
				lineU.SetLineStyle(kDashed);
                                lineU.SetLineColor(kRed);

				graph->GetXaxis()->SetTitle("M_{#gamma#gamma}");
				graph->GetYaxis()->SetTitle("Bias");
			
				graph->SetMarkerColor(4);
				graph->SetMarkerStyle(20);
				graph->SetMarkerSize(0.6);			

				graph->Draw("APL");
				lineD.Draw();
				lineU.Draw();
				graph->GetXaxis()->SetLimits(massLow - 5,massHigh + 5);

				latexLabel.SetTextSize(0.030);
                                latexLabel.SetNDC();
                                latexLabel.DrawLatex(0.13, 0.96, "CMS Preliminary 2012, #sqrt{s} = 8 TeV");
				latexLabel.DrawLatex(0.17, 0.88, Form("Category %d",i));
				latexLabel.DrawLatex(0.17, 0.83, Form("Truth : %s",(*truth).c_str()));
				latexLabel.DrawLatex(0.17, 0.78, Form("Test : %s",(*test).c_str()));
				latexLabel.DrawLatex(0.66, 0.88, Form("Bias max =  %.2f",maxBias));		
	
			
				graph->GetYaxis()->SetRangeUser(-0.5,0.5);
	
				c2->Print(Form("%s/CombinedPlots/cat%d_Truth_%s_Test_%s.pdf",output.c_str(),i,truth->c_str(),test->c_str()));
				c2->Print(Form("%s/CombinedPlots/cat%d_Truth_%s_Test_%s.png",output.c_str(),i,truth->c_str(),test->c_str()));	
				c2->Print(Form("%s/CombinedPlots/cat%d_Truth_%s_Test_%s.C",output.c_str(),i,truth->c_str(),test->c_str()));	

				c2->Clear();

				graph->Delete();
				graph = 0;
                                delete c2;
                                c2 = 0;

				//biasVector.erase(biasVector.begin(),biasVector.end());
				it_test++;
			}
			maxBiasVector.push_back(maxBiasVectorRow);
			cout << endl << "coucou 3" << endl;
		}

		fstream maxBiasFile;
		maxBiasFile.open(Form("%s/CombinedPlots/maxBias_cat%d.txt",output.c_str(),i),ios::out);

		cout << endl << "coucou 4" << endl;
		string maxBiasChain = "";
		for(int j = 0; j < nbTest; j++)
		{
			for(int k = 0; k < nbTruth; k++)
			{
				if(abs(maxBiasVector[k][j]) > 0.14) maxBiasChain = Form("\\textbf{%.2f}",maxBiasVector[k][j]);
				else maxBiasChain = Form("%.2f",maxBiasVector[k][j]);	
				if(k == 0) maxBiasFile << maxBiasChain;
				else maxBiasFile << " & " << maxBiasChain;
				if(nbTruth == 3 && k == 2) maxBiasFile << " & X ";
			}
			maxBiasFile << endl;
		}
		maxBiasFile.close();
		cout << endl << "coucou 5" << endl;

		truthFunctions.erase(truthFunctions.begin(),truthFunctions.end());
		testFunctions.erase(testFunctions.begin(),testFunctions.end());
		cout << endl << "coucou 6" << endl;
	}


	return 0;
}







