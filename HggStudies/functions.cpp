#include "functions.h"

int rowsNumberInFile(string filename)
{
        ifstream in(filename.c_str());

        string row = "";
        int nbRows = 0;

        while(std::getline(in, row)) nbRows++;

        in.close();

        return nbRows;
}

// --- Test the existence of a file --- //
bool FileExists(const char* FileName)
{
        FILE* fp = NULL;

        //will not work if you do not have read permissions

        //to the file, but if you don't have read, it

        //may as well not exist to begin with.

        fp = fopen( FileName, "rb" );
        if( fp != NULL )
        {
                fclose( fp );
                return true;
        }

        return false;
}


double Median(const TH1D * h1) {

   int n = h1->GetXaxis()->GetNbins(); 
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray();
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]);
}


//map<string,RooAbsPdf*> readDatafile(string datafile, RooRealVar *obs_var, int cat, bool truth, RooRealVar* &fraction)
map< pair<string, RooAbsPdf*>, pair<string*, double*> > readDatafile(string datafile, RooRealVar *obs_var, int cat, bool truth, RooRealVar* &fraction, string eVeto)
{
	string temp = "";
	int tempOrder = 1;
	//map<string,RooAbsPdf*> pdfs_map;
	//pair<RooAbsPdf*,RooRealVar*> pdf;
	map< pair<string, RooAbsPdf*>, pair<string*, double*> > pdfs_map; 
	pair< pair<RooAbsPdf*,RooRealVar*>, pair<string*, double*> > pdf;
	fstream file;
	file.open(datafile.c_str(),ios::in);

	if(file) cout << endl << "dat file open" << endl;
	else cout << endl << "A problem happened with your dat file" << endl; 
	
	while(temp != Form("cat%d",cat))
	{
		file >> temp;
	}
	
	file >> temp;

	if(truth == true)
	{
		while(true)
		{
			file >> temp;	
			if(temp == "TestFunctions") break;			
			file >> tempOrder;
			//cout << endl << "temp = " <<temp << ", tempOrder = " <<tempOrder <<endl;
			pdf = getPdf(temp,tempOrder,obs_var,fraction, cat, truth, eVeto);
			pdfs_map.insert( make_pair( make_pair(pdf.first.first->GetName(), pdf.first.first ), pdf.second ) );
			//pdfs_map.insert(pair<string,RooAbsPdf*>(pdf.first->GetName(),pdf.first));	
		}
	}
	else
	{
		while(temp != "TestFunctions")
		{
			file >> temp;
		}

		while(! file.eof())
		{
			file >> temp;  
			if(temp == Form("cat%d",cat+1)) break;  
                        file >> tempOrder;
                        pdf = getPdf(temp,tempOrder,obs_var,fraction, cat, truth, eVeto);
			pdfs_map.insert( make_pair( make_pair(pdf.first.first->GetName(), pdf.first.first ), pdf.second ) );
			//pdfs_map.insert(pair<string,RooAbsPdf*>(pdf.first->GetName(),pdf.first));   
		}
	}

	file.close();

	return pdfs_map;	

}

pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> >  getPdf(string type, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto)
{
	cout << endl << "coucou in getPdf : " << endl;
  
	if (type=="BernsteinPlusVoigtian") return getBernsteinPlusVoigtian(Form("BernsteinPlusVoigtian_%d",order),order,obs_var,fraction, cat, truth, eVeto); 
  	else if (type=="ChebychevPlusVoigtian") return getChebychevPlusVoigtian(Form("ChebychevPlusVoigtian_%d",order),order,obs_var,fraction, cat, truth, eVeto);
	else if (type=="ExponentialPlusVoigtian") return getExponentialPlusVoigtian(Form("ExponentialPlusVoigtian_%d",order),order,obs_var,fraction, cat, truth, eVeto);
	else if (type=="PowerLawPlusVoigtian") return getPowerLawPlusVoigtian(Form("PowerLawPlusVoigtian_%d",order),order,obs_var,fraction, cat, truth, eVeto);
	else if (type == "LaurentPlusVoigtian") return getLaurentPlusVoigtian(Form("LaurentPlusVoigtian_%d",order),order,obs_var,fraction, cat, truth, eVeto);
	else if (type == "PolynomialsPlusVoigtian") return getPolynomialsPlusVoigtian(Form("PolynomialsPlusVoigtian_%d",order),order,obs_var,fraction, cat, truth, eVeto);	

	else if (type=="Bernstein") return getBernstein(Form("Bernstein_%d",order),order,obs_var,fraction, cat, truth); 
        else if (type=="Chebychev") return getChebychev(Form("Chebychev_%d",order),order,obs_var,fraction, cat, truth);
        else if (type=="Exponential") return getExponential(Form("Exponential_%d",order),order,obs_var,fraction, cat, truth);
        else if (type=="PowerLaw") return getPowerLaw(Form("PowerLaw_%d",order),order,obs_var,fraction, cat, truth);
        else if (type == "Laurent") return getLaurent(Form("Laurent_%d",order),order,obs_var,fraction, cat, truth);

	else if (type=="BernsteinPlusCB") return getBernsteinCB(Form("BernsteinPlusCB_%d",order),order,obs_var,fraction, cat, truth);
        else if (type=="ChebychevPlusCB") return getChebychevCB(Form("ChebychevPlusCB_%d",order),order,obs_var,fraction, cat, truth);
        else if (type=="ExponentialPlusCB") return getExponentialCB(Form("ExponentialPlusCB_%d",order),order,obs_var,fraction, cat, truth);
        else if (type=="PowerLawPlusCB") return getPowerLawCB(Form("PowerLawPlusCB_%d",order),order,obs_var,fraction, cat, truth);
        else if (type == "LaurentPlusCB") return getLaurentCB(Form("LaurentPlusCB_%d",order),order,obs_var,fraction, cat, truth);	

	
	else if (type=="BernsteinPlusCBgaus") return getBernsteinCBgaus(Form("BernsteinPlusCBgaus_%d",order),order,obs_var,fraction, cat, truth);
	else if (type=="ChebychevPlusCBgaus") return getChebychevCBgaus(Form("ChebychevPlusCBgaus_%d",order),order,obs_var,fraction, cat, truth);
	else if (type=="ExponentialPlusCBgaus") return getExponentialCBgaus(Form("ExponentialPlusCBgaus_%d",order),order,obs_var,fraction, cat, truth);
	else if (type=="PowerLawPlusCBgaus") return getPowerLawCBgaus(Form("PowerLawPlusCBgaus_%d",order),order,obs_var,fraction, cat, truth);
	else if (type == "LaurentPlusCBgaus") return getLaurentCBgaus(Form("LaurentPlusCBgaus_%d",order),order,obs_var,fraction, cat, truth);	

	else if (type == "Voigtian") return getVoigtian("Voigtian",order,obs_var,fraction, cat, truth);
	else if (type == "CB") return getCB("CB",order,obs_var,fraction, cat, truth);
	else if (type == "CBgaus") return getCBgaus("CBgaus",order,obs_var,fraction, cat, truth);	
	else if (type == "CBpCB") return getCBpCB("CBpCB",order,obs_var,fraction, cat, truth);
	else if (type == "CBconvCB") return getCBconvCB("CBconvCB",order,obs_var,fraction, cat, truth);
	
	/*else if (type=="Chebychev") return getChebychevPlusVoigtian(Form("%s_cheb%d",ext,order),order); 
  	else if (type=="Exponential") return getExponentialSingle(Form("%s_exp%d",ext,order),order); 
  	else if (type=="PowerLaw") return getPowerLawSingle(Form("%s_pow%d",ext,order),order); 
  	else if (type=="Laurent") return getLaurentSeries(Form("%s_lau%d",ext,order),order);
  	*/else 
  	{
    		cerr << "ERROR -- getPdf() -- type " << type << " not recognised." << endl;
		pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > def;
		def.first.first = NULL;
		def.first.second = NULL;
		def.second.first = NULL;		
		def.second.second = NULL;
    		return def;
  	}
}


//RooAbsPdf* getBernsteinPlusVoigtian(string function_name, int order, RooRealVar *obs_var)
pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getBernsteinPlusVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto)
{
	RooArgList *coeffList1 = new RooArgList();
	
	RooRealVar *mean = NULL; 
        RooRealVar *sigma = NULL; 
        RooRealVar *width = NULL;

	double centerVal0, lowVal0, highVal0, lowVal0_test, highVal0_test;
	double centerVal1, lowVal1, highVal1, lowVal1_test, highVal1_test;
	double centerVal2, lowVal2, highVal2, lowVal2_test, highVal2_test;
	
	if(eVeto == "defaultE")
        {
                if(cat == 0) { centerVal0 = 90.40; centerVal1 = 2.49; centerVal2 = 1.74; }
                if(cat == 1) { centerVal0 = 90.33; centerVal1 = 3.45; centerVal2 = 1.00; }
                if(cat == 2) { centerVal0 = 90.34; centerVal1 = 4.38; centerVal2 = 0.48; }
                if(cat == 3) { centerVal0 = 89.87; centerVal1 = 4.92; centerVal2 = 0.93; }
                if(cat == 4) { centerVal0 = 89.62; centerVal1 = 2.49; centerVal2 = 2.22; }
                if(cat == 5) { centerVal0 = 89.99; centerVal1 = 3.81; centerVal2 = 1.04; }
                if(cat == 6) { centerVal0 = 90.19; centerVal1 = 3.04; centerVal2 = 1.38; }
                if(cat == 7) { centerVal0 = 90.13; centerVal1 = 3.75; centerVal2 = 1.42; }
                if(cat == 8) { centerVal0 = 89.75; centerVal1 = 2.49; centerVal2 = 0.0; } // ~ no entries
                if(cat == 9) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 1.41; }
                if(cat == 10) { centerVal0 = 90.28; centerVal1 = 2.49; centerVal2 = 1.57; }
                if(cat == 11) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 0.01; } //NO entries
                if(cat == 12) { centerVal0 = 90.24; centerVal1 = 3.21; centerVal2 = 1.70; }
                if(cat == 13) { centerVal0 = 90.11; centerVal1 = 5.32; centerVal2 = 0.40; }
        }
        if(eVeto == "notHasPixelSeed")
        {
                if(cat == 0) { centerVal0 = 90.40; centerVal1 = 2.49; centerVal2 = 1.74; }
                if(cat == 1) { centerVal0 = 90.32; centerVal1 = 3.48; centerVal2 = 1.00; }
                if(cat == 2) { centerVal0 = 90.33; centerVal1 = 4.42; centerVal2 = 0.48; }
                if(cat == 3) { centerVal0 = 89.84; centerVal1 = 4.80; centerVal2 = 1.01; }
                if(cat == 4) { centerVal0 = 89.60; centerVal1 = 2.49; centerVal2 = 2.24; }
                if(cat == 5) { centerVal0 = 89.99; centerVal1 = 3.09; centerVal2 = 1.42; }
                if(cat == 6) { centerVal0 = 90.14; centerVal1 = 3.37; centerVal2 = 1.22; }
                if(cat == 7) { centerVal0 = 90.10; centerVal1 = 3.82; centerVal2 = 1.38; }
                if(cat == 8) { centerVal0 = 89.75; centerVal1 = 2.49; centerVal2 = 0.0; } // ~ no entries
                if(cat == 9) { centerVal0 = 90.48; centerVal1 = 2.49; centerVal2 = 1.46; }
                if(cat == 10) { centerVal0 = 90.34; centerVal1 = 2.49; centerVal2 = 1.68; }
                if(cat == 11) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 0.01; } //NO entries
                if(cat == 12) { centerVal0 = 90.28; centerVal1 = 2.72; centerVal2 = 1.57; }
                if(cat == 13) { centerVal0 = 90.05; centerVal1 = 3.90; centerVal2 = 1.38; }
        }
	
	lowVal0 = centerVal0 - 1.5; highVal0 = centerVal0 + 1.5;
	lowVal1 = centerVal1 - 7.0; highVal1 = centerVal1 + 7.0;
	lowVal2 = centerVal2 - 7.0; highVal2 = centerVal2 + 7.0;

	lowVal1 = 2.49; //2.49 pour Z theo
	centerVal1 = 2.49;
	highVal1 = 2.49;

	//if(cat == 10) lowVal1 = 1.0;       
        //if(cat == 10) lowVal2 = 0.5;

	lowVal0_test = centerVal0 - 2.5; highVal0_test = centerVal0 + 2.5;
        lowVal1_test = centerVal1 - 8.0; highVal1_test = centerVal1 + 8.0;
        lowVal2_test = centerVal2 - 8.0; highVal2_test = centerVal2 + 8.0;;

        if(lowVal1 < 0) lowVal1 = 0.0;
        if(lowVal2 < 0) lowVal2 = 0.0;
        if(lowVal1_test < 0) lowVal1_test = 0.0;
        if(lowVal2_test < 0) lowVal2_test = 0.0;

	lowVal0 = lowVal0_test = highVal0 = highVal0_test = centerVal0;
        lowVal1 = lowVal1_test = highVal1 = highVal1_test = centerVal1;
        lowVal2 = lowVal2_test = highVal2 = highVal2_test = centerVal2;
	
	if(truth)
	{	/*
		mean = new RooRealVar("mean","mean",91.188,90.0,92.0);		
		sigma = new RooRealVar("sigma","sigma",2.5,2.5,5.0);
		width = new RooRealVar("width","width",1.5,0.0,9.0);
		*/
		mean = new RooRealVar("mean","mean",centerVal0,lowVal0,highVal0);          
                sigma = new RooRealVar("sigma","sigma",centerVal1,lowVal1,highVal1);
                width = new RooRealVar("width","width",centerVal2,lowVal2,highVal2);
	}
	else
	{
		mean = new RooRealVar("mean","mean",centerVal0,lowVal0_test,highVal0_test);
                sigma = new RooRealVar("sigma","sigma",centerVal1,lowVal1_test,highVal1_test);
                width = new RooRealVar("width","width",centerVal2,0.0,5.0);
	}

	pair<string*, double*> initPar;
        initPar.first = new string[order+3+1+1];
        initPar.second = new double[order+3+1+1];

	initPar.first[0] = "nParameters";
	initPar.second[0] = 3 + order + 1;

	initPar.first[1] = "mean";
	initPar.second[1] = mean->getVal();
	initPar.first[2] = "sigma";
        initPar.second[2] = sigma->getVal();
	initPar.first[3] = "width";
        initPar.second[3] = width->getVal();	


	RooVoigtian * voigt = new RooVoigtian("voigt","voigt",*obs_var,*mean,*sigma,*width);

  	for (int i = 0; i < order; i++)
	{
    		string name = Form("%s_p%d",function_name.c_str(),i);
    		//params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    		RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
    		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
    		//params.insert(pair<string,RooRealVar*>(name,param));
    		//prods.insert(pair<string,RooFormulaVar*>(name,form));
    		//coeffList->add(*params[name]);
  		coeffList1->add(*param);
		//coeffList1->add(*form);

		initPar.first[4 + i] = name;
        	initPar.second[4 + i] = param->getVal();

	}
	
	RooBernstein * bern = new RooBernstein("bern","bern",*obs_var,*coeffList1);

	cout << endl << "HERE !!!! fct name = " << function_name << endl;

	if(truth)
        {
                //if(cat == 0 || cat == 1 || cat == 2 || cat == 3 || cat == 4) fraction = new RooRealVar("fraction","fraction",0.99,0.5,1.0);
                //else if(cat == 8) fraction = new RooRealVar("fraction","fraction",0.89,0.0,0.99);
                //else fraction = new RooRealVar("fraction","fraction",0.90,0.0,1.0);
        	fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}
        else
        {
                //fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
        	fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}

	initPar.first[4 + order] = "fraction";
        initPar.second[4 + order] = fraction->getVal();


	RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*bern,*voigt,*fraction);

	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
	return_var.first.first = pdf;
	return_var.first.second = fraction;
	return_var.second = initPar;

  	return return_var;
}

pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getBernstein(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{
	RooArgList *coeffList1 = new RooArgList();

	pair<string*, double*> initPar;
        initPar.first = new string[order+2];
        initPar.second = new double[order+2];

	initPar.first[0] = "nParameters";
	initPar.second[0] = order;

  	for (int i = 0; i < order; i++)
	{
    		string name = Form("%s_p%d",function_name.c_str(),i);
    		//params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    		RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
    		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
    		//params.insert(pair<string,RooRealVar*>(name,param));
    		//prods.insert(pair<string,RooFormulaVar*>(name,form));
    		//coeffList->add(*params[name]);
  		coeffList1->add(*param);
		//coeffList1->add(*form);

		initPar.first[1 + i] = name;
        	initPar.second[1 + i] = param->getVal();

	}
	
	RooBernstein * bern = new RooBernstein("bern","bern",*obs_var,*coeffList1);

	cout << endl << "HERE !!!! fct name = " << function_name << endl;

	fraction = new RooRealVar("fraction","fraction",1.0,1.0,1.0);

	initPar.first[1 + order] = "fraction";
        initPar.second[1 + order] = fraction->getVal();

	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
	return_var.first.first = bern;
	return_var.first.second = fraction;
	return_var.second = initPar;

  	return return_var;
}



pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getBernsteinCBgaus(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{
	RooArgList *coeffList1 = new RooArgList();

	RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;
	
	RooRealVar * meanGaus = NULL;
	RooRealVar * sigmaGaus = NULL;

	if(truth)
	{
		meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);	
		sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
		alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
		nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);
	
		meanGaus = new RooRealVar("meanGaus","meanGaus",0.0,-0.5,0.5);
		sigmaGaus = new RooRealVar("sigmaGaus","sigmaGaus",1.5,0.0,5.0);
	}	
	else
	{
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

                meanGaus = new RooRealVar("meanGaus","meanGaus",0.0,-0.5,0.5);
                sigmaGaus = new RooRealVar("sigmaGaus","sigmaGaus",1.5,0.0,5.0);
	}

	pair<string*, double*> initPar;
        initPar.first = new string[order+6+1+1];
        initPar.second = new double[order+6+1+1];

	initPar.first[0] = "nParameters";
	initPar.second[0] = 6 + order + 1;

	initPar.first[1] = "meanCB";
	initPar.second[1] = meanCB->getVal();
	initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
	initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();	
	initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();
	initPar.first[5] = "meanGaus";
        initPar.second[5] = meanGaus->getVal();
	initPar.first[6] = "sigmaGaus";
        initPar.second[6] = sigmaGaus->getVal();

	RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);
	RooGaussian * gaus = new RooGaussian("gaus","gaus",*obs_var,*meanGaus,*sigmaGaus);	

	RooFFTConvPdf * conv  = new RooFFTConvPdf("conv","conv",*obs_var,*cb,*gaus);

  	for (int i = 0; i < order; i++)
	{
    		string name = Form("%s_p%d",function_name.c_str(),i);
    		//params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    		RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
    		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
    		//params.insert(pair<string,RooRealVar*>(name,param));
    		//prods.insert(pair<string,RooFormulaVar*>(name,form));
    		//coeffList->add(*params[name]);
  		coeffList1->add(*param);
		//coeffList1->add(*form);

		initPar.first[7 + i] = name;
        	initPar.second[7 + i] = param->getVal();

	}
	
	RooBernstein * bern = new RooBernstein("bern","bern",*obs_var,*coeffList1);

	cout << endl << "HERE !!!! fct name = " << function_name << endl;

	if(truth)
        {
                if(cat == 7) fraction = new RooRealVar("fraction","fraction",0.90,0.0,1.0);
                //else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
                else fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99);
        }
        else
        {
                fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
        }

	initPar.first[7 + order] = "fraction";
        initPar.second[7 + order] = fraction->getVal();


	RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*bern,*conv,*fraction);

	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
	return_var.first.first = pdf;
	return_var.first.second = fraction;
	return_var.second = initPar;

  	return return_var;
}

pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getBernsteinCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{
	RooArgList *coeffList1 = new RooArgList();

	RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;
	
	if(truth)
	{
		meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);	
		sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
		alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
		nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);
	
	}	
	else
	{
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

	}

	pair<string*, double*> initPar;
        initPar.first = new string[order+4+1+1];
        initPar.second = new double[order+4+1+1];

	initPar.first[0] = "nParameters";
	initPar.second[0] = 4 + order + 1;

	initPar.first[1] = "meanCB";
	initPar.second[1] = meanCB->getVal();
	initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
	initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();	
	initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();

	RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);

  	for (int i = 0; i < order; i++)
	{
    		string name = Form("%s_p%d",function_name.c_str(),i);
    		//params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    		RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
    		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
    		//params.insert(pair<string,RooRealVar*>(name,param));
    		//prods.insert(pair<string,RooFormulaVar*>(name,form));
    		//coeffList->add(*params[name]);
  		coeffList1->add(*param);
		//coeffList1->add(*form);

		initPar.first[5 + i] = name;
        	initPar.second[5 + i] = param->getVal();

	}
	
	RooBernstein * bern = new RooBernstein("bern","bern",*obs_var,*coeffList1);

	cout << endl << "HERE !!!! fct name = " << function_name << endl;

	if(truth)
        {
                //if(cat == 7) fraction = new RooRealVar("fraction","fraction",0.90,0.0,1.0);
                //else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
                //else fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99);
        	fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99);
	}
        else
        {
                fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
        }

	initPar.first[5 + order] = "fraction";
        initPar.second[5 + order] = fraction->getVal();


	RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*bern,*cb,*fraction);

	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
	return_var.first.first = pdf;
	return_var.first.second = fraction;
	return_var.second = initPar;

  	return return_var;
}

//RooAbsPdf* getChebychevPlusVoigtian(string function_name, int order, RooRealVar *obs_var)
pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getChebychevPlusVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto)
{
	RooArgList *coeffList1 = new RooArgList();

        RooRealVar *mean = NULL;           
        RooRealVar *sigma = NULL;          
        RooRealVar *width = NULL;

	double centerVal0, lowVal0, highVal0, lowVal0_test, highVal0_test;
	double centerVal1, lowVal1, highVal1, lowVal1_test, highVal1_test;
	double centerVal2, lowVal2, highVal2, lowVal2_test, highVal2_test;
	
	if(eVeto == "defaultE")
        {
                if(cat == 0) { centerVal0 = 90.40; centerVal1 = 2.49; centerVal2 = 1.74; }
                if(cat == 1) { centerVal0 = 90.33; centerVal1 = 3.45; centerVal2 = 1.00; }
                if(cat == 2) { centerVal0 = 90.34; centerVal1 = 4.38; centerVal2 = 0.48; }
                if(cat == 3) { centerVal0 = 89.87; centerVal1 = 4.92; centerVal2 = 0.93; }
                if(cat == 4) { centerVal0 = 89.62; centerVal1 = 2.49; centerVal2 = 2.22; }
                if(cat == 5) { centerVal0 = 89.99; centerVal1 = 3.81; centerVal2 = 1.04; }
                if(cat == 6) { centerVal0 = 90.19; centerVal1 = 3.04; centerVal2 = 1.38; }
                if(cat == 7) { centerVal0 = 90.13; centerVal1 = 3.75; centerVal2 = 1.42; }
                if(cat == 8) { centerVal0 = 89.75; centerVal1 = 2.49; centerVal2 = 0.0; } // ~ no entries
                if(cat == 9) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 1.41; }
                if(cat == 10) { centerVal0 = 90.28; centerVal1 = 2.49; centerVal2 = 1.57; }
                if(cat == 11) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 0.01; } //NO entries
                if(cat == 12) { centerVal0 = 90.24; centerVal1 = 3.21; centerVal2 = 1.70; }
                if(cat == 13) { centerVal0 = 90.11; centerVal1 = 5.32; centerVal2 = 0.40; }
        }
        if(eVeto == "notHasPixelSeed")
        {
                if(cat == 0) { centerVal0 = 90.40; centerVal1 = 2.49; centerVal2 = 1.74; }
                if(cat == 1) { centerVal0 = 90.32; centerVal1 = 3.48; centerVal2 = 1.00; }
                if(cat == 2) { centerVal0 = 90.33; centerVal1 = 4.42; centerVal2 = 0.48; }
                if(cat == 3) { centerVal0 = 89.84; centerVal1 = 4.80; centerVal2 = 1.01; }
                if(cat == 4) { centerVal0 = 89.60; centerVal1 = 2.49; centerVal2 = 2.24; }
                if(cat == 5) { centerVal0 = 89.99; centerVal1 = 3.09; centerVal2 = 1.42; }
                if(cat == 6) { centerVal0 = 90.14; centerVal1 = 3.37; centerVal2 = 1.22; }
                if(cat == 7) { centerVal0 = 90.10; centerVal1 = 3.82; centerVal2 = 1.38; }
                if(cat == 8) { centerVal0 = 89.75; centerVal1 = 2.49; centerVal2 = 0.0; } // ~ no entries
                if(cat == 9) { centerVal0 = 90.48; centerVal1 = 2.49; centerVal2 = 1.46; }
                if(cat == 10) { centerVal0 = 90.34; centerVal1 = 2.49; centerVal2 = 1.68; }
                if(cat == 11) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 0.01; } //NO entries
                if(cat == 12) { centerVal0 = 90.28; centerVal1 = 2.72; centerVal2 = 1.57; }
                if(cat == 13) { centerVal0 = 90.05; centerVal1 = 3.90; centerVal2 = 1.38; }
        }

	lowVal0 = centerVal0 - 1.5; highVal0 = centerVal0 + 1.5;
	lowVal1 = centerVal1 - 7.0; highVal1 = centerVal1 + 7.0;
	lowVal2 = centerVal2 - 7.0; highVal2 = centerVal2 + 7.0;

	//if(cat == 10) lowVal1 = 1.0; 
        //if(cat == 10) lowVal2 = 0.5;
	lowVal1 = 2.49; //2.49 pour Z theo
	centerVal1 = 2.49;
	highVal1 = 2.49;

	lowVal0_test = centerVal0 - 2.5; highVal0_test = centerVal0 + 2.5;
        lowVal1_test = centerVal1 - 8.0; highVal1_test = centerVal1 + 8.0;
        lowVal2_test = centerVal2 - 8.0; highVal2_test = centerVal2 + 8.0;

        if(lowVal1 < 0) lowVal1 = 0.0;
        if(lowVal2 < 0) lowVal2 = 0.0;
        if(lowVal1_test < 0) lowVal1_test = 0.0;
        if(lowVal2_test < 0) lowVal2_test = 0.0;

	lowVal0 = lowVal0_test = highVal0 = highVal0_test = centerVal0;
        lowVal1 = lowVal1_test = highVal1 = highVal1_test = centerVal1;
        lowVal2 = lowVal2_test = highVal2 = highVal2_test = centerVal2;
	
	if(truth)
	{	/*
		mean = new RooRealVar("mean","mean",91.188,90.0,92.0);		
		sigma = new RooRealVar("sigma","sigma",2.5,2.5,5.0);
		width = new RooRealVar("width","width",1.5,0.0,9.0);
		*/
		mean = new RooRealVar("mean","mean",centerVal0,lowVal0,highVal0);          
                sigma = new RooRealVar("sigma","sigma",centerVal1,lowVal1,highVal1);
                width = new RooRealVar("width","width",centerVal2,lowVal2,highVal2);
	}
	else
	{
		mean = new RooRealVar("mean","mean",centerVal0,lowVal0_test,highVal0_test);
                sigma = new RooRealVar("sigma","sigma",centerVal1,lowVal1_test,highVal1_test);
                width = new RooRealVar("width","width",centerVal2,0.0,5.0);
	}

	pair<string*, double*> initPar;
        initPar.first = new string[order+3+1+1];
        initPar.second = new double[order+3+1+1];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 3 + order + 1;

        initPar.first[1] = "mean";
        initPar.second[1] = mean->getVal();
        initPar.first[2] = "sigma";
        initPar.second[2] = sigma->getVal();
        initPar.first[3] = "width";
        initPar.second[3] = width->getVal();


	RooVoigtian * voigt = new RooVoigtian("voigt","voigt",*obs_var,*mean,*sigma,*width);

  	for (int i = 0; i < order; i++)
	{
    		string name = Form("%s_p%d",function_name.c_str(),i);
    		//params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    		RooRealVar *param = NULL;
		if(truth) param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
		else
		{
			param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
			//if(i == 0) param = new RooRealVar(name.c_str(),name.c_str(),-1.5,-3.0,3.0);
			//else param = new RooRealVar(name.c_str(),name.c_str(),0.01,-1.0,1.0);

		}
    		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
    		//params.insert(pair<string,RooRealVar*>(name,param));
    		//prods.insert(pair<string,RooFormulaVar*>(name,form));
    		//coeffList->add(*params[name]);
  		coeffList1->add(*param);
		//coeffList1->add(*form);
		
		initPar.first[4 + i] = name;
                initPar.second[4 + i] = param->getVal();
	}
	
	RooChebychev * che = new RooChebychev("che","che",*obs_var,*coeffList1);

	if(truth)
        {
                //if(cat == 0 || cat == 1 || cat == 2 || cat == 3 || cat == 4) fraction = new RooRealVar("fraction","fraction",0.99,0.5,1.0);
                //else if(cat == 8) fraction = new RooRealVar("fraction","fraction",0.89,0.0,0.99);
		//else fraction = new RooRealVar("fraction","fraction",1.0,0.0,1.0);
        	fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}
        else
        {
                //fraction = new RooRealVar("fraction","fraction",0.95,0.90,0.99);
        	fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}


	initPar.first[4 + order] = "fraction";
        initPar.second[4 + order] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*che,*voigt,*fraction);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;

        return return_var;

}

pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getChebychev(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{
	RooArgList *coeffList1 = new RooArgList();

	pair<string*, double*> initPar;
        initPar.first = new string[order+2];
        initPar.second = new double[order+2];

	initPar.first[0] = "nParameters";
        initPar.second[0] = order;

  	for (int i = 0; i < order; i++)
	{
    		string name = Form("%s_p%d",function_name.c_str(),i);
    		//params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    		RooRealVar *param = NULL;
		if(truth) param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
		else
		{
			param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
			//if(i == 0) param = new RooRealVar(name.c_str(),name.c_str(),-1.5,-3.0,3.0);
			//else param = new RooRealVar(name.c_str(),name.c_str(),0.01,-1.0,1.0);

		}
    		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
    		//params.insert(pair<string,RooRealVar*>(name,param));
    		//prods.insert(pair<string,RooFormulaVar*>(name,form));
    		//coeffList->add(*params[name]);
  		coeffList1->add(*param);
		//coeffList1->add(*form);
		
		initPar.first[1 + i] = name;
                initPar.second[1 + i] = param->getVal();
	}
	
	RooChebychev * che = new RooChebychev("che","che",*obs_var,*coeffList1);
	
	fraction = new RooRealVar("fraction","fraction",1.0,1.0,1.0);

	initPar.first[1 + order] = "fraction";
        initPar.second[1 + order] = fraction->getVal();


        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = che;
        return_var.first.second = fraction;
        return_var.second = initPar;

        return return_var;

}



pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getChebychevCBgaus(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{
	RooArgList *coeffList1 = new RooArgList();

        RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;

        RooRealVar * meanGaus = NULL;
        RooRealVar * sigmaGaus = NULL;

        if(truth)
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

                meanGaus = new RooRealVar("meanGaus","meanGaus",0.0,-0.5,0.5);
                sigmaGaus = new RooRealVar("sigmaGaus","sigmaGaus",1.5,0.0,5.0);
        }
        else
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

                meanGaus = new RooRealVar("meanGaus","meanGaus",0.0,-0.5,0.5);
                sigmaGaus = new RooRealVar("sigmaGaus","sigmaGaus",1.5,0.0,5.0);
        }


	pair<string*, double*> initPar;
        initPar.first = new string[order+6+1+1];
        initPar.second = new double[order+6+1+1];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 6 + order + 1;

        initPar.first[1] = "meanCB";
        initPar.second[1] = meanCB->getVal();
        initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
        initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();
        initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();
        initPar.first[5] = "meanGaus";
        initPar.second[5] = meanGaus->getVal();
        initPar.first[6] = "sigmaGaus";
        initPar.second[6] = sigmaGaus->getVal();

        RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);
        RooGaussian * gaus = new RooGaussian("gaus","gaus",*obs_var,*meanGaus,*sigmaGaus);

        RooFFTConvPdf * conv  = new RooFFTConvPdf("conv","conv",*obs_var,*cb,*gaus);


  	for (int i = 0; i < order; i++)
	{
    		string name = Form("%s_p%d",function_name.c_str(),i);
    		//params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    		RooRealVar *param = NULL;
		if(truth) param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
		else
		{
			param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
			//if(i == 0) param = new RooRealVar(name.c_str(),name.c_str(),-1.5,-3.0,3.0);
			//else param = new RooRealVar(name.c_str(),name.c_str(),0.01,-1.0,1.0);

		}
    		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
    		//params.insert(pair<string,RooRealVar*>(name,param));
    		//prods.insert(pair<string,RooFormulaVar*>(name,form));
    		//coeffList->add(*params[name]);
  		coeffList1->add(*param);
		//coeffList1->add(*form);
		
		initPar.first[7 + i] = name;
                initPar.second[7 + i] = param->getVal();
	}
	
	RooChebychev * che = new RooChebychev("che","che",*obs_var,*coeffList1);

	if(truth)
        {
                if(cat == 7) fraction = new RooRealVar("fraction","fraction",0.90,0.0,1.0);
                //else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
                else fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99);
        }
        else
        {
                fraction = new RooRealVar("fraction","fraction",0.95,0.90,0.99);
        }


	initPar.first[7 + order] = "fraction";
        initPar.second[7 + order] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*che,*conv,*fraction);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;

        return return_var;

}


pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getChebychevCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{
	RooArgList *coeffList1 = new RooArgList();

        RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;

        if(truth)
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

        }
        else
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

        }


	pair<string*, double*> initPar;
        initPar.first = new string[order+4+1+1];
        initPar.second = new double[order+4+1+1];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 4 + order + 1;

        initPar.first[1] = "meanCB";
        initPar.second[1] = meanCB->getVal();
        initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
        initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();
        initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();

        RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);


  	for (int i = 0; i < order; i++)
	{
    		string name = Form("%s_p%d",function_name.c_str(),i);
    		//params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    		RooRealVar *param = NULL;
		if(truth) param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
		else
		{
			param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
			//if(i == 0) param = new RooRealVar(name.c_str(),name.c_str(),-1.5,-3.0,3.0);
			//else param = new RooRealVar(name.c_str(),name.c_str(),0.01,-1.0,1.0);

		}
    		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
    		//params.insert(pair<string,RooRealVar*>(name,param));
    		//prods.insert(pair<string,RooFormulaVar*>(name,form));
    		//coeffList->add(*params[name]);
  		coeffList1->add(*param);
		//coeffList1->add(*form);
		
		initPar.first[5 + i] = name;
                initPar.second[5 + i] = param->getVal();
	}
	
	RooChebychev * che = new RooChebychev("che","che",*obs_var,*coeffList1);

	if(truth)
        {
                //if(cat == 7) fraction = new RooRealVar("fraction","fraction",0.90,0.0,1.0);
                //else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
                //else fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99);
        	fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99);
	}
        else
        {
                //fraction = new RooRealVar("fraction","fraction",0.95,0.90,0.99);
        	fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99);
	}


	initPar.first[5 + order] = "fraction";
        initPar.second[5 + order] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*che,*cb,*fraction);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;

        return return_var;

}


//RooAbsPdf* getExponentialPlusVoigtian(string function_name, int order, RooRealVar *obs_var)
pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getExponentialPlusVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto)
{
	RooMsgService::instance().setGlobalKillBelow(ERROR);
	RooArgList *coeffList1 = new RooArgList();

        RooRealVar *mean = NULL;           
        RooRealVar *sigma = NULL;          
        RooRealVar *width = NULL;

	
	double centerVal0, lowVal0, highVal0, lowVal0_test, highVal0_test;
	double centerVal1, lowVal1, highVal1, lowVal1_test, highVal1_test;
	double centerVal2, lowVal2, highVal2, lowVal2_test, highVal2_test;
	
	if(eVeto == "defaultE")
        {
                if(cat == 0) { centerVal0 = 90.40; centerVal1 = 2.49; centerVal2 = 1.74; }
                if(cat == 1) { centerVal0 = 90.33; centerVal1 = 3.45; centerVal2 = 1.00; }
                if(cat == 2) { centerVal0 = 90.34; centerVal1 = 4.38; centerVal2 = 0.48; }
                if(cat == 3) { centerVal0 = 89.87; centerVal1 = 4.92; centerVal2 = 0.93; }
                if(cat == 4) { centerVal0 = 89.62; centerVal1 = 2.49; centerVal2 = 2.22; }
                if(cat == 5) { centerVal0 = 89.99; centerVal1 = 3.81; centerVal2 = 1.04; }
                if(cat == 6) { centerVal0 = 90.19; centerVal1 = 3.04; centerVal2 = 1.38; }
                if(cat == 7) { centerVal0 = 90.13; centerVal1 = 3.75; centerVal2 = 1.42; }
                if(cat == 8) { centerVal0 = 89.75; centerVal1 = 2.49; centerVal2 = 0.0; } // ~ no entries
                if(cat == 9) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 1.41; }
                if(cat == 10) { centerVal0 = 90.28; centerVal1 = 2.49; centerVal2 = 1.57; }
                if(cat == 11) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 0.01; } //NO entries
                if(cat == 12) { centerVal0 = 90.24; centerVal1 = 3.21; centerVal2 = 1.70; }
                if(cat == 13) { centerVal0 = 90.11; centerVal1 = 5.32; centerVal2 = 0.40; }
        }
        if(eVeto == "notHasPixelSeed")
        {
                if(cat == 0) { centerVal0 = 90.40; centerVal1 = 2.49; centerVal2 = 1.74; }
                if(cat == 1) { centerVal0 = 90.32; centerVal1 = 3.48; centerVal2 = 1.00; }
                if(cat == 2) { centerVal0 = 90.33; centerVal1 = 4.42; centerVal2 = 0.48; }
                if(cat == 3) { centerVal0 = 89.84; centerVal1 = 4.80; centerVal2 = 1.01; }
                if(cat == 4) { centerVal0 = 89.60; centerVal1 = 2.49; centerVal2 = 2.24; }
                if(cat == 5) { centerVal0 = 89.99; centerVal1 = 3.09; centerVal2 = 1.42; }
                if(cat == 6) { centerVal0 = 90.14; centerVal1 = 3.37; centerVal2 = 1.22; }
                if(cat == 7) { centerVal0 = 90.10; centerVal1 = 3.82; centerVal2 = 1.38; }
                if(cat == 8) { centerVal0 = 89.75; centerVal1 = 2.49; centerVal2 = 0.0; } // ~ no entries
                if(cat == 9) { centerVal0 = 90.48; centerVal1 = 2.49; centerVal2 = 1.46; }
                if(cat == 10) { centerVal0 = 90.34; centerVal1 = 2.49; centerVal2 = 1.68; }
                if(cat == 11) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 0.01; } //NO entries
                if(cat == 12) { centerVal0 = 90.28; centerVal1 = 2.72; centerVal2 = 1.57; }
                if(cat == 13) { centerVal0 = 90.05; centerVal1 = 3.90; centerVal2 = 1.38; }
        }

	lowVal0 = centerVal0 - 1.5; highVal0 = centerVal0 + 1.5;
	lowVal1 = centerVal1 - 7.0; highVal1 = centerVal1 + 7.0;
	lowVal2 = centerVal2 - 7.0; highVal2 = centerVal2 + 7.0;

	//if(cat == 1) lowVal1 = 1.0;
        //if(cat == 1) lowVal2 = 0.5;

	lowVal1 = 2.49; //2.49 pour Z theo
	centerVal1 = 2.49;
	highVal1 = 2.49;

	lowVal0_test = centerVal0 - 2.5; highVal0_test = centerVal0 + 2.5;
        lowVal1_test = centerVal1 - 8.0; highVal1_test = centerVal1 + 8.0;
        lowVal2_test = centerVal2 - 8.0; highVal2_test = centerVal2 + 8.0;;

        if(lowVal1 < 0) lowVal1 = 0.0;
        if(lowVal2 < 0) lowVal2 = 0.0;
        if(lowVal1_test < 0) lowVal1_test = 0.0;
        if(lowVal2_test < 0) lowVal2_test = 0.0;

	lowVal0 = lowVal0_test = highVal0 = highVal0_test = centerVal0;
        lowVal1 = lowVal1_test = highVal1 = highVal1_test = centerVal1;
        lowVal2 = lowVal2_test = highVal2 = highVal2_test = centerVal2;
	
	if(truth)
	{	/*
		mean = new RooRealVar("mean","mean",91.188,90.0,92.0);		
		sigma = new RooRealVar("sigma","sigma",2.5,2.5,5.0);
		width = new RooRealVar("width","width",1.5,0.0,9.0);
		*/
		mean = new RooRealVar("mean","mean",centerVal0,lowVal0,highVal0);          
                sigma = new RooRealVar("sigma","sigma",centerVal1,lowVal1,highVal1);
                width = new RooRealVar("width","width",centerVal2,lowVal2,highVal2);
	}
	else
	{
		mean = new RooRealVar("mean","mean",centerVal0,lowVal0_test,highVal0_test);
                sigma = new RooRealVar("sigma","sigma",centerVal1,lowVal1_test,highVal1_test);
                width = new RooRealVar("width","width",centerVal2,0.0,5.0);
	}

	pair<string*, double*> initPar;
        initPar.first = new string[2*order+3+1+1];
        initPar.second = new double[2*order+3+1+1];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 3 + 2*order + 1;

        initPar.first[1] = "mean";
        initPar.second[1] = mean->getVal();
        initPar.first[2] = "sigma";
        initPar.second[2] = sigma->getVal();
        initPar.first[3] = "width";
        initPar.second[3] = width->getVal();


        RooVoigtian * voigt = new RooVoigtian("voigt","voigt",*obs_var,*mean,*sigma,*width);

	coeffList1->add(*obs_var);

	string formula = "";

        for (int i = 1; i <= order; i++)
        {
                string name1 = Form("%s_p1_%d",function_name.c_str(),i);
		string name2 = Form("%s_p2_%d",function_name.c_str(),i);
                //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
		RooRealVar *param1 = NULL;
		//if(i == 1) param1 = new RooRealVar(name1.c_str(),name1.c_str(),-2.0,-10.0,0.0);
		//else param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-1.0,1.0);
		param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-10.0,10.0);
		RooRealVar *param2 = new RooRealVar(name2.c_str(),name2.c_str(),1,0,100);
                //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),start,low,high);
		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
                //params.insert(pair<string,RooRealVar*>(name,param));
                //prods.insert(pair<string,RooFormulaVar*>(name,form));
                //coeffList->add(*params[name]);
                if(i == 1) formula += Form("@%d * exp(@0 * @%d)",2*i,2*i-1);
                else formula += Form(" + @%d * exp(@0 * @%d)",2*i,2*i-1);
                coeffList1->add(*param1);
                coeffList1->add(*param2);
		//coeffList1->add(*form);

		initPar.first[3 + 2*i-1] = name1;
                initPar.second[3 + 2*i-1] = param1->getVal();
		initPar.first[3 + 2*i] = name2;
                initPar.second[3 + 2*i] = param2->getVal();

        }


	cout << endl << "Formula = " << formula << endl; 

	RooGenericPdf * exp = new RooGenericPdf("exp","exp",formula.c_str(),*coeffList1);
	//RooExponential * exp  = new RooExponential("exp","exp",*obs_var,*coeffList1);
	//RooExponentialSum * exp  = new RooExponentialSum("exp","exp",*obs_var,*coeffList1);

	if(truth)
        {
                //if(cat == 0 || cat == 1 || cat == 2 || cat == 3 || cat == 4) fraction = new RooRealVar("fraction","fraction",0.99,0.5,1.0);
		//else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
                //else fraction = new RooRealVar("fraction","fraction",1.0,0.0,1.0);
        	fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}
        else
        {
                //fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
		fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
        }


	initPar.first[3 + 2*order + 1] = "fraction";
        initPar.second[3 + 2*order + 1] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*exp,*voigt,*fraction);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;

        return return_var;


}


pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getExponential(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{
	RooMsgService::instance().setGlobalKillBelow(ERROR);
	RooArgList *coeffList1 = new RooArgList();

	pair<string*, double*> initPar;
        initPar.first = new string[2*order+2];
        initPar.second = new double[2*order+2];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 2*order;

	coeffList1->add(*obs_var);

	string formula = "";

        for (int i = 1; i <= order; i++)
        {
                string name1 = Form("%s_p1_%d",function_name.c_str(),i);
		string name2 = Form("%s_p2_%d",function_name.c_str(),i);
                //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
		RooRealVar *param1 = NULL;
		//if(i == 1) param1 = new RooRealVar(name1.c_str(),name1.c_str(),-2.0,-10.0,0.0);
		//else param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-1.0,1.0);
		param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-10.0,10.0);
		RooRealVar *param2 = new RooRealVar(name2.c_str(),name2.c_str(),1,0,100);
                //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),start,low,high);
		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
                //params.insert(pair<string,RooRealVar*>(name,param));
                //prods.insert(pair<string,RooFormulaVar*>(name,form));
                //coeffList->add(*params[name]);
                if(i == 1) formula += Form("@%d * exp(@0 * @%d)",2*i,2*i-1);
                else formula += Form(" + @%d * exp(@0 * @%d)",2*i,2*i-1);
                coeffList1->add(*param1);
                coeffList1->add(*param2);
		//coeffList1->add(*form);

		initPar.first[2*i-1] = name1;
                initPar.second[2*i-1] = param1->getVal();
		initPar.first[2*i] = name2;
                initPar.second[2*i] = param2->getVal();

        }


	cout << endl << "Formula = " << formula << endl; 

	RooGenericPdf * exp = new RooGenericPdf("exp","exp",formula.c_str(),*coeffList1);
	//RooExponential * exp  = new RooExponential("exp","exp",*obs_var,*coeffList1);
	//RooExponentialSum * exp  = new RooExponentialSum("exp","exp",*obs_var,*coeffList1);

	fraction = new RooRealVar("fraction","fraction",1.0,1.0,1.0);

	initPar.first[2*order + 1] = "fraction";
        initPar.second[2*order + 1] = fraction->getVal();

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = exp;
        return_var.first.second = fraction;
        return_var.second = initPar;

        return return_var;


}






pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getExponentialCBgaus(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{
	RooMsgService::instance().setGlobalKillBelow(ERROR);
	RooArgList *coeffList1 = new RooArgList();

	RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;

        RooRealVar * meanGaus = NULL;
        RooRealVar * sigmaGaus = NULL;

        if(truth)
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

                meanGaus = new RooRealVar("meanGaus","meanGaus",0.0,-0.5,0.5);
                sigmaGaus = new RooRealVar("sigmaGaus","sigmaGaus",1.5,0.0,5.0);
        }
        else
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

                meanGaus = new RooRealVar("meanGaus","meanGaus",0.0,-0.5,0.5);
                sigmaGaus = new RooRealVar("sigmaGaus","sigmaGaus",1.5,0.0,5.0);
        }
	

	pair<string*, double*> initPar;
        initPar.first = new string[2*order+6+1+1];
        initPar.second = new double[2*order+6+1+1];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 6 + 2*order + 1;

        initPar.first[1] = "meanCB";
        initPar.second[1] = meanCB->getVal();
        initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
        initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();
        initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();
        initPar.first[5] = "meanGaus";
        initPar.second[5] = meanGaus->getVal();
        initPar.first[6] = "sigmaGaus";
        initPar.second[6] = sigmaGaus->getVal();

        RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);
        RooGaussian * gaus = new RooGaussian("gaus","gaus",*obs_var,*meanGaus,*sigmaGaus);

        RooFFTConvPdf * conv  = new RooFFTConvPdf("conv","conv",*obs_var,*cb,*gaus);	


	coeffList1->add(*obs_var);

	string formula = "";

        for (int i = 1; i <= order; i++)
        {
                string name1 = Form("%s_p1_%d",function_name.c_str(),i);
		string name2 = Form("%s_p2_%d",function_name.c_str(),i);
                //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
		RooRealVar *param1 = NULL;
		//if(i == 1) param1 = new RooRealVar(name1.c_str(),name1.c_str(),-2.0,-10.0,0.0);
		//else param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-1.0,1.0);
		param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-10.0,10.0);
		RooRealVar *param2 = new RooRealVar(name2.c_str(),name2.c_str(),1,0,100);
                //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),start,low,high);
		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
                //params.insert(pair<string,RooRealVar*>(name,param));
                //prods.insert(pair<string,RooFormulaVar*>(name,form));
                //coeffList->add(*params[name]);
                if(i == 1) formula += Form("@%d * exp(@0 * @%d)",2*i,2*i-1);
                else formula += Form(" + @%d * exp(@0 * @%d)",2*i,2*i-1);
                coeffList1->add(*param1);
                coeffList1->add(*param2);
		//coeffList1->add(*form);

		initPar.first[6 + 2*i-1] = name1;
                initPar.second[6 + 2*i-1] = param1->getVal();
		initPar.first[6 + 2*i] = name2;
                initPar.second[6 + 2*i] = param2->getVal();

        }


	cout << endl << "Formula = " << formula << endl; 

	RooGenericPdf * exp = new RooGenericPdf("exp","exp",formula.c_str(),*coeffList1);
	//RooExponential * exp  = new RooExponential("exp","exp",*obs_var,*coeffList1);
	//RooExponentialSum * exp  = new RooExponentialSum("exp","exp",*obs_var,*coeffList1);

	if(truth)
        {
                if(cat == 7) fraction = new RooRealVar("fraction","fraction",0.90,0.0,1.0);
                //else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
                else fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99);
        }
        else
        {
                fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
        }


	initPar.first[6 + 2*order + 1] = "fraction";
        initPar.second[6 + 2*order + 1] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*exp,*conv,*fraction);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;

        return return_var;


}


pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getExponentialCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{
	RooMsgService::instance().setGlobalKillBelow(ERROR);
	RooArgList *coeffList1 = new RooArgList();

	RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;

        if(truth)
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

        }
        else
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

        }
	

	pair<string*, double*> initPar;
        initPar.first = new string[2*order+4+1+1];
        initPar.second = new double[2*order+4+1+1];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 4 + 2*order + 1;

        initPar.first[1] = "meanCB";
        initPar.second[1] = meanCB->getVal();
        initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
        initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();
        initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();

        RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);

	coeffList1->add(*obs_var);

	string formula = "";

        for (int i = 1; i <= order; i++)
        {
                string name1 = Form("%s_p1_%d",function_name.c_str(),i);
		string name2 = Form("%s_p2_%d",function_name.c_str(),i);
                //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
		RooRealVar *param1 = NULL;
		//if(i == 1) param1 = new RooRealVar(name1.c_str(),name1.c_str(),-2.0,-10.0,0.0);
		//else param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-1.0,1.0);
		param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-10.0,10.0);
		RooRealVar *param2 = new RooRealVar(name2.c_str(),name2.c_str(),1,0,100);
                //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),start,low,high);
		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
                //params.insert(pair<string,RooRealVar*>(name,param));
                //prods.insert(pair<string,RooFormulaVar*>(name,form));
                //coeffList->add(*params[name]);
                if(i == 1) formula += Form("@%d * exp(@0 * @%d)",2*i,2*i-1);
                else formula += Form(" + @%d * exp(@0 * @%d)",2*i,2*i-1);
                coeffList1->add(*param1);
                coeffList1->add(*param2);
		//coeffList1->add(*form);

		initPar.first[4 + 2*i-1] = name1;
                initPar.second[4 + 2*i-1] = param1->getVal();
		initPar.first[4 + 2*i] = name2;
                initPar.second[4 + 2*i] = param2->getVal();

        }


	cout << endl << "Formula = " << formula << endl; 

	RooGenericPdf * exp = new RooGenericPdf("exp","exp",formula.c_str(),*coeffList1);
	//RooExponential * exp  = new RooExponential("exp","exp",*obs_var,*coeffList1);
	//RooExponentialSum * exp  = new RooExponentialSum("exp","exp",*obs_var,*coeffList1);

	if(truth)
        {
                //if(cat == 7) fraction = new RooRealVar("fraction","fraction",0.90,0.0,1.0);
                //else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
                //else fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99);
        	fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}
        else
        {
                //fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
        	fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}


	initPar.first[4 + 2*order + 1] = "fraction";
        initPar.second[4 + 2*order + 1] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*exp,*cb,*fraction);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;

        return return_var;


}



//RooAbsPdf* getPowerLawPlusVoigtian(string function_name, int order, RooRealVar *obs_var)
pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getPowerLawPlusVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto)
{

	RooArgList *coeffList1 = new RooArgList();

        RooRealVar *mean = NULL;           
        RooRealVar *sigma = NULL;          
        RooRealVar *width = NULL;

		
	double centerVal0, lowVal0, highVal0, lowVal0_test, highVal0_test;
	double centerVal1, lowVal1, highVal1, lowVal1_test, highVal1_test;
	double centerVal2, lowVal2, highVal2, lowVal2_test, highVal2_test;
	
	if(eVeto == "defaultE")
        {
                if(cat == 0) { centerVal0 = 90.40; centerVal1 = 2.49; centerVal2 = 1.74; }
                if(cat == 1) { centerVal0 = 90.33; centerVal1 = 3.45; centerVal2 = 1.00; }
                if(cat == 2) { centerVal0 = 90.34; centerVal1 = 4.38; centerVal2 = 0.48; }
                if(cat == 3) { centerVal0 = 89.87; centerVal1 = 4.92; centerVal2 = 0.93; }
                if(cat == 4) { centerVal0 = 89.62; centerVal1 = 2.49; centerVal2 = 2.22; }
                if(cat == 5) { centerVal0 = 89.99; centerVal1 = 3.81; centerVal2 = 1.04; }
                if(cat == 6) { centerVal0 = 90.19; centerVal1 = 3.04; centerVal2 = 1.38; }
                if(cat == 7) { centerVal0 = 90.13; centerVal1 = 3.75; centerVal2 = 1.42; }
                if(cat == 8) { centerVal0 = 89.75; centerVal1 = 2.49; centerVal2 = 0.0; } // ~ no entries
                if(cat == 9) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 1.41; }
                if(cat == 10) { centerVal0 = 90.28; centerVal1 = 2.49; centerVal2 = 1.57; }
                if(cat == 11) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 0.01; } //NO entries
                if(cat == 12) { centerVal0 = 90.24; centerVal1 = 3.21; centerVal2 = 1.70; }
                if(cat == 13) { centerVal0 = 90.11; centerVal1 = 5.32; centerVal2 = 0.40; }
        }
        if(eVeto == "notHasPixelSeed")
        {
                if(cat == 0) { centerVal0 = 90.40; centerVal1 = 2.49; centerVal2 = 1.74; }
                if(cat == 1) { centerVal0 = 90.32; centerVal1 = 3.48; centerVal2 = 1.00; }
                if(cat == 2) { centerVal0 = 90.33; centerVal1 = 4.42; centerVal2 = 0.48; }
                if(cat == 3) { centerVal0 = 89.84; centerVal1 = 4.80; centerVal2 = 1.01; }
                if(cat == 4) { centerVal0 = 89.60; centerVal1 = 2.49; centerVal2 = 2.24; }
                if(cat == 5) { centerVal0 = 89.99; centerVal1 = 3.09; centerVal2 = 1.42; }
                if(cat == 6) { centerVal0 = 90.14; centerVal1 = 3.37; centerVal2 = 1.22; }
                if(cat == 7) { centerVal0 = 90.10; centerVal1 = 3.82; centerVal2 = 1.38; }
                if(cat == 8) { centerVal0 = 89.75; centerVal1 = 2.49; centerVal2 = 0.0; } // ~ no entries
                if(cat == 9) { centerVal0 = 90.48; centerVal1 = 2.49; centerVal2 = 1.46; }
                if(cat == 10) { centerVal0 = 90.34; centerVal1 = 2.49; centerVal2 = 1.68; }
                if(cat == 11) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 0.01; } //NO entries
                if(cat == 12) { centerVal0 = 90.28; centerVal1 = 2.72; centerVal2 = 1.57; }
                if(cat == 13) { centerVal0 = 90.05; centerVal1 = 3.90; centerVal2 = 1.38; }
        }

	lowVal0 = centerVal0 - 1.5; highVal0 = centerVal0 + 1.5;
	lowVal1 = centerVal1 - 7.0; highVal1 = centerVal1 + 7.0;
	lowVal2 = centerVal2 - 7.0; highVal2 = centerVal2 + 7.0;

	lowVal0_test = centerVal0 - 2.5; highVal0_test = centerVal0 + 2.5;
        lowVal1_test = centerVal1 - 8.0; highVal1_test = centerVal1 + 8.0;
        lowVal2_test = centerVal2 - 8.0; highVal2_test = centerVal2 + 8.0;;


	lowVal1 = 2.49; //2.49 pour Z theo
	centerVal1 = 2.49;
	highVal1 = 2.49;

        if(lowVal1 < 0) lowVal1 = 0.0;
        if(lowVal2 < 0) lowVal2 = 0.0;
        if(lowVal1_test < 0) lowVal1_test = 0.0;
        if(lowVal2_test < 0) lowVal2_test = 0.0;

	lowVal0 = lowVal0_test = highVal0 = highVal0_test = centerVal0;
        lowVal1 = lowVal1_test = highVal1 = highVal1_test = centerVal1;
        lowVal2 = lowVal2_test = highVal2 = highVal2_test = centerVal2;
	
	if(truth)
	{	/*
		mean = new RooRealVar("mean","mean",91.188,90.0,92.0);		
		sigma = new RooRealVar("sigma","sigma",2.5,2.5,5.0);
		width = new RooRealVar("width","width",1.5,0.0,9.0);
		*/
		mean = new RooRealVar("mean","mean",centerVal0,lowVal0,highVal0);          
                sigma = new RooRealVar("sigma","sigma",centerVal1,lowVal1,highVal1);
                width = new RooRealVar("width","width",centerVal2,lowVal2,highVal2);
	}
	else
	{
		mean = new RooRealVar("mean","mean",centerVal0,lowVal0_test,highVal0_test);
                sigma = new RooRealVar("sigma","sigma",centerVal1,lowVal1_test,highVal1_test);
                width = new RooRealVar("width","width",centerVal2,0.0,5.0);
	}

	pair<string*, double*> initPar;
        initPar.first = new string[2*order+3+1+1];
        initPar.second = new double[2*order+3+1+1];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 3 + 2* order +1;

        initPar.first[1] = "mean";
        initPar.second[1] = mean->getVal();
        initPar.first[2] = "sigma";
        initPar.second[2] = sigma->getVal();
        initPar.first[3] = "width";
        initPar.second[3] = width->getVal();



        RooVoigtian * voigt = new RooVoigtian("voigt","voigt",*obs_var,*mean,*sigma,*width);

	coeffList1->add(*obs_var);

	string formula = "";

        for (int i = 1; i <= order; i++)
        {
                string name1 = Form("%s_p1_%d",function_name.c_str(),i);
		string name2 = Form("%s_p2_%d",function_name.c_str(),i);
                //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
                RooRealVar *param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-1.0,1.0);
		RooRealVar *param2 = new RooRealVar(name2.c_str(),name2.c_str(),1,0,100);
                //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),start,low,high);
		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
                //params.insert(pair<string,RooRealVar*>(name,param));
                //prods.insert(pair<string,RooFormulaVar*>(name,form));
                //coeffList->add(*params[name]);
		if(i == 1) formula += Form("@%d * pow(@0,@%d)",2*i,2*i-1);
		else formula += Form(" + @%d * pow(@0,@%d)",2*i,2*i-1);
                coeffList1->add(*param1);
		coeffList1->add(*param2);
                //coeffList1->add(*form);
	
		initPar.first[3 + 2*i-1] = name1;
                initPar.second[3 + 2*i-1] = param1->getVal();
                initPar.first[3 + 2*i] = name2;
                initPar.second[3 + 2*i] = param2->getVal();

        }


	//cout << endl << "Formula = " << formula << endl; 

	RooGenericPdf * pow = new RooGenericPdf("pow","pow",formula.c_str(),*coeffList1);

	if(truth)
        {
                //if(cat == 0 || cat == 1 || cat == 2 || cat == 3 || cat == 4) fraction = new RooRealVar("fraction","fraction",0.99,0.5,1.0);
		//else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
                //else fraction = new RooRealVar("fraction","fraction",1.0,0.0,1.0);
        	fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}
        else
        {
                //fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
        	fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}

	initPar.first[3 + 2*order + 1] = "fraction";
        initPar.second[3 + 2*order + 1] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*pow,*voigt,*fraction);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;

        return return_var;

}

pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getPowerLaw(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{

	RooArgList *coeffList1 = new RooArgList();

	pair<string*, double*> initPar;
        initPar.first = new string[2*order+2];
        initPar.second = new double[2*order+2];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 2* order;

	coeffList1->add(*obs_var);

	string formula = "";

        for (int i = 1; i <= order; i++)
        {
                string name1 = Form("%s_p1_%d",function_name.c_str(),i);
		string name2 = Form("%s_p2_%d",function_name.c_str(),i);
                //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
                RooRealVar *param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-1.0,1.0);
		RooRealVar *param2 = new RooRealVar(name2.c_str(),name2.c_str(),1,0,100);
                //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),start,low,high);
		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
                //params.insert(pair<string,RooRealVar*>(name,param));
                //prods.insert(pair<string,RooFormulaVar*>(name,form));
                //coeffList->add(*params[name]);
		if(i == 1) formula += Form("@%d * pow(@0,@%d)",2*i,2*i-1);
		else formula += Form(" + @%d * pow(@0,@%d)",2*i,2*i-1);
                coeffList1->add(*param1);
		coeffList1->add(*param2);
                //coeffList1->add(*form);
	
		initPar.first[2*i-1] = name1;
                initPar.second[2*i-1] = param1->getVal();
                initPar.first[2*i] = name2;
                initPar.second[2*i] = param2->getVal();

        }


	//cout << endl << "Formula = " << formula << endl; 

	RooGenericPdf * pow = new RooGenericPdf("pow","pow",formula.c_str(),*coeffList1);

	fraction = new RooRealVar("fraction","fraction",1.0,1.0,1.0);

	initPar.first[2*order + 1] = "fraction";
        initPar.second[2*order + 1] = fraction->getVal();

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pow;
        return_var.first.second = fraction;
        return_var.second = initPar;

        return return_var;

}





pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getPowerLawCBgaus(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{

	RooArgList *coeffList1 = new RooArgList();

        RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;

        RooRealVar * meanGaus = NULL;
        RooRealVar * sigmaGaus = NULL;

        if(truth)
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

                meanGaus = new RooRealVar("meanGaus","meanGaus",0.0,-0.5,0.5);
                sigmaGaus = new RooRealVar("sigmaGaus","sigmaGaus",1.5,0.0,5.0);
        }
        else
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

                meanGaus = new RooRealVar("meanGaus","meanGaus",0.0,-0.5,0.5);
                sigmaGaus = new RooRealVar("sigmaGaus","sigmaGaus",1.5,0.0,5.0);
        }


	pair<string*, double*> initPar;
        initPar.first = new string[2*order+6+1+1];
        initPar.second = new double[2*order+6+1+1];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 6 + 2* order +1;

        initPar.first[1] = "meanCB";
        initPar.second[1] = meanCB->getVal();
        initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
        initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();
        initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();
        initPar.first[5] = "meanGaus";
        initPar.second[5] = meanGaus->getVal();
        initPar.first[6] = "sigmaGaus";
        initPar.second[6] = sigmaGaus->getVal();

        RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);
        RooGaussian * gaus = new RooGaussian("gaus","gaus",*obs_var,*meanGaus,*sigmaGaus);

        RooFFTConvPdf * conv  = new RooFFTConvPdf("conv","conv",*obs_var,*cb,*gaus);

	coeffList1->add(*obs_var);

	string formula = "";

        for (int i = 1; i <= order; i++)
        {
                string name1 = Form("%s_p1_%d",function_name.c_str(),i);
		string name2 = Form("%s_p2_%d",function_name.c_str(),i);
                //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
                RooRealVar *param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-1.0,1.0);
		RooRealVar *param2 = new RooRealVar(name2.c_str(),name2.c_str(),1,0,100);
                //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),start,low,high);
		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
                //params.insert(pair<string,RooRealVar*>(name,param));
                //prods.insert(pair<string,RooFormulaVar*>(name,form));
                //coeffList->add(*params[name]);
		if(i == 1) formula += Form("@%d * pow(@0,@%d)",2*i,2*i-1);
		else formula += Form(" + @%d * pow(@0,@%d)",2*i,2*i-1);
                coeffList1->add(*param1);
		coeffList1->add(*param2);
                //coeffList1->add(*form);
	
		initPar.first[6 + 2*i-1] = name1;
                initPar.second[6 + 2*i-1] = param1->getVal();
                initPar.first[6 + 2*i] = name2;
                initPar.second[6 + 2*i] = param2->getVal();

        }


	//cout << endl << "Formula = " << formula << endl; 

	RooGenericPdf * pow = new RooGenericPdf("pow","pow",formula.c_str(),*coeffList1);

	if(truth)
        {
                if(cat == 7) fraction = new RooRealVar("fraction","fraction",0.90,0.0,1.0);
                //else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
                else fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99);
        }
        else
        {
                fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
        }

	initPar.first[6 + 2*order + 1] = "fraction";
        initPar.second[6 + 2*order + 1] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*pow,*conv,*fraction);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;

        return return_var;

}


pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getPowerLawCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{

	RooArgList *coeffList1 = new RooArgList();

        RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;

        if(truth)
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

        }
        else
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

        }


	pair<string*, double*> initPar;
        initPar.first = new string[2*order+4+1+1];
        initPar.second = new double[2*order+4+1+1];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 4 + 2* order +1;

        initPar.first[1] = "meanCB";
        initPar.second[1] = meanCB->getVal();
        initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
        initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();
        initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();

        RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);

	coeffList1->add(*obs_var);

	string formula = "";

        for (int i = 1; i <= order; i++)
        {
                string name1 = Form("%s_p1_%d",function_name.c_str(),i);
		string name2 = Form("%s_p2_%d",function_name.c_str(),i);
                //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
                RooRealVar *param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-1.0,1.0);
		RooRealVar *param2 = new RooRealVar(name2.c_str(),name2.c_str(),1,0,100);
                //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),start,low,high);
		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
                //params.insert(pair<string,RooRealVar*>(name,param));
                //prods.insert(pair<string,RooFormulaVar*>(name,form));
                //coeffList->add(*params[name]);
		if(i == 1) formula += Form("@%d * pow(@0,@%d)",2*i,2*i-1);
		else formula += Form(" + @%d * pow(@0,@%d)",2*i,2*i-1);
                coeffList1->add(*param1);
		coeffList1->add(*param2);
                //coeffList1->add(*form);
	
		initPar.first[4 + 2*i-1] = name1;
                initPar.second[4 + 2*i-1] = param1->getVal();
                initPar.first[4 + 2*i] = name2;
                initPar.second[4 + 2*i] = param2->getVal();

        }


	//cout << endl << "Formula = " << formula << endl; 

	RooGenericPdf * pow = new RooGenericPdf("pow","pow",formula.c_str(),*coeffList1);

	if(truth)
        {
                //if(cat == 7) fraction = new RooRealVar("fraction","fraction",0.90,0.0,1.0);
                //else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
                //else fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99);
        	fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}
        else
        {
                //fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
        	fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}

	initPar.first[4 + 2*order + 1] = "fraction";
        initPar.second[4 + 2*order + 1] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*pow,*cb,*fraction);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;

        return return_var;

}


//RooAbsPdf* getLaurentPlusVoigtian(string function_name, int order, RooRealVar *obs_var)
pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getLaurentPlusVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto)
{

	RooArgList *coeffList1 = new RooArgList();

        RooRealVar *mean = NULL;           
        RooRealVar *sigma = NULL;          
        RooRealVar *width = NULL;
		
	double centerVal0, lowVal0, highVal0, lowVal0_test, highVal0_test;
	double centerVal1, lowVal1, highVal1, lowVal1_test, highVal1_test;
	double centerVal2, lowVal2, highVal2, lowVal2_test, highVal2_test;
	
	if(eVeto == "defaultE")
        {
                if(cat == 0) { centerVal0 = 90.40; centerVal1 = 2.49; centerVal2 = 1.74; }
                if(cat == 1) { centerVal0 = 90.33; centerVal1 = 3.45; centerVal2 = 1.00; }
                if(cat == 2) { centerVal0 = 90.34; centerVal1 = 4.38; centerVal2 = 0.48; }
                if(cat == 3) { centerVal0 = 89.87; centerVal1 = 4.92; centerVal2 = 0.93; }
                if(cat == 4) { centerVal0 = 89.62; centerVal1 = 2.49; centerVal2 = 2.22; }
                if(cat == 5) { centerVal0 = 89.99; centerVal1 = 3.81; centerVal2 = 1.04; }
                if(cat == 6) { centerVal0 = 90.19; centerVal1 = 3.04; centerVal2 = 1.38; }
                if(cat == 7) { centerVal0 = 90.13; centerVal1 = 3.75; centerVal2 = 1.42; }
                if(cat == 8) { centerVal0 = 89.75; centerVal1 = 2.49; centerVal2 = 0.0; } // ~ no entries
                if(cat == 9) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 1.41; }
                if(cat == 10) { centerVal0 = 90.28; centerVal1 = 2.49; centerVal2 = 1.57; }
                if(cat == 11) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 0.01; } //NO entries
                if(cat == 12) { centerVal0 = 90.24; centerVal1 = 3.21; centerVal2 = 1.70; }
                if(cat == 13) { centerVal0 = 90.11; centerVal1 = 5.32; centerVal2 = 0.40; }
        }
        if(eVeto == "notHasPixelSeed")
        {
                if(cat == 0) { centerVal0 = 90.40; centerVal1 = 2.49; centerVal2 = 1.74; }
                if(cat == 1) { centerVal0 = 90.32; centerVal1 = 3.48; centerVal2 = 1.00; }
                if(cat == 2) { centerVal0 = 90.33; centerVal1 = 4.42; centerVal2 = 0.48; }
                if(cat == 3) { centerVal0 = 89.84; centerVal1 = 4.80; centerVal2 = 1.01; }
                if(cat == 4) { centerVal0 = 89.60; centerVal1 = 2.49; centerVal2 = 2.24; }
                if(cat == 5) { centerVal0 = 89.99; centerVal1 = 3.09; centerVal2 = 1.42; }
                if(cat == 6) { centerVal0 = 90.14; centerVal1 = 3.37; centerVal2 = 1.22; }
                if(cat == 7) { centerVal0 = 90.10; centerVal1 = 3.82; centerVal2 = 1.38; }
                if(cat == 8) { centerVal0 = 89.75; centerVal1 = 2.49; centerVal2 = 0.0; } // ~ no entries
                if(cat == 9) { centerVal0 = 90.48; centerVal1 = 2.49; centerVal2 = 1.46; }
                if(cat == 10) { centerVal0 = 90.34; centerVal1 = 2.49; centerVal2 = 1.68; }
                if(cat == 11) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 0.01; } //NO entries
                if(cat == 12) { centerVal0 = 90.28; centerVal1 = 2.72; centerVal2 = 1.57; }
                if(cat == 13) { centerVal0 = 90.05; centerVal1 = 3.90; centerVal2 = 1.38; }
        }

	lowVal0 = centerVal0 - 1.5; highVal0 = centerVal0 + 1.5;
	lowVal1 = centerVal1 - 7.0; highVal1 = centerVal1 + 7.0;
	lowVal2 = centerVal2 - 7.0; highVal2 = centerVal2 + 7.0;

	//if(cat == 1) lowVal1 = 1.0;
        //if(cat == 1) lowVal2 = 0.5;
        //if(cat == 10) lowVal1 = 1.0;       
        //if(cat == 10) lowVal2 = 0.5;


	lowVal0_test = centerVal0 - 2.5; highVal0_test = centerVal0 + 2.5;
        lowVal1_test = centerVal1 - 8.0; highVal1_test = centerVal1 + 8.0;
        lowVal2_test = centerVal2 - 8.0; highVal2_test = centerVal2 + 8.0;;	

	lowVal1 = 2.49; //2.49 pour Z theo
	centerVal1 = 2.49;
	highVal1 = 2.49;

	if(lowVal1 < 0) lowVal1 = 0.0;
	if(lowVal2 < 0) lowVal2 = 0.0;
	if(lowVal1_test < 0) lowVal1_test = 0.0;
        if(lowVal2_test < 0) lowVal2_test = 0.0;

	lowVal0 = lowVal0_test = highVal0 = highVal0_test = centerVal0;
        lowVal1 = lowVal1_test = highVal1 = highVal1_test = centerVal1;
        lowVal2 = lowVal2_test = highVal2 = highVal2_test = centerVal2;
	
	if(truth)
	{	/*
		mean = new RooRealVar("mean","mean",91.188,90.0,92.0);		
		sigma = new RooRealVar("sigma","sigma",2.5,2.5,5.0);
		width = new RooRealVar("width","width",1.5,0.0,9.0);
		*/
		mean = new RooRealVar("mean","mean",centerVal0,lowVal0,highVal0);          
                sigma = new RooRealVar("sigma","sigma",centerVal1,lowVal1,highVal1);
                width = new RooRealVar("width","width",centerVal2,lowVal2,highVal2);
	}
	else
	{
		mean = new RooRealVar("mean","mean",centerVal0,lowVal0_test,highVal0_test);
                sigma = new RooRealVar("sigma","sigma",centerVal1,lowVal1_test,highVal1_test);
                width = new RooRealVar("width","width",centerVal2,0.0,5.0);
	}

	pair<string*, double*> initPar;
        initPar.first = new string[2*order+3+1+1];
        initPar.second = new double[2*order+3+1+1];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 3 + order + 1;

        initPar.first[1] = "mean";
        initPar.second[1] = mean->getVal();
        initPar.first[2] = "sigma";
        initPar.second[2] = sigma->getVal();
        initPar.first[3] = "width";
        initPar.second[3] = width->getVal();


        RooVoigtian * voigt = new RooVoigtian("voigt","voigt",*obs_var,*mean,*sigma,*width);

	coeffList1->add(*obs_var);

	string formula = "";

        for (int i = 1; i <= order; i++)
        {
                //string name1 = Form("%s_p1_%d",function_name.c_str(),i);
		string name2 = Form("%s_p2_%d",function_name.c_str(),i);
                //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
                //RooRealVar *param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-1.0,1.0);
		RooRealVar *param2 = new RooRealVar(name2.c_str(),name2.c_str(),1,0,100);
                //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),start,low,high);
		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
                //params.insert(pair<string,RooRealVar*>(name,param));
                //prods.insert(pair<string,RooFormulaVar*>(name,form));
                //coeffList->add(*params[name]);
		
		int expo = -4;
		for(int j = 1; j <= i; j++)
		{
			expo += pow(-1,j) * (j-1); 
		}

		if(i == 1) formula += Form("@%d * pow(@0,%d)",i,expo);
		else formula += Form(" + @%d * pow(@0,%d)",i,expo);
                //coeffList1->add(*param1);
		coeffList1->add(*param2);
                //coeffList1->add(*form);

		initPar.first[3 + i] = name2;
                initPar.second[3 + i] = param2->getVal();

        }


	//cout << endl << "Formula = " << formula << endl; 

	RooGenericPdf * lau = new RooGenericPdf("lau","lau",formula.c_str(),*coeffList1);


	if(truth)
	{
        	//if(cat == 0 || cat == 1 || cat == 2 || cat == 3 || cat == 4) fraction = new RooRealVar("fraction","fraction",0.98,0.0,1.0);
		//else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
		//else fraction = new RooRealVar("fraction","fraction",1.0,0.0,1.0); 
		fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}
	else
	{
		//fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
		fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}


	initPar.first[4 + order] = "fraction";
        initPar.second[4 + order] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*lau,*voigt,*fraction);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;
	
        return return_var;

}

pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getLaurent(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{

	RooArgList *coeffList1 = new RooArgList();

	pair<string*, double*> initPar;
        initPar.first = new string[order+2];
        initPar.second = new double[order+2];

	initPar.first[0] = "nParameters";
        initPar.second[0] = order;

	coeffList1->add(*obs_var);

	string formula = "";

        for (int i = 1; i <= order; i++)
        {
                //string name1 = Form("%s_p1_%d",function_name.c_str(),i);
		string name2 = Form("%s_p2_%d",function_name.c_str(),i);
                //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
                //RooRealVar *param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-1.0,1.0);
		RooRealVar *param2 = new RooRealVar(name2.c_str(),name2.c_str(),1,0,100);
                //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),start,low,high);
		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
                //params.insert(pair<string,RooRealVar*>(name,param));
                //prods.insert(pair<string,RooFormulaVar*>(name,form));
                //coeffList->add(*params[name]);
		
		int expo = -4;
		for(int j = 1; j <= i; j++)
		{
			expo += pow(-1,j) * (j-1); 
		}

		if(i == 1) formula += Form("@%d * pow(@0,%d)",i,expo);
		else formula += Form(" + @%d * pow(@0,%d)",i,expo);
                //coeffList1->add(*param1);
		coeffList1->add(*param2);
                //coeffList1->add(*form);

		initPar.first[i] = name2;
                initPar.second[i] = param2->getVal();

        }


	//cout << endl << "Formula = " << formula << endl; 

	RooGenericPdf * lau = new RooGenericPdf("lau","lau",formula.c_str(),*coeffList1);

	fraction = new RooRealVar("fraction","fraction",1.0,1.0,1.0);

	initPar.first[1 + order] = "fraction";
        initPar.second[1 + order] = fraction->getVal();

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = lau;
        return_var.first.second = fraction;
        return_var.second = initPar;
	
        return return_var;

}




pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getLaurentCBgaus(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{

	RooArgList *coeffList1 = new RooArgList();

        RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;

        RooRealVar * meanGaus = NULL;
        RooRealVar * sigmaGaus = NULL;

        if(truth)
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

                meanGaus = new RooRealVar("meanGaus","meanGaus",0.0,-0.5,0.5);
                sigmaGaus = new RooRealVar("sigmaGaus","sigmaGaus",1.5,0.0,5.0);
        }
        else
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

                meanGaus = new RooRealVar("meanGaus","meanGaus",0.0,-0.5,0.5);
                sigmaGaus = new RooRealVar("sigmaGaus","sigmaGaus",1.5,0.0,5.0);
        }


	pair<string*, double*> initPar;
        initPar.first = new string[2*order+6+1+1];
        initPar.second = new double[2*order+6+1+1];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 6 + order + 1;

        initPar.first[1] = "meanCB";
        initPar.second[1] = meanCB->getVal();
        initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
        initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();
        initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();
        initPar.first[5] = "meanGaus";
        initPar.second[5] = meanGaus->getVal();
        initPar.first[6] = "sigmaGaus";
        initPar.second[6] = sigmaGaus->getVal();

        RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);
        RooGaussian * gaus = new RooGaussian("gaus","gaus",*obs_var,*meanGaus,*sigmaGaus);

        RooFFTConvPdf * conv  = new RooFFTConvPdf("conv","conv",*obs_var,*cb,*gaus);

	coeffList1->add(*obs_var);

	string formula = "";

        for (int i = 1; i <= order; i++)
        {
                //string name1 = Form("%s_p1_%d",function_name.c_str(),i);
		string name2 = Form("%s_p2_%d",function_name.c_str(),i);
                //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
                //RooRealVar *param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-1.0,1.0);
		RooRealVar *param2 = new RooRealVar(name2.c_str(),name2.c_str(),1,0,100);
                //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),start,low,high);
		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
                //params.insert(pair<string,RooRealVar*>(name,param));
                //prods.insert(pair<string,RooFormulaVar*>(name,form));
                //coeffList->add(*params[name]);
		
		int expo = -4;
		for(int j = 1; j <= i; j++)
		{
			expo += pow(-1,j) * (j-1); 
		}

		if(i == 1) formula += Form("@%d * pow(@0,%d)",i,expo);
		else formula += Form(" + @%d * pow(@0,%d)",i,expo);
                //coeffList1->add(*param1);
		coeffList1->add(*param2);
                //coeffList1->add(*form);

		initPar.first[6 + i] = name2;
                initPar.second[6 + i] = param2->getVal();

        }


	//cout << endl << "Formula = " << formula << endl; 

	RooGenericPdf * lau = new RooGenericPdf("lau","lau",formula.c_str(),*coeffList1);


	if(truth)
	{
		if(cat == 7) fraction = new RooRealVar("fraction","fraction",0.90,0.0,1.0);
        	//else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
		else fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99); 
	}
	else
	{
		fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
	}


	initPar.first[7 + order] = "fraction";
        initPar.second[7 + order] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*lau,*conv,*fraction);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;
	
        return return_var;

}



pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getLaurentCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{

	RooArgList *coeffList1 = new RooArgList();

        RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;

        if(truth)
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

        }
        else
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.5,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",10.0,0.0,100.0);

        }


	pair<string*, double*> initPar;
        initPar.first = new string[2*order+4+1+1];
        initPar.second = new double[2*order+4+1+1];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 4 + order + 1;

        initPar.first[1] = "meanCB";
        initPar.second[1] = meanCB->getVal();
        initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
        initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();
        initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();

        RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);

	coeffList1->add(*obs_var);

	string formula = "";

        for (int i = 1; i <= order; i++)
        {
                //string name1 = Form("%s_p1_%d",function_name.c_str(),i);
		string name2 = Form("%s_p2_%d",function_name.c_str(),i);
                //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
                //RooRealVar *param1 = new RooRealVar(name1.c_str(),name1.c_str(),-0.1,-1.0,1.0);
		RooRealVar *param2 = new RooRealVar(name2.c_str(),name2.c_str(),1,0,100);
                //RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),start,low,high);
		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
                //params.insert(pair<string,RooRealVar*>(name,param));
                //prods.insert(pair<string,RooFormulaVar*>(name,form));
                //coeffList->add(*params[name]);
		
		int expo = -4;
		for(int j = 1; j <= i; j++)
		{
			expo += pow(-1,j) * (j-1); 
		}

		if(i == 1) formula += Form("@%d * pow(@0,%d)",i,expo);
		else formula += Form(" + @%d * pow(@0,%d)",i,expo);
                //coeffList1->add(*param1);
		coeffList1->add(*param2);
                //coeffList1->add(*form);

		initPar.first[4 + i] = name2;
                initPar.second[4 + i] = param2->getVal();

        }


	//cout << endl << "Formula = " << formula << endl; 

	RooGenericPdf * lau = new RooGenericPdf("lau","lau",formula.c_str(),*coeffList1);


	if(truth)
	{
		//if(cat == 7) fraction = new RooRealVar("fraction","fraction",0.90,0.0,1.0);
        	//else if(cat == 3) fraction = new RooRealVar("fraction","fraction",0.99,0.99,0.995);
		//else fraction = new RooRealVar("fraction","fraction",0.95,0.0,0.99); 
		fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}
	else
	{
		//fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
		fraction = new RooRealVar("fraction","fraction",0.99,0.0,1.0);
	}


	initPar.first[5 + order] = "fraction";
        initPar.second[5 + order] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*lau,*cb,*fraction);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;
	
        return return_var;

}



pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getPolynomialsPlusVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth, string eVeto)
{
	RooArgList *coeffList1 = new RooArgList();
	
	RooRealVar *mean = NULL; 
        RooRealVar *sigma = NULL; 
        RooRealVar *width = NULL;

	double centerVal0, lowVal0, highVal0, lowVal0_test, highVal0_test;
	double centerVal1, lowVal1, highVal1, lowVal1_test, highVal1_test;
	double centerVal2, lowVal2, highVal2, lowVal2_test, highVal2_test;
	
	if(eVeto == "defaultE")
	{
		if(cat == 0) { centerVal0 = 90.40; centerVal1 = 2.49; centerVal2 = 1.74; } 
		if(cat == 1) { centerVal0 = 90.33; centerVal1 = 3.45; centerVal2 = 1.00; }
		if(cat == 2) { centerVal0 = 90.34; centerVal1 = 4.38; centerVal2 = 0.48; }
		if(cat == 3) { centerVal0 = 89.87; centerVal1 = 4.92; centerVal2 = 0.93; }
		if(cat == 4) { centerVal0 = 89.62; centerVal1 = 2.49; centerVal2 = 2.22; }
		if(cat == 5) { centerVal0 = 89.99; centerVal1 = 3.81; centerVal2 = 1.04; }
		if(cat == 6) { centerVal0 = 90.19; centerVal1 = 3.04; centerVal2 = 1.38; }
		if(cat == 7) { centerVal0 = 90.13; centerVal1 = 3.75; centerVal2 = 1.42; }
		if(cat == 8) { centerVal0 = 89.75; centerVal1 = 2.49; centerVal2 = 0.0; } // ~ no entries
		if(cat == 9) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 1.41; }
		if(cat == 10) { centerVal0 = 90.28; centerVal1 = 2.49; centerVal2 = 1.57; }
		if(cat == 11) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 0.01; } //NO entries
		if(cat == 12) { centerVal0 = 90.24; centerVal1 = 3.21; centerVal2 = 1.70; }
		if(cat == 13) { centerVal0 = 90.11; centerVal1 = 5.32; centerVal2 = 0.40; }	
	}
	if(eVeto == "notHasPixelSeed")
	{
		if(cat == 0) { centerVal0 = 90.40; centerVal1 = 2.49; centerVal2 = 1.74; } 
	        if(cat == 1) { centerVal0 = 90.32; centerVal1 = 3.48; centerVal2 = 1.00; }
	        if(cat == 2) { centerVal0 = 90.33; centerVal1 = 4.42; centerVal2 = 0.48; }
	        if(cat == 3) { centerVal0 = 89.84; centerVal1 = 4.80; centerVal2 = 1.01; }
	        if(cat == 4) { centerVal0 = 89.60; centerVal1 = 2.49; centerVal2 = 2.24; }
	        if(cat == 5) { centerVal0 = 89.99; centerVal1 = 3.09; centerVal2 = 1.42; }
	        if(cat == 6) { centerVal0 = 90.14; centerVal1 = 3.37; centerVal2 = 1.22; }
	        if(cat == 7) { centerVal0 = 90.10; centerVal1 = 3.82; centerVal2 = 1.38; }
	        if(cat == 8) { centerVal0 = 89.75; centerVal1 = 2.49; centerVal2 = 0.0; } // ~ no entries
	        if(cat == 9) { centerVal0 = 90.48; centerVal1 = 2.49; centerVal2 = 1.46; }
	        if(cat == 10) { centerVal0 = 90.34; centerVal1 = 2.49; centerVal2 = 1.68; }
	        if(cat == 11) { centerVal0 = 90.25; centerVal1 = 2.49; centerVal2 = 0.01; } //NO entries
	        if(cat == 12) { centerVal0 = 90.28; centerVal1 = 2.72; centerVal2 = 1.57; }
	        if(cat == 13) { centerVal0 = 90.05; centerVal1 = 3.90; centerVal2 = 1.38; }
	}

	lowVal0 = centerVal0 - 1.5; highVal0 = centerVal0 + 1.5;
	lowVal1 = centerVal1 - 7.0; highVal1 = centerVal1 + 7.0;
	lowVal2 = centerVal2 - 7.0; highVal2 = centerVal2 + 7.0;

	if(cat == 1) lowVal1 = 1.0;
	if(cat == 1) lowVal2 = 0.5;

	lowVal0_test = centerVal0 - 2.5; highVal0_test = centerVal0 + 2.5;
        lowVal1_test = centerVal1 - 8.0; highVal1_test = centerVal1 + 8.0;
        lowVal2_test = centerVal2 - 8.0; highVal2_test = centerVal2 + 8.0;

        if(lowVal1 < 0) lowVal1 = 0.0;
        if(lowVal2 < 0) lowVal2 = 0.0;
        if(lowVal1_test < 0) lowVal1_test = 0.0;
        if(lowVal2_test < 0) lowVal2_test = 0.0;

	lowVal0 = lowVal0_test = highVal0 = highVal0_test = centerVal0;
	lowVal1 = lowVal1_test = highVal1 = highVal1_test = centerVal1;
	lowVal2 = lowVal2_test = highVal2 = highVal2_test = centerVal2;

	
	if(truth)
	{	/*
		mean = new RooRealVar("mean","mean",91.188,90.0,92.0);		
		sigma = new RooRealVar("sigma","sigma",2.5,2.5,5.0);
		width = new RooRealVar("width","width",1.5,0.0,9.0);
		*/
		mean = new RooRealVar("mean","mean",centerVal0,lowVal0,highVal0);          
                sigma = new RooRealVar("sigma","sigma",centerVal1,lowVal1,highVal1);
                width = new RooRealVar("width","width",centerVal2,lowVal2,highVal2);
	}
	else
	{
		mean = new RooRealVar("mean","mean",centerVal0,lowVal0_test,highVal0_test);
                sigma = new RooRealVar("sigma","sigma",centerVal1,lowVal1_test,highVal1_test);
                width = new RooRealVar("width","width",centerVal2,0.0,5.0);
	}

	pair<string*, double*> initPar;
        initPar.first = new string[order+3+1+1];
        initPar.second = new double[order+3+1+1];

	initPar.first[0] = "nParameters";
	initPar.second[0] = 3 + order + 1;

	initPar.first[1] = "mean";
	initPar.second[1] = mean->getVal();
	initPar.first[2] = "sigma";
        initPar.second[2] = sigma->getVal();
	initPar.first[3] = "width";
        initPar.second[3] = width->getVal();	


	RooVoigtian * voigt = new RooVoigtian("voigt","voigt",*obs_var,*mean,*sigma,*width);

  	for (int i = 0; i < order; i++)
	{
    		string name = Form("%s_p%d",function_name.c_str(),i);
    		//params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    		RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01,-100.,100.);
    		//RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0",RooArgList(*param));
    		//params.insert(pair<string,RooRealVar*>(name,param));
    		//prods.insert(pair<string,RooFormulaVar*>(name,form));
    		//coeffList->add(*params[name]);
  		coeffList1->add(*param);
		//coeffList1->add(*form);

		initPar.first[4 + i] = name;
        	initPar.second[4 + i] = param->getVal();

	}
	
	RooPolynomial * poly = new RooPolynomial("poly", "poly", *obs_var, *coeffList1, 1);

	cout << endl << "HERE !!!! fct name = " << function_name << endl;

	if(truth)
        {
                if(cat == 0 || cat == 1 || cat == 2 || cat == 3 || cat == 4) fraction = new RooRealVar("fraction","fraction",0.99,0.5,1.0);
                else if(cat == 8) fraction = new RooRealVar("fraction","fraction",0.89,0.0,0.99);
                else fraction = new RooRealVar("fraction","fraction",0.90,0.0,1.0);
        }
        else
        {
                fraction = new RooRealVar("fraction","fraction",0.95,0.0,1.0);
        }

	initPar.first[4 + order] = "fraction";
        initPar.second[4 + order] = fraction->getVal();


	RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*poly,*voigt,*fraction);

	pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
	return_var.first.first = pdf;
	return_var.first.second = fraction;
	return_var.second = initPar;

  	return return_var;
}


pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getVoigtian(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{
	cout << endl << "coucou in getVoigtian : " << endl;
        RooRealVar *mean = NULL;           
        RooRealVar *sigma = NULL;          
        RooRealVar *width = NULL;

        if(truth)
        {
                mean = new RooRealVar("mean","mean",91.188,88.0,94.0);          
                sigma = new RooRealVar("sigma","sigma",2.49,2.49,9.0);
                width = new RooRealVar("width","width",1.5,0.0,9.0);
        }
        else
        {
                mean = new RooRealVar("mean","mean",91.188,88.0,94.0);
                sigma = new RooRealVar("sigma","sigma",2.49,2.49,9.0);
                width = new RooRealVar("width","width",1.5,0.1,9.0);

        }

	pair<string*, double*> initPar;
        initPar.first = new string[4];
        initPar.second = new double[4];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 3 + 1;

        initPar.first[1] = "mean";
        initPar.second[1] = mean->getVal();
        initPar.first[2] = "sigma";
        initPar.second[2] = sigma->getVal();
        initPar.first[3] = "width";
        initPar.second[3] = width->getVal();

	cout << endl << "coucou in getVoigtian 2 : " << endl;

        RooVoigtian * voigt = new RooVoigtian("voigt","voigt",*obs_var,*mean,*sigma,*width);

	cout << endl << "coucou in getVoigtian 3 : " << endl;
        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = voigt;
        return_var.first.second = fraction;
        return_var.second = initPar;
	
	cout << endl << "coucou in getVoigtian 4: " << endl;
        return return_var;

}

pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{


        RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;

        if(truth)
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.0,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",2.7,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",2.0,0.0,100.0);
        }
        else
        {
		meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.0,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",2.7,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",2.0,0.0,100.0);
        }


	pair<string*, double*> initPar;
        initPar.first = new string[5];
        initPar.second = new double[5];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 4 + 1;

        initPar.first[1] = "meanCB";
        initPar.second[1] = meanCB->getVal();
        initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
        initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();
        initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();

        RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = cb;
        return_var.first.second = fraction;
        return_var.second = initPar;
	
        return return_var;

}




pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getCBgaus(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{


        RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;

        RooRealVar * meanGaus = NULL;
        RooRealVar * sigmaGaus = NULL;

        if(truth)
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.0,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",2.7,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",2.0,0.0,100.0);

                meanGaus = new RooRealVar("meanGaus","meanGaus",0.0,-0.5,0.5);
                sigmaGaus = new RooRealVar("sigmaGaus","sigmaGaus",0.35,0.0,5.0);
        }
        else
        {
		meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.0,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",2.7,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",2.0,0.0,100.0);

                meanGaus = new RooRealVar("meanGaus","meanGaus",0.0,-0.5,0.5);
                sigmaGaus = new RooRealVar("sigmaGaus","sigmaGaus",0.35,0.0,5.0);	

        }


	pair<string*, double*> initPar;
        initPar.first = new string[7];
        initPar.second = new double[7];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 6 + 1;

        initPar.first[1] = "meanCB";
        initPar.second[1] = meanCB->getVal();
        initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
        initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();
        initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();
        initPar.first[5] = "meanGaus";
        initPar.second[5] = meanGaus->getVal();
        initPar.first[6] = "sigmaGaus";
        initPar.second[6] = sigmaGaus->getVal();

        RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);
        RooGaussian * gaus = new RooGaussian("gaus","gaus",*obs_var,*meanGaus,*sigmaGaus);

        RooFFTConvPdf * conv  = new RooFFTConvPdf("conv","conv",*obs_var,*cb,*gaus);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = conv;
        return_var.first.second = fraction;
        return_var.second = initPar;
	
        return return_var;

}

pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getCBpCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{
	cout << endl << "coucou in getCBpCB 1" << endl;

        RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;

	RooRealVar * meanCB2 = NULL;
        RooRealVar * sigmaCB2 = NULL;
        RooRealVar * alphaCB2 = NULL;
        RooRealVar * nCB2  = NULL;
	cout << endl << "coucou in getCBpCB 2" << endl;

        if(truth)
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.0,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",5.0,0.0,1000.0);

		meanCB2 = new RooRealVar("meanCB2","meanCB2",91.188,84.0,98.0);
                sigmaCB2 = new RooRealVar("sigmaCB2","sigmaCB2",3.0,0.0,5.0);
                alphaCB2 = new RooRealVar("alphaCB2","alphaCB2",-0.35,-3.0,-0.1);
                nCB2 = new RooRealVar("nCB2","nCB2",6.0,0.0,1000.0);



        }
        else
        {
		meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.0,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.35,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",5.0,0.0,1000.0);

                meanCB2 = new RooRealVar("meanCB2","meanCB2",91.188,84.0,98.0);
                sigmaCB2 = new RooRealVar("sigmaCB2","sigmaCB2",3.0,0.0,5.0);
                alphaCB2 = new RooRealVar("alphaCB2","alphaCB2",-0.35,-3.0,-0.1);
                nCB2 = new RooRealVar("nCB2","nCB2",6.0,0.0,1000.0);		
        }

	cout << endl << "coucou in getCBpCB 3" << endl;

	pair<string*, double*> initPar;
        initPar.first = new string[10];
        initPar.second = new double[10];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 8 + 1 + 1;

        initPar.first[1] = "meanCB";
        initPar.second[1] = meanCB->getVal();
        initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
        initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();
        initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();
	initPar.first[5] = "meanCB2";
        initPar.second[5] = meanCB2->getVal();
        initPar.first[6] = "sigmaCB2";
        initPar.second[6] = sigmaCB2->getVal();
        initPar.first[7] = "alphaCB2";
        initPar.second[7] = alphaCB2->getVal();
        initPar.first[8] = "nCB2";
        initPar.second[8] = nCB2->getVal();

	cout << endl << "coucou in getCBpCB 4" << endl;



        RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);
	RooCBShape * cb2 = new RooCBShape("cb2","cb2",*obs_var,*meanCB2,*sigmaCB2,*alphaCB2,*nCB2); //FIXME

	fraction = new RooRealVar("fraction","fraction",0.97,0.80,1.0);

        initPar.first[9] = "fraction";
        initPar.second[9] = fraction->getVal();


        RooAddPdf *pdf = new RooAddPdf(function_name.c_str(),function_name.c_str(),*cb,*cb2,*fraction);



        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = pdf;
        return_var.first.second = fraction;
        return_var.second = initPar;
	
        return return_var;

}

pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > getCBconvCB(string function_name, int order, RooRealVar *obs_var, RooRealVar* &fraction, int cat, bool truth)
{

        RooRealVar * meanCB = NULL;
        RooRealVar * sigmaCB = NULL;
        RooRealVar * alphaCB = NULL;
        RooRealVar * nCB  = NULL;

	RooRealVar * meanCB2 = NULL;
        RooRealVar * sigmaCB2 = NULL;
        RooRealVar * alphaCB2 = NULL;
        RooRealVar * nCB2  = NULL;

        if(truth)
        {
                meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.44,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.40,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",2.0,0.0,100.0);

		meanCB2 = new RooRealVar("meanCB2","meanCB2",0.0,-3.0,3.0);
                sigmaCB2 = new RooRealVar("sigmaCB2","sigmaCB2",0.11,0.0,5.0);
                alphaCB2 = new RooRealVar("alphaCB2","alphaCB2",-0.85,-3.0,0.0);
                nCB2 = new RooRealVar("nCB2","nCB2",1.35,0.0,100.0);



        }
        else
        {
		meanCB = new RooRealVar("meanCB","meanCB",90.0,88.0,92.0);
                sigmaCB = new RooRealVar("sigmaCB","sigmaCB",2.44,0.0,5.0);
                alphaCB = new RooRealVar("alphaCB","alphaCB",1.40,0.0,3.0);
                nCB = new RooRealVar("nCB","nCB",2.0,0.0,100.0);

                meanCB2 = new RooRealVar("meanCB2","meanCB2",0.0,-3.0,3.0);
                sigmaCB2 = new RooRealVar("sigmaCB2","sigmaCB2",0.11,0.0,5.0);
                alphaCB2 = new RooRealVar("alphaCB2","alphaCB2",-0.85,-3.0,0.0);
                nCB2 = new RooRealVar("nCB2","nCB2",1.35,0.0,100.0);	

        }


	pair<string*, double*> initPar;
        initPar.first = new string[10];
        initPar.second = new double[10];

	initPar.first[0] = "nParameters";
        initPar.second[0] = 8 + 1;

        initPar.first[1] = "meanCB";
        initPar.second[1] = meanCB->getVal();
        initPar.first[2] = "sigmaCB";
        initPar.second[2] = sigmaCB->getVal();
        initPar.first[3] = "alphaCB";
        initPar.second[3] = alphaCB->getVal();
        initPar.first[4] = "nCB";
        initPar.second[4] = nCB->getVal();
	initPar.first[5] = "meanCB2";
        initPar.second[5] = meanCB2->getVal();
        initPar.first[6] = "sigmaCB2";
        initPar.second[6] = sigmaCB2->getVal();
        initPar.first[7] = "alphaCB2";
        initPar.second[7] = alphaCB2->getVal();
        initPar.first[8] = "nCB2";
        initPar.second[8] = nCB2->getVal();




        RooCBShape * cb = new RooCBShape("cb","cb",*obs_var,*meanCB,*sigmaCB,*alphaCB,*nCB);
	RooCBShape * cb2 = new RooCBShape("cb2","cb2",*obs_var,*meanCB2,*sigmaCB2,*alphaCB2,*nCB2);
	
	RooFFTConvPdf * conv  = new RooFFTConvPdf("conv","conv",*obs_var,*cb,*cb2);

        pair< pair<RooAbsPdf*,RooRealVar*>,pair<string*, double*> > return_var;
        return_var.first.first = conv;
        return_var.first.second = fraction;
        return_var.second = initPar;
	
        return return_var;

}





Double_t chiSquare(RooPlot* plot, char* pdfname, char* histname) //int* fewBins 
{
  // Calculate the chi^2/NDOF of this curve with respect to the histogram
  // 'hist' accounting nFitParam floating parameters in case the curve
  // was the result of a fit

  // Find curve object
  RooCurve* curve = (RooCurve*) plot->findObject(pdfname, RooCurve::Class());
  //RooCurve* curve = plot->getCurve(pdfname);  
  //curve->Print();

  if (!curve) {
    cout<<endl << "cit::RooChi2Calculator(plotname=" << plot->GetName()
         << ")::chiSquare(..) cannot find curve" << endl ;
    return 0 ;
  }

  // Find histogram object
  RooHist* hist = (RooHist*) plot->findObject(histname, RooHist::Class()) ;
  //RooHist* hist = plot->getHist(histname);

  if (!hist) {
    cout<<endl << "cit::RooChi2Calculator(plotname=" << plot->GetName()
         << ")::chiSquare(..) cannot find histogram" << endl ;
    return 0 ;
  }


  Int_t i,np = hist->GetN() ;
  Double_t x,y,/*eyl,eyh,*/ xl,xh ;

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

  Double_t chisq(0) ;
  for (i=0 ; i<np ; i++) {   

    // Retrieve histogram contents
    hist->GetPoint(i,x,y) ;
    xl = x - hist->GetEXlow()[i] ;
    xh = x + hist->GetEXhigh()[i] ;
    // eyl = hist->GetEYlow()[i] ;
    // eyh = hist->GetEYhigh()[i] ;

    // Check if the whole bin is in range of curve
    if (xl < xstart || xstop < xh) continue ;

    //if(y != 0 && y < 35.0)
    //if(y == 0)
    if(y < 5)
    {
    	cout<<endl<<"Too few entries : "<<y<<" in the bin : "<<i<<"  >> Need to reduce the binning for the p-value calculation!"<<endl;
	//*fewBins = 1;
	//break;
	continue;
	
    }
    //else *fewBins = 0;

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

    if (avg < 5.) {
      cout << "cit::RooChi2Calculator(plotname=" << plot->GetName()
			    << ")::chiSquare(..) expectation in bin "
			    << i << " is " << avg << " < 5!" << endl ;
    }

    // JV: Use the expected number of events for the y uncertainty,
    // See (33.34) of http://pdg.lbl.gov/2011/reviews/rpp2011-rev-statistics.pdf

	
    // Add pull^2 to chisq
    if (avg != 0) {      
      Double_t resid = y - avg;
      chisq += (resid * resid / avg) ;
    }
  }
/*
  int nFitParam = fitParameters.size() / 2;
  fitParameters.push_back(chisq / (nbin - nFitParam));
  fitParameters.push_back(nbin - nFitParam);
  fitParameters.push_back(TMath::Prob(chisq, nbin - nFitParam));
*/

  //return chisq / (nbin - nFitParam) ;
   return chisq ;

}

