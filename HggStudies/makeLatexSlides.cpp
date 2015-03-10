// !!!!! the order of list of polynomials in your dat file should match with the one in TruthChoice.cpp !!!! 

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
	

	string configFile = "config.dat";
	string truthDirectory = "truth"; 
	string biasDirectory = "bias";
	int ncats = 14;
	int catStart = 0;
	string output = "Latex.txt";

	if( argc > 1 ) configFile = argv[1];
	if( argc > 2 ) truthDirectory = argv[2];
	if( argc > 3 ) biasDirectory = argv[3];
	if( argc > 4 ) 
        {
                std::stringstream ss ( argv[4] );
                ss >> ncats;
        }
	if( argc > 5 ) 
        {
                std::stringstream ss ( argv[5] );
                ss >> catStart;
        }	
	if( argc > 6 ) output = argv[6];
  
	string * labels;

	bool cic = false;
        if(ncats > 15 ) cic = true;

	string ciclabels[16] = { 
                "Untagged 0 (High P_{T})"
                ,"Untagged 1 (High P_{T})"
                ,"Untagged 2 (High P_{T})"
                ,"Untagged 3 (High P_{T})"
                ,"Untagged 4 (Low P_{T})"
                ,"Untagged 5 (Low P_{T})"
                ,"Untagged 6 (Low P_{T})"
                ,"Untagged 7 (Low P_{T})"
                ,"Dijet Tag 0"
                ,"Dijet Tag 1"
                ,"VH Lepton Tight"
                ,"VH Lepton Loose"
                ,"VH MET Tag"
                ,"t$\\bar{t}$H Lepton Tag"
                ,"t$\\bar{t}$H Multijet Tag" 
                ,"VH Dijet Tag"
        };

	string massfactlabels[14] = { 
                "Untagged 0 btd > 0.76"
                ,"Untagged 1 btd > 0.36"
                ,"Untagged 2 btd > 0.0"
                ,"Untagged 3 btd > -0.42"
                ,"Untagged 4 bdt > -0.78"
                ,"Dijet Tag 0"
                ,"Dijet Tag 1"
                ,"Dijet Tag 2"
                ,"VH Lepton Tight"
                ,"VH Lepton Loose"
                ,"VH MET Tag"
                ,"t$\\bar{t}$H Tag"
                ,"t$\\bar{t}$H Multijet Tag"
                ,"VH Dijet Tag"
        };

	string analysisName = "";

        if( cic ) 
	{ 
		labels = ciclabels; 
		analysisName = "Cutbased";
	}
        else 
	{ 
		labels = massfactlabels; 
		analysisName = "MassFactorizedMVA";
	}
 
	pair< vector< vector<string> >,vector< vector<int> > > truthVector;
	pair< vector< vector<string> >,vector< vector<int> > > testVector; 

	string tempChain = "";
	int tempOrder = 0;
	string tempLine;

	fstream fconfigFile;
	fconfigFile.open(configFile.c_str(),ios::in);

	fconfigFile >> tempChain;

	for(int i = catStart; i < ncats; i++)
	{
		cout << endl << "i = " << i;
		string chainNext = Form("cat%d",i+1);
		pair< vector<string>, vector<int> > tempTruthVector;
		pair< vector<string>, vector<int> > tempTestVector;

		fconfigFile >> tempChain;
		while(true)
                {
			fconfigFile >> tempChain;
			if(tempChain == "TestFunctions") break;		 	
			fconfigFile >> tempOrder;
			tempTruthVector.first.push_back(tempChain);
			tempTruthVector.second.push_back(tempOrder);
			cout << endl << "tempChain = " <<tempChain;
			cout << endl << "tempOrder = " <<tempOrder;

		}
		while(! fconfigFile.eof())
                {
                        fconfigFile >> tempChain;
                        if(tempChain == chainNext) break;
                        fconfigFile >> tempOrder;
                        tempTestVector.first.push_back(tempChain);
                        tempTestVector.second.push_back(tempOrder);
			cout << endl << "tempChain = " <<tempChain;
                        cout << endl << "tempOrder = " <<tempOrder;

                }

		cout << endl << "tempTruthVector.first.size() = " << tempTruthVector.first.size();
		cout << endl << "tempTestVector.first.size() = " << tempTestVector.first.size();

		truthVector.first.push_back(tempTruthVector.first);
		truthVector.second.push_back(tempTruthVector.second);
		testVector.first.push_back(tempTestVector.first);
                testVector.second.push_back(tempTestVector.second);

		cout << endl << "truthVector.first[i-catStart-catStart].size() = " << truthVector.first[i-catStart-catStart].size();
		cout << endl << "testVector.first[i-catStart-catStart].size() = " << testVector.first[i-catStart-catStart].size();


	}	


	fconfigFile.close();



	system(Form("rm %s",output.c_str()));
	fstream outFile;
	outFile.open(output.c_str(),ios::out);

	for(int i = catStart; i < ncats; i++)
        {
		cout << endl << "i = " << i;
		outFile << "\\begin{frame}{Truth model determination and Bias study : Results}" << endl;
		outFile << "\\begin{center} \\scriptsize" << endl;
		outFile << "\\renewcommand{\\arraystretch}{1.2}" << endl;
		outFile << "\\begin{tabular}{|c|c|c|c|c|}" << endl;
		outFile << "   \\hline" << endl;
		outFile << Form("   \\multicolumn{5}{|c|}{\\textbf{%s Class %d (%s)}} \\\\",analysisName.c_str(),i,labels[i-catStart].c_str()) << endl; 
		outFile << "   \\hline" << endl;
		outFile << "   \\multicolumn{5}{|c|}{Test for suitable Truth Models} \\\\" << endl;
		outFile << "   \\hline" << endl;
		outFile << "   Truth Model & df & $NLL_{N}$ & $\\chi^{2}(\\Delta NLL_{N+1})$ & $p( \\chi^{2} > \\chi^{2}(\\Delta NLL_{N+1}))$ \\\\" << endl;
		outFile << "   \\hline" << endl;
		fstream truthFile;
		truthFile.open(Form("%s/Chosen/Chosen_cat%d.txt",truthDirectory.c_str(),i),ios::in);
		cout << endl << Form("%s/Chosen/Chosen_cat%d.txt",truthDirectory.c_str(),i);
		for(int j = 0; j < truthVector.first[i-catStart].size(); j++)
		{
			cout << endl << "j = " << j;
			cout << endl << "truthVector.first[i-catStart].size() = " <<truthVector.first[i-catStart].size();
			getline(truthFile,tempLine);
			cout << endl << "tempLine" << tempLine;
			getline(truthFile,tempLine);
			cout << endl << "tempLine" << tempLine;
			if(truthVector.first[i-catStart][j] == "BernsteinPlusVoigtian") outFile << truthVector.second[i-catStart][j] << "Ber + Voi &"<< truthVector.second[i-catStart][j] + 1 <<" + 3 & "<< tempLine <<"\\\\" <<endl;
			if(truthVector.first[i-catStart][j] == "ChebychevPlusVoigtian") outFile << truthVector.second[i-catStart][j] << "Che + Voi &"<< truthVector.second[i-catStart][j] + 1 <<" + 3 & "<< tempLine <<"\\\\" <<endl;
			if(truthVector.first[i-catStart][j] == "LaurentPlusVoigtian") outFile << truthVector.second[i-catStart][j] << "Lau + Voi &"<< truthVector.second[i-catStart][j] + 1 <<" + 3 & "<< tempLine <<"\\\\" <<endl;
			if(truthVector.first[i-catStart][j] == "ExponentialPlusVoigtian") outFile << truthVector.second[i-catStart][j] << "Exp + Voi &"<< truthVector.second[i-catStart][j] * 2  <<" + 3 & "<< tempLine <<"\\\\" <<endl;
			if(truthVector.first[i-catStart][j] == "PowerLawPlusVoigtian") outFile << truthVector.second[i-catStart][j] << "Pow + Voi &"<< truthVector.second[i-catStart][j] * 2 <<" + 3 & "<< tempLine <<"\\\\" <<endl;	
			
			if(truthVector.first[i-catStart][j] == "Bernstein") outFile << truthVector.second[i-catStart][j] << "Ber &"<< truthVector.second[i-catStart][j] + 1 <<" & "<< tempLine <<"\\\\" <<endl; 
                        if(truthVector.first[i-catStart][j] == "Chebychev") outFile << truthVector.second[i-catStart][j] << "Che &"<< truthVector.second[i-catStart][j] + 1 <<" & "<< tempLine <<"\\\\" <<endl; 
                        if(truthVector.first[i-catStart][j] == "Laurent") outFile << truthVector.second[i-catStart][j] << "Lau &"<< truthVector.second[i-catStart][j] + 1 <<" & "<< tempLine <<"\\\\" <<endl; 
                        if(truthVector.first[i-catStart][j] == "Exponential") outFile << truthVector.second[i-catStart][j] << "Exp &"<< truthVector.second[i-catStart][j] * 2  <<" & "<< tempLine <<"\\\\" <<endl; 
                        if(truthVector.first[i-catStart][j] == "PowerLaw") outFile << truthVector.second[i-catStart][j] << "Pow &"<< truthVector.second[i-catStart][j] * 2 <<" & "<< tempLine <<"\\\\" <<endl;
	
		}
		truthFile.close();
		outFile << "   \\hline" << endl;
		outFile << "\\end{tabular}" << endl;
		if(truthVector.first[i-catStart].size() > 4) outFile << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
		else outFile << "\\begin{tabular}{|c|c|c|c|c|}" << endl;
                outFile << "   \\hline" << endl;	
		if(truthVector.first[i-catStart].size() > 4) outFile << "   \\multicolumn{6}{|c|}{Test for maximum Bias} \\\\" << endl;
		else outFile << "   \\multicolumn{5}{|c|}{Test for maximum Bias} \\\\" << endl;
	        outFile << "   \\hline" << endl;
		outFile << "Fit Model";
		for(int j = 0; j < truthVector.first[i-catStart].size(); j++)
                {
			if(truthVector.first[i-catStart][j] == "BernsteinPlusVoigtian") outFile << " & " <<truthVector.second[i-catStart][j] <<"Ber + Voi ";
			if(truthVector.first[i-catStart][j] == "ChebychevPlusVoigtian") outFile << " & " <<truthVector.second[i-catStart][j] <<"Che + Voi ";
			if(truthVector.first[i-catStart][j] == "LaurentPlusVoigtian") outFile << " & " <<truthVector.second[i-catStart][j] <<"Lau + Voi ";
			if(truthVector.first[i-catStart][j] == "ExponentialPlusVoigtian") outFile << " & " <<truthVector.second[i-catStart][j] <<"Exp + Voi ";
			if(truthVector.first[i-catStart][j] == "PowerLawPlusVoigtian") outFile << " & " <<truthVector.second[i-catStart][j] <<"Pow + Voi ";
		
			if(truthVector.first[i-catStart][j] == "Bernstein") outFile << " & " <<truthVector.second[i-catStart][j] <<"Ber ";
                        if(truthVector.first[i-catStart][j] == "Chebychev") outFile << " & " <<truthVector.second[i-catStart][j] <<"Che ";
                        if(truthVector.first[i-catStart][j] == "Laurent") outFile << " & " <<truthVector.second[i-catStart][j] <<"Lau ";
                        if(truthVector.first[i-catStart][j] == "Exponential") outFile << " & " <<truthVector.second[i-catStart][j] <<"Exp ";
                        if(truthVector.first[i-catStart][j] == "PowerLaw") outFile << " & " <<truthVector.second[i-catStart][j] <<"Pow ";
		}
		if(truthVector.first[i-catStart].size() == 3) outFile << " & X";
		outFile << " \\\\" << endl;
		outFile << "   \\hline" << endl;	
		fstream biasFile;
                biasFile.open(Form("%s/CombinedPlots/maxBias_cat%d.txt",biasDirectory.c_str(),i),ios::in);
                cout << endl << Form("%s/CombinedPlots/maxBias_cat%d.txt",biasDirectory.c_str(),i);
		for(int j = 0; j < testVector.first[i-catStart].size(); j++)
	        {
			getline(biasFile,tempLine);
			if(testVector.first[i-catStart][j] == "BernsteinPlusVoigtian") outFile << testVector.second[i-catStart][j] <<"Ber + Voi & ";
			if(testVector.first[i-catStart][j] == "ChebychevPlusVoigtian") outFile << testVector.second[i-catStart][j] <<"Che + Voi & ";
			if(testVector.first[i-catStart][j] == "LaurentPlusVoigtian") outFile << testVector.second[i-catStart][j] <<"Lau + Voi & ";
			if(testVector.first[i-catStart][j] == "ExponentialPlusVoigtian") outFile << testVector.second[i-catStart][j] <<"Exp + Voi & ";
			if(testVector.first[i-catStart][j] == "PowerLawPlusVoigtian") outFile << testVector.second[i-catStart][j] <<"Pow + Voi & "; 

			if(testVector.first[i-catStart][j] == "Bernstein") outFile << testVector.second[i-catStart][j] <<"Ber & ";
                        if(testVector.first[i-catStart][j] == "Chebychev") outFile << testVector.second[i-catStart][j] <<"Che & ";
                        if(testVector.first[i-catStart][j] == "Laurent") outFile << testVector.second[i-catStart][j] <<"Lau & ";
                        if(testVector.first[i-catStart][j] == "Exponential") outFile << testVector.second[i-catStart][j] <<"Exp & ";
                        if(testVector.first[i-catStart][j] == "PowerLaw") outFile << testVector.second[i-catStart][j] <<"Pow & ";

			outFile << tempLine << " \\\\" << endl;
			cout << endl << tempLine;	
	
	        }
		biasFile.close();

		outFile << "   \\hline" << endl;
		outFile << "\\end{tabular}" << endl;
		outFile << "\\end{center}" << endl;
		outFile << "\\end{frame}" << endl;
		outFile << "" << endl;

	}
/*	
    4Ber + Voi & -0.12 & \textbf{0.42} & \textbf{-0.18} & X \\ 
    3Che + Voi & 0.13 & 0.07 & 0.08 & X  \\
    2Exp + Voi & \textbf{0.69} & \textbf{0.78} & \textbf{0.42} & X \\
    2Exp + Voi & \textbf{0.69} & \textbf{0.78} & \textbf{0.42} & X \\
    \hline
\end{tabular}
\end{center}
\end{frame}
*/

	outFile.close();

	
	return 0;
}
