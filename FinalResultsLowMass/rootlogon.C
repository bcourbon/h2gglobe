{
        gSystem->AddIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include");
	gSystem->Load("../../BackgroundProfileFitting/lib/libBackgroundProfileFitting.so");
	gSystem->Load("../../../../lib/slc6_amd64_gcc472/libHiggsAnalysisGBRLikelihood.so");


}
