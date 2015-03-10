#include "TFile.h"

int rootLsFile()
{

	TFile *fh1 = TFile::Open("histo_massfacmva_legacy_defaultE_60_180_Pt036.root");
	fh1->ls();	



	return 0;
}
