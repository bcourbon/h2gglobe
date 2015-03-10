#!/bin/bash 

queue=1nd
outputDirectory=BiasJobs_massfacmva_legacy_defaultE_80_120_PTl036_PTs03_DynPtCut45_NoVoi_noChi2GF_thesis_v3_5000_cat11
bkdWorspace=CMS-HGG_massfacmva_legacy_defaultE_80_120_PTl036_PTs03_DynPtCut45.root
sigWorspace=CMS-HGG_massfacmva_legacy_defaultE_80_120_PTl036_PTs03_DynPtCut45.root
datafile=config_massfacmva_legacy_defaultE_80_120_PTl036_PTs03_DynPtCut45_NoVoi_noChi2GF_v3.dat
expectedSignal=0.0
muLow=-1000.0
muHigh=1000.0
muStep=0.01
blind=false
toysPerJob=2000
binning=40

##for cat in 0
##for cat in `seq 5 13`
##for cat in `seq 6 13`
for cat in 11
do

	for mass in '95' '100' '105' '110' 
	do
		##for job in `seq 0 50`
		##for job in `seq 0 4`
		for job in '0'
		do
	
			bsub -q ${queue} -o ${outputDirectory}/cat_${cat}_expectSignalMass_${mass}_job_${job}.log sub_biasJobs.sh ${outputDirectory} ${bkdWorspace} ${sigWorspace} ${datafile} ${cat} ${expectedSignal} ${mass} ${muLow} ${muHigh} ${muStep} ${blind} ${toysPerJob} ${job} ${binning}	

		done
	done
done



##bsub -q 8nh -o sub_biasJobs.sh.log sub_biasJobs.sh test_subBias2 CMS-HGG_refit_BernPlusVoi_NewAddNoExt.root CMS-HGG_refit_BernPlusVoi_NewAddNoExt.root config.dat 0 0.0 125 -3.0 3.0 0.01 false 10 0
