#!/bin/bash 

queue=1nd
outputDirectory=BiasJobs_massfacmva_legacy_notHasPixelSeed_80_120_lastStaticCutsP
bkdWorspace=CMS-HGG_massfacmva_legacy_notHasPixelSeed_80_120_lastStaticCuts.root
sigWorspace=CMS-HGG_massfacmva_legacy_notHasPixelSeed_80_120_lastStaticCuts.root
datafile=config_massfacmva_legacy_notHasPixelSeed_80_120_lastStaticCutsP.dat
expectedSignal=0.0
muLow=-1000.0
muHigh=1000.0
muStep=0.01
blind=false
toysPerJob=2000
binning=40
eVeto=notHasPixelSeed

for cat in `seq 0 13`
do

	for mass in '95' '100' '105' '110' 
	do
		for job in '0'
		do
	
			bsub -q ${queue} -o ${outputDirectory}/cat_${cat}_expectSignalMass_${mass}_job_${job}.log sub_biasJobs.sh ${outputDirectory} ${bkdWorspace} ${sigWorspace} ${datafile} ${cat} ${expectedSignal} ${mass} ${muLow} ${muHigh} ${muStep} ${blind} ${toysPerJob} ${job} ${binning} ${eVeto}	
			##bsub -q ${queue} -o ${outputDirectory}/cat_${cat}_expectSignalMass_${mass}_job_${job}.log sub_biasJobs2.sh ${outputDirectory} ${bkdWorspace} ${sigWorspace} ${datafile} ${cat} ${expectedSignal} ${mass} ${muLow} ${muHigh} ${muStep} ${blind} ${toysPerJob} ${job} ${binning}	

		done
	done
done



##bsub -q 8nh -o sub_biasJobs.sh.log sub_biasJobs.sh test_subBias2 CMS-HGG_refit_BernPlusVoi_NewAddNoExt.root CMS-HGG_refit_BernPlusVoi_NewAddNoExt.root config.dat 0 0.0 125 -3.0 3.0 0.01 false 10 0
