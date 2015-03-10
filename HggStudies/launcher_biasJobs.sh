#!/bin/bash 

queue=1nd 
outputDirectory=BiasJobs_massFact_80_180_nothaspixseed_5_m05_02_02_14
bkdWorspace=CMS-HGG_massfacmva_legacy_80_180_notHasPixelSeed.root
sigWorspace=CMS-HGG_massfacmva_legacy_80_180_notHasPixelSeed.root
datafile=config_massfacmva_legacy_80_180_notHasPixelSeed_5_m05.dat
expectedSignal=0.0
muLow=-3.0
muHigh=3.0
muStep=0.01
blind=false
toysPerJob=1000
binning=100

##for cat in 0
for cat in `seq 0 13`
do

	for mass in '90' '100' '110' '120' '130' '140' '150'
	##for mass in '90'
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
