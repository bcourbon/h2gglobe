#!/bin/bash 

queue=8nh 
outputDirectory=BiasJobs_massFact_09_01_14_v2_defaultEveto
bkdWorspace=CMS-HGG_massfit_mva_80_180.root
sigWorspace=CMS-HGG_massfit_mva_80_180.root
datafile=config_v3_defaultEveto.dat
expectedSignal=0.0
muLow=-3.0
muHigh=3.0
muStep=0.01
blind=false
toysPerJob=1000


##for cat in 0
for cat in `seq 0 8`
do

	for mass in '90' '100' '110' '120' '130' '140' '150'
	do
		##for job in `seq 0 50`
		##for job in `seq 0 4`
		for job in '0'
		do
	
			bsub -q ${queue} -o ${outputDirectory}/cat_${cat}_expectSignalMass_${mass}_job_${job}.log sub_biasJobs.sh ${outputDirectory} ${bkdWorspace} ${sigWorspace} ${datafile} ${cat} ${expectedSignal} ${mass} ${muLow} ${muHigh} ${muStep} ${blind} ${toysPerJob} ${job}	

		done
	done
done



##bsub -q 8nh -o sub_biasJobs.sh.log sub_biasJobs.sh test_subBias2 CMS-HGG_refit_BernPlusVoi_NewAddNoExt.root CMS-HGG_refit_BernPlusVoi_NewAddNoExt.root config.dat 0 0.0 125 -3.0 3.0 0.01 false 10 0
