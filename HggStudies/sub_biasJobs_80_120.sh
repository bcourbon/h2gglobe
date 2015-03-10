#!/bin/bash 
##set -x

echo "----START----"

cd /afs/cern.ch/work/l/lsgandur/CMSSW_6_1_1_legacy/src/BiasStudy_15_01_14/ 
eval `scramv1 runtime -sh`

##source $HOME/loadRoot.sh

ls -ltrh

mkdir -p ${1}

touch ${1}/cat_${5}_expectSignalMass_${7}_job_${13}.run

ls 

if ( ./BiasStudy1.exe ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14}) then
	touch ${1}/cat_${5}_expectSignalMass_${7}_job_${13}.done
else
	touch ${1}/cat_${5}_expectSignalMass_${7}_job_${13}.fail
fi

rm ${1}/cat_${5}_expectSignalMass_${7}_job_${13}.run

echo "----END----"


##bsub -q 8nh -o sub_biasJobs.sh.log sub_biasJobs.sh test_subBias2 CMS-HGG_refit_BernPlusVoi_NewAddNoExt.root CMS-HGG_refit_BernPlusVoi_NewAddNoExt.root config.dat 0 0.0 125 -3.0 3.0 0.01 false 10 0
