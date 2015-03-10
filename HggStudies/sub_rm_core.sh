#!/bin/bash 
##set -x

echo "----START----"

cd /afs/cern.ch/work/l/lsgandur/CMSSW_6_2_11_sl6_BiasStudy/src/BiasStudy 
eval `scramv1 runtime -sh`

for i in `seq 1 100`
do
	rm core.*
	sleep 5m
	echo "$i * 5 minutes"
done


echo "----END----"


