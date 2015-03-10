#!/bin/bash 

inputDirectory=${1}
datafile=${2}
massLow=${3}
massHigh=${4}
ncats=${5}

##for cat in `seq 5 ${5}`
for cat in '11'
do
	./BiasPlotter_2.exe ${inputDirectory} ${datafile} ${cat} ${massLow} ${massHigh}
done



