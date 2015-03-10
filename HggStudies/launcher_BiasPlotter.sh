#!/bin/bash 

inputDirectory=${1}
datafile=${2}
massLow=${3}
massHigh=${4}
ncats=${5}

for cat in `seq 0 ${5}`
##for cat in '4' '8'
do
	./BiasPlotter.exe ${inputDirectory} ${datafile} ${cat} ${massLow} ${massHigh}
	##./BiasPlotter_2.exe ${inputDirectory} ${datafile} ${cat} ${massLow} ${massHigh}
done



