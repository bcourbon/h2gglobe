#!/bin/bash 

outDirectory=${1}
rootFile=${2}
ncats=${3}

for cat in `seq 0 ${3}`
do
	./massPlot.exe ${outDirectory} ${rootFile} ${cat} 

done



