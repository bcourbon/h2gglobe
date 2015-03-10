#!/bin/bash 

echo "----START----"

for i in `seq 1 4`
do
	rm core.*
	sleep 5s
	ls
done


echo "----END----"


