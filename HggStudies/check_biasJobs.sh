#! /usr/local/bin/bash -l

cd ${1}

echo "Failed jobs : "
ls *.fail 2>/dev/null 

echo "Jobs in run : "
ls *.run 2>/dev/null

echo "Jobs succedeed : "
ls *.done 2>/dev/null

cd -

