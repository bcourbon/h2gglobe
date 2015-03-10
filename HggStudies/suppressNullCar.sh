#! /usr/local/bin/bash -l

tr -d '\000' < ${1} > myfile.out
cp myfile.out ${1}
