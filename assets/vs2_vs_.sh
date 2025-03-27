#!/usr/bin/env bash


mkdir others
mv  *.txt others
cd others
mv *complete* ..
cd ..

for i in *complete*; do mv $i `echo $i | sed 's/vs/_vs_/g'` ; done
