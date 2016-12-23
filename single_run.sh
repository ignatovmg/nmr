#!/bin/bash

FILE_NAME=./sandbox/refined/CIL_00001.dat

./build/main ./sandbox/groups/eq_groups ${FILE_NAME} ./sandbox/chi_scores `find ./sandbox/matrix/*`

python src/assignfreqs.py

#TEXT=`head -n-1 ./sandbox/chi_scores`
#echo $TEXT > ./sandbox/chi_scores
