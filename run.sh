#!/bin/bash

FILE_LIST=`find ./sandbox/refined/*`

# clean file
> ./sandbox/chi_scores

# calc chi_scores
for FILE in ${FILE_LIST}
do
./build/main "./sandbox/groups/eq_groups" ${FILE} "./sandbox/chi_scores" `find ./sandbox/matrix/*`
done


