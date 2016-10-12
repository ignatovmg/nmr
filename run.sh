#!/bin/bash

# experiment
NOESY=data/noesy.csv

# preprocess data
python src/dataprep.py
python src/preprocess.py ${NOESY}

FILE_LIST=`find ./sandbox/refined/*`

# clean file
> ./sandbox/chi_scores

# calc scores
for FILE in ${FILE_LIST}
do
./build/main "./sandbox/groups/eq_groups" ${FILE} "./sandbox/chi_scores" `find ./sandbox/matrix/*`
done

# sort results
python src/sortres.py


