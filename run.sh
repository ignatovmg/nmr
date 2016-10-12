#!/bin/bash

# experiment
NOESY=data/noesy.csv

# path to complexes
PDBS=data/RUN2

# preprocess data
python src/dataprep.py ${PDBS}
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

echo
echo "Scores can be found in sandbox"
