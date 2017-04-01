#!/bin/bash

RESTRAINTS=data/cilengitit_rosi/correct_distance_restraints_reind.csv #data/cilengitit_rosi/correct_rest_schdrop_reind.csv #data/cilengitit_rosi/correct_distance_restraints_reind.csv
PDB_FOLDER=data/bst #data/cilengitit_rosi #data/RUN2_aligned
OUT_FILE=sandbox/719715rotamers10k/distance_filtered #/md_dist #/distance_filtered_60k
THREADS=6

mkdir sandbox/719715rotamers10k

for ID in $(seq 0 $((THREADS-1))) 
do
	> ${OUT_FILE}_${ID}
done

python src/distance_filter.py ${RESTRAINTS} ${PDB_FOLDER} ${OUT_FILE} ${THREADS}

> ${OUT_FILE}

for ID in $(seq 0 $((THREADS-1)))
do
	cat ${OUT_FILE}_${ID} >> ${OUT_FILE}
	rm -f ${OUT_FILE}_${ID} 
done

LC_ALL=C sort -k2 -n ${OUT_FILE} > bufer
cat bufer > ${OUT_FILE} 
rm -f bufer

