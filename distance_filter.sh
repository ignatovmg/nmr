#!/bin/bash

RESTRAINTS=data/formisha/dph/restraints.csv #data/cilengitit_rosi/correct_rest_schdrop_reind.csv #data/cilengitit_rosi/correct_distance_restraints_reind.csv
PDB_FOLDER=data/runs/100D_min_k7 #data/cilengitit_rosi #data/RUN2_aligned
OUT_FILE=sandbox/100D_min_k7/distance_filtered #/md_dist #/distance_filtered_60k
THREADS=8

mkdir sandbox/1M_min_k7

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

