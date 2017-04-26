#!/bin/bash

DIR=data/runs/100k
OUT=data/runs/100k_cor
mkdir ${OUT}

CUR=0
for SUFFIX in {0..9}; do
	for PDB in ${DIR}/*${SUFFIX}.pdb; do
		cat ${PDB} | sed -e "64,82 s/MVA/VAL/g" -e "44,63 s/DPN/DPH/g" -re "s/([1-3])(H[HG][1-2])/\2\1/g" > ${OUT}/$(basename ${PDB})

		CUR=$((CUR+1))
		printf "\t%.0s ${SUFFIX}: $CUR\r" {0..${SUFFIX}}
	done &
done;
