#!/bin/bash

SPRK=10
PDB_DIR=100k
CONF_DIR=data/runs/raw/${PDB_DIR}
RAW_DIR=${CONF_DIR}/raw
COR_DIR=${CONF_DIR}/corrected

echo "convert pdbs?"
read ANS

if [[ ${ANS} == "y" ]]; then

	if [ -d ${COR_DIR} ]; then
		echo "${COR_DIR} exists, are you sure?"
		read ANS
	fi

	mkdir ${COR_DIR}

	CUR=0
	for SUFFIX in {0..9}; do
		for PDB in ${RAW_DIR}/*${SUFFIX}.pdb; do
			cat ${PDB} | sed -e "64,82 s/MVA/VAL/g" -e "44,63 s/DPN/DPH/g" -re "s/([1-3])(H[HG][1-2])/\2\1/g" > ${COR_DIR}/$(basename ${PDB})

			CUR=$((CUR+1))
			printf "\t%.0s ${SUFFIX}: ${CUR}\r" {0..${SUFFIX}}
		done &
	done

	wait
fi

MAIN_DIR=data/good

echo "proceed to min? SPRK=${SPRK}"
read

RESTRAINTS=${MAIN_DIR}/restraints.csv
SPRING_DIR=${MAIN_DIR}

SPR_COORD=${SPRING_DIR}/springs.pdb
MIN_INPUT=${SPRING_DIR}/springs_input


> ${SPR_COORD}
> ${MIN_INPUT}

python -c '

import sys

out1 = open(sys.argv[2], "w")
out2 = open(sys.argv[3], "w")

fk = float(int(sys.argv[4]));

with open(sys.argv[1], "r") as f:
	lines = f.readlines()
	out1.write("%i\n%s\n" % (len(lines) - 1, sys.argv[3]))
	
	for line in lines[1:]:
		spt = line.split()
		atom1 = int(spt[0])
		atom2 = int(spt[1])
		dst12 = (float(spt[2]) + float(spt[3]))/2.0
		
		out1.write("%.3f\n%.3f\n" % (dst12, fk))
		out2.write("ATOM%7i\n" % atom1)
		out2.write("ATOM%7i\n" % atom2)

out2.write("END\n")

out1.close()
out2.close()	
			
' ${RESTRAINTS} ${MIN_INPUT} ${SPR_COORD} ${SPRK} || exit






echo "springs generated, start min?"
read

PSF=${MAIN_DIR}/md.psf #data/cilengitit_rosi/mutated3_mod.psf #data/cilengitit_rosi/CilMD_mod.psf
RTF=${MAIN_DIR}/prm/top_all36_prot_met.rtf
PRM=${MAIN_DIR}/prm/par_all36_prot_met.prm
FIX=dummy
SPR=${MAIN_DIR}/springs_input

RAW_PDB=${COR_DIR}
RES_DIR=data/runs/springs/${PDB_DIR}/fk${SPRK}
MIN_PDB=${RES_DIR}/pdb

if [ -d ${RES_DIR} ]; then
		echo "${RES_DIR} exists, are you sure?"
		read
fi

mkdir -p ${MIN_PDB}

OUTPUT_PREFIX=${RES_DIR}/min


CUR=0
for SUFFIX in {0..9}; do
	OUTNAME=${OUTPUT_PREFIX}${SUFFIX}.log
	> ${OUTNAME}
	for PDB in ${RAW_PDB}/*${SUFFIX}.pdb
	do
		OUT=${MIN_PDB}/$(basename ${PDB})
		CUR=$((CUR+1))
		echo "${PDB}: $CUR"
	
		ENERGY=$(./build/spring_pair_min ${PDB} ${PSF} ${PRM} ${RTF} ${FIX} ${SPR} ${OUT} | tail -n1) || ! echo "Error" || exit 1
		echo "${PDB} ${ENERGY}" >> ${OUTNAME}
	done &
done

wait

> ${OUTPUT_PREFIX}.log
for FILE in ${OUTPUT_PREFIX}?.log; do cat ${FILE} >> ${OUTPUT_PREFIX}.log && rm -f ${FILE}; done;

echo "min finished, start dist filter?"
read


OUT_FILE=${RES_DIR}/distance_filtered
THREADS=8

for ID in $(seq 0 $((THREADS-1))) 
do
	> ${OUT_FILE}_${ID}
done

python src/distance_filter.py ${RESTRAINTS} ${MIN_PDB} ${OUT_FILE} ${THREADS}

> ${OUT_FILE}

for ID in $(seq 0 $((THREADS-1)))
do
	cat ${OUT_FILE}_${ID} >> ${OUT_FILE}
	rm -f ${OUT_FILE}_${ID} 
done

LC_ALL=C sort -k2 -n ${OUT_FILE} > bufer
cat bufer > ${OUT_FILE} 
rm -f bufer

echo "alles gut"
cd ${RES_DIR}

paste <(sort -k1 distance_filtered) <(sort -k1 ${OUTPUT_PREFIX}.log | cut -d' ' -f2) | LC_ALL=C sort -k2 -n > combo

cd -





