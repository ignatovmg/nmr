#!/bin/bash

NOEK=50
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

echo "proceed to min? NOEK=${NOEK}"

PSF=${MAIN_DIR}/md.psf
RTF=${MAIN_DIR}/prm/top_all36_prot_met.rtf
PRM=${MAIN_DIR}/prm/par_all36_prot_met.prm
FIX=dummy

RAW_PDB=${COR_DIR}
RES_DIR=data/runs/noe/${PDB_DIR}/fk${NOEK}
MIN_PDB=${RES_DIR}/pdb

NOESY=data/good/tvol.csv
python src/preprocess.py ${NOESY}

GRP="sandbox/groups/eq_groups"
EXP="sandbox/matrix/0"

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
	
		ENERGY=$(./build/noe_min ${PDB} ${PSF} ${PRM} ${RTF} dummy dummy ${GRP} ${EXP} ${OUT} | tail -n1) || ! echo "Error" || exit 1
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





