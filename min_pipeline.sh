#!/bin/bash

SPRK=14
DIR=data/runs/100k
OUT=data/runs/100k_cor
mkdir ${OUT}

echo "convert pdbs?"
read ANS

if [[ ${ANS} == "y" ]]; then
	CUR=0
	for SUFFIX in {0..9}; do
		for PDB in ${DIR}/*${SUFFIX}.pdb; do
			cat ${PDB} | sed -e "64,82 s/MVA/VAL/g" -e "44,63 s/DPN/DPH/g" -re "s/([1-3])(H[HG][1-2])/\2\1/g" > ${OUT}/$(basename ${PDB})

			CUR=$((CUR+1))
			printf "\t%.0s ${SUFFIX}: $CUR\r" {0..${SUFFIX}}
		done &
	done

	mv ${DIR} ${DIR}_raw || exit
	mv ${OUT} ${DIR} || exit

	wait
fi

echo "proceed to min? SPRK=${SPRK}"
read

RESTRAINTS=data/formisha/dph/restraints.csv
DIRRR=data/formisha/dph

SPR_COORD=${DIRRR}/springs.pdb
MIN_INPUT=${DIRRR}/springs_input


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

PSF=data/formisha/dph/some_cmin.psf #data/cilengitit_rosi/mutated3_mod.psf #data/cilengitit_rosi/CilMD_mod.psf
RTF=data/formisha/dph/top_all36_prot_met.rtf #data/toppar/top_all22_prot_num.rtf #data/toppar/top_all22_prot_num.rtf #data/toppar/top_all22_prot_num.rtf #data/toppar/top_all22_prot.rtf #data/toppar/top_all36_prot.rtf #data/cilengitit_rosi/charmm_param.rtf #data/top_all22_prot_modified.rtf #data/cilengitit_rosi/charmm_param.rtf
PRM=data/formisha/dph/par_all36_prot_met.prm #data/toppar/par_all22_prot_mod.prm #data/cilengitit_rosi/charmm_param.prm #data/toppar/top_all22_prot.prm #data/toppar/top_all22_prot.prm #data/cilengitit_rosi/charmm_param.prm #data/toppar/par_all36_prot.prm   #data/cilengitit_rosi/charmm_param.prm #data/toppar/par_all22_prot.prm   #data/cilengitit_rosi/charmm_param.prm
FIX=dummy
SPR=data/formisha/dph/springs_input

INDIR=${DIR}
OUTDIR=${INDIR}_min_k${SPRK}
OUTPUT_PREFIX=sandbox/$(basename ${OUTDIR})/spring_min
FIN_DIR=$(dirname ${OUTPUT_PREFIX})
mkdir ${OUTDIR}
mkdir ${FIN_DIR}

CUR=0
for SUFFIX in {0..9}; do
	OUTNAME=${OUTPUT_PREFIX}${SUFFIX}.log
	> ${OUTNAME}
	for PDB in ${INDIR}/*${SUFFIX}.pdb
	do
		OUT=${OUTDIR}/$(basename ${PDB})
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


#data/cilengitit_rosi/correct_rest_schdrop_reind.csv #data/cilengitit_rosi/correct_distance_restraints_reind.csv
PDB_FOLDER=${OUTDIR} #data/cilengitit_rosi #data/RUN2_aligned
OUT_FILE=${FIN_DIR}/distance_filtered #/md_dist #/distance_filtered_60k
THREADS=8

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

echo "alles gut"
cd ${FIN_DIR}

paste <(sort -k1 distance_filtered) <(sort -k1 spring_min.log | cut -d' ' -f2) | LC_ALL=C sort -k2 -n > combo

cd -





