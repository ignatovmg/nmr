#!/bin/bash

PDB=data/formisha/dph/some_cmin.pdb #data/cilengitit_rosi/mutated3_mod.pdb #data/cilengitit_rosi/CilMD_mod.pdb
PSF=data/formisha/dph/some_cmin.psf #data/cilengitit_rosi/mutated3_mod.psf #data/cilengitit_rosi/CilMD_mod.psf
RTF=data/formisha/dph/top_all36_prot_met.rtf #data/toppar/top_all22_prot_num.rtf #data/toppar/top_all22_prot_num.rtf #data/toppar/top_all22_prot_num.rtf #data/toppar/top_all22_prot.rtf #data/toppar/top_all36_prot.rtf #data/cilengitit_rosi/charmm_param.rtf #data/top_all22_prot_modified.rtf #data/cilengitit_rosi/charmm_param.rtf
PRM=data/formisha/dph/par_all36_prot_met.prm #data/toppar/par_all22_prot_mod.prm #data/cilengitit_rosi/charmm_param.prm #data/toppar/top_all22_prot.prm #data/toppar/top_all22_prot.prm #data/cilengitit_rosi/charmm_param.prm #data/toppar/par_all36_prot.prm   #data/cilengitit_rosi/charmm_param.prm #data/toppar/par_all22_prot.prm   #data/cilengitit_rosi/charmm_param.prm
FIX=dummy
SPR=data/formisha/dph/springs_input
OUT=data/formisha/dph/some_cmin_minimized.pdb

#./build/spring_pair_min ${PDB} ${PSF} ${PRM} ${RTF} ${FIX} ${SPR} ${OUT}
#exit

PSF=data/formisha/dph/some_cmin.psf #data/cilengitit_rosi/mutated3_mod.psf #data/cilengitit_rosi/CilMD_mod.psf
RTF=data/formisha/dph/top_all36_prot_met.rtf #data/toppar/top_all22_prot_num.rtf #data/toppar/top_all22_prot_num.rtf #data/toppar/top_all22_prot_num.rtf #data/toppar/top_all22_prot.rtf #data/toppar/top_all36_prot.rtf #data/cilengitit_rosi/charmm_param.rtf #data/top_all22_prot_modified.rtf #data/cilengitit_rosi/charmm_param.rtf
PRM=data/formisha/dph/par_all36_prot_met.prm #data/toppar/par_all22_prot_mod.prm #data/cilengitit_rosi/charmm_param.prm #data/toppar/top_all22_prot.prm #data/toppar/top_all22_prot.prm #data/cilengitit_rosi/charmm_param.prm #data/toppar/par_all36_prot.prm   #data/cilengitit_rosi/charmm_param.prm #data/toppar/par_all22_prot.prm   #data/cilengitit_rosi/charmm_param.prm
FIX=dummy
SPR=data/formisha/dph/springs_input

INDIR=data/runs/100D
OUTDIR=data/runs/100D_min_k14
OUTPUT_PREFIX=sandbox/100D_min_k14/spring_min
mkdir ${OUTDIR}

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

for FILE in ${OUTPUT_PREFIX}*.log; do cat ${FILE} >> ${OUTPUT_PREFIX}.log && rm -f ${FILE}; done;
