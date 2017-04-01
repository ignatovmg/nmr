#!/bin/bash

REF_PDB=data/cilengitit_rosi/CilMD.pdb
IN_FOLDER=data/RUN2
OUT_FOLDER=data/RUN2_aligned
RMSD_FILE=${OUT_FOLDER}/rmsd_backbone

python src/rmsd.py ${REF_PDB} ${IN_FOLDER} ${OUT_FOLDER} ${RMSD_FILE}
