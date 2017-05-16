#!/bin/bash

PDB=data/good/md.pdb
PSF=data/good/md.psf
RTF=data/good/prm/top_all36_prot_met.rtf
PRM=data/good/prm/par_all36_prot_met.prm
FIX=dummy
SPR=data/good/springs_input
OUT=./min.pdb

# preprocess data
NOESY=data/good/tvol.csv
python src/preprocess.py ${NOESY}

GRP="sandbox/groups/eq_groups"
EXP="sandbox/matrix/0"

./build/noe_min ${PDB} ${PSF} ${PRM} ${RTF} ${FIX} ${SPR} ${GRP} ${EXP} ${OUT}
