#!/bin/bash

RESTRAINTS=data/formisha/dph/restraints.csv
PDB_FILE=sandbox/100k_noemin_k5/some_cmin_noemin.pdb
PDB_FILE=data/formisha/dph/some_cmin.pdb
PDB_FILE=data/runs/100k_noemin_k5/CIL_006151.pdb
OUT_FOLDER=sandbox/100k_noemin_k5

NUM_PAIR=9

python src/distance_filter_verbose.py ${RESTRAINTS} ${PDB_FILE} ${OUT_FOLDER}

PEN_FILE=$(basename ${PDB_FILE})
PEN_FILE=${OUT_FOLDER}/penalties_${PEN_FILE%.pdb}
#PEN_FILE=${RESTRAINTS}

echo ${PEN_FILE}

python src/pymol_draw_dists.py ${PDB_FILE} ${PEN_FILE} ${NUM_PAIR}
