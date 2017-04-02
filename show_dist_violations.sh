#!/bin/bash

RESTRAINTS=data/cilengitit_rosi/correct_distance_restraints.csv #data/cilengitit_rosi/correct_rest_schdrop_reind.csv
PDB_FILE=data/1000k/CIL_904308.pdb  #data/1000k/CIL_466811.pdb #data/4_of_60k/C35_000130#001.pdb #data/4_of_60k/C52_000198#002.pdb #data/1000k/CIL_719715.pdb # 466811 965811 543475
PDB_FILE=data/cilengitit_rosi/CilMD.pdb
OUT_FOLDER=sandbox/1000k

NUM_PAIR=100

python src/distance_filter_verbose.py ${RESTRAINTS} ${PDB_FILE} ${OUT_FOLDER}

PEN_FILE=$(basename ${PDB_FILE})
PEN_FILE=${OUT_FOLDER}/penalties_${PEN_FILE%.pdb}
PEN_FILE=${RESTRAINTS}

echo ${PEN_FILE}

python src/pymol_draw_dists.py ${PDB_FILE} ${PEN_FILE} ${NUM_PAIR}
