#!/bin/bash

RESTRAINTS=data/cilengitit_rosi/correct_distance_restraints.csv  #data/formisha/dph/restraints.csv  #correct_distance_restraints.csv # test_restraints.csv #data/cilengitit_rosi/correct_rest_schdrop_reind.csv
#data/1000k/CIL_466811.pdb #data/4_of_60k/C35_000130#001.pdb #data/4_of_60k/C52_000198#002.pdb #data/1000k/CIL_719715.pdb # 466811 965811 543475
PDB_FILE=data/60k_min/CIL_033989.pdb
PDB_FILE=data/cilengitit_rosi/CilMD.pdb  #sandbox/60k_min_k14/some_cmin_minimized.pdb
OUT_FOLDER=./

NUM_PAIR=9

python src/distance_filter_verbose.py ${RESTRAINTS} ${PDB_FILE} ${OUT_FOLDER}

PEN_FILE=$(basename ${PDB_FILE})
PEN_FILE=${OUT_FOLDER}/penalties_${PEN_FILE%.pdb}
#PEN_FILE=${RESTRAINTS}

echo ${PEN_FILE}

python src/pymol_draw_dists.py ${PDB_FILE} ${PEN_FILE} ${NUM_PAIR}
