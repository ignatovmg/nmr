#!/bin/bash

# experiment
NOESY=data/good/tvol.csv

# path to complexes
PDB=data/good/md.pdb
PSF=data/good/md.psf
RTF=data/good/prm/top_all36_prot_met.rtf
PRM=data/good/prm/par_all36_prot_met.prm
GRP="sandbox/groups/eq_groups"
EXP="sandbox/matrix/0"

# preprocess data
python src/dataprep.py   ${PDB}
python src/preprocess.py ${NOESY}

./build/main ${PDB} ${PSF} ${PRM} ${RTF} ${GRP} ${EXP} || exit 1

echo
echo "Scores can be found in sandbox"
