#!/bin/bash

# experiment
NOESY=data/noe_testing/some_cmin_xplor_volume

# path to complexes
PDB=data/formisha/dph/some_cmin.pdb
PSF=data/formisha/dph/some_cmin.psf
PRM=data/formisha/dph/par_all36_prot_met.prm
RTF=data/formisha/dph/top_all36_prot_met.rtf
GRP="sandbox/groups/eq_groups"
EXP="sandbox/matrix/0"

# preprocess data
python src/dataprep.py   ${PDB}
python src/preprocess.py ${NOESY}

./build/main ${PDB} ${PSF} ${PRM} ${RTF} ${GRP} ${EXP} || exit 1

echo
echo "Scores can be found in sandbox"
