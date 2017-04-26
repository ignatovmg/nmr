#!/bin/bash

function progress_bar {
	CUR=$1; LEN=50; MAX=100;

	PROG=$(( ${CUR}*${LEN}/${MAX} ))
	FILLED_PART=$( printf '⬛%.0s' $(seq 0 ${PROG}) )
	EMPTY__PART=$( printf ' %.0s' $(seq 0 $(( ${LEN}-${PROG} ))) )
	 PERCENTAGE=$( printf '%3i%%' ${CUR} )

	echo -ne "\r\033[1;33m[${FILLED_PART}${EMPTY__PART}${PERCENTAGE}]\033[0m"
}

# experiment
NOESY=data/cilengitit_rosi/rosi_vol.csv
NOESY=data/noe_testing/some_cmin_xplor_volume

# path to complexes
PDBS=data/cilengitit_rosi/CilMD.pdb
PDBS=data/noe_testing/some_cmin_xplor.pdb

# preprocess data
python src/dataprep.py ${PDBS}
python src/preprocess.py ${NOESY}

FILE_LIST=`find ./sandbox/refined/*`

# clean file
> ./sandbox/chi_scores_md

# calc scores
ARR=(${FILE_LIST})
MAX_NUM=${#ARR[@]}
COUNTER=0
echo

# main cycle
time for FILE in ${FILE_LIST}
do

SPECTRUM_FILES=$(find ./sandbox/matrix/*) || exit 1

./build/main "./sandbox/groups/eq_groups" ${FILE} "./sandbox/chi_scores_md" ${SPECTRUM_FILES} || exit 1

# for fun
COUNTER=$((COUNTER+100))
progress_bar $((COUNTER/MAX_NUM))

done


# sort results
#python src/sortres.py

echo
echo "Scores can be found in sandbox"
