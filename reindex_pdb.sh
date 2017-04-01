#!/bin/bash

REFPDB=CilMD.pdb;
TARPDB=CilMD_mod.pdb;

#awk '
#	FNR == NR { 
#		if ($1 == "ATOM" || $1 == "HETATM") {
#			print substr($0, 6, 8);
#			a[substr($0, 22, 26)][substr($0, 17, 20)][substr($0, 12, 16)] = substr($0, 6, 11);
#		}
#		next;
#	}
	
#	($1 == "ATOM") || ($1 == "HETATM") {
#		if ((substr($0, 22, 26), substr($0, 17, 20), substr($0, 12, 16)) in a) {
#			print a[substr($0, 22, 26)][substr($0, 17, 20)][substr($0, 12, 16)];
#		}
#	}' ${REFPDB} ${TARPDB}

#awk '$1=="ATOM" { 
#	print; 
#	}' ${REFPDB} ${TARPDB}
	

	
#if ((substr($0, 23, 26), substr($0, 18, 20), substr($0, 13, 16)) in a)	

while IFS='' read -r line || [[ -n "$line" ]]; do
	if [[ ${line:0:6} == "ATOM  " || ${line:0:6} == "HETATM" ]]; then
		echo ${line:17:3}${line:12:4}
		#arr["${line:17:3}""${line:12:4}"]="${line:6:5}"
	fi
done < ${REFPDB}

echo "arr[*]"
	
