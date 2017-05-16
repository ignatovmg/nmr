#!/bin/bash

DIST_FILE=data/good/restraints.csv
DIR=data/good/

SPR_COORD=${DIR}/springs.pdb
MIN_INPUT=${DIR}/springs_input


> ${SPR_COORD}
> ${MIN_INPUT}

python -c '

import sys

out1 = open(sys.argv[2], "w")
out2 = open(sys.argv[3], "w")

fk = 14.0;

with open(sys.argv[1], "r") as f:
	lines = f.readlines()
	out1.write("%i\n%s\n" % (len(lines) - 1, sys.argv[3]))
	
	for line in lines[1:]:
		spt = line.split()
		atom1 = int(spt[0])
		atom2 = int(spt[1])
		dst12 = (float(spt[2]) + float(spt[3]))/2.0
		
		out1.write("%.3f\n%.3f\n" % (dst12, fk))
		out2.write("ATOM%7i\n" % atom1)
		out2.write("ATOM%7i\n" % atom2)

out2.write("END\n")

out1.close()
out2.close()	
			
' ${DIST_FILE} ${MIN_INPUT} ${SPR_COORD} ${DIR}
	
	

