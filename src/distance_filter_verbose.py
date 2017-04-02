import pandas as pd
import numpy as np
import glob
import sys
import traceback
import os

k = 1.0

restraints_file = sys.argv[1]
pdb_file = sys.argv[2]
out_dir = sys.argv[3]

restraints = pd.read_csv(restraints_file, sep='\t')
restraints = restraints.dropna()

#print("%7i: %s" % (counter, pdb_file))
atoms = {}

with open(pdb_file, 'r') as f:
	for line in f.readlines():
		if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
			atoms.update( { int(line[6:11]): np.array([ float(line[30:38]), float(line[38:46]), float(line[46:54]) ]) } )

penalty = 0.0
penalties = []

for row in range(restraints.shape[0]):
	i = restraints.iloc[row, 0]
	j = restraints.iloc[row, 1]
	min_d = restraints.iloc[row, 2]
	max_d = restraints.iloc[row, 3]

	distance = np.sqrt(((atoms[i] - atoms[j]) * (atoms[i] - atoms[j])).sum())
	if min_d > distance:
		penalty += k*(min_d - distance)*(min_d - distance);
		penalties.append([abs(min_d - distance), i, j])
		
	if max_d < distance:
		penalty += k*(max_d - distance)*(max_d - distance);
		penalties.append([abs(max_d - distance), i, j])

name = out_dir+"/penalties_"+os.path.basename(pdb_file)[:-4]
penalties_file = open(name, "w") 
penalties = sorted(penalties, reverse=True)

for item in penalties:
	penalties_file.write("%i\t%i\t%.3f\n" % (item[1], item[2], item[0]))

penalties_file.close()

print("Penalty: %f\n" % penalty)

		
