# run while in nmr directory

import pandas as pd
import numpy as np
import os
import glob
import sys

pdbs = sys.argv[1]

if not os.path.isdir('./sandbox'): os.mkdir('./sandbox')
if not os.path.isdir('./sandbox/refined'): os.mkdir('./sandbox/refined')
else:
	filelist = glob.glob('./sandbox/refined/*')
	for f in filelist:
		os.remove(f)
    
files = []
if pdbs[-3:] == 'pdb':
	files = [pdbs]
else:
	if not os.path.isdir(pdbs): exit('There is no %s dir in current dir' % pdbs)
	files = glob.glob('%s/CIL_??????.pdb' % pdbs)

#os.chdir(pdbs)

# find all files with .pdb extension, get H from them and put into another file
for fname in files:
	out = open('./sandbox/refined/'+fname.split('/')[-1][:-4], 'w')
	with open(fname, 'r') as f:
		for line in f.readlines():
			if line.startswith('ATOM') or line.startswith('HETATM'):
				if 'H' == line[13] or 'H' == line[12]:
					idx = line[6:11].strip()
					name = line[12:16].strip()
					x = float(line[30:38].strip())
					y = float(line[38:46].strip())
					z = float(line[46:54].strip())
					out.write('%s\t%s\t%.3f\t%.3f\t%.3f\n' % (idx, name, x, y, z))
	out.close()
	
f = open(files[0], 'r')
proton_list = open('./sandbox/protons', 'w')
for line in f.readlines():
	if line.startswith('ATOM') or line.startswith('HETATM'):
		if 'H' == line[13] or 'H' == line[12]:
			idx = line[6:11].strip()
			proton_list.write(idx+"\n")
			
proton_list.close()
f.close()
			

	
	






