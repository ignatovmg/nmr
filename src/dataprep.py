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
    
if not os.path.isdir(pdbs): exit('There is no %s dir in current dir' % pdbs)
#os.chdir(pdbs)

# find all files with .pdb extension, get H from them and put into another file
for fname in glob.glob('%s/*.pdb' % pdbs):
	'''df = pd.read_csv(fname, sep='\s+', skiprows=1, header=None)
	df = df.ix[(df[0] == 'ATOM') & (df[2].str.contains('^[0-9]?H', regex=True)), [1,2,6,7,8]]
	df.columns = np.arange(5)
	fun = lambda col: df[col].map(lambda x: '%.3f' % x)
	df[2] = fun(2); df[3] = fun(3); df[4] = fun(4)
	df[0] = df[0].map(int)
	df.to_csv('./sandbox/refined/'+fname.split('/')[-1][:-3]+'dat', header=False, sep='\t', index=False)'''
	
	out = open('./sandbox/refined/'+fname.split('/')[-1][:-3]+'dat', 'w')
	with open(fname, 'r') as f:
		for line in f.readlines():
			if line.startswith('ATOM'):
				if 'H' in line[12:16]:
					idx = line[6:11].strip()
					name = line[12:16].strip()
					x = float(line[30:38].strip())
					y = float(line[38:46].strip())
					z = float(line[46:54].strip())
					out.write('%s\t%s\t%.3f\t%.3f\t%.3f\n' % (idx, name, x, y, z))
	out.close()
	
	






