# run while in nmr directory

import pandas as pd
import numpy as np
import os
import glob

if not os.path.isdir('./sandbox'): os.mkdir('sandbox')
if not os.path.isdir('./sandbox/refined'): os.mkdir('./sandbox/refined')
if not os.path.isdir('./data'): exit('There is no data dir in current dir')
os.chdir('data')

# find all files with .pdb extension, get H from them and put into another file
for fname in glob.glob('*.pdb'):
	df = pd.read_csv(fname, sep='\s+', skiprows=1, header=None)
	df = df.ix[(df[0] == 'ATOM') & df[2].str.startswith('H'), [1,2,6,7,8]]
	df.columns = np.arange(5)
	fun = lambda col: df[col].map(lambda x: '%.3f' % x)
	df[2] = fun(2); df[3] = fun(3); df[4] = fun(4)
	df[0] = df[0].map(int)
	df.to_csv('../sandbox/refined/'+fname[:-3]+'dat', header=False, sep='\t', index=False)
	
	
	






