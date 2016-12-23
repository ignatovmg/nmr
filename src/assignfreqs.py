import pandas as pd
import numpy as np

comp = pd.read_csv("sandbox/computed_matrix", sep = '\t', header=None)
for i in range(comp.shape[0]):
    comp.iloc[i, i] = 0.0
    
noesy=pd.read_csv("data/noesy2.csv", sep = ',')

groups = {}
groups_file = open("sandbox/groups/eq_groups", 'r')
for line in groups_file.readlines():
    nums = line.split()
    groups.update({';'.join(nums[1:]) : int(nums[0])})
groups_file.close()
groups

freqs = list(noesy.iloc[:, 1])+list(noesy.iloc[:, 2])
names = list(noesy.iloc[:, 4])+list(noesy.iloc[:, 5])

group2freq = [-1]*len(groups)

for freq, name in zip(freqs, names):
    group = name.split('|')[0]
    group2freq[groups[group]] = freq
    
def sym(m):
    for i in range(m.shape[0]):
        for j in range(i, m.shape[0]):
            if np.isnan(m.iloc[i, j]):
                m.iloc[i, j] = m.iloc[j, i] 
            if np.isnan(m.iloc[j, i]):
                m.iloc[j, i] = m.iloc[i, j]
    return m

def assign_freqs(matrix, freqs):
	matrix.columns = freqs
	matrix.index = freqs
	matrix = matrix.sort_index(0)
	matrix = matrix.sort_index(1)
	
	matrix.columns = ['{:.3f}'.format(x) for x in matrix.columns]
	matrix.index = ['{:.3f}'.format(x) for x in matrix.index]
	
	for i in range(matrix.shape[0]):
		for j in range(0, i+1):
			matrix.iloc[i, j] = np.NaN
			
	matrix = pd.DataFrame(matrix.stack())
	#matrix = matrix.dropna()
	return matrix
     
comp = assign_freqs(comp, group2freq)

exp = pd.read_csv("sandbox/matrix/00000", sep = '\t', header=None)
exp = exp.pivot(index=0, columns=1, values=2)
empty = pd.DataFrame(0, columns=range(len(groups)), index=range(len(groups)))
exp = empty+exp
sym(exp)
exp = assign_freqs(exp, group2freq)

comp[1] = exp
comp = comp.dropna()
comp = comp.applymap(lambda x: "{:.0f}".format(x))
comp = comp.reset_index()
comp.columns = ['f1', 'f2', 'computed', 'experiment']
comp.to_csv('sandbox/volume2freqs', sep='\t', index=None)

