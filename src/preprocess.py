
# coding: utf-8

# In[6]:

import pandas as pd
import numpy as np
import seaborn as sns
import os
import glob
import sys


# In[7]:

perm_path = "./sandbox/matrix/"
summary = "./sandbox/groups/"
protons = "./sandbox/protons"

if not os.path.isdir('./sandbox/matrix/'):
	os.mkdir('./sandbox/matrix/')
else:
	filelist = glob.glob('./sandbox/matrix/*')
	for f in filelist:
		os.remove(f)
   
if not os.path.isdir('./sandbox/groups/'): 
	os.mkdir('./sandbox/groups/')
else:
	filelist = glob.glob('./sandbox/groups/*')
	for f in filelist:
		os.remove(f)

# In[8]:

noesy = pd.read_csv(sys.argv[1], sep='\t', header=None)
noesy = noesy.iloc[:, 1:4]
noesy.columns = list(range(noesy.shape[1]))
noesy.iloc[:, 1:3] = noesy.iloc[:, 1:3].applymap(str)


# In[9]:

# find groups of equivalent protons

groups = set()
for item in list(noesy[1])+list(noesy[2]):
    for chunk in item.split('|'):
        groups.add(chunk)
        
f = open(summary+'eq_groups', 'w')
counter = 0
groups_dict = {}
for item in sorted(groups):
    if ';' in item:
        f.write("%i\t%s\n" % (counter, '\t'.join(item.split(';'))))
        #groups_dict.update({counter: tuple(map(int, item.split(';')))})
    else:
        f.write("%i\t%s\n" % (counter, item))
        #groups_dict.update({counter: tuple([int(item)])})
    groups_dict.update({item: counter})
    counter+= 1

f2 = open(protons, "r")

current_proton_list = ";".join(groups_dict.keys()).split(";")
for proton in f2.readlines():
	proton = proton[:-1]
	if proton not in current_proton_list:
		groups_dict.update({proton: counter})
		f.write("%i\t%s\n" % (counter, proton))
		counter += 1
		
f2.close()



# In[10]:

noesy.iloc[:, 1:] = noesy.iloc[:, 1:].applymap(lambda x: '|'.join([str(groups_dict[i]) for i in x.split('|')]))


# In[11]:

swap_groups = {}

for i in range(noesy.shape[0]):
    for j in range(1, 3):
        if '|' in noesy.iloc[i, j]:
            key = tuple(noesy.iloc[i, j].split('|'))
            if (key[1], key[0]) in swap_groups.keys():
                key = (key[1], key[0])
            value = [(i, j)]
            if key in swap_groups.keys():  
                value = swap_groups[key]+[(i, j)]
            swap_groups.update({key: value})


# In[12]:

def get_perms(n):
    if n == 0:
        return [[0]]
    if n == 1:
        return [[0], [1]]
    a = get_perms(int(n/2.))
    b = get_perms(n - int(n/2.))
    c = []
    for i in a:
        for j in b:
            c.append(i+j)
    return c

def do_swap(df, cells):
    for i, j in cells:
        val = df.iloc[i, j].split('|')
        df.iloc[i, j] = val[1]+'|'+val[0]
    return df
        
def print_matrix(df):
    return pd.concat([df.iloc[:, 1:].applymap(lambda x: x.split('|')[0]), df.iloc[:, 0]], axis=1)

perms = get_perms(len(swap_groups))
pairs = swap_groups.keys()
for perm in perms:
    noesy_tmp = noesy.copy()
    for flag, key in zip(perm, pairs):
        if flag == 1:
            noesy_tmp = do_swap(noesy_tmp, swap_groups[key])
    print_matrix(noesy_tmp).to_csv(perm_path+''.join([str(x) for x in perm]), sep='\t', header=None, index=False)
            


# In[13]:

f = open(summary+'swap_grps', 'w')
for i, j in pairs:
    f.write('%i\t%i\n' % (int(i), int(j)))
f.close()

