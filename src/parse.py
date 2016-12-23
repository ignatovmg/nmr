# designed specifically for parsing the output of the Mark program

# coding: utf-8

# In[86]:

import re


# In[87]:

pdbfile = '/home/mikhail/Projects/nmr/data/Cyclo-dA/useful/c_dA_final_dnm.pdb'
subject = '/home/mikhail/Projects/nmr/data/Cyclo-dA/useful/C_dA_300_pks.txt'
outpath = '/home/mikhail/Projects/nmr/data/Cyclo-dA/useful/out'


# In[95]:

atoms = dict()
with open(pdbfile, 'r') as f:
    for line in f.readlines():
        if line.startswith('ATOM'):
            segment = line[66:76].strip()
            residue = line[22:26].strip()
            atomidx = line[12:16].strip()
            atomnum = line[6 :11].strip()
            if 'H' in atomidx:
                atoms.update({(segment, residue, atomidx):int(atomnum)})


# In[96]:

out = open(outpath, 'w')

pattern = re.compile('.+resid ([0-9]+) .+name (.+) .+ segid (.+)\)')

with open(subject, 'r') as f:
    a = b = c = -1
    for line in f.readlines():
        match = re.match(pattern, line)
        if match:
            if line.startswith('assign'):
                if a != -1:
                    out.write('\t%s\t%s\t%s\n' % (c, a, b))
                a = atoms[(match.group(3), match.group(1), match.group(2))]
            else:
                b = atoms[(match.group(3), match.group(1), match.group(2))]
        else:
            volume = re.match('\s+(\S+)\s+(\S+)', line)
            if volume:
                c = volume.group(1)
        
out.close()

