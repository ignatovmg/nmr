
# coding: utf-8

# In[38]:

import Bio.PDB
import sys
import glob
import os


# In[37]:

ref_filename = sys.argv[1]
alt_list = glob.glob(sys.argv[2]+'/*.pdb')
out_path = sys.argv[3]
rmsd_path = sys.argv[4]

rmsd_file = open(rmsd_path, 'w')

io=Bio.PDB.PDBIO()

print("Loading PDB file %s" % ref_filename)
structure = Bio.PDB.PDBParser(QUIET=True).get_structure('ref', ref_filename)[0]

print("Everything aligned to first model...")
ref_model = list(structure.get_chains())[0]
ref_atoms = []
for res in structure.get_residues():
    ref_atoms += [res['N'], res['CA'], res['C']]

for alt_filename in alt_list :
    #Build paired lists of c-alpha atoms
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure('alt', alt_filename)[0]
    alt_model = list(structure.get_chains())[0]
    
    alt_atoms = []
    for i in [3, 4, 5, 1, 2]:
        alt_atoms += [alt_model[i]['N'], alt_model[i]['CA'], alt_model[i]['C']]

    #Align these paired atom lists:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, alt_atoms)
        
    super_imposer.apply(structure.get_atoms())
    io.set_structure(structure)
    io.save(out_path+"/"+os.path.basename(alt_filename)[:-4]+'_aln.pdb')

    rmsd_file.write("%s %.2f\n" % (os.path.basename(alt_filename), super_imposer.rms))
    
rmsd_file.close()


# In[ ]:



