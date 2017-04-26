import sys

refpdb = sys.argv[1]
tarpdb = sys.argv[2]
result = sys.argv[3]

arr = {}
count1 = 0
with open(refpdb, 'r') as f:
	for line in f.readlines():
		if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
			key = (line[17:20], line[12:16])
			resi_id = line[22:26]
			atom_id = line[6:11]
			
			if key in arr.keys():
				arr[key].append((resi_id, atom_id))
			else:
				arr.update({key : [(resi_id, atom_id)]})
			
			count1 += 1
		
conv = {}
with open(tarpdb, 'r') as f:		
	for line in f.readlines():
		if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
			init_key = (line[17:20], line[12:16])
			resi_id = line[22:26]
			atom_id = line[6:11]
			
			if (resi_id, atom_id) in conv.keys():
				continue
			
			key = ();
			tmp_key = (init_key[0], init_key[1][3]+init_key[1][1:3]+init_key[1][0])
			if tmp_key in arr.keys():
				key = tmp_key
			else:
				key = init_key

			if key in arr.keys():
				if len(arr[key]) == 1:
					conv.update({(resi_id, atom_id): ''.join([atom_id, key[1], key[0], resi_id, arr[key][0][1], key[1], key[0], arr[key][0][0]])})
				else:
					string = raw_input('Enter line from reference pdb for atom %s %s %s %s in target pdb: ' % (atom_id, key[1], key[0], resi_id))
					
					corr_atom_id = string[6:11]
					corr_atom_nm = string[12:16]
					corr_resi_nm = string[17:20]
					corr_resi_id = string[22:26]
					
					conv.update({(resi_id, atom_id): ([atom_id, key[1], key[0], resi_id], [corr_atom_id,corr_atom_nm, corr_resi_nm, corr_resi_id])})
			else:
					string = raw_input('Enter line from reference pdb for atom %s %s %s %s in target pdb: ' % (atom_id, key[1], key[0], resi_id))
					
					corr_atom_id = string[6:11]
					corr_atom_nm = string[12:16]
					corr_resi_nm = string[17:20]
					corr_resi_id = string[22:26]
					
					conv.update({(resi_id, atom_id): ''.join([atom_id, key[1], key[0], resi_id, corr_atom_id, corr_atom_nm, corr_resi_nm, corr_resi_id])})
			count1 -= 1

if count1 != 0:
	print('Atom number doesnt match! Difference: %i' % count1)
									
with open(result, 'w') as f:
	for line in conv.values():
		f.write(line+'\n')
	f.close()
