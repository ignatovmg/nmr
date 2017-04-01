import pandas as pd
import numpy as np
import glob
import sys
import traceback

from multiprocessing import Process
import time

k = 1.0

def main(restraints, filelist, out_file, norm, thread_id):

	print("%s thread started\n" % thread_id)
	out = open(out_file, 'a')
	counter = 0

	for pdb_file in filelist:
		penalty=np.nan
		
		try:
			#print("%7i: %s" % (counter, pdb_file), end='\r')
			counter += 1
			atoms = {}

			with open(pdb_file, 'r') as f:
				for line in f.readlines():
					if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
						atoms.update( { int(line[6:11]): np.array([ float(line[30:38]), float(line[38:46]), float(line[46:54]) ]) } )

			penalty = 0.0
			for row in range(restraints.shape[0]):
				i = restraints[row, 0]
				j = restraints[row, 1]
				min_d = restraints[row, 2]
				max_d = restraints[row, 3]
				
	
				distance = np.sqrt(((atoms[i] - atoms[j]) * (atoms[i] - atoms[j])).sum())
				if min_d > distance:
					penalty += (min_d - distance)*(min_d - distance);
				
				if max_d < distance:
					penalty += (max_d - distance)*(max_d - distance);
			
			#penalty /= norm
			penalty *= k
				

		except Exception as e:
			print("%s: %s" % ( thread_id, time.ctime(time.time()) ) )
			print("Exception with %s:\n" % pdb_file)
			print("Traceback:\n")
			print(traceback.format_exc())
			print(e)
			continue
		
		else:
			out.write("%s\t" % pdb_file)
			out.write("%.3f\n" % penalty)

	out.close()


restraints_file = sys.argv[1]
pdb_folder = sys.argv[2]
out_file   = sys.argv[3] # appends penalty to the end of it
thread_num = int(sys.argv[4]);

restraints = pd.read_csv(restraints_file, sep='\t')
restraints = np.array(restraints.dropna())
norm = (restraints[:, 2] + restraints[:, 3]) / 2.0
norm = sum(norm * norm)

filelist = glob.glob(pdb_folder+'/*.pdb')
assert(len(filelist) >= thread_num)

if __name__ == '__main__':	
	for thread_id in range(thread_num):	
		try:
			args = (restraints.copy(), filelist[ thread_id : : thread_num ], out_file+"_"+str(thread_id), norm, thread_id, )
			p = Process(target=main, args=args)
			p.start()
			
		except:
	   		print("Error: unable to start thread %i" % thread_id)
	   
	p.join()
		
