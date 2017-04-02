import __main__
__main__.pymol_argv = [ 'pymol' ]

import sys, time, os
import pymol

pymol.finish_launching()

from pymol import cmd

cmd.load(sys.argv[1])
penalties = sys.argv[2]
max_num = int(sys.argv[3])

#cmd.set("dash_width", 4)
#cmd.set("dash_width", 4)

count = 0
with open(penalties, "r") as f:
	for line in f.readlines()[:max_num]:
		split = line.split()
		
		name = split[0]+"_"+split[1]
		cmd.distance(name, "id "+split[0], "id "+split[1])
		
		if max_num <= 10:
			cmd.color("br"+str(9-count), name)
		
		count += 1
