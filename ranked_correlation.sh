#!/bin/bash

LIST1=${1}
LIST2=${2}

time python -c '

import pandas as pd
import numpy as np
import sys
import os

basename = lambda x: os.path.basename(os.path.splitext(x)[0])[:9]

list1 = list(map(basename, pd.read_csv(sys.argv[1], sep="\t", header=None)[0]))
list2 = list(map(basename, pd.read_csv(sys.argv[2], sep="\t", header=None)[0]))

list2 = dict(zip(list2, range(len(list2))))

# https://ragrawal.wordpress.com/2013/01/18/comparing-ranked-list/

measure = 0
pairs = 0

for i in range(len(list1)-1):
	print(i, end="\r")
	for j in range(i+1, len(list1)):
		#pairs += 1
		if (list2[list1[i]] - list2[list1[j]] < 0):
			measure += 1
		else:
			measure -= 1
			
n = len(list1)			
measure /= n*(n-1)/2.0
print(measure)

' ${LIST1} ${LIST2}
		
