import sys
import numpy as np
inp = np.loadtxt("out.txt", dtype='i')
#print(inp)
t = np.transpose(inp)
for r in t:
	for c in r:
		print(c, end = ' ')
	print()
