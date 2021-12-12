import numpy as np
inp = np.loadtxt("matrix3.txt", dtype='i', delimiter=' ')
print(inp)
print(np.transpose(inp))
