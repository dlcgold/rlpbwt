import numpy as np
import pandas as pd

filein = "../../Documenti/macs/11k.macs"
# filein = "../input/sample_new.txt"
print("reading file")
dft = pd.read_csv(filein,
                  skiprows=2,
                  sep="\t",
                  names=["col"],
                  usecols=[4])
dft.dropna()
print(dft)
# df = pd.read_csv(filein,
#                  skiprows=2,
#                  sep="\t",
#                  names=["col"],
#                  usecols=[4],
#                  converters={'col': lambda x: np.array(list(x),
#                                                        dtype=np.byte)})[:-1]
# df.drop(len(df.index) - 1)
# print("traspose")
# m = np.flip(np.stack(df['col']).transpose(), 1)
# # print(m[0])
# print("writing")
# np.savetxt("../output/test.txt", m, fmt='%i', delimiter='', newline='')
