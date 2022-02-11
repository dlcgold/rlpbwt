# import pandas as pd
#
# df = pd.read_csv('../input/sample_new.txt',
#                  sep='\t',
#                  skiprows=2,
#                  header=None,
#                  names=["col"],
#                  usecols=[4])
# # df.drop(df.head(1).index, inplace=True)
# df.drop(df.tail(1).index,
#         inplace=True)
# #df.drop(df.columns[[0, 1, 2, 3]], axis = 1, inplace = True)
# print(df)

with open('../input/sample_new.txt') as f:
    count = 0
    with open('../output/sample.fa', 'w') as out:
        #out.write("header:\n")
        tra = []
        for line in f:
            if line.strip().split(" ")[0] == 'TOTAL_SAMPLES:':
                break
            if count > 1:
                #out.write(line.strip().split("\t")[4])
                tra.append(line.strip().split("\t")[4])
            count += 1
        for i in range(0, len(tra[0])):
            for j in range(0, len(tra)):
                #print(tra[j][i], end=" ")
                out.write(tra[j][i])
            #print()