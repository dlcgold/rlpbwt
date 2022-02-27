import sys
import getopt
import subprocess


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('tsv_to_single_row.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('tsv_to_single_row.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    with open(inputfile) as f:
        count = 0
        with open(outputfile, 'w') as out:
            col = []
            print("read file")
            for line in f:
                c = str(count)
                print(c, end="\r")
                if line.strip().split()[0] == 'TOTAL_SAMPLES:':
                    break
                if count > 1:
                    # out.write(line.strip().split("\t")[4])
                    col.append(line.strip().split()[4])
                count += 1
            print("end read file")
            print("inverse matrix")
            col = col[::-1]
            print("end inverse matrix of measures: ")
            print(len(col[0]))
            print(len(col))
            print("begin \"stretch\" matrix")
            ## TODO too slow
            for c in range(0, len(col[0])):
                print(c, end="\r")
                for r in range(0, len(col)):
                    out.write(col[r][c])
            print("end \"stretch\" matrix")


if __name__ == "__main__":
    main(sys.argv[1:])
