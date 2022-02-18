import sys
import getopt
import subprocess


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('matrix_to_slp.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('matrix_to_slp.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    with open(inputfile) as f:
        count = 0
        with open(outputfile, 'w') as out:
            col = []
            for line in f:
                if line.strip().split(" ")[0] == 'TOTAL_SAMPLES:':
                    break
                if count > 1:
                    # out.write(line.strip().split("\t")[4])
                    col.append(line.strip().split("\t")[4])
                count += 1
            col = col[::-1]
            print(len(col[0]))
            print(len(col))
            for c in range(0, len(col[0])):
                for r in range(0, len(col)):
                    out.write(col[r][c])
    ## shapedslp


if __name__ == "__main__":
    main(sys.argv[1:])
