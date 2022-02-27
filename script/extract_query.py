import sys
import getopt
import subprocess


def main(argv):
    inputfile = ''
    outputfile = ''
    queryfile = ''
    querynumber = 0
    try:
        opts, args = getopt.getopt(argv, "hi:o:q:n:",
                                   ["ifile=", "ofile=", "qfile=", "num="])
    except getopt.GetoptError:
        print('extract_query.py -i <inputfile> -o <outputfile> -q <queryfile> '
              '-n <querynumber>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(
                'extract_query.py -i <inputfile> -o <outputfile> -q '
                '<queryfile> -n <querynumber>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-q", "--qfile"):
            queryfile = arg
        elif opt in ("-n", "--qnum"):
            querynumber = int(arg)
    with open(inputfile) as f:
        count = 0
        with open(outputfile, "w") as out, open(queryfile, "w") as outq:
            for line in f:
                if line.strip().split()[0] == 'TOTAL_SAMPLES:':
                    out.write(line)
                    outq.write(line)
                    break
                if count == 0 or count == 1:
                    out.write(line)
                    outq.write(line)
                if count > 1:
                    out.write(line[:-querynumber])
                    restline = "\t".join(line.strip().split()[:-1]) + "\t"
                    query = restline + line.strip().split()[4][-querynumber:]
                    outq.write(query)
                    out.write("\n")
                    outq.write("\n")
                count += 1


if __name__ == "__main__":
    main(sys.argv[1:])
