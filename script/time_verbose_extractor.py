import sys
import getopt


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('time_verbose_extractor.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('time_verbose_extractor.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    with open(inputfile, "r") as f:
        with open(outputfile, "w") as out:
            for line in f:
                line = line[1:-1]
                tokens = line.split(sep=":")
                if tokens[0] == "Command being timed":
                    # out.write(tokens[1].lstrip())
                    # out.write("\n")
                    if './rlpbwt' in tokens[1]:
                        out.write("rlpbwt\n")
                        if '-N' in tokens[1]:
                            out.write("naive\n")
                        if '-B' in tokens[1]:
                            out.write("bitvectors\n")
                        if '-P' in tokens[1]:
                            if '-e' in tokens[1]:
                                out.write("panel extended\n")
                            else:
                                out.write("panel\n")
                        if '-S' in tokens[1]:
                            if '-e' in tokens[1]:
                                out.write("slp extended\n")
                            else:
                                out.write("slp\n")
                    if './pbwt' in tokens[1]:
                        out.write("pbwt\n")
                if tokens[0] == "User time (seconds)":
                    out.write(tokens[1].lstrip())
                    out.write("\n")
                if tokens[0] == "System time (seconds)":
                    out.write(tokens[1].lstrip())
                    out.write("\n")
                if tokens[0] == "Maximum resident set size (kbytes)":
                    out.write(tokens[1].lstrip())
                    out.write("\n")


if __name__ == "__main__":
    main(sys.argv[1:])
