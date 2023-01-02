import sys
import getopt
from pysam import VariantFile


def main(argv):
    inputfile = ''
    outputfile = ''
    vcf = False
    try:
        opts, args = getopt.getopt(argv, "hi:o:v", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('tsv_to_single_row.py -i <inputfile> -o <outputfile> -v/--vcf')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('tsv_to_single_row.py -i <inputfile> -o <outputfile> -v/--vcf')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-v", "--vcf"):
            vcf = True
    if not vcf:
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
                for c in range(0, len(col[0])):
                    print(c, end="\r")
                    for r in range(0, len(col)):
                        out.write(col[r][c])
                print("end \"stretch\" matrix")
    else:
        bcf_in = VariantFile(inputfile)  
        col = []    
        count = 0
        for rec in bcf_in.fetch(): 
            c = str(count)
            print(c, end="\r")
            col.append("".join([f"{s['GT'][0]}{s['GT'][1]}" for s in rec.samples.values()]))
            count += 1
        print("end read file")
        print("inverse matrix")
        col = col[::-1]
        print("end inverse matrix of measures: ")
        print(len(col[0]))
        print(len(col))
        print("begin \"stretch\" matrix")
        with open(outputfile, 'w') as out:
            for c in range(0, len(col[0])):
                print(c, end="\r")
                tmp = ""
                for r in range(0, len(col)):
                    tmp += col[r][c]
                out.write(tmp)
        print("end \"stretch\" matrix")

if __name__ == "__main__":
    main(sys.argv[1:])
