import sys, time, argparse, subprocess, os.path, threading, tempfile, getopt


def execute_command(command):
    try:
        subprocess.check_call(command.split())
    except subprocess.CalledProcessError:
        print("Error executing command line:", flush=True)
        return False
    return True


dirname = '../build/' + os.path.dirname(os.path.relpath(__file__))
currdir = 'tmp/'
bigrepair_dirname = os.path.join(dirname, "_deps/bigrepair-src")
shaped_slp_dirname = os.path.join(dirname, "_deps/shaped_slp-build")


class extract:
    def __init__(self, filename, output):
        self.filename = filename
        self.output = output

    def run(self):
        filename = self.filename
        output = self.output
        exe = "tsv_to_single_row.py"
        command = "python {exe} -i {file} -o {currdir}{out}.txt".format(
            exe=os.path.join('.', exe), file=filename, currdir=currdir,
            out=output)
        if not execute_command(command):
            return


class bigrepair:
    def __init__(self, filename):
        self.filename = filename

    def run(self):
        filename = self.filename
        exe = "bigrepair"
        command = "{exe} {currdir}{file}.txt".format(
            exe=os.path.join(bigrepair_dirname,
                             exe), currdir=currdir, file=filename)
        if not execute_command(command):
            return


class shapedslp:
    def __init__(self, filename, output):
        self.filename = filename
        self.output = output

    def run(self):
        filename = self.filename
        outfile = self.output
        grammar = "SelfShapedSlp_SdSd_Sd"

        exe = "SlpEncBuild"
        command = "{exe} -i {curr}{file}.txt -o {out} -e {g} -f Bigrepair".format(
            exe=os.path.join(shaped_slp_dirname, exe), curr=currdir,
            file=filename, out=outfile, g=grammar)

        print("==== ShapedSLP construction.\nCommand:", command, flush=True)
        if not execute_command(command):
            return


def main(argv):
    inputfile = ''
    outputfile = ''
    if not os.path.isdir("tmp"):
        os.mkdir("tmp")
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('build_slp.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('build_slp.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    outdir = outputfile.split('/')[0] + "/" + "".join(
        outputfile.split('/')[1::-2]) + "/"
    outname = outputfile.split('/')[-1].split('.')[0]

    print("extract matrix")
    matr = extract(inputfile, outname)
    matr.run()

    print("end extract matrix \nbegin bigrepair")
    bigr = bigrepair(outname)
    bigr.run()

    print("end bigrepair \nbegin shapedslp")
    slp = shapedslp(outname, outputfile)
    slp.run()
    print("end shapedslp \nclearing")
    if os.path.isfile("rs_temp_output"):
        os.remove("rs_temp_output")

    if os.path.isdir("tmp"):
        directory = "tmp"
        for filename in os.listdir(directory):
            f = os.path.join(directory, filename)
            if os.path.isfile(f):
                os.remove(f)

        os.rmdir("tmp")


if __name__ == "__main__":
    main(sys.argv[1:])
