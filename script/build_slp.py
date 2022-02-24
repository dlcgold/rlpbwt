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
repair_dirname = os.path.join(bigrepair_dirname, "repair")
largeb_repair_dirname = os.path.join(bigrepair_dirname, "largeb_repair")

repair_exe = os.path.join(repair_dirname, "irepair")
largerepair_exe = os.path.join(largeb_repair_dirname, "largeb_irepair")
bigrepair_exe = os.path.join(bigrepair_dirname, "bigrepair")
despair_exe = os.path.join(repair_dirname, "despair")
integer_despair_exe = os.path.join(repair_dirname, "idespair")
preprocess_exe = os.path.join(bigrepair_dirname, "procdic")
integer_preprocess_exe = os.path.join(bigrepair_dirname, "iprocdic")
postprocess_exe = os.path.join(bigrepair_dirname, "postproc")
integer_postprocess_exe = os.path.join(bigrepair_dirname, "ipostproc")
shaped_slp = os.path.join(shaped_slp_dirname, "SlpEncBuild")

print(dirname)
print(bigrepair_dirname)
print(shaped_slp_dirname)
print(repair_dirname)
print(largeb_repair_dirname)


class extract():
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


class bigrepair():
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


class shapedslp():
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
        if execute_command(command) != True:
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

    matr = extract(inputfile, outname)
    matr.run()

    bigr = bigrepair(outname)
    bigr.run()

    slp = shapedslp(outname, outputfile)
    slp.run()

    os.remove("rs_temp_output")
    os.rmdir("tmp")

if __name__ == "__main__":
    main(sys.argv[1:])
