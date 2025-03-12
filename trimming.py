import sys
import os
import gzip


def print_usage():
    print(f"Usage: python {sys.argv[0]} <inputfile> <bpToKeep | max> "
          "[-trim5 bp] [-flowcellID flowcell] [-addEnd 1 | 2] "
          "[-replace old new | blank] [-renameIDs prefix] [-stdout]")
    sys.exit(1)


def open_file(filename):
    """Opens a file, supporting gzipped formats."""
    if filename.endswith(".gz"):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


def main(argv):
    if len(argv) < 3:
        print_usage()

    input_filename = argv[1]
    trim = "max" if argv[2] == "max" else int(argv[2])
    output_filename = f"{os.path.splitext(os.path.basename(input_filename))[0]}.{trim}mers.fastq"

    # Option flags
    options = {
        "stdout": "-stdout" in argv,
        "flowcellID": argv[argv.index("-flowcellID") + 1] if "-flowcellID" in argv else None,
        "renameIDs": argv[argv.index("-renameIDs") + 1] if "-renameIDs" in argv else None,
        "trim5": int(argv[argv.index("-trim5") + 1]) if "-trim5" in argv else 0,
        "addEnd": argv[argv.index("-addEnd") + 1] if "-addEnd" in argv else None,
        "replace": (
            argv[argv.index("-replace") + 1],
            "" if argv[argv.index("-replace") + 2] == "blank" else argv[argv.index("-replace") + 2]
        ) if "-replace" in argv else None
    }

    # Open output file if not using stdout
    outfile = sys.stdout if options["stdout"] else open(output_filename, "w")

    # Open input file
    input_file = sys.stdin if input_filename == "-" else open_file(input_filename)

    i, j, shorter = 0, 0, 0
    for line in input_file:
        line = line.strip()
        if i == 0 and line.startswith("@"):
            j += 1
            ID = line.replace(" ", "_")
            if options["flowcellID"] and options["flowcellID"] not in line:
                ID = f"@{options['flowcellID']}_{ID[1:]}"
            if options["replace"]:
                ID = ID.replace(*options["replace"])
            if options["renameIDs"]:
                ID = f"@{options['renameIDs']}{j}"
            if options["addEnd"]:
                ID = f"{ID}/{options['addEnd']}"
            i = 1
        elif i == 1:
            sequence = line[options["trim5"]:] if options["trim5"] else line
            sequence = sequence.replace(".", "N")
            sequence = sequence[:trim] if trim != "max" else sequence
            i = 2
        elif i == 2 and line.startswith("+"):
            plus = "+"
            i = 3
        elif i == 3:
            scores = line[options["trim5"]:] if options["trim5"] else line
            scores = scores[:trim] if trim != "max" else scores
            if len(scores) < trim and trim != "max":
                shorter += 1
                continue
            # Write output
            outfile.write(f"{ID}\n{sequence}\n{plus}\n{scores}\n")
            i = 0

    # Close output file if not stdout
    if not options["stdout"]:
        outfile.close()

    if shorter:
        print(f"{shorter} sequences were shorter than the desired length.")


if __name__ == "__main__":
    main(sys.argv)