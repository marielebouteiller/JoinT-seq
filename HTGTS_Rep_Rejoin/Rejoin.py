import sys, os
import argparse


def parse_args():

    # parses input arguments with defined flags
    parser = argparse.ArgumentParser(
        description="Rep seq complete post-processing script"
    )
    parser.add_argument(
        "-od",
        dest="output",
        type=str,
        required=True,
        help="output directory path: where the post-processed data will be stored",
    )
    parser.add_argument(
        "-id",
        dest="input",
        type=str,
        required=True,
        help="input directory path: directory of all individual runs by default",
    )
    parser.add_argument(
        "-r",
        dest="results",
        type=str,
        required=True,
        help="results path: directory where the master and results tlx are",
    )
    parser.add_argument(
        "-m",
        dest="meta",
        type=str,
        required=True,
        help="metadata file associated with igblast run",
    )
    parser.add_argument(
        "-g",
        dest="genome",
        type=str,
        required=True,
        help="species that will be used to generate bedgraphs (mouse/human)",
    )
    parser.add_argument(
        "--full",
        dest="full",
        action="store_const",
        const=1,
        default=0,
        help="set flag to inlcude MH and Insert data in final output",
    )

    args = parser.parse_args()
    return args


def first(args):
    print("Starting the Rejoin process!")

    # pull the full path of all the scripts
    pathname = os.path.dirname(sys.argv[0])
    fullpath = os.path.abspath(pathname)

    library_dictionary = {}
    for line in open(args.meta):
        if line.startswith("Library"):
            continue
        else:
            l = line.strip().split("\t")
            library = l[0]
            sequencing = l[1]
            key = library + "_" + sequencing
            library_dictionary[key] = ",".join(l[2:])
    print("Metadata file was successfully read in.")

    # runs the initial processing of the joinged reads from igblast
    os.system(
        "python {2}/RepSeqPP.py -id {0} -od {1} -m {3}".format(
            args.input, args.output, fullpath, args.meta
        )
    )
    directory_list = os.listdir("%s" % (args.input))
    for item in directory_list:
        # ensures only relevant items are looked at
        if "." in item or item == "unmatched":
            continue
        elif item in library_dictionary:

            head, tail = os.path.split(item)

            directory = args.output
            filename = item + "_PP.xls"
            master = args.results + "/" + tail + "/" + tail + ".tlx"
            result_TLX = args.results + "/" + tail + "/" + tail + "_result.tlx"
            full = args.full
            meta_info = library_dictionary[item]
            genome = args.genome

            os.system(
                "python {5}/BreaksitePP.py {0} {1} {2} {3} {4} {6} {7}".format(
                    directory,
                    filename,
                    master,
                    result_TLX,
                    full,
                    fullpath,
                    meta_info,
                    genome,
                )
            )
    os.system(
        "(head -n 1 {0}/values.xls && tail -n +2 {0}/values.xls | sort) > {0}/sorted_values.xls".format(
            args.output
        )
    )
    os.system("rm -f {0}/values.xls".format(args.output))
    print("Breaksite analysis and rejoining completed for all libraries!")


def main():
    args = parse_args()
    first(args)


main()

