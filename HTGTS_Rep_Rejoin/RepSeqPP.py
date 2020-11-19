import sys, os
import argparse

PYTHON = "python3"


def parse_args():
    # parses input arguments with defined flags
    parser = argparse.ArgumentParser(description="Rep seq post-processing script")
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
        "-m",
        dest="meta",
        type=str,
        required=True,
        help="metadata information passed through as comma seperated string",
    )
    args = parser.parse_args()
    return args


def process(args):
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
    # creates output directory if it does not exist
    os.system("mkdir -p %s" % (args.output))

    # takes all directories from input directory and does second step of post-processing with gen_post_proc script
    directory_list = os.listdir("%s" % (args.input))
    for item in directory_list:
        # ensures only relevant items are looked at
        if "." in item or item == "unmatched":
            continue
        elif item in library_dictionary:
            filename = args.input + "/" + item
            meta_info = library_dictionary[item]
            print(
                "{4} {0}/GenPP.py {1} {2} {3}".format(
                    fullpath, filename, args.output, meta_info, PYTHON
                )
            )
            os.system(
                "{4} {0}/GenPP.py {1} {2} {3}".format(
                    fullpath, filename, args.output, meta_info, PYTHON
                )
            )


def main():
    args = parse_args()
    process(args)


main()
