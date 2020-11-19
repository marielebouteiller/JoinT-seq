import os, sys
import os.path
from Bio import SeqIO

filename = sys.argv[1]
output_folder = sys.argv[2]
meta_info = sys.argv[3].split(",")
head, tail = os.path.split(filename)

USE_UNJOIN_READS = True

print("Starting the general post-proccessing step for {0}".format(filename))
# open the gzipped joined fasta file and place sequence as a value in a dictionary with qname as a key
link_dict = {}


# i removed this line if it is already unzipped
os.system("gunzip %s/reads_fasta/%s_join.fa.gz" % (filename, tail))


if USE_UNJOIN_READS:
    os.system("gunzip %s/reads_fasta/%s_unjoinR1.fa.gz" % (filename, tail))

    # Merge the join and unjoinR1 fasta files
    # open the join and unjoin files, create a new file which concatenates both
    with open(f"{filename}/reads_fasta/{tail}_join.fa", "r") as join_fasta_file, open(
        f"{filename}/reads_fasta/{tail}_unjoinR1.fa", "r"
    ) as unjoinR1_fasta_file, open(
        f"{filename}/reads_fasta/{tail}_join_unjoinR1.fa", "w"
    ) as merge_fasta_file:

        # copy the join.fasta file as is
        for l_join in join_fasta_file:
            merge_fasta_file.write(l_join)

        # copy the unjoinR1 file as is
        for l_unjoinR1 in unjoinR1_fasta_file:
            merge_fasta_file.write(l_unjoinR1)

    # Merge the join and unjoinR1 result XLS files
    root_path = filename + "/igblast_results/" + tail
    with open(root_path + "_join.all.xls", "r") as join_xls, open(
        root_path + "_unjoinR1.all.xls", "r"
    ) as unjoin_xls, open(root_path + "_join_unjoinR1.all.xls", "w") as merge_xls:

        # copy the join file as is
        for l_join in join_xls:
            merge_xls.write(l_join)

        # skip the first line of the unjoin file, and copy the rest as is
        next(unjoin_xls)
        for l_unjoinR1 in unjoin_xls:
            merge_xls.write(l_unjoinR1)

fasta_file = (
    "%s/reads_fasta/%s_join_unjoinR1.fa" % (filename, tail)
    if USE_UNJOIN_READS
    else "%s/reads_fasta/%s_join.fa" % (filename, tail)
)


result_file = (
    filename + "/igblast_results/" + tail + "_join_unjoinR1.all.xls"
    if USE_UNJOIN_READS
    else filename + "/igblast_results/" + tail + "_join.all.xls"
)


for seq_record in SeqIO.parse(fasta_file, "fasta"):
    link_dict[str(seq_record.id)] = str(seq_record.seq)

# define the "breaksite" sequence and preset the number of unmodified and total sequences
bs_seq = meta_info[12] + meta_info[13]
amplicon = meta_info[14].upper()
plus_count = 0.0
total_count = 0.0
minus_count = 0.0

# open the output file to write to
with open(output_folder + "/" + tail + "_PP.xls", "w") as output_file:
    # open the joined results file and begins the processing and filtering on each read
    for line in open(result_file, "r"):
        # if header line, append new field and write to file
        if line.startswith("Qname"):
            output_file.write(
                "%s\t%s\t%s\t%s\t%s\t%s\n"
                % (
                    line.strip(),
                    "sequence",
                    "V_end_loc",
                    "J_start_loc",
                    "VDJ_concat",
                    "Pass(+/-)",
                )
            )
            continue

        # split the line reads and assign to variables
        l = line.strip().split("\t")
        a = line.strip()
        Qname = l[0]
        V_end = l[7]
        V_D_Junc = l[8]
        D_Junc = l[9]
        D_J_Junc = l[10]
        V_J_Junc = l[11]
        J_start = l[12]
        sequence = "NOPE"

        # concatenate everything from the "vend" to the "jstart" columns
        # then check if the 10nt breaksite sequence we are interested can be found in it
        # set as modified or unmodified accordingly, incrementing the unmodified counter as needed
        middle = ""
        if V_D_Junc != "N/A" and V_D_Junc != "-":
            if V_D_Junc[0] == "(":
                middle = middle + V_D_Junc[1:-1]
            else:
                middle = middle + V_D_Junc
        if D_Junc != "N/A" and D_Junc != "-":
            if D_Junc[0] == "(":
                middle = middle + D_Junc[1:-1]
            else:
                middle = middle + D_Junc
        if D_J_Junc != "N/A" and D_J_Junc != "-":
            if D_J_Junc[0] == "(":
                middle = middle + D_J_Junc[1:-1]
            else:
                middle = middle + D_J_Junc
        if V_J_Junc != "N/A" and V_J_Junc != "-":
            if V_J_Junc[0] == "(":
                middle = middle + V_J_Junc[1:-1]
            else:
                middle = middle + V_J_Junc
        if J_start != "N/A":
            concat = V_end + middle + J_start
        else:
            continue

        # pull the sequence from the dictionary with the qname
        sequence = link_dict[Qname]
        # find location of the concatenated sequence and its length
        concatInd = sequence.find(concat)
        concatLength = len(concat)
        # extend the bait and prey sequences to 20 nt from the breaksite and remove read if it doesn't reach that length
        vSeq = sequence[concatInd - 5 : concatInd + 5]
        jSeq = sequence[concatInd + concatLength - 5 : concatInd + concatLength + 5]

        if len(vSeq) < 10 or len(jSeq) < 10:
            continue

        if vSeq in amplicon and jSeq in amplicon:
            # checks if modified
            if bs_seq in concat:
                Pass = "+"
                plus_count += 1.0
            else:
                Pass = "-"
                minus_count += 1.0

            # find rough location "vend" and "jstart"
            v_loc = 0
            j_loc = 0

            # increment total count
            total_count += 1.0
            # write the file line with new information appended to the ended
            output_file.write(
                "%s\t%s\t%s\t%s\t%s\t%s\n" % (a, sequence, v_loc, j_loc, concat, Pass)
            )

# creates the qname and sequence file for the unodified reads AFTER all filtering
with open(output_folder + "/" + tail + "_unmod.xls", "w") as unmod_file:
    for line in open(output_folder + "/" + tail + "_PP.xls", "r"):
        if not line.startswith("Qname"):
            l = line.strip().split("\t")
            if l[38] == "+":
                qname = l[0]
                readSequence = l[34]
                unmod_file.write("%s\t%s\n" % (qname, readSequence))

# create a values file and do calculations of unmodfied and total count, and percent modified
percent = (1 - (float(plus_count) / float(total_count))) * 100
if not (os.path.isfile("%s/values.xls" % output_folder)):
    with open("%s/values.xls" % output_folder, "w") as values_file:
        # values_file.write("library\tnon-modified\ttotal\tpercent_modified\n")
        values_file.write("library\tmodified\tunmodified\ttotal\tpercent_unmodified\n")

if os.path.isfile("%s/values.xls" % output_folder):
    with open("%s/values.xls" % output_folder, "a") as values_file:
        values_file.write(
            "%s\t%d\t%d\t%d\t%f\n"
            % (tail, minus_count, plus_count, total_count, percent)
        )

print("Succesfully completed the post-proccessing step for {0}".format(filename))
