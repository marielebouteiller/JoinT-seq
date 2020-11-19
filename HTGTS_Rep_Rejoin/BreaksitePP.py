# test
import os, sys

directory = sys.argv[1]
filename = sys.argv[2]
master = sys.argv[3]
result_TLX = sys.argv[4]
full = sys.argv[5]
meta_info = sys.argv[6]
genome = sys.argv[-1]
head, tail = os.path.split(filename)
tail = tail[:-4]
cwd = os.getcwd()

# pull the full path of all the scripts
pathname = os.path.dirname(sys.argv[0])
fullpath = os.path.abspath(pathname)

# runs the CutExtend script to gather the alignments using a >=20 nt string on both sides of breaksite
print(
    "python {3}/CutExtend.py {2}/{0} {4}> {2}/{1}_CE.xls".format(
        filename, tail, directory, fullpath, meta_info
    )
)
os.system(
    "python {3}/CutExtend.py {2}/{0} {4}> {2}/{1}_CE.xls".format(
        filename, tail, directory, fullpath, meta_info
    )
)

# os.system("python {fullpath}/CutExtend.py {directory}/{filename} {metainfo} > {directory}/{tail}_CE.xls".format(
# 	filename=filename,
# 	tail=tail,
# 	directory=directory,
# 	fullpath=fullpath,
# 	meta_info=meta_info))

# os.system(f"python {fullpath}/CutExtend.py {directory}/{filename} {meta_info} > {directory}/{tail}_CE.xls")


# takes all the qnames and place into a temp file
os.system("awk '{{print $1}}' {1}/{0}_CE.xls > Qnames_CE.txt".format(tail, directory))
# find all lines in master.tlx that have a qname in the qnames file
os.system(
    "grep -Fwf Qnames_CE.txt {0} > {2}/{1}_CE.tlx".format(master, tail, directory)
)
# take first of each of the repeated qnames
os.system(
    "awk '!array[$1]++' {1}/{0}_CE.tlx > {1}/{0}_uniqueQ.tlx".format(tail, directory)
)

# sort the alignment info and the extracted from the master tlx and sitck them together side by side
# os.system("(head -n 1 {1}/{0}_CE.xls && tail -n +2 {1}/{0}_CE.xls | sort) > {1}/{0}_CE_sorted.xls".format(tail,directory))
# os.system("(head -n 1 {1}/{0}_uniqueQ.tlx && tail -n +2 {1}/{0}_uniqueQ.tlx | sort) > {1}/{0}_uniqueQ_sorted.tlx".format(tail,directory))
# os.system('paste --delimiters="\t" {1}/{0}_uniqueQ_sorted.tlx {1}/{0}_CE_sorted.xls > {1}/{0}_uniqueQwVJ.tlx'.format(tail,directory))
os.system(
    "python {2}/SortPaste.py {1}/{0}_uniqueQ.tlx {1}/{0}_CE.xls > {1}/{0}_uniqueQwVJ.tlx".format(
        tail, directory, fullpath
    )
)

# fill all the blank columns with "NA" and edit the different columns in the tlx file
# must go back and make the 61818750 an input
os.system(
    'awk \'BEGIN {{ FS = OFS = "\t" }} {{ for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }}; 1\' {1}/{0}_uniqueQwVJ.tlx > {1}/{0}_filled_uniqueQwVJ.tlx '.format(
        tail, directory
    )
)

# parse the relevant meta data info for this particular library and places in list
meta_info = sys.argv[6].split(",")

# appropirately restructure to fit the tlx file format with the corrected values pulled from the read sequences
# either inlcudes the data with the MH and Insert if prompted, or does not
if full == "1":
    if meta_info[5] == "+":
        os.system(
            f'awk -v OFS=\'\t\' \'{{print $1,$2,"{meta_info[2]}",{meta_info[3]} + $37 - 1 ,"1",{meta_info[3]} + $37 - 1 ," ","{meta_info[2]}",{meta_info[3]}, {meta_info[3]}+ $34 -1,$11,$12,$35 + {len(meta_info[8])},$38 + {len(meta_info[8])},$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$40,$41,$34,$35,$36,$37,$38,$39}}\' {directory}/{tail}_filled_uniqueQwVJ.tlx > {directory}/{tail}_proc_uniqueQwVJ.tlx'
        )
    else:  # in - orientation
        # MLB :I change meta_info[3] into meta_info[4] and I remove -1 for B_Rstart, i.e. {2}-$34 instead of {2}-1-$34, and the length of the barcode to B_Qend and Qstart ($35 and $38)
        # check if the Junction is ok for negative orientation - when I modified Amplicon prey start, I had to adjust the junction for positive orientation -> now this should be ok with {meta_info[4]} - $37
        os.system(
            f'awk -v OFS=\'\t\' \'{{print $1,$2,"{meta_info[2]}",{meta_info[4]} - $37,"-1"," ",{meta_info[4]} - $37,"{meta_info[2]}",{meta_info[4]} - $34,{meta_info[4]} - 1,"-1",$12,$35 + {len(meta_info[8])},$38 + {len(meta_info[8])},$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$40,$41,$34,$35,$36,$37,$38,$39}}\' {directory}/{tail}_filled_uniqueQwVJ.tlx > {directory}/{tail}_proc_uniqueQwVJ.tlx'
        )
    os.system(
        "{{ printf 'Qname\tJuncID\tRname\tJunction\tStrand\tRstart\tRend\tB_Rname\tB_Rstart\tB_Rend\tB_Strand\tB_Qstart\tB_Qend\tQstart\tQend\tQlen\tB_Cigar\tCigar\tSeq\tJ_Seq\tBarcode\tunaligned\tbaitonly\tuncut\tmisprimed\tfreqcut\tlargegap\tmapqual\tbreaksite\tsequential\trepeatseq\tduplicate\tInsert\tMH\tAmplicon_bait_end\tReadSeq_bait_end\tbait_string\tAmplicon_prey_start\tReadSeq_prey_start\tprey_string\n'; cat {1}/{0}_proc_uniqueQwVJ.tlx | tail -n +2; }} > {1}/{0}_proc2_uniqueQwVJ.tlx".format(
            tail, directory
        )
    )
else:
    if meta_info[5] == "+":
        os.system(
            'awk -v OFS=\'\t\' \'{{print $1,$2,"{3}",{2} + $37,"1",{2} + $37," ","{3}",{2}, {2} + $34,$11,$12,$35,$38,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$34,$35,$36,$37,$38,$39}}\' {1}/{0}_filled_uniqueQwVJ.tlx > {1}/{0}_proc_uniqueQwVJ.tlx'.format(
                tail, directory, meta_info[3], meta_info[2]
            )
        )
    else:
        os.system(
            'awk -v OFS=\'\t\' \'{{print $1,$2,"{3}",{2} - 1 - $37,"-1"," ",{2} - 1 - $37,"{3}",{2} - 1 - $34,{2} - 1,"-1",$12,$35,$38,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$34,$35,$36,$37,$38,$39}}\' {1}/{0}_filled_uniqueQwVJ.tlx > {1}/{0}_proc_uniqueQwVJ.tlx'.format(
                tail, directory, meta_info[3], meta_info[2]
            )
        )
    os.system(
        "{{ printf 'Qname\tJuncID\tRname\tJunction\tStrand\tRstart\tRend\tB_Rname\tB_Rstart\tB_Rend\tB_Strand\tB_Qstart\tB_Qend\tQstart\tQend\tQlen\tB_Cigar\tCigar\tSeq\tJ_Seq\tBarcode\tunaligned\tbaitonly\tuncut\tmisprimed\tfreqcut\tlargegap\tmapqual\tbreaksite\tsequential\trepeatseq\tduplicate\tAmplicon_bait_end\tReadSeq_bait_end\tbait_string\tAmplicon_prey_start\tReadSeq_prey_start\tprey_string\n'; cat {1}/{0}_proc_uniqueQwVJ.tlx | tail -n +2; }} > {1}/{0}_proc2_uniqueQwVJ.tlx".format(
            tail, directory
        )
    )


# append the result tlx to the end of our processed file and then take first of each of the repeated qnames
os.system(
    "cat {2}/{0}_proc2_uniqueQwVJ.tlx {1} > {2}/{0}_combined.tlx".format(
        tail, result_TLX, directory
    )
)
os.system(
    "awk '!array[$1]++' {1}/{0}_combined.tlx > {1}/{0}_plotting.tlx".format(
        tail, directory
    )
)

# appends the MH and Insert data for all the translocated reads
if full == "1":
    os.system(
        "python {2}/transloc_fill.py {1}/{0}_plotting.tlx".format(
            tail, directory, fullpath
        )
    )
    os.system("mv {1}/{0}_plotting2.tlx {1}/{0}_plotting.tlx".format(tail, directory))

# create the limit conditionals for where the b_cigar data will be used to effectively parse
# lower_conditional = (
#     meta_info[14].upper().find(meta_info[12] + meta_info[13]) + len(meta_info[8]) + 1
# )
lower_conditional = meta_info[14].upper().find(meta_info[12] + meta_info[13]) + 1
upper_conditional = lower_conditional + 10
print(
    "python {2}/b_cigar.py {1}/{0}_plotting.tlx {3} {4} > {1}/{0}_bcigar.tlx ".format(
        tail, directory, fullpath, lower_conditional, upper_conditional
    )
)
os.system(
    "python {2}/b_cigar.py {1}/{0}_plotting.tlx {3} {4} > {1}/{0}_bcigar.tlx ".format(
        tail, directory, fullpath, lower_conditional, upper_conditional
    )
)

# create bedgraphs from all the plotting tlxs
# copy over TLX2BED
print(
    "python {3}/tlx2bed.py -f {1}/{0}_bcigar.tlx -g {2} --v3".format(
        tail, directory, genome, fullpath
    )
)
os.system(
    "python {3}/tlx2bed.py -f {1}/{0}_bcigar.tlx -g {2} --v3".format(
        tail, directory, genome, fullpath
    )
)

# remove intermediary files
os.system(
    "rm -f {1}/Qnames_CE.txt {1}/{0}_CE.tlx {1}/{0}_CE.xls {1}/{0}_uniqueQ.tlx {1}/{0}_CE_sorted.xls {1}/{0}_uniqueQ_sorted.tlx {1}/{0}_uniqueQwVJ.tlx {1}/{0}_filled_uniqueQwVJ.tlx {1}/{0}_proc_uniqueQwVJ.tlx {1}/{0}_proc2_uniqueQwVJ.tlx {1}/{0}_combined.tlx".format(
        tail, directory
    )
)
os.system(
    "rm -f {1}/{0}_plotting.tlx {1}/{0}.xls {1}/{0}_unmod.xls".format(
        tail, directory
    )
)
shorttail = tail[:-3]
os.system(
    "rm -f {1}/{0}_unmod.xls".format(
        shorttail, directory
    )
)
