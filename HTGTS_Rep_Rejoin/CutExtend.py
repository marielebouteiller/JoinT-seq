import os, sys
import re

filename = sys.argv[1]
meta_info = sys.argv[2].split(",")

# dictionary that will have the key of a qname
bait_prey_readSeq = {}

# defines the variables of the "v", "j", and the amplicon
# later must be accepted as an input
amplicon = meta_info[14].upper()
bs_seq = meta_info[12] + meta_info[13]
baitAmp = amplicon[: amplicon.find(bs_seq) + 5]
preyAmp = amplicon[amplicon.find(bs_seq) + 5 :]

# opens the output file from the gen_post_proc run
# pulls out only modified reads, reads with >= 10 nucleotides from the breaksite (both ends)
# uses the first/last 5 nt of the concatenated "break" sequence and aligns to bait/prey respectively
# increases depth into concatenated sequence (+1 nt inward) as long as it matches alignment
for line in open(filename, "r"):
    if not line.startswith("Qname"):
        l = line.strip().split("\t")
        if l[38] == "-":
            qname = l[0]
            readSequence = l[34]
            concat = l[37]

            if len(concat) < 10:
                continue

            # finds location of concatenated sequence in read sequence and length of concatenated
            concatInd = readSequence.find(concat)
            concatLength = len(concat)

            # pulls the initial 10 nt string from the "v" and "j" of the read sequence
            baitSeq = readSequence[concatInd - 5 : concatInd + 5]
            preySeq = readSequence[
                concatInd + concatLength - 5 : concatInd + concatLength + 5
            ]

            # define flags that will sequence increment into the breaksite and adding to 10 nt string
            baitFlag = 1
            preyFlag = 1

            # keep iterating +1 nt inward and match to respective half of amplicon
            # does iteration for "v"
            while baitFlag > 0:
                tempBait = readSequence[concatInd - 5 : concatInd + 5 + baitFlag]
                if tempBait in baitAmp:
                    baitSeq = tempBait
                    baitFlag += 1
                else:
                    baitFlag = 0

            # does iteration for "j"
            while preyFlag > 0:
                tempPrey = readSequence[
                    concatInd
                    + concatLength
                    - 5
                    - preyFlag : concatInd
                    + concatLength
                    + 5
                ]
                if tempPrey in preyAmp:
                    preySeq = tempPrey
                    preyFlag += 1
                else:
                    preyFlag = 0
            # gives the value inserted as well

            thing1 = readSequence.find(baitSeq) + len(baitSeq)
            thing2 = readSequence.find(preySeq) + 1

            if (thing1) < (thing2):
                include_insert = readSequence[thing1 : thing2 - 1]
                MH = "NA"
            elif (thing1) == (thing2):
                include_insert = "NA"
                MH = readSequence[thing1 - 1 : thing2]
            else:
                include_insert = "NA"
                difference = thing1 - thing2
                MH = readSequence[thing1 - difference - 1 : thing2 + difference]
            full = readSequence[concatInd - 5 : concatInd + concatLength + 5]
            # put values of the >= 10 nt sequence for both sides of breaksite, the concatenated sequence, and the read sequence
            bait_prey_readSeq[qname] = [
                baitSeq,
                preySeq,
                concat,
                readSequence,
                include_insert,
                MH,
                full,
                thing1,
                thing2,
            ]

# prints header of the output file
print(
    "Qname\tAmplicon_bait_end\tReadSeq_bait_end\tbait_string\tAmplicon_prey_start\tReadSeq_prey_start\tprey_string\tInsert\tMH\tfull\t"
)

# for each read that passed the prior filters
# print qname
# print location of "v end" to amplicon and read sequence, the >=10 nt string of used for alignment of the "v"
# print same but for "j start"
for key in bait_prey_readSeq:

    Amplicon_bait_end = amplicon.find(bait_prey_readSeq[key][0]) + len(
        bait_prey_readSeq[key][0]
    )
    ReadSeq_bait_end = bait_prey_readSeq[key][3].find(bait_prey_readSeq[key][0]) + len(
        bait_prey_readSeq[key][0]
    )
    bait_string = bait_prey_readSeq[key][0]

    Amplicon_prey_start = (
        amplicon.find(bait_prey_readSeq[key][1]) + 1
    )  # I add +1 to be consistent with the other counts. This counts starts at 1, as the others MLB
    ReadSeq_prey_start = bait_prey_readSeq[key][3].find(bait_prey_readSeq[key][1]) + 1
    prey_string = bait_prey_readSeq[key][1]

    insert_seq = bait_prey_readSeq[key][4]
    MH = bait_prey_readSeq[key][5]
    full = bait_prey_readSeq[key][6]

    print(
        "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t"
        % (
            key,
            Amplicon_bait_end,
            ReadSeq_bait_end,
            bait_string,
            Amplicon_prey_start,
            ReadSeq_prey_start,
            prey_string,
            insert_seq,
            MH,
            full,
        )
    )
