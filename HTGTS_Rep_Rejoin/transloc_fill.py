import os, sys

#take in inputs: the tlx file which we are appending the MH and Insert to
filename = sys.argv[1]

#open the outfile to write to
outfile = filename.replace("plotting.tlx","plotting2.tlx")
f = open(outfile,"w")

#opens original plotting tlx file and calculates the MH and Insert
#done based on the B_Qend and Qstart coordinates on the sequence
for line in open(filename):
	if line.startswith("Qname"):
		f.write(line)
	else:
		l = line.strip().split("\t")
		if len(l) < 18:
			continue
		sequence = l[18]
		if len(l) < 38:
			thing1 = int(l[12])
			thing2 = int(l[13])

			if (thing1) < (thing2):
				include_insert = sequence[thing1:thing2-1]
				if len(include_insert) == 0:
					include_insert = "NA"
				MH = "NA"
			elif (thing1) == (thing2):
				include_insert = "NA"
				MH = sequence[thing1-1:thing2]
			else:
				include_insert = "NA"
				difference = thing1-thing2
				MH = sequence[thing1-difference-1:thing2 + difference]
			f.write("%s\t%s\t%s\n" %(line[:-1],include_insert,MH) )
		else:
			f.write(line)
