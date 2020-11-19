import os, sys

master_tlx = sys.argv[1]
aligned_tlx = sys.argv[2]

alignment_dict = {}

for line in open(aligned_tlx):
	if not line.startswith("Qname"):
		l = line.strip().split('\t')
		alignment_dict[l[0]] = '\t'.join(l)
	else:
		header_append = line[:-1]


to_print = []
for line in open(master_tlx):
	if not line.startswith("Qname"):
		l = line.strip().split('\t')
		if l[0] in alignment_dict:
			to_print.append('\t'.join(l) + "\t" + alignment_dict[l[0]])
	else:
		to_print.append(line[:-1]+"\t"+header_append)

for item in to_print:
	print(item)
