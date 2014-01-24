import gzip,pdb
from Bio import SeqIO

infile = open("data/keep.txt")
infile.next()
components = {}
for line in infile:
	components[line.split()[1].replace('"',"")] = ""

#go through the fasta file and find the longest component
for rec in SeqIO.parse(gzip.open("data/trinity_out_dir/Trinity.fasta.gz"),"fasta"):
	comp = "_".join(rec.id.split("_")[:2])
	if comp in components:
		if not components[comp] or components[comp][1] < len(rec):
			components[comp] = [rec.id,len(rec)]

for rec in SeqIO.parse(gzip.open("data/trinity_out_dir/Trinity.fasta.gz"),"fasta"):
	comp = "_".join(rec.id.split("_")[:2])
	if comp in components and rec.id == components[comp][0]:
		rec.description = rec.id
		rec.id = comp
		print rec.format("fasta"),
