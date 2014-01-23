from Bio import SeqIO

for rec in SeqIO.parse(open("data/trinity_out_dir/Trinity_noERCC.fasta"),"fasta"):
	if rec.id[0:4] == "ERCC":
		gene = rec.id
	else:
		gene = "_".join(rec.id.split("_")[0:2])
		
	print "%s\t%s" % (gene,rec.id)

