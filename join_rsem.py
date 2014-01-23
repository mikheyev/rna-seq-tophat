from glob import glob
from collections import defaultdict
import pdb
#combine results of RSEM analysis
files = glob("data/rsem/*genes.results")
names = [f.split("/")[-1].split(".")[0] for f in files]

outfile = open("data/rsem/rsem.csv","w")
outfile.write("id,gene,count,fpkm\n")
for f,name in zip(files,names):
	infile=open(f)	
	infile.next()
	for line in infile:
		line=line.rstrip().split()
		outfile.write("%s,%s,%s,%s\n" % (name,line[0],line[4],line[-1]))
	infile.close()
outfile.close()

# counts = defaultdict(list)
# fpkm = defaultdict(list)

# genes = []
# first = 1
# for idx,f in enumerate(files):
# 	infile=open(f)
# 	infile.next()
# 	for line in infile:
# 		line = line.rstrip().split()
# 		if first:
# 			genes.append(line[0])
# 		counts[line[0]].append(line[4])
# 		fpkm[line[0]].append(line[-1])
# 	first = 0	
# 	infile.close()


# counts_file = open("data/rsem/counts.csv","w")
# counts_file.write(",".join(["gene"] + names)+"\n")
# fpkm_file = open("data/rsem/fpkm.csv","w")
# fpkm_file.write(",".join(["gene"] + names)+"\n")
# for i in genes:
# 	counts_file.write(",".join([i]+counts[i])+"\n")
# 	fpkm_file.write(",".join([i]+fpkm[i])+"\n")
# counts_file.close()
# fpkm_file.close()

