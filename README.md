# Gene expression associated with haplo- and pleiometrotic founding in *Pogonomyrmex californicus*

## Preparing reads

move from SQC folder (**move.sh**), then **trim.sh** to perform quality trimming

## Assemble transcriptome

**trinity.sh** gives us 311,726 transcripts

## Mapping reads back with rsem

First, prepare reference
```
python  component_to_trans_map.py > data/trinity_out_dir/gene2isoform.txt
rsem-prepare-reference --transcript-to-gene-map data/trinity_out_dir/gene2isoform.txt --no-polyA data/trinity_out_dir/Trinity_noERCC.fasta data/trinity_out_dir/Trinity_noERCC
```
Next, use **rsem.sh** to map reads. The reads are aggregated by **join_rsem.py**.

## BLAST search of transcripts vs nr

first, use **filter_assembly.py** to sub-select only sufficiently abundant components (see **californicus.R** for details)

```
blastx -num_threads 1 -max_target_seqs 10 -db /genefs/MikheyevU/NCBI/nr -evalue 0.0001  -outfmt 5 -query ${a["$SGE_TASK_ID"-1]} -out  $(basename ${a["$SGE_TASK_ID"-1]} "fa")"xml"
```

