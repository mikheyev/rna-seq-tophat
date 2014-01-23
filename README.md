# Gene expression associated with haplo- and pleiometrotic founding in *Pogonomyrmex californicus*

## Preparing reads

move from SQC folder (**move.sh**), then **trim.sh** to perform quality trimming

## Assemble transcriptome

**trinity.sh** gives us 311,726 transcripts

## Mapping reads back with rsem

First, prepare reference
```
/apps/MikheyevU/sasha/trinity/trinity/util/RSEM_util/run_RSEM_align_n_estimate.pl  --transcripts trinity_out_dir/Trinity.fasta --just_prep_reference 
```
Next, use **rsem.sh** to map reads.


## BLAST search of transcripts vs nr

### Determine ERCC transcripts

```
makeblastdb -in ERCC92.fa -out ERCC -dbtype nucl
blastn -db ERCC -outfmt 6 -ungapped -word_size 20 -max_target_seqs 1 -evalue 1e-100 -query ../../trinity_out_dir/Trinity.fasta -out hits.txt

```

### Find hits to nr

