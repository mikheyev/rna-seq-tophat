#!/bin/bash
#$ -q long
#$ -j y
#$ -cwd
#$ -N cutadapt
#$ -l h_vmem=4G
#$ -l virtual_free=4G

. $HOME/.bashrc

export TEMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

a=(raw/*/*.fastq)
base=$(basename ${a["SGE_TASK_ID"-1]} "_1.fastq")
f=${a["SGE_TASK_ID"-1]}

cutadapt --match-read-wildcards -a CTGTCTCTTATACACATCT -g AGATGTGTATAAGAGACAG -g GGTATCAACGCAGAGTACATGGGNNNN -a NNNNCCCATGTACTCTGCGTTGATACCACTG -g CGGCCGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTN $f -o trimmed/$(basename $f)

#SGE_TASK_ID=1
#a=(raw/*/*1.fastq)
#b=(raw/*/*2.fastq)
#base=$(basename ${a["SGE_TASK_ID"-1]} "_1.fastq")
#f=${a["SGE_TASK_ID"-1]}
#r=${b["SGE_TASK_ID"-1]}

#cutadapt --minimum-length 50 --match-read-wildcards -a CTGTCTCTTATACACATCT -g AGATGTGTATAAGAGACAG -g GGTATCAACGCAGAGTACATGGGNNNN -a NNNNCCCATGTACTCTGCGTTGATACCACTG -g CGGCCGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTN --paired-output trimmed/$(basename $r).tmp -o timmed/$(basename $f).tmp $f $r
