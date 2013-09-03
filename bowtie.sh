#!/bin/bash
#$ -q long
#$ -j y
#$ -cwd
#$ -N bowtie
#$ -l h_vmem=20G
#$ -l virtual_free=20G

#re-map raw reads to cufflinks transcripts

. $HOME/.bashrc

export TEMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

a=(*/*/*1.fastq)
b=(*/*/*2.fastq)
base=$(basename ${a["SGE_TASK_ID"-1]} "_1.fastq")
f=${a["SGE_TASK_ID"-1]}
r=${b["SGE_TASK_ID"-1]}

bowtie2  --very-sensitive-local -p 4 -x merged_asm/merged_spike -1 $f -2 $r |samtools view -Su -F4 - | novosort -t /genefs/MikheyevU/temp --compression 0 -c 2 -i -m 19G -o merged_asm/$base.bam -
