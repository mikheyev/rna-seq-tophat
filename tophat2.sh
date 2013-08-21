#!/bin/bash
#$ -q long
#$ -j y
#$ -cwd
#$ -N tophat
#$ -l h_vmem=20G
#$ -l virtual_free=20G

. $HOME/.bashrc

export TEMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

a=(rnaseq2/*1.fastq)
b=(rnaseq2/*2.fastq)
c=$(rnaseq2/*.fastq | awk 'NR%3==0')
base=$(basename ${a["SGE_TASK_ID"-1]} "_1.fastq")
f=${a["SGE_TASK_ID"-1]}
r=${b["SGE_TASK_ID"-1]}
u=${c["SGE_TASK_ID"-1]}

# try using the pbar reference gene set annotations
tophat2 -z pigz --G ../pogo/final/pbar_OGSv1.2.gff3 -j junctions.merged.junc --no-novel-junc --b2-very-sensitive -p 2 -o rnaseq2/$base"_paired" ../pogo/Pbar $f $r 
tophat2 -z pigz --G ../pogo/final/pbar_OGSv1.2.gff3 -j junctions.merged.junc --no-novel-junc --b2-very-sensitive -p 2 -o rnaseq2/$base"_single" ../pogo/Pbar u 

##extract unmapped reads

if [[ ! -e unmapped ]]; then
    mkdir unmapped
fi

novosort -t /genefs/MikheyevU/temp -c 2 --compression 0 --namesort rnaseq/$base/unmapped.bam > unmapped/$base.temp.bam
bam bam2FastQ --in unmapped/$base.temp.bam  --readname --outBase unmapped/$base 2> /dev/null
rm unmapped/$base.temp.bam