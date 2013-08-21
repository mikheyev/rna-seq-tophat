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

a=(*/*/*1.fastq)
b=(*/*/*2.fastq)
base=$(basename ${a["SGE_TASK_ID"-1]} "_1.fastq")
f=${a["SGE_TASK_ID"-1]}
r=${b["SGE_TASK_ID"-1]}

# try using the pbar reference gene set annotations
#tophat2 -z pigz --G ../pogo/final/pbar_OGSv1.2.gff3 --b2-very-sensitive -p 2 -o rnaseq/$base ../pogo/Pbar $f $r 

##extract unmapped reads

if [[ ! -e rnaseq2 ]]; then
    mkdir rnaseq2
fi

novosort -t /genefs/MikheyevU/temp -c 2 --compression 0 --namesort rnaseq/$base/unmapped.bam > rnaseq2/$base.temp.bam
bam bam2FastQ --in rnaseq2/$base.temp.bam  --readname --outBase rnaseq2/$base 2> /dev/null
rm rnaseq2/$base.temp.bam