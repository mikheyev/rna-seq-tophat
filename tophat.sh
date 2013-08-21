#!/bin/bash
#$ -q long
#$ -j y
#$ -cwd
#$ -N tophat
#$ -l h_vmem=4G
#$ -l virtual_free=4G

. $HOME/.bashrc
a=(*/*/*1.fastq)
b=(*/*/*2.fastq)
#SGE_TASK_ID=1
base=$(basename ${a["SGE_TASK_ID"-1]} "_1.fastq")
f=${a["SGE_TASK_ID"-1]}
r=${b["SGE_TASK_ID"-1]}

tophat2 -z pigz --b2-very-sensitive -p 2 -o rnaseq/$base ../pogo/Pbar $f $r 