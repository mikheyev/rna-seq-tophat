#!/bin/bash
#$ -q short
#$ -j y
#$ -cwd
#$ -l h_vmem=4G
#$ -l virtual_free=4G
#$ -N trim
. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

a=(rnaseq2/*1.fastq)
b=(rnaseq2/*2.fastq)
c=(`ls rnaseq2/*.fastq | awk 'NR%3==0' |tr "\n" " "`)
base=$(basename ${a["SGE_TASK_ID"-1]} "_1.fastq")
f=${a["SGE_TASK_ID"-1]}
r=${b["SGE_TASK_ID"-1]}
u=${c["SGE_TASK_ID"-1]}

trim_galore --dont_gzip --length 50 --retain_unpaired --output_dir ./trimmed_unpaired --paired $f $r
trim_galore --dont_gzip --length 50 --output_dir ./trimmed_unpaired $u
