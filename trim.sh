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

a=(noadapt/*1.fastq)
b=(noadapt/*2.fastq)
base=$(basename ${a["SGE_TASK_ID"-1]} "_1.fastq")
f=${a["SGE_TASK_ID"-1]}
r=${b["SGE_TASK_ID"-1]}

trim_galore --length 50 --retain_unpaired --output_dir ./trimmed --paired $f $r

