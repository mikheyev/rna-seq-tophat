#!/bin/bash
#$ -q long
#$ -j y
#$ -cwd
#$ -l h_vmem=10G
#$ -l virtual_free=10G
#$ -N rs
. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

a=(data/trimmed/*val_1.fq) #44
b=(data/trimmed/*val_2.fq)
f=${a["SGE_TASK_ID"-1]}
r=${b["SGE_TASK_ID"-1]}
base=`basename $f _1_val_1.fq`

rsem-calculate-expression -p 4 --paired-end $f $r  data/trinity_out_dir/Trinity_noERCC $base

