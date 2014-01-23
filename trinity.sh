#!/bin/bash
#$ -q genomics
#$ -j y
#$ -cwd
#$ -l h_vmem=200G
#$ -l virtual_free=200G
#$ -N assemble_californicus
. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

Trinity.pl --seqType fq --JM 190G --left ./data/trimmed_1.fq  --right ./data/trimmed_2.fq  --single ./data/trimmed_u.fq --CPU 12