#!/bin/bash
#$ -q short
#$ -j y
#$ -cwd
#$ -l h_vmem=10G
#$ -l virtual_free=10G
#$ -N merge
. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

#novosort -o merged.bam -t /genefs/MikheyevU/temp -c 2 rnaseq*/*/accepted_hits.bam

for i in rnaseq/*/accepted_hits.bam
do
	samtools view -c $i
done

for i in rnaseq/*/unmapped.bam
do
	samtools view -c $i
done