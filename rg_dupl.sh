#!/bin/bash
#$ -q long
#$ -j y
#$ -cwd
#$ -l h_vmem=10G
#$ -l virtual_free=10G
#$ -N rg
. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

a=(tophat/*/accepted_hits.bam)
ref=ref_spike/ref_spike.fa
f=${a["SGE_TASK_ID"-1]} 
base=$(echo $f | cut -d "/" -f2)

alias GA="java -Xmx5g -jar /apps/MikheyevU/sasha/GATK/GenomeAnalysisTK.jar"

java -Xmx9g -Djava.io.tmpdir=/genefs/MikheyevU/temp -jar /apps/MikheyevU/picard-tools-1.66/AddOrReplaceReadGroups.jar INPUT=$f OUTPUT=$(dirname $f)/rg.bam RGID=$base RGLB=$base RGPL=ILLUMINA RGSM=$base RGPU=$base

java -Xmx9g -Djava.io.tmpdir=/genefs/MikheyevU/temp -jar /apps/MikheyevU/picard-tools-1.66/MarkDuplicates.jar METRICS_FILE=$base"_duplicates.txt" ASSUME_SORTED=1 INPUT=$(dirname $f)/rg.bam OUTPUT=$(dirname $f)/nodup.bam MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=512 TMP_DIR=/genefs/MikheyevU/temp VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=1 COMPRESSION_LEVEL=0


