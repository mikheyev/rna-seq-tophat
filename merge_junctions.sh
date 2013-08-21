#!/bin/bash
#$ -q short
#$ -j y
#$ -cwd
#$ -l h_vmem=10G
#$ -l virtual_free=10G
#$ -N junctions
. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp


cat rnaseq/*/junctions.bed |grep -v "track name=" | bed_to_juncs | sort -k 1,4 -u | sort -k 1,1 > junctions.merged.junc 
