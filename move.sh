#!/bin/bash
#$ -q short
#$ -j y
#$ -cwd
#$ -l h_vmem=4G
#$ -l virtual_free=4G
#$ -N move
. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

for i in /genefs/SQC/Illumina/130809_SN1077_0126_AD2A38ACXX/Unaligned_Index/Project_MikheyevU/Sample_??_*
do 
    dname=$(echo $(basename $i) | sed 's/Sample_//')
    mkdir raw/$dname
    gunzip -c $i/*R1*gz > raw/$dname/$dname"_1.fastq"
    gunzip -c $i/*R2*gz > raw/$dname/$dname"_2.fastq"
done

