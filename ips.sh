#!/bin/bash
#$ -N scan
#$ -q long
#$ -j y
#$ -cwd
#$ -l h_vmem=5G
#$ -l virtual_free=5G
. $HOME/.bashrc
a=(`pwd`/data/blast/*.fa) #673
file=${a["$SGE_TASK_ID"-1]}
export TMP=/genefs/MikheyevU/temp

file_loc=`pwd`/data/blast
cd /apps/MikheyevU/iprscan/interproscan-5-44.0/
java -Xmx4G -jar interproscan-5.jar -t n -dp --goterms --iprlookup -f XML -T $TMP -o $file_loc/$(basename $file ".fa")"_ips.xml" -i $file
