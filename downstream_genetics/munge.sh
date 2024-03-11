#!/bin/bash
#$ -R y
#$ -cwd
#$ -j y
#$ -o logs/munge.log
#$ -l h_vmem=64g
#$ -l h_rt=12:00:00
#$ -hold_jid 28291556

source /broad/software/scripts/useuse
use Python-3.6
use Anaconda

source activate ldsc

files=$(ls *.tsv)
for file in $files; do
        python ~/ldsc/munge_sumstats.py --out ${file%%.*} --merge-alleles /humgen/atgu1/fs03/shared_resources/ldsc_reference/w_hm3.snplist --sumstats $file
done


