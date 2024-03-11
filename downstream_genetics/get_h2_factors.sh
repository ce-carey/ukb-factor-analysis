#!/bin/bash
#$ -R y
#$ -cwd
#$ -j y
#$ -o logs/h2.log
#$ -l h_vmem=64g
#$ -l h_rt=12:00:00
#$ -hold_jid 28291556,28291611

source /broad/software/scripts/useuse
use Python-3.6
use Anaconda

source activate ldsc

files=$(ls *.sumstats.gz)
for file in $files; do
        python ~/ldsc/ldsc.py --chisq-max 999999.0 --h2 ${file%%.*}.sumstats.gz --ref-ld-chr /humgen/atgu1/fs03/shared_resources/ldsc_reference/baselineLD_v1.1/baselineLD. --out h2_${file%%.*} --overlap-annot --frqfile-chr /humgen/atgu1/fs03/shared_resources/ldsc_reference/1000G_Phase3_frq/1000G.EUR.QC. --w-ld-chr /humgen/atgu1/fs03/shared_resources/ldsc_reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --print-coefficients
done

