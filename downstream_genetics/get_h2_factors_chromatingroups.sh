#!/bin/bash
#$ -cwd
#$ -j y
#$ -o logs/h2_chromatingroups.log
#$ -l h_vmem=16g
#$ -l h_rt=48:00:00

source /broad/software/scripts/useuse
use Anaconda

source activate ldsc

files=$(ls *.sumstats.gz)
for file in $files; do
	python ~/ldsc/ldsc.py --h2-cts $file --ref-ld-chr /stanley/robinson/ccarey/ldsc_reference/1000G_EUR_Phase3_baseline/baseline. --out h2_${file%%.*}_chromatingroups --ref-ld-chr-cts /stanley/robinson/ccarey/roadmap_cell_type_groups/Roadmap_ctg.ldcts --w-ld-chr /stanley/robinson/ccarey/ldsc_reference/weights_hm3_no_hla/weights.
done
