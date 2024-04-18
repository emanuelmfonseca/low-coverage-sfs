#!/bin/bash

#SBATCH --job-name=make_genomics_map
#SBATCH --output=outfiles/%x-%j.out
#SBATCH --error=outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:05:00

mkdir -p temp
for depth in 3x 5x 10x 30x
do
    # make genomic map
    module load bcftools
    for i in gatk_called_vcf_${depth}/NA*.vcf.gz
    do 
        echo `bcftools query -l $i`
        echo $i
    done | paste - - > "all_sample_map_${depth}"
done