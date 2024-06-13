#!/bin/bash

#SBATCH --job-name=vcftools_process_aa_ref
#SBATCH --output=outfiles/%x-%j.out
#SBATCH --error=outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:05:00

# load required modules
module load samtools

# variables for file paths
aa_ref_dir=homo_sapiens_ancestor_GRCh38
aa_ref_file=homo_sapiens_ancestor_20.fa

### Prep reference ancestral allele file for vcftools
# change the first line in ref fasta to ">20" instead of "ANCESTOR_for_chromosome..."
cat ${aa_ref_dir}/${aa_ref_file} | sed 's,^>.*,>20,' | bgzip > ${aa_ref_dir}/${aa_ref_file}.gz
# generate index for the processed reference file
samtools faidx ${aa_ref_dir}/${aa_ref_file}.gz
# test to make sure the above works as intended, shoulve give a TGAC string to output
samtools faidx ${aa_ref_dir}/${aa_ref_file}.gz 20:1000000-1000020