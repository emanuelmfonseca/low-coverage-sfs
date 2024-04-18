#!/bin/bash

#SBATCH --job-name=vcftools_fill_aa
#SBATCH --output=outfiles/%x-%j.out
#SBATCH --error=outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:05:00

# load required modules
module load samtools
module load vcftools

# need to have this path locally for perl to work properly for vcftools
# otherwise get this error: Can't locate Vcf.pm in @INC...
export PERL5LIB=/xdisk/rgutenk/lnt/projects/low_coverage/software/vcftools-vcftools-d511f46/src/perl/

# variables for file paths
aa_ref_dir="../reference/homo_sapiens_ancestor_GRCh38"
aa_ref_file="homo_sapiens_ancestor_20.fa"

for depth in 3x 5x 10x 30x
do
    data_vcf_dir="gatk_combined_vcf_${depth}"
    for file in $(find ${data_vcf_dir} -type f -name "all_sample_*_only_SNP_filtered.vcf"); do
        # get filename
        data_vcf_file="$(basename "$file")"
        # fill-aa for each file
        cat $file | fill-aa -a "${aa_ref_dir}/${aa_ref_file}.gz" | gzip -c > "${data_vcf_dir}/aa_annotated_${data_vcf_file}.gz"
    done
done
