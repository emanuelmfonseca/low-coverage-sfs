#!/bin/bash

#SBATCH --job-name=combine_gvcf_genomicsDB_chr20
#SBATCH --output=outfiles/combine_gvcf/%x-%A-%a.out
#SBATCH --error=outfiles/combine_gvcf/%x-%A-%a.err
#SBATCH --account=rgutenk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00
#SBATCH --array=0-3

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load bcftools
module load gatk

mkdir -p temp

i=${SLURM_ARRAY_TASK_ID}
depths=("3x" "5x" "10x" "30x")
depth=${depths[${i}]}

# combine gvcf
gatk --java-options "-Xmx4g -Xms4g" \
   GenomicsDBImport \
   --genomicsdb-workspace-path my_database_chr20_${depth} \
   --sample-name-map all_sample_map_${depth} \
   -L chr20 \
   --tmp-dir temp \
   --reader-threads 20
   
# filter final gvcf

# make dir for final outputs
res_dir="gatk_combined_vcf_${depth}"
mkdir -p $res_dir

# reference genome
ref_path="../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"

gatk --java-options "-Xmx4g" GenotypeGVCFs \
    -R $ref_path \
    -V gendb://my_database_chr20_${depth} \
    -O ${res_dir}/all_sample_gatk_chr20.vcf.gz \
    --tmp-dir temp

# select only SNP
gatk SelectVariants \
    -R $ref_path \
    -V ${res_dir}/all_sample_gatk_chr20.vcf.gz \
    --select-type-to-include SNP \
    -O ${res_dir}/all_sample_gatk_chr20_only_SNP.vcf
 
# filter snps
gatk --java-options "-Xmx4G" VariantFiltration \
    -R $ref_path \
    -V ${res_dir}/all_sample_gatk_chr20_only_SNP.vcf \
    --filter-name 'QD_filter' --filter-expression 'QD<2.0' \
    --filter-name 'FS_filter' --filter-expression 'FS>60.0' \
    --filter-name 'MQ_filter' --filter-expression 'MQ<40.0' \
    --filter-name 'SOR_filter' --filter-expression 'SOR>10.0' \
    -O ${res_dir}/all_sample_gatk_chr20_only_SNP_filtered.vcf
